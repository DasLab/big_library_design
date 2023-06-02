import numpy as np
import pandas as pd
import os
from math import floor
from itertools import product, chain
from random import shuffle, sample, choices, random
from tqdm import tqdm
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
from arnie.bpps import bpps
from arnie.utils import convert_dotbracket_to_bp_list
from Levenshtein import distance as edit_distance
from pybktree import BKTree


'''
Example construct:
     26 nt.   5' fixed sequence
    240 nt.   Region of interest
     20 nt.   8 bp barcode with UUCG tetraloop (65,536 possible barcodes)
      1 nt.   Single A spacer
     20 nt.   3' fixed sequence
'''
BASES = ['A', 'C', 'G', 'T']
SEQ5 = 'GGGAACGACTCGAGTAGAGTCGAAAA'
SEQ3 = 'AAAAGAAACAACAACAACAAC'  # with the single A spacer added already
TETRALOOP = 'TTCG'

'''
Example procedure:
    for each WT, get single mutants and souble mutants of interaction
    generate pad to make all equal length
    each sequence add 5', random 8bp barcode, 3'
Example structure checks:

'''

MAXPROB_NONINTERACT = 0.09
MINPROB_PAIRED = 0.75
MINAVGPROB_PAIRED = 0.85
MINPROB_UNPAIRED = 0.75
MINAVGPROB_UNPAIRED = 0.85


###############################################################################
# fasta file manipulations
###############################################################################


def combine_fastas(fastas, out_fasta):
    '''
    Given list of fastas, appends all and saves new fasta.
    '''

    all_seqs = []
    for fasta in fastas:
        all_seqs.extend(list(SeqIO.parse(fasta, "fasta")))

    print(f'Combined and saved {fastas} to {out_fasta}')

    SeqIO.write(all_seqs, out_fasta, "fasta")


def format_fasta_for_submission(fasta, out_file, file_format='twist'):
    '''
    From a fasta, save sequences in a format for library submission
    For agilent and twist: csv
    For custom_array: txt with each line a new sequence
    '''

    if file_format == 'twist' or file_format == 'agilent':

        if out_file[-4:] != '.csv':
            print('twist requires .csv files, change out_file')

        names = []
        seqs = []
        for seq in SeqIO.parse(fasta, "fasta"):
            names.append(seq.name)
            seqs.append(_get_dna_from_SeqRecord(seq))
        df = pd.DataFrame(np.array([names, seqs]).T,
                          columns=["name", "sequence"])

        df.to_csv(out_file, index=False)
        print(f'Written {out_file} for submission to twist or agilent.')

    elif file_format == 'custom_array':

        if out_file[-4:] != '.txt':
            print('custom_array requires .txt files, change out_file')

        with open(out_file, 'w') as f:
            for seq in SeqIO.parse(fasta, "fasta"):
                f.write(_get_dna_from_SeqRecord(seq))
                f.write('\n')
            print(f'Written {out_file} for submission to custom_array.')

    else:
        print('file_format not supported, available: custom_array twist agilent')


def randomly_select_seqs(fasta, reject_file, prop):
    '''
    Given a fasta randomly select N sequences and save to new fasta
    '''

    all_seqs = list(SeqIO.parse(fasta, "fasta"))
    N = round(len(all_seqs)*prop)

    if len(all_seqs) > N:
        seq_indices = list(range(len(all_seqs)))
        shuffle(seq_indices)
        pass_ind = seq_indices[:N]
        rejected_ind = seq_indices[N:]
        pass_ind.sort()
        rejected_ind.sort()
        pass_seqs = [seq for i, seq in enumerate(all_seqs) if i in pass_ind]
        rejected_seqs = [seq for i, seq in enumerate(
            all_seqs) if i in rejected_ind]
        SeqIO.write(rejected_seqs, reject_file, "fasta")
        SeqIO.write(pass_seqs, fasta, "fasta")
        print(f'Written {fasta} for with {N} selected sequences and {reject_file} with {len(all_seqs)-N}.')

    else:

        SeqIO.write(all_seqs, out_fasta, "fasta")
        print(f'WARNING: Written {fasta} but had less than {N} sequences so none removed.')


def get_same_length(fasta):
    '''
    Checks if all sequences in a given fasta are the same length
    '''

    seqs = list(SeqIO.parse(fasta, "fasta"))
    return _get_same_length(seqs)[0]


def check_sequences_contents(fasta,seq5=SEQ5, seq3=SEQ3,bases=BASES):
    problems = []
    names = []
    sequences = []
    seqs = list(SeqIO.parse(fasta, "fasta"))
    for seq_rec in seqs:
        seq = _get_dna_from_SeqRecord(seq_rec)
        if seq in sequences:
            problems.append(f"{seq_rec.name} is a repeated sequence.")
        if seq_rec.name in names:
            problems.append(f"{seq_rec.name} is a repeated name.")
        if seq[:len(seq5)] != seq5:
            problems.append(f"{seq_rec.name} has incorrect 5' sequence.")
        if seq[-len(seq3):] != seq3:
            problems.append(f"{seq_rec.name} has incorrect 3' sequence.")
        for n in seq:
            if n not in bases:
                problems.append(f"{seq_rec.name} has non {bases} base.")
        sequences.append(seq)
        names.append(seq_rec.name)
    return len(problems)==0, problems

def remove_seqs_already_in_other_file(fasta, other_fasta, out_file):
    '''
    Given 2 fasta files, removes any sequence in the first that
    has same name as a sequence in the second, save non-duplicated
    sequences of the first.
    DOES NOT check seq itself and does not change the second fasta!
    '''

    all_seqs = list(SeqIO.parse(fasta, "fasta"))
    other_seqs = list(SeqIO.parse(other_fasta, "fasta"))
    good_seqs = _remove_seqs_in_other_list(all_seqs, other_seqs)
    SeqIO.write(good_seqs, out_file, "fasta")


def get_used_barcodes(fasta, start, end):
    '''
    from a fasta file return all sequences between start and end inclusive
    '''

    all_seqs = list(SeqIO.parse(fasta, "fasta"))
    barcodes = []
    for record in all_seqs:
        seq = _get_dna_from_SeqRecord(record)
        # end is inclusive
        barcode = seq[start:end+1]
        barcodes.append(str(barcode))
    return barcodes


###############################################################################
# other utils
###############################################################################


def get_bp_set_from_dotbracket(dotbracket):
    '''
    Given a dotbracket structure, return a list of base-pairs
    IGNORES pseudoknots
    '''

    return convert_dotbracket_to_bp_list(dotbracket)


###############################################################################
# get desired sequences
###############################################################################


def _fill_in_any_incomplete(seq,seqs):
    incomplete_seq = {'N':['A','C','T','G'],
                        'R':['A','G'],
                        'Y':['C','T'],
                        'S':['C','G'],
                        'W':['A','T'],
                        'K':['T','G'],
                        'M':['A','C'],
                        'B':['C','T','G'],
                        'D':['A','T','G'],
                        'H':['A','C','T'],
                        'V':['A','C','G'],}
    if seq == '':
        return seqs
    elif seq[0] in incomplete_seq:
        new_seqs = []
        for n in incomplete_seq[seq[0]]:
            potential_seq = n+seq[1:]
            potential_seqs = [s[:-len(potential_seq)]+potential_seq for s in seqs]
            new_seqs.extend(potential_seqs)
        return _fill_in_any_incomplete(seq[1:],new_seqs)
    else:
        return _fill_in_any_incomplete(seq[1:],seqs)


def get_windows(fasta, window_length, window_slide, out_fasta=None,
                circularize=False, fraction_use=1, reverse_complement=False):
    '''
    Get sliding windows from an inputted fasta file.

    Args:
        fasta (str): fasta file containing sequence to get windows of
        window_length (int): length (in number nt) of desired windows
        window_slide (int): number nucleotide to slide between windows
        out_fasta (str): if specified, save windows to fasta (default None)
        circularize (bool): whether to circularize the sequence thus once it
            reaches the end, create windows with 3'-5' connected (default False)
        fraction_use (float): proportion (from start) of genome to use, 
            defaults all (default 1)

    Returns:
        list of SeqRecord with the windows
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_start-end with start and
        end inclusive and indexing from 0.
    '''

    print(f'Getting all sliding windows, {window_length}nt every {window_slide}nt.')

    # get sequences and initialize
    seqs = list(SeqIO.parse(fasta, "fasta"))
    windows = []
    unused_windows = []

    for seq_rec in seqs:

        # loop through sequence every window_slide nucleotides
        seq = _get_dna_from_SeqRecord(seq_rec)
        window_limit = floor(fraction_use*len(seq))

        for i in range(0, len(seq), window_slide):

            # when we hit the end of sequence
            if i+window_length > len(seq):

                # add last window
                if not circularize:
                    a, b = len(seq)-window_length, len(seq)
                    new_seq = seq[a:b]
                    namenum = f'{a}-{b-1}'
                    new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                              f'{seq_rec.name}_{namenum}',
                                              '', '')
                    windows.append(new_rec)
                    break

                # or circularize and add from 5' until very end
                else:
                    a, b, c = i, len(seq), window_length-len(seq)+i
                    new_seq = seq[a:b]+seq[:c]
                    namenum = f'{a}-{c-1}'
            # otherwise just add window as normal
            else:
                a, b = i, i+window_length
                new_seq = seq[a:b]
                namenum = f'{a}-{b-1}'
          

            new_seqs = _fill_in_any_incomplete(new_seq,[new_seq])
            for j,new_seq in enumerate(new_seqs):
                if len(new_seqs) == 1:
                    name = f'{seq_rec.name}_{namenum}'
                else:
                    name = f'{seq_rec.name}_amb{j}_{namenum}'
                # save with name inclusive!
                new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                          name, '', '')
                if fraction_use != 1 and ((i + window_length-1) > window_limit):
                    unused_windows.append(new_rec)
                else:
                    windows.append(new_rec)
                if reverse_complement:
                    rc_rec = SeqIO.SeqRecord(Seq.Seq(_get_reverse_complement(new_seq)),
                                          f'{name}_rc', '', '')
                    if fraction_use != 1 and ((i + window_length-1) > window_limit):
                        unused_windows.append(rc_rec)
                    else:
                        windows.append(rc_rec)


    # remove and save fraction unused
    if fraction_use != 1:
        unused_file = f'{out_fasta.rsplit(".",1)[0]}_unused.{out_fasta.rsplit(".",1)[1]}'
        SeqIO.write(unused_windows, unused_file, "fasta")
        print(f'Saved unused windows to {unused_file}.')

    # save file
    if out_fasta is not None:
        SeqIO.write(windows, out_fasta, "fasta")
        print(f'Saved windows to {out_fasta}.')

    # return list of windows
    return windows


def get_all_single_mutants(fasta, out_fasta=None, bases=BASES):
    '''
    Get all single mutants from sequences in a fasta file.

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        out_fasta (str): if specified, save mutants to fasta (default None)
        bases (list): what bases can be selected as mutants (default ['A', 'C', 'G', 'T'])

    Returns:
        list of SeqRecord with the mutants
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_#wt-mutant eg for sequence P4P6
        mutant with 138 mutated from G to A is: P4P6_138G-A
        Indexing from 0.
        Generally, total is 3*length seq
    '''

    print("Getting all single mutants.")

    # get all sequences and initialize
    all_WT = SeqIO.parse(fasta, "fasta")
    all_single_mutants = []

    for record in all_WT:
        seq = _get_dna_from_SeqRecord(record)

        # at each position, get single mutants
        for i in range(len(seq)):
            for mut in bases:
                if mut != seq[i]:
                    name = f' {record.id}_{i}{seq[i]}-{mut}'
                    new_seq = seq[:i]+mut+seq[i+1:]
                    new_mut = SeqIO.SeqRecord(Seq.Seq(new_seq), name, '', '')
                    all_single_mutants.append(new_mut)

    # save file
    if out_fasta is not None:
        SeqIO.write(all_single_mutants, out_fasta, "fasta")
        print(f'Saved single mutants to {out_fasta}.')
    return all_single_mutants


def get_all_double_mutants(fasta, regionAs, regionBs,
                           out_fasta=None, bases=BASES):
    '''
    Get all double mutants between specified region for all sequences in a fasta file.

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        regionAs (list of lists of int): for each sequence, a list of indices specifying first region
        regionBs (list of lists of int): for each sequence, a list of indices specifying second region
        out_fasta (str): if specified, save mutants to fasta (default None)
        bases (list): what bases can be selected as mutants (default ['A', 'C', 'G', 'T'])

    Returns:
        list of SeqRecord with the mutants
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_#wt-mutant eg for sequence P4P6
        mutant with 138 mutated from G to A and 150 C to T is: P4P6_138G-A_150C-T
        Indexing from 0 on mutants are ordered by nucleotide number.
        Generally, total is 9*lengthA*lengthB
    '''

    print("Getting all double mutants.")

    # get all sequences and initialize
    all_WT = list(SeqIO.parse(fasta, "fasta"))
    all_double_mutants = []
    print(regionAs,regionBs)

    # check user specified regions for every sequence in the fasta
    if len(all_WT) != len(regionAs) or len(all_WT) != len(regionBs):
        print(f'WARNING: you must list regions (as a list of lists) for all sequences in the fasta, number sequences in fasta {len(all_WT)} and number of regions specified {len(regionAs)} {len(regionBs)}')

    for record, regionA_s, regionB_s in zip(all_WT, regionAs, regionBs):

        seq = _get_dna_from_SeqRecord(record)

        for regionA,regionB in zip(regionA_s,regionB_s):
            # at each position pair, get double mutants
            for i in regionA:
                for mutA in bases:
                    if mutA != seq[i]:
                        for j in regionB:
                            for mutB in bases:
                                if mutB != seq[j]:

                                    # get mutant
                                    # name according to convention, in order, index at 1
                                    if i == j:
                                        continue
                                    elif i < j:
                                        new_seq = seq[:i]+mutA + \
                                            seq[i+1:j]+mutB+seq[j+1:]
                                        name = f' {record.id}_{i}{seq[i]}-{mutA}_{j}{seq[j]}-{mutB}'
                                    else:
                                        new_seq = seq[:j]+mutB + \
                                            seq[j+1:i]+mutA+seq[i+1:]
                                        name = f' {record.id}_{j}{seq[j]}-{mutB}_{i}{seq[i]}-{mutA}'
                                    new_mut = SeqIO.SeqRecord(
                                        Seq.Seq(new_seq), name, '', '')
                                    all_double_mutants.append(new_mut)

    # save file
    if out_fasta is not None:
        SeqIO.write(all_double_mutants, out_fasta, "fasta")
        print(f'Saved all double mutants between the 2 regions to {out_fasta}.')

    return all_double_mutants


def get_wcf_rescue_mutants(fasta, bp_sets, out_fasta=None,
                           wfc_base_pairs=['AT', 'TA', 'CG', 'GC']):
    '''
    Get all Watson-Crick-Franklin rescue mutants between specified base-pairs for all sequences in a fasta file.

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        bp_sets (list of lists of pairs): for each sequence, a list of base-pairs as a tuple or list of int
        out_fasta (str): if specified, save mutants to fasta (default None)
        wfc_base_pairs (list): list of potential base pairs to rescue with
            (default ['AT', 'TA', 'CG', 'GC'])

    Returns:
        list of SeqRecord with the mutants
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_#wt-mutant eg for sequence P4P6
        mutant with 138 mutated from G to A and 150 C to T is: P4P6_138G-A_150C-T
        Indexing from 0 on mutants are ordered by nucleotide number.
        Generally, total is 3*numbps
    '''

    print("Getting all rescue mutants.")

    # get all sequences and initialize
    all_WT = list(SeqIO.parse(fasta, "fasta"))
    all_rescue_mutants = []

    # check user specified set of basepairs for every sequence in the fasta
    if len(all_WT) != len(bp_sets):
        print(f'WARNING: bps must be a list, one for each sequence in fasta, of lists of basepairs to rescue for that sequence. You have {len(all_WT)} inputted sequences and {len(bp_sets)} base-pair sets.')

    for record, bps in zip(all_WT, bp_sets):

        seq = _get_dna_from_SeqRecord(record)

        # at each base pair, get rescue mutants
        for bp in bps:
            current_bp = seq[bp[0]]+seq[bp[1]]
            for new_bp in wfc_base_pairs:
                if new_bp != current_bp:
                    # get mutant
                    # name according to convention, in order, index at 1
                    if bp[0] < bp[1]:
                        name = f' {record.id}_{bp[0]}{seq[bp[0]]}-{new_bp[0]}_{bp[1]}{seq[bp[1]]}-{new_bp[1]}'
                        new_seq = seq[:bp[0]]+new_bp[0] + \
                            seq[bp[0]+1:bp[1]]+new_bp[1]+seq[bp[1]+1:]
                    else:
                        name = f' {record.id}_{bp[1]}{seq[bp[1]]}-{new_bp[1]}_{bp[0]}{seq[bp[0]]}-{new_bp[0]}'
                        new_seq = seq[:bp[1]]+new_bp[1] + \
                            seq[bp[1]+1:bp[0]]+new_bp[0]+seq[bp[0]+1:]
                    new_mut = SeqIO.SeqRecord(
                        Seq.Seq(new_seq), name, '', '')
                    all_rescue_mutants.append(new_mut)

    # save file
    if out_fasta is not None:
        SeqIO.write(all_rescue_mutants, out_fasta, "fasta")
        print(f'Saved all rescue mutants to {out_fasta}.')

    return all_rescue_mutants


###############################################################################
# add library parts
###############################################################################


def add_known_pads(fasta, pad5_dict, pad3_dict, out_fasta=None):
    '''
    From fasta prepend and append given sequences

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        pad5_dict (dict int:str): sequence length : pad to add on 5'
        pad3_dict (dict int:str): sequence length : pad to add on 3'
        out_fasta (str): if specified, save mutants to fasta (default None)

    Returns:
        list of SeqRecord with pads added
        # pad# where numbers are length of pad
        if out_fasta specified, also saves these to fasta file naming
    '''

    # get all sequences and initialize
    all_WT = list(SeqIO.parse(fasta, "fasta"))
    all_seqs = []

    for record in all_WT:
        seq = _get_dna_from_SeqRecord(record)

        # get pad for this length
        pad5 = pad5_dict[len(seq)]
        pad3 = pad3_dict[len(seq)]
        name = f' {record.id}_{len(pad5)}pad{len(pad3)}'
        new_seq = SeqIO.SeqRecord(Seq.Seq(pad5+seq+pad3), name, '', '')
        all_seqs.append(new_seq)

    # save file
    if out_fasta is not None:
        SeqIO.write(all_seqs, out_fasta, "fasta")
        print(f'Saved all with correct constant pad added to {out_fasta}.')

    return all_seqs


def get_all_barcodes(out_fasta=None, num_bp=8, num5hang=0, num3hang=0,
                     polyA5=0, polyA3=0,
                     loop=TETRALOOP, bases=BASES):
    '''
    Return all barcodes of specified structure

    Args:
        out_fasta (str): if specified, save mutants to fasta (default None)
        num_bp (int): length of stem (default 8)
        num5hang (int): length of random 5' hang (default 0)
        num3hang (int): length of random 3' hang (default 0)
        polyA5 (int): length of polyA 5' hang (placed before random) (default 0)
        polyA3 (int): length of polyA 3' hang (placed after random) (default 0)
        loop (str): sequence of loop
        bases (list): list of bases that can be used

    Returns:
        list of SeqRecord of all possible barcodes
        if out_fasta specified, also saves these to fasta file

    '''

    print("Getting all possible barcodes.")
    all_barcodes = []

    # get all possible combinations of bases for random/barcode regions
    for x in product(bases, repeat=num5hang+num_bp+num3hang):
        uid = ''.join(x)

        # split barcode in stem and hang regions
        hang5 = uid[:num5hang]
        if num3hang == 0:
            hang3 = ''
            stemA = uid[num5hang:]
        else:
            hang3 = uid[-num3hang:]
            stemA = uid[num5hang:-num3hang]
        stemB = _get_reverse_complement(stemA)

        # put all barcode parts together
        seq = ("A"*polyA5)+hang5+stemA+loop+stemB+hang3+("A"*polyA3)
        name = f' stem{stemA}_{hang5}hang{hang3}_{polyA5}polyA{polyA3}'
        seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), name, '', '')
        all_barcodes.append(seq_rec)

    # save
    if out_fasta is not None:
        SeqIO.write(all_barcodes, out_fasta, "fasta")
        print(f'Saved all barcodes to {out_fasta}.')

    return all_barcodes


def add_pad(fasta, out_fasta, bases=BASES, share_pad='same_length',
            epsilon_punpaired=MINPROB_UNPAIRED,
            epsilon_interaction=MAXPROB_NONINTERACT,
            epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
            epsilon_paired=MINPROB_PAIRED,
            epsilon_avg_paired=MINAVGPROB_PAIRED,
            loop=TETRALOOP, hang=3, polyAhang=0, min_num_samples=30,
            max_prop_bad=0.05, pad_side='both',
            min_length_stem=4, max_length_stem=12,
            num_pads_reduce=100, percent_reduce_prob=10,pad_to_length=None):
    '''
    Given a fasta of sequence pad all sequence to the same length

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        out_fasta (str): if specified, save mutants to fasta (default None)
        share_pad (str): either all different pads (none),
            all same for sequences of same length (same_length),
            all same where sequence of different length are just truncated (all)
        epsilon_interaction (float): Maximum base-pair-probability for 2 regions to be considered non-interacting
        epsilon_punpaired (float): Minimum probability unpaired for region to be unpaired
        epsilon_avg_punpaired (float): Average probability unpaired for region to be unpaired
        epsilon_paired (float): Minimum base-pair-probability for 2 regions to be considered paired
        epsilon_avg_paired (float): Average base-pair-probability for 2 regions to be considered paired
        loop (str): constant sequence for loop of stem-loop
        hang (int): distance between sequence and the stem will be this or this+1, random sequence (default 3)
        polyAhang (int): distance between sequence and the stem will be this or this+1, polyA (default 3)
        min_length_stem (int): minimum number base-pairs to form stem
        max_length_stem (int): maximum number of base-pairs to form stem
        min_num_samples (int): minimum number of samples to test pad structure (default 30)
        max_prop_bad (float): of the samples sequences, proportion that can be bad (default 0.05)
        pad_side (str): whether split pad ('both') or force to be exclusively 5' or 3'
        num_pads_reduce (int): number of pads to try before reducing probability (default 100)
        percent_reduce_prob (float): percent to reduce probabilities each time (default 10)

    Returns:
        list of SeqRecord with pads added
        # pad# where numbers are length of pad
        if out_fasta specified, also saves these to fasta file naming added _
    '''

    # get sequences and sort by length
    seqs = list(SeqIO.parse(fasta, "fasta"))
    seq_by_length = {}
    for seq_rec in seqs:
        len_seq = len(seq_rec.seq)
        if len_seq in seq_by_length:
            seq_by_length[len_seq].append(seq_rec)
        else:
            seq_by_length[len_seq] = [seq_rec]

    # for each length, randomly select 1% of sequences or min_num_samples
    selected_sec = []
    max_bad_structs = {}
    for len_seq, seq_group in seq_by_length.items():
        number_to_select = min(len(seq_group), max(
            min_num_samples, len(seq_group)*0.01))
        number_to_select = round(number_to_select)
        selected_sec.append(sample(seq_group, k=number_to_select))
        max_bad_structs[len_seq] = floor(max_prop_bad*number_to_select)
        print(f'Pad for length {len_seq} search using {number_to_select} sequences for structure check for each length.')

    # get max length of pad needed
    if pad_to_length is None:
        desired_len = max(seq_by_length.keys())
    else:
        desired_len = pad_to_length
    porp_reduce = (1-(percent_reduce_prob/100))

    # if want all pads to be random
    if share_pad == 'none':
        print("Getting random pads for each sequence, not guaranteed to be unique.")
        padded_seqs = []

        for seq_length, seqs in seq_by_length.items():
            # for each length get the structure of the 3' and 5' pad
            pad_length = desired_len-seq_length
            structs, regions = _get_5_3_split(pad_length, hang, polyAhang,
                                              min_length_stem, max_length_stem,
                                              pad_side, loop, seq_length)
            # for each sequence search for good_pad
            for seq in list(seqs):
                good_pad = False
                while not good_pad:
                    # get random pad
                    pad5 = _get_random_barcode(num_bp=structs["5"]['bp'],
                                               num3hang=structs["5"]['hang_rand'],
                                               polyA3=structs["5"]['hang_polyA'],
                                               loop=structs["5"]['loop'])
                    pad3 = _get_random_barcode(num_bp=structs["3"]['bp'],
                                               num5hang=structs["3"]['hang_rand'],
                                               polyA5=structs["3"]['hang_polyA'],
                                               loop=structs["3"]['loop'])
                    pad5 = _get_dna_from_SeqRecord(pad5)
                    pad3 = _get_dna_from_SeqRecord(pad3)

                    # get full sequence and check structure
                    full_seq = pad5+_get_dna_from_SeqRecord(seq)+pad3
                    full_seq_name = f'{seq.name}_{len(pad5)}pad{len(pad3)}'
                    struct_results = check_struct_bpp(full_seq,
                                                      regions['unpaired'], regions['pairedA'],
                                                      regions['pairedB'], regions['noninteractA'], regions['noninteractB'],
                                                      epsilon_interaction=epsilon_interaction,
                                                      epsilon_punpaired=epsilon_punpaired,
                                                      epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                      epsilon_paired=epsilon_paired,
                                                      epsilon_avg_paired=epsilon_avg_paired)
                    good_pad = struct_results["pass"]

                # if pass add final sequence
                padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
                                             full_seq_name, '', '')
                padded_seqs.append(padded_seq)

    # if want sequences of same length to share same pad
    elif share_pad == 'same_length':
        pads_by_len = {}
        for seqs in selected_sec:
            # get length of pad for this group, if 0 done
            len_group = len(str(seqs[0].seq))
            length_to_add = desired_len-len_group
            if length_to_add == 0:
                pads_by_len[len_group] = ['', '']
                continue

            print(f"Searching for pad for group len {len_group}")

            # split length of pad to either end of the sequence and get regions
            structs, regions = _get_5_3_split(length_to_add, hang, polyAhang,
                                              min_length_stem, max_length_stem,
                                              pad_side, loop, len_group)

            print(f"Finding a {structs['5']['N']}nt 5' pad with {structs['5']['bp']}bp stem {structs['5']['hang_rand']}nt random hang {structs['5']['hang_polyA']}nt polyA hang.")
            print(f"Finding a {structs['3']['N']}nt 3' pad with {structs['3']['bp']}bp stem {structs['3']['hang_rand']}nt random hang {structs['3']['hang_polyA']}nt polyA hang.")

            # loop through to find pad that works
            good_pad = False
            seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
            # current_pad = 0
            while not good_pad:

                # if tried enough relax the probability contstraints
                prob_factor = {}
                for type_error, count in seq_count.items():
                    mutiplier = (count // num_pads_reduce)
                    prob_factor[type_error] = porp_reduce**mutiplier
                    if (count-mutiplier) % num_pads_reduce == 0 and count != 0:
                        seq_count[type_error] += 1
                        print(f'For {len_group}, failed to find pad from {count-mutiplier} pads because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')

                # get random pad
                pad5 = _get_random_barcode(num_bp=structs["5"]['bp'],
                                           num3hang=structs["5"]['hang_rand'],
                                           polyA3=structs["5"]['hang_polyA'],
                                           loop=structs["5"]['loop'])
                pad3 = _get_random_barcode(num_bp=structs["3"]['bp'],
                                           num5hang=structs["3"]['hang_rand'],
                                           polyA5=structs["3"]['hang_polyA'],
                                           loop=structs["3"]['loop'])

                # chek all samples sequences
                bad_count = 0
                for i, seq in enumerate(seqs):
                    # if i % 50 == 0 and i != 0:
                    #    print(i)

                    # get full sequence and check its structure
                    full_seq = pad5 + _get_dna_from_SeqRecord(seq) + pad3
                    struct_results = check_struct_bpp(full_seq,
                                                      regions['unpaired'], regions['pairedA'],
                                                      regions['pairedB'], regions['noninteractA'], regions['noninteractB'],
                                                      epsilon_interaction=epsilon_interaction,
                                                      epsilon_punpaired=epsilon_punpaired,
                                                      epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                      epsilon_paired=epsilon_paired,
                                                      epsilon_avg_paired=epsilon_avg_paired,
                                                      prob_factor=prob_factor)

                    # if not good structure add to count
                    # stop if reached maximal bad
                    good_pad = struct_results["pass"]
                    if not good_pad:
                        seq_count = _update_seq_count(seq_count,struct_results)
                        bad_count += 1
                        if bad_count >= max_bad_structs[len_group]:
                            break

            pads_by_len[len_group] = [_get_dna_from_SeqRecord(
                pad5), _get_dna_from_SeqRecord(pad3)]

        print('Found pads, now adding pads.')
        padded_seqs = []
        for seq_length, seqs in seq_by_length.items():
            pad5, pad3 = pads_by_len[seq_length]
            for seq in list(seqs):
                full_seq = pad5+_get_dna_from_SeqRecord(seq)+pad3
                full_seq_name = f'{seq.name}_{len(pad5)}pad{len(pad3)}'
                padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
                                             full_seq_name, '', '')
                padded_seqs.append(padded_seq)

    # if want all sequences to share same pad (truncated as needed)
    elif share_pad == 'all':
        print(f"Searching for pad for group all sequences")
        print('WARNING depending on your set of sequences and structure desired this may be a near impossible random search space!')

        # get all pad structures for the pad to be shared across all
        pad_lengths = [desired_len-x for x in seq_by_length.keys()]
        all_structs, cutoff_by_len, regions_by_len = _get_5_3_split_multi(pad_lengths, hang, polyAhang,
                                                                          min_length_stem, max_length_stem,
                                                                          pad_side, loop, list(seq_by_length.keys()))

        # loop until find a good pad
        good_pad = False
        seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
        while not good_pad:

            # get random pad
            pad5_full = ''
            pad3_full = ''
            for structs in all_structs:
                pad5_full = pad5_full+_get_random_barcode(num_bp=structs["5"]['bp'],
                                                          num3hang=structs["5"]['hang_rand'],
                                                          polyA3=structs["5"]['hang_polyA'],
                                                          loop=structs["5"]['loop'])
                pad3_full = _get_random_barcode(num_bp=structs["3"]['bp'],
                                                num5hang=structs["3"]['hang_rand'],
                                                polyA5=structs["3"]['hang_polyA'],
                                                loop=structs["3"]['loop']) + pad3_full

            # for all sequences
            for seqs in selected_sec:
                bad_count = 0

                # if tried enough relax the probability constraints
                prob_factor = {}
                for type_error, count in seq_count.items():
                    mutiplier = (count // num_pads_reduce)
                    prob_factor[type_error] = porp_reduce**mutiplier
                    if (count-mutiplier) % num_pads_reduce == 0 and count != 0:
                        seq_count[type_error] += 1
                        print(f'For {len_group}, failed to find pad from {count-mutiplier} pads because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')

                # for each length get regions and cutoff of full pads
                len_group = len(seqs[0].seq)
                regions = regions_by_len[len_group]
                cutoff = cutoff_by_len[len_group]
                pad5 = _get_dna_from_SeqRecord(pad5_full[:cutoff[0]])
                if cutoff[1] == 0:
                    pad3 = ''
                else:
                    pad3 = _get_dna_from_SeqRecord(pad3_full[-cutoff[1]:])

                for i, seq in enumerate(seqs):
                    # if i%10==0 and i!=0:
                    #    print(i)

                    if len_group != desired_len:

                        # get full sequence and check its structure
                        full_seq = pad5 + _get_dna_from_SeqRecord(seq) + pad3
                        struct_results = check_struct_bpp(full_seq, regions['unpaired'], regions['pairedA'],
                                                          regions['pairedB'], regions['noninteractA'], regions['noninteractB'],
                                                          epsilon_interaction=epsilon_interaction,
                                                          epsilon_punpaired=epsilon_punpaired,
                                                          epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                          epsilon_paired=epsilon_paired,
                                                          epsilon_avg_paired=epsilon_avg_paired,
                                                          prob_factor=prob_factor)

                        # if not good structure add to count
                        # stop if reached maximal bad
                        good_pad = struct_results["pass"]
                        if not good_pad:
                            seq_count = _update_seq_count(seq_count,struct_results)
                            bad_count += 1
                            if bad_count >= max_bad_structs[len_group]:
                                break

        # add pads to the sequences
        padded_seqs = []
        for seq_rec in seqs:
            len_group = len(seq_rec.seq)
            cutoff = cutoff_by_len[len_group]
            if length_to_add != 0:
                pad5 = _get_dna_from_SeqRecord(pad5_full[:cutoff[0]])
                if cutoff[1] == 0:
                    pad3 = ''
                else:
                    pad3 = _get_dna_from_SeqRecord(pad3_full[-cutoff[1]:])
                seq = _get_dna_from_SeqRecord(seq_rec)
                full_seq = pad5 + seq + pad3
                name = f'{seq_rec.name}_{len(pad5)}pad{len(pad3)}'
                padded_seqs.append(SeqIO.SeqRecord(
                    Seq.Seq(full_seq), name, '', ''))
            else:
                padded_seqs.append(seq_rec)

    else:
        print("ERROR share_pad option not recognized.")

    # save
    if out_fasta is not None:
        SeqIO.write(padded_seqs, out_fasta, "fasta")
        print(f'Saved all padded sequences to {out_fasta}.')

    return padded_seqs


def get_edit_distance(barcodeA, barcodeB, num_bp=0, len_loop=0):
    if num_bp == 0:
        dist = edit_distance(barcodeA, barcodeB)
    else:
        start_stem = -len_loop-(2*num_bp)
        end_stem = -len_loop-num_bp
        stemA = barcodeA[start_stem:end_stem]
        stemB = barcodeB[start_stem:end_stem]
        otherA = barcodeA[:start_stem]+barcodeB[end_stem:-num_bp]
        otherB = barcodeB[:start_stem]+barcodeB[end_stem:-num_bp]
        dist = 2*edit_distance(stemA, stemB)
        dist += edit_distance(otherA, otherB)
    return dist


def add_fixed_seq_and_barcode(fasta, out_fasta=None, seq5=SEQ5, seq3=SEQ3,
                              num_bp=8, num5hang=0, num5polyA=4,
                              loop=TETRALOOP,
                              epsilon_interaction=MAXPROB_NONINTERACT,
                              epsilon_punpaired=MINPROB_UNPAIRED,
                              epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                              epsilon_paired=MINPROB_PAIRED,
                              epsilon_avg_paired=MINAVGPROB_PAIRED,
                              save_image_folder=None, save_bpp_fig=0,
                              punpaired_chunk_size=500, used_barcodes=[],
                              num_barcode_before_reduce=100,
                              percent_reduce_prob=10, min_edit=2):
    '''
    From a fasta of sequences, add constant regions and barcodes

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        out_fasta (str): if specified, save mutants to fasta (default None)
        seq5 (str): the constant sequence to place at 5' end of all sequences
        seq3 (str): the constant sequence to place at 3' end of all sequences
        num_bp (int): number of base pairs to have in barcode stem (default 8)
        num5hang (int): number of random nucleotides to put before the stem (default 0)
        num5polyA (int): number of A to put before the barcode (stem and 5 random hang if applicable) (default 4)
        loop (str): sequence of loop to have in the hairpin
        epsilon_interaction (float): Maximum base-pair-probability for 2 regions to be considered non-interacting
        epsilon_punpaired (float): Minimum probability unpaired for region to be unpaired
        epsilon_avg_punpaired (float): Average probability unpaired for region to be unpaired
        epsilon_paired (float): Minimum base-pair-probability for 2 regions to be considered paired
        epsilon_avg_paired (float): Average base-pair-probability for 2 regions to be considered paired
        save_image_folder (str): folder to save images to
        save_bpp_fig (float): proportion of sequences to save base-pair-probability matrix figure, 0 is none, 1 is all
        punpaired_chunk_size (int): max number of sequences to plot on each p unpaired plot (default 500)
        used_barcodes (list): list of sequences not to use as barcodes (default [])
        num_barcode_before_reduce (int): for each sequence, the number of barcodes to try
            before reducing the probability thresholds by 10% (default 100)
        percent_reduce_prob (float): percent to reduce probability threshold each time (default 10)
        min_edit (int): minimum edit distance of barcodes, base-pairs and num5hang included (default 2)

    Returns:
        list of SeqRecord which are library ready
        if out_fasta specified, also saves these to fasta file, _libraryready appended to names

    '''

    # if image folder does not exist create it
    if save_image_folder is not None:
        if not os.path.exists(save_image_folder):
            print(f'{save_image_folder} did not exists, creating.')
            os.makedirs(save_image_folder)

    # get and randomly shuffle all potential barcodes
    all_uids = get_all_barcodes(num_bp=num_bp, loop=loop, num5hang=num5hang,
                                polyA5=num5polyA)
    if len(used_barcodes) != 0:
         if len(all_uids[0]) != len(used_barcodes[0]):
             print('ERROR: usd barcodes are not the correct length')
    
    shuffle(all_uids)

    # read sequences, check all same length
    seqs = list(SeqIO.parse(fasta, "fasta"))
    same_length, seq_len = _get_same_length(seqs)
    if not same_length:
        print("ERROR: sequences not all same length.")

    # initialize variables
    all_full_seqs, p_unpaireds, rejected_uids = [], {}, []
    pad_lines, new_lines, seqs_for_labeling, muts = [], [], [], []
    current_uid, chunk_count, image_count = 0, 0, 0
    seq5 = _get_dna_from_SeqRecord(seq5)
    seq3 = _get_dna_from_SeqRecord(seq3)
    loop = _get_dna_from_SeqRecord(loop)
    pun_xlabel = [i if i % 10 == 0 else '' for i in range(seq_len+len(seq5)+len(seq3)+len(all_uids[0]))]
    porp_reduce = (1-(percent_reduce_prob/100))

    # find structural regions
    regions = [len(seq5), len(seq5)+seq_len,
               len(seq5)+seq_len+num5hang+num5polyA+(2*num_bp)+len(loop)]
    regionA = list(range(regions[0], regions[1]))
    regionB = list(range(regions[1], regions[2]))
    region_unpaired = [list(range(regions[1],
                                  regions[1]+num5hang+num5polyA)),
                       list(range(regions[2]-num_bp-len(loop),
                                  regions[2]-num_bp))]
    region_paired_A = list(range(regions[1]+num5hang+num5polyA,
                                 regions[2]-num_bp-len(loop)))
    region_paired_B = list(range(regions[2]-num_bp, regions[2]))[::-1]

    # initialize way to keep track of barcodes used
    if num_bp == 0 and min_edit > 1:
        added_barcodes = BKTree(edit_distance)
        for barcode in used_barcodes:
            added_barcodes.add(barcode)
    elif min_edit > 2:
        def edit_dist_stem(barcodeA, barcodeB): return get_edit_distance(
            barcodeA, barcodeB, num_bp, len(loop))
        added_barcodes = BKTree(edit_dist_stem)
        for barcode in used_barcodes:
            added_barcodes.add(barcode)
    elif (min_edit == 2 and num5hang != 0):
        added_barcodes = {}
        start_stem = -len(loop)-(2*num_bp)
        end_stem = -len(loop)-num_bp
        for barcode in used_barcodes:
            stemA = barcode[start_stem:end_stem]
            otherA = str(barcode[:start_stem]+barcode[end_stem:-num_bp])
            if stemA in added_barcodes:
                added_barcodes[stemA].append(otherA)
            else:
                added_barcodes[stemA] = [otherA]

    print("Adding 5', barcode, 3'.")

    for seq_rec in tqdm(list(seqs)):

        # get sequence and name
        seq = _get_dna_from_SeqRecord(seq_rec)
        name = seq_rec.name+'_libraryready'

        # initialize values
        uid_good, mutate_polyA = False, False
        seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
        chunk_count += 1
        lines = regions.copy()
        num_barcode = num_barcode_before_reduce

        # if has pad, add lines for padded region
        if 'pad' in name:
            pad5_length,pad3_length = _get_pad_from_name(name)
            new_lines = [len(seq5)+pad5_length,
                 len(seq5)+len(seq)-pad3_length]
            lines.extend(new_lines)
        else:
            pad5_length = 0
            new_lines = []
        pad_lines.append(new_lines)

        # renumber mutations because of 5' additions from the desired sequence
        mutations = _get_mutations_from_name(seq_rec.name)
        mutations = [x+len(seq5)+pad5_length for x in mutations]
        muts.append(mutations)

        # search barcodes until find one with correct structure predicted
        while not uid_good:

            # get barcode, get new if barcode already used
            uid = all_uids[current_uid].seq
            if (num_bp == 0 and min_edit > 1) or (min_edit > 2):
                close_barcodes = added_barcodes.find(uid, min_edit)
                while close_barcodes != []:
                    rejected_uids.append(all_uids[current_uid])
                    current_uid += 1
                    uid = all_uids[current_uid].seq
                    close_barcodes = added_barcodes.find(uid, min_edit)
            elif (min_edit == 2 and num5hang != 0):
                stem = uid[start_stem:end_stem]
                other = str(uid[:start_stem]+uid[end_stem:-num_bp])
                if stem in added_barcodes:
                    dists = [edit_distance(other, otherB)
                             for otherB in added_barcodes[stem]]
                    while min(dists) < min_edit:
                        rejected_uids.append(all_uids[current_uid])
                        current_uid += 1
                        uid = all_uids[current_uid].seq
                        stem = uid[start_stem:end_stem]
                        other = str(uid[:start_stem]+uid[end_stem:-num_bp])
                        if stem in added_barcodes:
                            dists = [edit_distance(other, otherB)
                                     for otherB in added_barcodes[stem]]
                        else:
                            dists = [2]
            else:
                while str(uid) in used_barcodes:
                    current_uid += 1
                    uid = all_uids[current_uid].seq

            # if have looped through enough times reduce probability thresholds
            prob_factor = {}
            for type_error, count in seq_count.items():
                mutiplier = (count // num_barcode)
                prob_factor[type_error] = porp_reduce**mutiplier

                if (count-mutiplier) % num_barcode == 0 and count != 0:
                    seq_count[type_error] += 1
                    print(f'{name}, failed to find barcode from {count-mutiplier} barcodes because of {type_error}, reducing probabilities by {prob_factor[type_error]}.')
                # when this is 3x, if there is a ployA allow this to mutate
                if mutiplier >= 3 and num5polyA != 0 and not mutate_polyA:
                    mutate_polyA = True
                    seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
                    num_barcode = num_barcode_before_reduce*3
                    print(f'{name}, failed to find barcodes now allowing polyA to be random sequence.')

            # create full sequence
            if mutate_polyA:
                new_hang = _get_random_barcode(
                    num_bp=0, num5hang=num5polyA, loop='').seq
                # new_loop = _get_random_barcode(num_bp=0, num5hang=len(loop),loop='').seq
                # uid = new_hang+uid[num5polyA:num5polyA+num_bp]+new_loop+uid[num5polyA+num_bp+len(loop):]
                uid = new_hang+uid[num5polyA:]
            full_seq = f'{seq5}{seq}{uid}{seq3}'

            # check if structure correct, save picture if specified and chance has it
            if (save_image_folder is not None) and (random() < save_bpp_fig):
                save_image = f'{save_image_folder}/{name}.png'
                plot_lines = lines
            else:
                save_image, plot_lines = None, None

            struct_results = check_struct_bpp(full_seq, region_unpaired,
                                              region_paired_A, region_paired_B,
                                              regionA, regionB,
                                              epsilon_interaction=epsilon_interaction,
                                              epsilon_punpaired=epsilon_punpaired,
                                              epsilon_avg_punpaired=epsilon_avg_punpaired,
                                              epsilon_paired=epsilon_paired,
                                              epsilon_avg_paired=epsilon_avg_paired,
                                              prob_factor=prob_factor,
                                              lines=lines,
                                              save_image=f'{save_image_folder}/{name}.png', mutants=mutations)
            uid_good = struct_results["pass"]

            # if barcode is incorrect structure loop through again
            if not uid_good:
                rejected_uids.append(all_uids[current_uid])
                seq_count = _update_seq_count(seq_count,struct_results)

            current_uid += 1

            # if looped through all barcodes, go back to previously rejected barcodes and continue loop
            if current_uid == len(all_uids):
                all_uids = rejected_uids
                current_uid, rejected_uids = 0, []

        # once found barcode add its sequence and probability unpaired to list
        seqs_for_labeling.append(full_seq)
        all_full_seqs.append(SeqIO.SeqRecord(Seq.Seq(full_seq), name, '', ''))
        p_unpaireds[name] = struct_results["p_unpaired"]
        if (num_bp == 0 and min_edit > 1) or (min_edit > 2):
            added_barcodes.add(uid)
        elif (min_edit == 2 and num5hang != 0):
            if stem not in added_barcodes:
                added_barcodes[stem] = [other]
            else:
                added_barcodes[stem].append(other)
        else:
            used_barcodes.append(uid)

        # when processed enough sequences
        # save probability unpaired figures and reset parameters
        if ((save_image_folder is not None) and
                (chunk_count == punpaired_chunk_size)):
            plot_punpaired(p_unpaireds, pun_xlabel,
                           seqs_for_labeling, muts,
                           regions, pad_lines,
                           f'{save_image_folder}/all_p_unpaired{image_count}.png')
            chunk_count, p_unpaireds = 0, {}
            pad_lines, seqs_for_labeling, muts = [], [], []
            image_count += 1

    # save left over probability unpaired figures
    if save_image_folder is not None and len(p_unpaireds) != 0:
        plot_punpaired(p_unpaireds, pun_xlabel,
                       seqs_for_labeling, muts,
                       regions, pad_lines,
                       f'{save_image_folder}/all_p_unpaired{image_count}.png')
    # save
    if out_fasta is not None:
        SeqIO.write(all_full_seqs, out_fasta, "fasta")
        print(f'Saved all full sequences to {out_fasta}.')

    return all_full_seqs


###############################################################################
# check structures
###############################################################################


def check_struct_bpp(seq, regions_unpaired=None, region_paired_A=None,
                     region_paired_B=None,
                     regionA=None, regionB=None,
                     epsilon_interaction=MAXPROB_NONINTERACT,
                     epsilon_punpaired=MINPROB_UNPAIRED,
                     epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                     epsilon_paired=MINPROB_PAIRED,
                     epsilon_avg_paired=MINAVGPROB_PAIRED,
                     prob_factor={'unpaired': 1,
                                  'paired': 1, 'interaction': 1},
                     lines=None, save_image=None, mutants=[]):
    '''
    Check if sequence as a base-pair-probaility matrix of desired structure

    Args:
        seq (str or SeqRecord): sequence
        region_unpaired (list ints): list of indices for positions in the sequence that need to be unpaired
        region_paired_A (list ints): list of indices for positions in the sequence that need to be paired,
            specifically ordered to pair with indices in region_paired_B
        region_paired_B (list ints): list of indices for positions in the sequence that need to be paired,
            specifically ordered to pair with indices in region_paired_A
        regionA (list ints): list of indices for positions in the sequence that need to not interact with regionB
        regionB (list ints): list of indices for positions in the sequence that need to not interact with regionA
        epsilon_interaction (float): Maximum base-pair-probability for 2 regions to be considered non-interacting
        epsilon_punpaired (float): Minimum probability unpaired for region to be unpaired
        epsilon_avg_punpaired (float): Average probability unpaired for region to be unpaired
        epsilon_paired (float): Minimum base-pair-probability for 2 regions to be considered paired
        epsilon_avg_paired (float): Average base-pair-probability for 2 regions to be considered paired
        prob_factor (float): the factor to multiply probabilities by (divide for epsilon_interaction) (default 1.0)
        save_image (str): if specified save an image of bpp to this file (default None)
        lines (list ints): list of sequences positions to draw lines in bpp figure (default None)

    Returns:
        bool if structure is ok
        array of the probability each nucleotide is unpaired

    Verbose possibilities:
        print(interaction_check.max(),punpaired_check.min(),paired_check.min())
    '''

    # get base-pair probability matrix and probability unpaired
    bpp = bpps(_get_dna_from_SeqRecord(seq),
               package='eternafold')
    p_unpaired = 1-bpp.sum(axis=0)

    # for regions specified, check if structure is ok
    results = {"unpaired_fail": [], "paired_fail": [],
               "interaction_fail": [], "p_unpaired": p_unpaired}
    if regions_unpaired is not None:
        for region_unpaired in regions_unpaired:
            if region_unpaired != []:
                punpaired_check = p_unpaired[region_unpaired]
                if ((punpaired_check.mean() < epsilon_avg_punpaired*prob_factor['unpaired']) or
                        (punpaired_check.min() < epsilon_punpaired*prob_factor['unpaired'])):
                    results["unpaired_fail"].extend(region_unpaired)

    if regionA is not None and regionA != [] and regionB != []:
        interaction_check = bpp[regionA][:, regionB]
        if (interaction_check.max() > epsilon_interaction/prob_factor['interaction']):
            results["interaction_fail"].append([regionA, regionB])

    if region_paired_A is not None and region_paired_A != []:
        paired_check = bpp[region_paired_A, region_paired_B]
        if ((paired_check.mean() < epsilon_avg_paired*prob_factor['paired']) or
                (paired_check.min() < epsilon_paired*prob_factor['paired'])):
            results["paired_fail"].append([region_paired_A, region_paired_B])

    # if good save image if needed
    if results["unpaired_fail"] + results["interaction_fail"] + results["paired_fail"] == []:
        results['pass'] = True
        if save_image is not None:
            plot_bpp(bpp, seq, save_image, lines, mutants=mutants)
    else:
        results['pass'] = False

    return results


###############################################################################
# visualize structures
###############################################################################


def plot_bpp(bpp, seq, save_image, lines=[], cmap='gist_heat_r',
             scale_factor=0.12, line_color='grey', xyticks_size=8, dpi=80,
             freq_report_nuc_number=10, mutant_color='cyan', mutants=[]):
    '''
    Plot base pair probability matrix and save image

    Args:
        bpp (array): square array with probability each base paired
        seq (str): sequence
        lines (list of ints): list of positions to draw vertical and horizontal lines
        save_image (str): location to save image to
        scale_factor (float): size of figure relative to length of sequence (default 0.12)
        cmap (str): cmap for bpp (default 'gist_heat_r')
        line_color (str): color to make lines (default 'grey')
        xyticks_size (int): size of labels for y and x (seq) (default 8)
        dpi (int): dpi to save image in (default 100)
        freq_report_nuc_number(int): on x and y axis add nucleotide number
            indexing from 0, every freq_report_nuc_number nucleotides (default 10)
    '''

    fig = plt.figure(figsize=(len(seq)*scale_factor, len(seq)*scale_factor))

    # plot base-pair probability as heatmap
    pos = plt.imshow(bpp, origin='upper', cmap=cmap)
    ax = fig.gca()
    cax = ax.inset_axes([1.04, 0.2, 0.05, 0.6])
    fig.colorbar(pos, shrink=0.5, ax=ax, cax=cax)

    # add nucleotide numbers to seq and plot on x and y axis
    xlabels = []
    ylabels = []
    for i, s in enumerate(seq):
        if i % freq_report_nuc_number == 0:
            ylabels.append(f'{i} {s}')
            xlabels.append(f'{s}\n{i}')
        else:
            xlabels.append(s)
            ylabels.append(s)
    plt.xticks(range(len(seq)), xlabels, size=xyticks_size)
    plt.yticks(range(len(seq)), ylabels, size=xyticks_size)
    # plot vertical and horizantal lines to mark sequence regions
    for line in lines:
        plt.hlines(line, 0, len(seq), color=line_color)
        plt.vlines(line, 0, len(seq), color=line_color)
    for mut in mutants:
        ax.get_xticklabels()[mut].set_color(mutant_color)
        ax.get_yticklabels()[mut].set_color(mutant_color)
        plt.hlines(mut, 0, len(seq), color=mutant_color)
        plt.vlines(mut, 0, len(seq), color=mutant_color)

    # formatting
    plt.xlim(0, len(seq))
    plt.ylim(len(seq), 0)

    # save
    plt.savefig(save_image, bbox_inches='tight', dpi=dpi)
    plt.close()


def plot_punpaired(p_unpaired, xlabels, seqs, muts, lines, pad_lines, save_image,
                   cmap='gist_heat_r', linewidth=8, line_color='cyan', pad_line_color='lime',
                   scale_factor=0.2, seq_color='grey', mutant_color='cyan',
                   xyticks_size=8, dpi=80):
    '''
    Plot probability probability unpaired for sequences and save image

    Args:
        p_unpaired (dict): dictionary with sequence_name:array of p unpaired, all p unpaired must be same length
        xlabels (list): list of xlabels the length of sequences
        seqs (list of str): list of sequences to plot
        muts (list of lists of int): for each sequence a list of mutation locations
        lines (list of ints): list of sequence positions to draw vertical lines over all sequences
        pad_lines (list of list of ints); for each sequence a list of location to draw lines
        save_image (str): location to save image to
        scale_factor (float): size of figure relative to length of sequence and number of sequences (default 0.3)
        cmap (str): cmap for bpp (default 'gist_heat_r')
        linewidth (int): width of vertical lines to draw (lines and pad_line) (default 8)
        line_color (str): color of lines (default 'cyan')
        pad_line_color (str): color of pad_lines (default 'lime')
        seq_color (str): color to make sequence (default 'grey')
        mutant_color (str): color to make mutations (default 'cyan')
        xyticks_size (int): size of labels for y (sequence name) and x (xlabels) (default 8)
        dpi (int): dpi to save image in (default 100)
    '''

    fig = plt.figure(figsize=(len(seqs[0])*scale_factor,
                              len(p_unpaired)*scale_factor))

    # plot p unpaired as heatmap
    df = pd.DataFrame(p_unpaired).T
    pos = plt.imshow(df, cmap=cmap)
    fig.colorbar(pos, shrink=0.5, location='top')

    # plot sequence ontop of heatmap, coloring the mutations
    ax = plt.gca()
    for i, (seq, mut) in enumerate(zip(seqs, muts)):
        for j, nuc in enumerate(seq):
            if j in mut:
                text = ax.text(j, i, nuc, ha="center",
                               va="center", color=mutant_color, weight='bold')
            else:
                text = ax.text(j, i, nuc, ha="center",
                               va="center", color=seq_color)
    plt.yticks(range(len(df)), df.index, size=xyticks_size)
    plt.xticks(range(len(xlabels)), xlabels, size=xyticks_size)

    # plot vertical lines to deliniate regions of sequence
    y1, y2 = ax.get_ylim()
    for line in lines:
        plt.vlines(line-0.5, y1+1, y2-1, color=line_color, linewidth=linewidth)
    for i, line in enumerate(pad_lines):
        if line != []:
            plt.vlines(line[0]-0.5, i+0.5, i-0.5,
                       color=pad_line_color, linewidth=linewidth)
            plt.vlines(line[1]-0.5, i+0.5, i-0.5,
                       color=pad_line_color, linewidth=linewidth)

    # formatting
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_ylim((y1, y2))
    ax.tick_params(axis='both', which='both', length=0)

    # save
    plt.savefig(save_image, bbox_inches='tight', dpi=dpi)
    plt.close()


def plot_all_bpp_from_fasta(fasta, save_image_folder):
    '''
    from a fasta file, just plot all bpps
    '''
    seqs = list(SeqIO.parse(fasta, "fasta"))
    for seq_rec in tqdm(seqs):
        seq = _get_dna_from_SeqRecord(seq_rec)
        bpp = bpps(seq, package='eternafold')
        save_image = f'{save_image_folder}/{seq_rec.name}.png'
        plot_bpp(bpp, seq, save_image)


###############################################################################
# helpers
###############################################################################


def _get_reverse_complement(seq):
    '''
    Return reverse complement of sequence, converts to DNA
    '''

    dna = _get_dna_from_SeqRecord(seq)
    reverse = dna[::-1]
    complements = {'T': 'A', 'A': 'T', 'C': 'G', 'G': 'C'}
    reverse_complement = ''.join([complements.get(s, s) for s in reverse])
    return reverse_complement


def _remove_seqs_in_other_list(seqsA, seqsB):
    '''
    Given 2 lists of SeqRecord, remove any SeqRecord in A
    that has same name as SeqRecord in B.
    DOES NOT check seq itself!
    '''

    names = [n.name for n in seqsB]
    good_seqs = []
    for seq_rec in seqsA:
        if seq_rec.name not in names:
            good_seqs.append(seq_rec)
    return good_seqs


def _get_dna_from_SeqRecord(seqrecord):
    '''
    From a sequence, string or SeqRecord, return dna str version
    '''

    if type(seqrecord) == str:
        dna = seqrecord.upper().replace("U", "T")
    else:
        dna = str(seqrecord.seq).upper().replace("U", "T")
    return dna


def _get_all_rand_seq(length, bases=BASES):
    '''
    for a given length, return all random sequences (SeqRecord) using bases
    '''

    all_seq = []
    for x in product(bases, repeat=length):
        seq = ''.join(x)
        seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), f' {seq}', '', '')
        all_seq.append(seq_rec)
    return all_seq


def _get_random_barcode(num_bp=8, num5hang=0, num3hang=0,
                        polyA5=0, polyA3=0,
                        loop=TETRALOOP, bases=BASES):
    '''
    Return random barcodes of specified structure

    Args:
        num_bp (int): length of stem (default 8)
        num5hang (int): length of random 5' hang (default 0)
        num3hang (int): length of random 3' hang (default 0)
        polyA5 (int): length of polyA 5' hang (placed before random) (default 0)
        polyA3 (int): length of polyA 3' hang (placed after random) (default 0)
        loop (str): sequence of loop
        bases (list): list of bases that can be used

    Returns:
        SeqRecord of a random barcode

    '''

    # get random combinations of bases for random/barcode regions
    uid = ''.join(choices(bases, k=num5hang+num_bp+num3hang))

    # split barcode in stem and hang regions
    hang5 = uid[:num5hang]
    if num3hang == 0:
        hang3 = ''
        stemA = uid[num5hang:]
    else:
        hang3 = uid[-num3hang:]
        stemA = uid[num5hang:-num3hang]
    stemB = _get_reverse_complement(stemA)

    # put all barcode parts together
    seq = ("A"*polyA5)+hang5+stemA+loop+stemB+hang3+("A"*polyA3)
    name = f' stem{stemA}_{hang5}hang{hang3}_{polyA5}polyA{polyA3}'
    seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), name, '', '')

    return seq_rec


def _get_mutations_from_name(name):
    '''
    using the naming convention _#N-N_ where # is nucleotide number
    N are the nucleotides native and mutant, find these and return the number
    '''

    nucs = []
    for name_part in name.split('_'):
        if len(name_part) > 2:
            if name_part[-2] == '-':
                nucnum_A, nuc_B = name_part.split('-')
                num, nuc_A = nucnum_A[:-1], nucnum_A[-1]
                if nuc_A in BASES and nuc_B in BASES:
                    nucs.append(int(num))
    return nucs


def _get_5_3_split(length, hang, polyAhang, min_length_stem, max_length_stem, pad_side, loop, seq_len):

    # initialize variables
    loop_len = len(loop)
    initial_structure = {"N": 0, 'bp': 0, 'hang_rand': hang,
                         'hang_polyA': polyAhang, 'loop': loop}
    structs = {"5": initial_structure.copy(),
               "3": initial_structure.copy()}
    min_pad_for_stem = min_length_stem*2 + loop_len + hang + polyAhang
    max_pad_for_stem = max_length_stem*2 + loop_len + hang + polyAhang
    unpaired_length = loop_len+hang+polyAhang
    # if only pad one side
    if pad_side == "3'" or pad_side == "5'":
        # make a helix if long enough
        if length >= min_pad_for_stem:
            num_bp = (length-unpaired_length)//2
            if num_bp > max_length_stem:
                print("WARNING: stem too long, fix not implemented.")
        else:
            num_bp = 0

        # assign to correct side
        if pad_side == "3'":
            structs["3"]['N'] = length
            structs["3"]['bp'] = num_bp
        elif pad_side == "5'":
            structs["5"]['N'] = length
            structs["5"]['bp'] = num_bp

    # if pad both, make the split
    elif pad_side == 'both':
        # if too short for 1 stem, split into 2 unstructured regions
        if length < min_pad_for_stem:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = length//2
            if length % 2 != 0:
                structs["3"]['N'] += 1

        # if enough for 1 stem but not too much for 1 stem, make just 1 stem
        elif length < max_pad_for_stem:
            structs["3"]['N'] = length
            structs["3"]['bp'] = (length-unpaired_length)//2

        # if too much for 1 stem but not enough for 2 stems, make 1 stem as long as possible
        elif length < 2*min_pad_for_stem:
            structs["5"]['N'] = length-max_pad_for_stem
            structs["3"]['N'] = max_pad_for_stem
            structs["3"]['bp'] = max_length_stem

        # if enough for 2 stems, make 2
        else:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = length//2
            if length % 2 != 0:
                structs["3'"]['N'] += 1
            structs["5"]['bp'] = (structs["5"]['N']-unpaired_length)//2
            structs["3"]['bp'] = (structs["3"]['N']-unpaired_length)//2
            if (num_bp5 > max_length_stem) or (num_bp3 > max_length_stem):
                print("WARNING: stem too long, fix not implemented.")
    else:
        print("ERROR: pad_side not recognized.")

    # get hang and loop lengths
    for side, struct in structs.items():
        if struct['bp'] == 0:
            struct['loop'] = ''
        if struct['N'] == 0:
            struct['hang_polyA'] = 0
            struct['hang_rand'] = 0
        elif struct['N'] != (struct['bp']*2)+struct['hang_polyA']+struct['hang_rand']+len(struct['loop']):
            if struct['hang_rand'] == 0 and struct['hang_polyA'] != 0:
                struct['hang_polyA'] = struct['N'] - \
                    ((struct['bp']*2)+len(struct['loop']))
            else:
                struct['hang_rand'] = struct['N'] - \
                    ((struct['bp']*2)+struct['hang_polyA']+len(struct['loop']))

    # get parts that should and shouldn't be paired
    regions = {}
    parts = [structs["5"]['bp'], len(structs["5"]['loop']), structs["5"]['bp'],
             structs["5"]['hang_polyA'] + structs["5"]['hang_rand'],
             seq_len, structs["3"]['hang_polyA'] + structs["3"]['hang_rand'],
             structs["3"]['bp'], len(structs["3"]['loop']), structs["3"]['bp']]
    regions['unpaired'] = [list(range(sum(parts[:1]),
                                      sum(parts[:2]))),
                           list(range(sum(parts[:3]),
                                      sum(parts[:4]))),
                           list(range(sum(parts[:5]),
                                      sum(parts[:6]))),
                           list(range(sum(parts[:7]),
                                      sum(parts[:8])))]
    regions['pairedA'] = list(range(sum(parts[:0]),
                                    sum(parts[:1])))
    regions['pairedB'] = list(range(sum(parts[:2]),
                                    sum(parts[:3])))[::-1]
    regions['noninteractA'] = list(range(sum(parts[:0]),
                                         sum(parts[:4])))
    regions['noninteractB'] = list(range(sum(parts[:4]),
                                         sum(parts[:5])))
    regions['pairedA'].extend(list(range(sum(parts[:6]),
                                         sum(parts[:7]))))
    regions['pairedB'].extend(list(range(sum(parts[:8]),
                                         sum(parts[:9])))[::-1])
    regions['noninteractA'].extend(list(range(sum(parts[:5]),
                                              sum(parts[:9]))))

    return structs, regions


def _get_pad_from_name(name):
    name.split('pad')
    pad5_length = int(name.split('pad')[0].split("_")[-1])
    pad3_length = int(name.split('pad')[1].split("_")[0])
    return pad5_length,pad3_length

def _update_seq_count(seq_count,struct_results):
    if struct_results["unpaired_fail"] != []:
        seq_count['unpaired'] += 1
    if struct_results["interaction_fail"] != []:
        seq_count['interaction'] += 1
    if struct_results["paired_fail"] != []:
        seq_count['paired'] += 1
    return seq_count

def _get_5_3_split_multi(lengths, hang, polyAhang, min_length_stem,
                         max_length_stem, pad_side, loop, seq_lens):
    # length are lengths of pads
    # for use in all same pad code, should work if we keep track of cut points
    lengths.sort()
    seq_lens.sort(reverse=True)
    current_length = 0
    all_structs = []
    current_regions = {'unpaired': [], 'pairedA': [], 'pairedB': [],
                       'noninteractA': [], 'noninteractB': [], }
    cutoff_by_len = {}
    regions_by_len = {}
    current_cutoffs = [0, 0]
    for full_length, seq_len in zip(lengths, seq_lens):
        length = full_length - current_length
        structs, regions = _get_5_3_split(length, hang, polyAhang,
                                          min_length_stem, max_length_stem,
                                          pad_side, loop, seq_len)
        print(f"Finding a {structs['5']['N']}nt 5' pad with {structs['5']['bp']}bp stem {structs['5']['hang_rand']}nt random hang {structs['5']['hang_polyA']}nt polyA hang.")
        print(f"Finding a {structs['3']['N']}nt 3' pad with {structs['3']['bp']}bp stem {structs['3']['hang_rand']}nt random hang {structs['3']['hang_polyA']}nt polyA hang.")

        for part, region in regions.items():
            if part == 'noninteractB':
                current_regions[part] = region
            else:
                current_regions[part].extend(region)
        current_length += length
        all_structs.append(structs)
        current_cutoffs = [current_cutoffs[0]+structs['5']
                           ['N'], current_cutoffs[1]+structs['3']['N']]
        cutoff_by_len[seq_len] = current_cutoffs
        regions_by_len[seq_len] = current_regions.copy()

    return all_structs, cutoff_by_len, regions_by_len


def _get_same_length(seqs):
    '''
    Checks if all sequences in a given fasta are the same length
    '''

    length = len(seqs[0].seq)

    for seq_rec in seqs[1:]:
        len_seq = len(seq_rec)
        if length != len_seq:
            return False, np.nan

    return True, length


def get_regions_for_doublemut(doublemuts):
    '''
    # TODO probably should be input double mutants and this goes to helper
    '''

    regionAss, regionBss = [], []
    for mutstr in doublemuts:
        regionAs, regionBs = [], []
        for mutstrreg in mutstr.split('..'):
            strA, strB = mutstrreg.split('.')
            regionA = []
            for nucrange in strA.split(','):
                nucrange = [int(x) for x in nucrange.split('-')]
                if len(nucrange) == 2:
                    regionA.extend(list(range(nucrange[0], nucrange[1]+1)))
                else:
                    regionA.extend(nucrange)
            regionAs.append(regionA)
            regionB = []
            for nucrange in strB.split(','):
                nucrange = [int(x) for x in nucrange.split('-')]
                if len(nucrange) == 2:
                    regionB.extend(list(range(nucrange[0], nucrange[1]+1)))
                else:
                    regionB.extend(nucrange)
            regionBs.append(regionB)
        regionAss.append(regionAs)
        regionBss.append(regionBs)
    return regionAss, regionBss


###############################################################################
# UNDER CONSTRUCTION
###############################################################################


def plot_punpaired_from_fasta(fasta, save_image):
    # NOT well tested
    seqs = list(SeqIO.parse(fasta, "fasta"))
    p_unpaireds = {}
    seqs_list = []
    muts = []
    for seq_rec in seqs:
        seq = str(seq_rec.seq).upper().replace("T", "U")
        bpp = bpps(seq,
                   package='eternafold')
        p_unpaired = 1-bpp.sum(axis=0)
        p_unpaireds[seq_rec.name] = p_unpaired
        seqs_list.append(seq)
        muts.append([5, 15])
    labels = []
    for i in range(len(seqs[0])):
        if i % 10 == 0:
            labels.append(i)
        else:
            labels.append('')

    plot_punpaired(p_unpaireds, labels, seqs, muts, [], [], save_image)

