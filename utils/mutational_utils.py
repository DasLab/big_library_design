import numpy as np
import pandas as pd
import os
from itertools import product, chain
from random import shuffle, sample, choices, random
from tqdm import tqdm
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
from arnie.bpps import bpps
from arnie.utils import convert_dotbracket_to_bp_list


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


def randomly_select_seqs(fasta, out_file, N):
    '''
    Given a fasta randomly select N sequences and save to new fasta
    '''

    all_seqs = list(SeqIO.parse(fasta, "fasta"))
    if len(seqs) > N:
        all_seqs = sample(all_seqs, N)

        SeqIO.write(all_seqs, out_fasta, "fasta")

        print(f'Written {out_file} for with {N} selected sequences.')

    else:

        SeqIO.write(all_seqs, out_fasta, "fasta")

        print(f'WARNING: Written {out_file} but {fasta} had less than {N} sequences so none removed.')


def get_same_length(fasta):
    '''
    Checks if all sequences in a given fasta are the same length
    '''

    seqs = list(SeqIO.parse(fasta, "fasta"))

    length = len(seqs[0].seq)

    for seq_rec in seqs[1:]:
        if length != len_seq:
            return False

    return True


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


###############################################################################
# other utils
###############################################################################


def get_bp_set_from_dotbracket(dotbracket):
    '''
    Given a dotbracket structure, return a list of base-pairs
    IGNORES pseudoknots TODO
    '''

    return convert_dotbracket_to_bp_list(dotbracket)


###############################################################################
# get desired sequences
###############################################################################


def get_windows(fasta, window_length, window_slide, out_fasta=None,
                circularize=False):
    '''
    Get sliding windows from an inputted fasta file.

    Args:
        fasta (str): fasta file containing sequence to get windows of
        window_length (int): length (in number nt) of desired windows
        window_slide (int): number nucleotide to slide between windows
        out_fasta (str): if specified, save windows to fasta (default None)
        circularize (bool): whether to circularize the sequence thus once it 
            reaches the end, create windows with 3'-5' connected (default False)

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

    for seq_rec in seqs:

        # loop through sequence every window_slide nucleotides
        seq = _get_dna_from_SeqRecord(seq_rec)
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

            # save with name inclusive!
            new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                      f'{seq_rec.name}_{namenum}', '', '')
            windows.append(new_rec)

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
        mutant with 138 mutanted from G to A is: P4P6_138G-A
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


def get_all_double_mutants(fasta, regionAs, regionBs, out_fasta=None, bases=BASES):
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
        mutant with 138 mutanted from G to A and 150 C to T is: P4P6_138G-A_150C-T
        Indexing from 0 on mutants are ordered by nucleotide number.
        Generally, total is 9*lengthA*lengthB
    '''

    print("Getting all double mutants.")

    # get all sequences and initialize
    all_WT = list(SeqIO.parse(fasta, "fasta"))
    all_double_mutants = []

    # check user specified regions for every sequence in the fasta
    if len(all_WT) != len(regionAs) or len(all_WT) != len(regionBs):
        print(f'WARNING: you must list regions (as a list of lists) for all sequences in the fasta, number sequences in fasta {len(all_WT)} and number of regions specified {len(regionAs)} {len(regionBs)}')

    for record, regionA, regionB in zip(all_WT, regionAs, regionBs):

        seq = _get_dna_from_SeqRecord(record)

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


def get_wcf_rescue_mutants(fasta, bp_sets, out_fasta=None):
    '''
    Get all watson-crick-franklin rescue mutants between specified basepairs for all sequences in a fasta file.

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        bp_sets (list of lists of pairs): for each sequence, a list of basepairs as a tuple or list of int
        out_fasta (str): if specified, save mutants to fasta (default None)

    Returns:
        list of SeqRecord with the mutants
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_#wt-mutant eg for sequence P4P6
        mutant with 138 mutanted from G to A and 150 C to T is: P4P6_138G-A_150C-T
        Indexing from 0 on mutants are ordered by nucleotide number.
        Generally, total is 3*numbps
    '''

    # TODO can allow people to specify
    wfc_base_pairs = ['AT', 'TA', 'CG', 'GC']

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
        if out_fasta specified, also saves these to fasta file naming added _#pad# where numbers are length of pad
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

    '''

    # TODO
    # probably should add ability to randomly generate but this
    # is fast enough for these small barcode aka get_random_barcode
    print("Getting all possible barcodes.")
    all_barcodes = []
    for x in product(bases, repeat=num_bp+num5hang+num3hang):
        uid = ''.join(x)

        if num5hang+num3hang != 0:
            hang5 = uid[:num5hang]
            if num3hang == 0:
                hang3 = ''
                uid = uid[num5hang:]
            else:
                hang3 = uid[-num3hang:]
                uid = uid[num5hang:-num3hang]

            seq = ("A"*polyA5)+hang5+uid+loop + \
                _get_reverse_complement(uid)+hang3+("A"*polyA3)
            seq_rec = SeqIO.SeqRecord(Seq.Seq(seq),
                                      f' {uid}_{hang5}hang{hang3}', '', '')
        else:
            seq = ("A"*polyA5)+uid+loop + \
                _get_reverse_complement(uid)+("A"*polyA3)
            seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), f' {uid}', '', '')
        all_barcodes.append(seq_rec)
    if out_fasta is not None:
        SeqIO.write(all_barcodes, out_fasta, "fasta")
        print(f'Saved all barcodes to {out_fasta}.')
    return all_barcodes


def add_pad(fasta, out_fasta, bases=BASES, padding_type='SL_same_per_length',
            epsilon_punpaired=MINPROB_UNPAIRED,
            epsilon_interaction=MAXPROB_NONINTERACT,
            epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
            epsilon_paired=MINPROB_PAIRED,
            epsilon_avg_paired=MINAVGPROB_PAIRED,
            loop=TETRALOOP, min_hang=3, min_num_samples=30):
    '''

    '''

    # crurent options: rand_same_all SL_same_per_length

    seqs = list(SeqIO.parse(fasta, "fasta"))
    seq_by_length = {}
    for seq_rec in seqs:
        len_seq = len(seq_rec.seq)
        if len_seq in seq_by_length:
            seq_by_length[len_seq].append(seq_rec)
        else:
            seq_by_length[len_seq] = [seq_rec]

    # for each len randomly select 1% of sequences or 10
    selected_sec = []
    for len_seq, seq_group in seq_by_length.items():
        number_to_select = min(len(seq_group), max(
            min_num_samples, len(seq_group)*0.01))
        selected_sec.append(sample(seq_group, k=number_to_select))
        print(f'Pad for length {len_seq} search using {number_to_select} sequences for structure check for each length.')
    # get length of pad needed
    desired_len = max(seq_by_length.keys())

    pads_by_len = {}
    if padding_type == 'SL_same_per_length':
        for seqs in selected_sec:
            length_to_add = desired_len-len(str(seqs[0].seq))
            if length_to_add == 0:
                pads_by_len[len(str(seqs[0].seq))] = ['', '']
            else:

                pad_5n, pad_3n = _get_5_3_split(length_to_add)

                pad_5s, num_bp5, num_hang5, loop_len5 = _get_stem_pads(
                    pad_5n, "5'", loop, min_hang, bases)
                pad_3s, num_bp3, num_hang3, loop_len3 = _get_stem_pads(
                    pad_3n, "3'", loop, min_hang, bases)

                shuffle(pad_5s)
                shuffle(pad_3s)

                # numbp_loop_numbp_hang_seq_hang_numbp_loop_numbp
                part_lengths = [num_bp5, loop_len5, num_bp5, num_hang5, len(
                    str(seqs[0].seq)), num_hang3, num_bp3, loop_len3, num_bp3]

                print("Searching for 5' pad")
                region_unpaired = list(
                    range(sum(part_lengths[:1]), sum(part_lengths[:2])))
                region_unpaired.extend(
                    list(range(sum(part_lengths[:3]), sum(part_lengths[:4]))))

                region_paired_A = list(
                    range(sum(part_lengths[:0]), sum(part_lengths[:1])))
                region_paired_B = list(
                    range(sum(part_lengths[:2]), sum(part_lengths[:3])))[::-1]

                regionA = list(
                    range(sum(part_lengths[:0]), sum(part_lengths[:4])))

                regionB = list(
                    range(sum(part_lengths[:4]), sum(part_lengths[:5])))

                # loop through to find 5' that works
                if pad_5n == 0:
                    good5 = ''
                else:
                    good_pad = False
                    current_pad = 0
                    while not good_pad:

                        for i, seq in enumerate(seqs):
                            if i % 50 == 0 and i != 0:
                                print(i)
                            full_seq = pad_5s[current_pad].seq + \
                                _get_dna_from_SeqRecord(seq)

                            good_pad, p_unpaired = check_struct_bpp(
                                full_seq, region_unpaired, region_paired_A,
                                region_paired_B, regionA, regionB,
                                epsilon_interaction=epsilon_interaction,
                                epsilon_punpaired=epsilon_punpaired,
                                epsilon_avg_punpaired=epsilon_avg_punpaired,
                                epsilon_paired=epsilon_paired,
                                epsilon_avg_paired=epsilon_avg_paired)
                            if not good_pad:
                                break
                        current_pad += 1
                        if current_pad == len(pad_5s):
                            print("no pad found")

                    good5 = pad_5s[current_pad-1].seq

                print("Searching for 3' pad")
                region_unpaired = list(
                    range(sum(part_lengths[:5]), sum(part_lengths[:6])))
                region_unpaired.extend(
                    list(range(sum(part_lengths[:7]), sum(part_lengths[:8]))))
                region_paired_A = list(
                    range(sum(part_lengths[:6]), sum(part_lengths[:7])))
                region_paired_B = list(
                    range(sum(part_lengths[:8]), sum(part_lengths[:9])))[::-1]
                regionA = list(
                    range(sum(part_lengths[:5]), sum(part_lengths[:9])))

                # loop through to find 3' that works
                if pad_3n == 0:
                    pads_by_len[len(str(seqs[0].seq))] = [
                        good5, '']
                else:
                    good_pad = False
                    current_pad = 0
                    while not good_pad:

                        for i, seq in enumerate(seqs):
                            if i % 50 == 0 and i != 0:
                                print("b", i)

                            full_seq = good5 + \
                                _get_dna_from_SeqRecord(
                                    seq) + pad_3s[current_pad].seq
                            good_pad, p_unpaired = check_struct_bpp(
                                full_seq, region_unpaired, region_paired_A,
                                region_paired_B, regionA, regionB,
                                epsilon_interaction=epsilon_interaction,
                                epsilon_punpaired=epsilon_punpaired,
                                epsilon_avg_punpaired=epsilon_avg_punpaired,
                                epsilon_paired=epsilon_paired,
                                epsilon_avg_paired=epsilon_avg_paired)
                            if not good_pad:
                                break
                        current_pad += 1
                        if current_pad == len(pad_3s):
                            print("no pad found")

                    pads_by_len[len(str(seqs[0].seq))] = [
                        good5, pad_3s[current_pad-1].seq]

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

    elif padding_type == 'rand_same_all':
        selected_sec = list(chain(*selected_sec))

        max_length_to_add = desired_len-min(seq_by_length.keys())
        pad_5_len, pad_3_len = _get_5_3_split(max_length_to_add)
        good_pad = False
        while not good_pad:
            pad5 = ''.join(choices(bases, k=pad_5_len))
            pad3 = ''.join(choices(bases, k=pad_3_len))
            any_bad = False
            for i, seq in enumerate(selected_sec):
                # if i%10==0 and i!=0:
                #    print(i)
                length_to_add = desired_len-len(seq.seq)
                if length_to_add != 0:
                    pad_5n, pad_3n = _get_5_3_split(length_to_add)
                    full_seq = pad5[-pad_5n:] + \
                        _get_dna_from_SeqRecord(seq) + pad3[:pad_3n]
                    regionA = list(range(pad_5n))
                    regionA += list(range(pad_5n + len(seq.seq),
                                          pad_5n+len(seq.seq)+pad_3n))
                    regionB = list(range(pad_5n, pad_5n+len(seq.seq)))

                    this_good, p_unpaired = check_struct_bpp(full_seq, region_unpaired=regionA,
                                                             regionA=regionA, regionB=regionB,
                                                             epsilon_interaction=epsilon_interaction, epsilon_punpaired=epsilon_punpaired,
                                                             epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                             epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired)
                    if not this_good:
                        any_bad = True
                        break
            if not any_bad:
                good_pad = True

        # add pads to the sequences
        padded_seqs = []
        for seq_rec in seqs:
            length_to_add = desired_len-len(seq_rec.seq)
            if length_to_add != 0:
                pad_5n, pad_3n = _get_5_3_split(length_to_add)
                seq = _get_dna_from_SeqRecord(seq_rec)
                full_seq = pad5[-pad_5n:] + seq + pad3[:pad_3n]
                name = f'{seq_rec.name}_{pad_5n}pad{pad_3n}'
                padded_seqs.append(SeqIO.SeqRecord(
                    Seq.Seq(full_seq), name, '', ''))
            else:
                padded_seqs.append(seq_rec)

    if out_fasta is not None:
        SeqIO.write(padded_seqs, out_fasta, "fasta")
        print(f'Saved all padded sequences to {out_fasta}.')
    return padded_seqs


def add_fixed_seq_and_barcode(fasta, out_fasta=None, seq5=SEQ5, seq3=SEQ3,
                              num_bp=8, num5hang=0, num5polyA=4,
                              loop=TETRALOOP,
                              epsilon_interaction=MAXPROB_NONINTERACT,
                              epsilon_punpaired=MINPROB_UNPAIRED,
                              epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                              epsilon_paired=MINPROB_PAIRED,
                              epsilon_avg_paired=MINAVGPROB_PAIRED,
                              save_image_folder=None, save_bpp_fig=0,
                              punpaired_chunk_size=500, used_barcodes=None,
                              num_barcode_before_reduce=100):
    '''
    '''

    # if image folder does not exist create it
    if save_image_folder is not None:
        if not os.path.exists(save_image_folder):
            print(f'{save_image_folder} did not exists, creating.')
            os.makedirs(save_image_folder)

    # get and randomly shuffle all potential barcodes
    all_uids = get_all_barcodes(num_bp=num_bp, loop=loop, num5hang=num5hang,
                                polyA5=num5polyA)
    shuffle(all_uids)

    # read sequences and initialize variables
    seqs = SeqIO.parse(fasta, "fasta")
    all_full_seqs, p_unpaireds, rejected_uids = 0, [], {}, []
    pad_lines, new_lines, seqs_for_labeling, muts = [], [], [], []
    current_uid, chunk_count, image_count = 0, 0, 0

    print("Adding 5', barcode, 3'.")

    for seq_rec in tqdm(list(seqs)):
        chunk_count += 1

        # get sequence and name
        seq = _get_dna_from_SeqRecord(seq_rec)
        name = seq_rec.name+'_libraryready'
        mutations = _get_mutations_from_name(seq_rec.name)

        # 5_seq_num5polyA_num5hang_numbp_loop_numbp_3
        regionA = list(range(len(seq5), len(seq5)+len(seq)))
        regionB = list(range(len(seq5)+len(seq), len(seq5) +
                             len(seq)+len(loop)+num5hang+num5polyA+(2*num_bp)))
        region_unpaired = list(range(len(seq5)+len(seq),
                                     len(seq5) + len(seq)+num5hang+num5polyA))
        region_unpaired.extend(list(range(len(seq5)+len(seq)+num5hang+num5polyA+num_bp,
                                          len(seq5) + len(seq)+num5hang+num5polyA+num_bp+len(loop))))
        region_paired_A = list(range(len(seq5)+len(seq)+num5hang+num5polyA,
                                     len(seq5) + len(seq)+num5hang+num5polyA+num_bp))
        region_paired_B = list(range(len(seq5)+len(seq)+num5hang+num5polyA+num_bp+len(loop),
                                     len(seq5) + len(seq)+num5hang+num5polyA+(2*num_bp)+len(loop)))[::-1]
        lines = [len(seq5), len(seq5)+len(seq),
                 len(seq5)+len(seq)+len(loop)+num5hang+num5polyA+(2*num_bp)]

        # initialize values
        uid_good = False
        seq_count = 0

        if 'pad' in name:
            name.split('pad')
            pad5_length = int(name.split('pad')[0].split("_")[-1])
            pad3_length = int(name.split('pad')[1].split("_")[0])
            new_lines = [len(seq5)+pad5_length,
                         len(seq5)+len(seq)-pad3_length]
            lines.extend(new_lines)
        else:
            new_lines = []
        while not uid_good:
            uid = all_uids[current_uid].seq
            if used_barcodes is not None:
                while str(uid) in used_barcodes:
                    print('is true sometimes TEST')
                    current_uid += 1
                    uid = all_uids[current_uid].seq

            full_seq = f'{seq5}{seq}{uid}{seq3}'
            current_uid += 1

            prob_factor = 0.9**(1+(seq_count//num_barcode_before_reduce))
            if seq_count//num_barcode_before_reduce == 0 and seq_count != 0:
                print(f'For {seq}, failed to find barcode from {seq_count} barcodes, reducing (increasing for max prob) probabilities needed by a factor of {prob_factor}.')

            if save_image_folder is not None:

                if random() < save_bpp_fig:
                    uid_good, p_unpaired = check_struct_bpp(full_seq, region_unpaired, region_paired_A, region_paired_B, regionA, regionB,
                                                            epsilon_interaction=epsilon_interaction,
                                                            epsilon_punpaired=epsilon_punpaired, epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                            epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired,
                                                            lines=lines, prob_factor=prob_factor,
                                                            save_image=f'{save_image_folder}/{name}.png')
                else:
                    uid_good, p_unpaired = check_struct_bpp(full_seq, region_unpaired, region_paired_A, region_paired_B, regionA, regionB,
                                                            epsilon_interaction=epsilon_interaction, epsilon_punpaired=epsilon_punpaired,
                                                            epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                            epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired,
                                                            prob_factor=prob_factor)
            else:
                uid_good, p_unpaired = check_struct_bpp(full_seq, region_unpaired, region_paired_A, region_paired_B, regionA, regionB,
                                                        epsilon_interaction=epsilon_interaction, epsilon_punpaired=epsilon_punpaired,
                                                        epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                        epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired,
                                                        prob_factor=prob_factor)
            if not uid_good:
                rejected_uids.append(all_uids[current_uid-1])
                seq_count += 1
            if current_uid == len(all_uids):
                all_uids = rejected_uids
                current_uid = 0
                rejected_uids = []

        mutations = [x+len(seq5)+pad5_length for x in mutations]
        muts.append(mutations)
        pad_lines.append(new_lines)
        seqs_for_labeling.append(full_seq)
        all_full_seqs.append(SeqIO.SeqRecord(Seq.Seq(full_seq), name, '', ''))
        p_unpaireds[name] = p_unpaired
        if ((save_image_folder is not None) and
                (chunk_count == punpaired_chunk_size)):
            plot_punpaired(p_unpaireds,
                           [i if i % 10 == 0 else '' for i in range(
                               len(seqs_for_labeling[0]))],
                           seqs_for_labeling, muts,
                           [len(seq5), len(seq5)+len(seq), len(seq5) +
                            len(seq)+len(loop)+num5hang+num5polyA+(2*num_bp)],
                           pad_lines,
                           f'{save_image_folder}/all_p_unpaired{image_count}.png')

            chunk_count = 0
            image_count += 1
            p_unpaireds = {}
            pad_lines, seqs_for_labeling, muts = [], [], []
    if out_fasta is not None:
        SeqIO.write(all_full_seqs, out_fasta, "fasta")
        print(f'Saved all full sequences to {out_fasta}.')
    if save_image_folder is not None and len(p_unpaired) != 0:
        plot_punpaired(p_unpaireds,
                       [i if i % 10 == 0 else '' for i in range(
                           len(seqs_for_labeling[0]))],
                       seqs_for_labeling, muts,
                       [len(seq5), len(seq5)+len(seq), len(seq5) +
                        len(seq)+len(loop)+num5hang+num5polyA+(2*num_bp)],
                       pad_lines,
                       f'{save_image_folder}/all_p_unpaired.png')
    return all_full_seqs


###############################################################################
# check structures
###############################################################################


def check_struct_bpp(seq, region_unpaired=None, region_paired_A=None,
                     region_paired_B=None,
                     regionA=None, regionB=None,
                     epsilon_interaction=MAXPROB_NONINTERACT,
                     epsilon_punpaired=MINPROB_UNPAIRED,
                     epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                     epsilon_paired=MINPROB_PAIRED,
                     epsilon_avg_paired=MINAVGPROB_PAIRED,
                     prob_factor=1.0,
                     lines=None, save_image=None):
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
        epsilon_interaction (float): Maximum base-pair-probability for 2 regions to be considered noninteracting
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
        print(interaction_check.max(),punpaired_check.min())#paired_check.min(),
    '''

    # get base-pair probability matrix and probability unpaired
    bpp = bpps(_get_dna_from_SeqRecord(seq),
               package='eternafold')
    p_unpaired = 1-bpp.sum(axis=0)

    # for regions specified, check if structure is ok
    bad_struct = False
    if region_unpaired is not None:
        punpaired_check = p_unpaired[region_unpaired]
        bad_struct = bad_struct or (
            punpaired_check.min() < epsilon_punpaired*prob_factor)
        bad_struct = bad_struct or (
            punpaired_check.mean() < epsilon_avg_punpaired*prob_factor)

    if regionA is not None:
        interaction_check = bpp[regionA][:, regionB]
        bad_struct = bad_struct or (
            interaction_check.max() > epsilon_interaction/prob_factor)

    if region_paired_A is not None and region_paired_A != []:
        paired_check = bpp[region_paired_A, region_paired_B]
        bad_struct = bad_struct or (
            paired_check.mean() < epsilon_avg_paired*prob_factor)
        bad_struct = bad_struct or (
            paired_check.min() < epsilon_paired*prob_factor)

    # if bad structure return False
    if bad_struct:
        return False, p_unpaired

    # otherwise, save image if needed and return True
    else:
        if save_image is not None:
            plot_bpp(bpp, seq, lines, save_image)
        return True, p_unpaired


###############################################################################
# visualize structures
###############################################################################


def plot_bpp(bpp, seq, lines, save_image, cmap='gist_heat_r',
             scale_factor=0.12, line_color='grey', xyticks_size=8, dpi=100,
             freq_report_nuc_number=10):
    '''
    Plot base pair probability matrix and save image

    Args:
        bpp (array): square array with probability each base paired
        seq (str): sequence
        lines (list of ints): list of positions to draw vertical and horizantal lines
        save_image (str): location to save image to
        scale_factor (float): size of figure relative to length of sequence (default 0.12)
        cmap (str): cmap for bpp (default 'gist_heat_r')
        line_color (str): color to make lines (default 'grey')
        xyticks_size (int): size of labels for y and x (seq) (default 8)
        dpi (int): dpi to save image in (default 100)
        freq_report_nuc_number(int): on x and y axis add nucleotide number
            indexing from 0, every freq_report_nuc_number nucleotides (default 10)
    '''

    plt.figure(figsize=(len(seq)*scale_factor, len(seq)*scale_factor))

    # plot base-pair probability as heatmap
    plt.imshow(bpp, origin='upper', cmap=cmap)

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

    # formatting
    plt.xlim(0, len(seq))
    plt.ylim(len(seq), 0)

    # save
    plt.savefig(save_image, bbox_inches='tight', dpi=dpi)
    plt.close()


def plot_punpaired(p_unpaired, xlabels, seqs, muts, lines, pad_lines, save_image,
                   cmap='gist_heat_r', linewidth=8, line_color='cyan', pad_line_color='lime',
                   scale_factor=0.3, seq_color='grey', mutant_color='cyan',
                   xyticks_size=8, dpi=100):
    '''
    Plot probability probability unpaired for sequences and save image

    Args:
        p_unpaired (dict): dictionary with sequence_name:array of punpaired, all punpaired must be same length
        xlabels (list): list of xlabels the length of seuqences 
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

    plt.figure(figsize=(len(seqs[0])*scale_factor,
                        len(p_unpaired)*scale_factor))

    # plot p unpaired as heatmap
    df = pd.DataFrame(p_unpaired).T
    plt.imshow(df, cmap=cmap)

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


def _get_5_3_split(length):
    # 4+6+4+6 20
    if length < 15 and length % 2 == 0:
        pad_5_len, pad_3_len = length//2, length//2
    elif length < 15:
        pad_5_len, pad_3_len = length//2, 1+(length//2)
    elif length < 20:
        pad_5_len, pad_3_len = 0, length
    elif length < 32:
        pad_5_len, pad_3_len = length-20, 20
    elif length % 2 == 0:
        pad_5_len, pad_3_len = length//2, length//2
    else:
        pad_5_len, pad_3_len = length//2, 1+(length//2)
    return pad_5_len, pad_3_len


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


def _get_stem_pads(pad_length, side="5'", loop=TETRALOOP,
                   min_hang=3, bases=BASES, min_stem_length=4):
    '''
    Get all possible stems for pad

    Args:
        pad_length (int): the total length of desired pad
        side (str): either 5' or 3' pad
        loop (str): sequence of the loop in stem
        min_hang (int): the minimum length of sequence between pad and sequence of interest, 
            will only min_hang or min_hang+1 (default 3)
        bases (list): list of possible bases we can use
        min_stem_length (int): minimum stem length to allow 
            if not long enough to form this stem, return random sequence (default 4)

    Returns:
        list of barcodes
        number of base-pairs in stem
        number of hanging nucleotides between stem and sequence of interest
        number of nucleotides in loop
    '''

    # TODO probably let specifiy hang has to be polyA

    # for things too short for stem return random sequence
    min_pad_for_stem = min_stem_length*2 + len(loop) + min_hang
    if pad_length < min_pad_for_stem:
        barcodes = _get_all_rand_seq(pad_length, bases)
        return barcodes, 0, 0, pad_length

    else:
        # get number base pairs and length of hang
        num_bp = (pad_length-len(loop)-min_hang)//2
        num_hang = pad_length-len(loop)-(2*num_bp)
        
        # for 5' should be on 3' ((((....))))....
        if side == "5'":
            barcodes = get_all_barcodes(num_bp=num_bp, num3hang=num_hang,
                                        loop=loop, bases=bases)

        # for 3' should be on 5' ....((((....))))
        elif side == "3'":
            barcodes = get_all_barcodes(num_bp=num_bp, num5hang=num_hang,
                                        loop=loop, bases=bases)

        else:
            print("ERROR side must be 5' or 3'")
        
        return barcodes, num_bp, num_hang, len(loop)

    
def get_used_barcodes(fasta, start, end):
    '''
    from a fasta file return all sequences between start and end inclusive
    '''

    # TODO the add barcodes code should have an  option of
    # fasta file with barcodes to not use, start, end
    # and this code should move to helper?
    # inclusive

    all_seqs = list(SeqIO.parse(fasta, "fasta"))
    barcodes = []
    for record in all_seqs:
        seq = _get_dna_from_SeqRecord(record)

        barcode = seq[start:end+1]
        barcodes.append(str(seq))
    return barcodes


def get_regions_for_doublemut(doublemuts):
    '''
    # TODO probably should be input double mutants and this goes to helper
    '''

    regionAs, regionBs = [], []
    for mutstr in doublemuts:
        strA, strB = mutstr.split('.')
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
    return regionAs, regionBs


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
