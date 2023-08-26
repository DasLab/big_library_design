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
from glob import glob


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

def split_fasta_file(fasta,N,folder,name,offset=0):
    if not os.path.isdir(folder):
        os.mkdir(folder)

    seqs = list(SeqIO.parse(fasta, "fasta"))
    len_seqs = (len(seqs)//(N))+1


    for n in range(N):
        if not os.path.isdir(f'{folder}/{n+offset}'):
            os.mkdir(f'{folder}/{n+offset}')

        seqset = seqs[len_seqs*n:min(len_seqs*(n+1),len(seqs))]
        SeqIO.write(seqset,f'{folder}/{n+offset}/{name}',"fasta")

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
            names.append(get_name_SeqRecord(seq))
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


def remove_library_elements(fasta,out_fasta,lenbarcode,len5=len(SEQ5),len3=len(SEQ3)):
    seqs = parse_input_sequences(fasta)
    strip_seqs = []
    for seq_rec in seqs:
        seq = seq_rec.seq
        pad5_length,pad3_length = _get_pad_length_from_name(get_name_SeqRecord(seq_rec))
        new_seq = seq[len5+pad5_length:len(seq)-(len3+lenbarcode+pad3_length)]
        name = seq_rec.name
        new_name = []
        for x in name.split('_'):
            if 'pad' not in x and 'libraryready' not in x:
                new_name.append(x)
        strip_seqs.append(SeqIO.SeqRecord(new_seq,'_'.join(new_name),'',''))
    SeqIO.write(strip_seqs, out_fasta, "fasta")




def parse_input_sequences(seqs):
    if type(seqs)==str:
        if os.path.isfile(seqs):
            extension = seqs.rsplit('.',1)[1]
            if extension == 'fasta' or extension == 'fa':
                return list(SeqIO.parse(seqs, "fasta"))
            elif extension == 'tsv':
                df = pd.read_table(seqs)
                cols = df.columns
                name_seq_columns = [['id','sequence']]
                for name_col, seq_col in name_seq_columns:
                    if name_col in cols and seq_col in cols:
                        df[seq_col] = df[seq_col].str.strip()
                        return df.apply(lambda row: SeqIO.SeqRecord(Seq.Seq(_get_dna_from_SeqRecord(row[seq_col].strip())),
                                                  name=str(row[name_col])), axis=1).to_list()
                print("ERROR recognized tsv but did not recognize column nammes, please submmit and issue to request new input format.")
            else:
                print("ERROR unrecognized input file extension, please submit an issue to have the input format recognized.")
        else:
            print("ERROR recognized input as string, but file does not extist.")
    elif type(seqs)==list and isinstance(seqs[0],SeqIO.SeqRecord):
        return seqs
    else:
        print("ERROR unrecognized input, please submit an issue to have the input format recognized.")


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

    seqs = parse_input_sequences(fasta)
    return _get_same_length(seqs)


def check_sequences_contents(fasta,seq5=SEQ5, seq3=SEQ3,bases=BASES,content_only=False):
    problems = []
    names = []
    sequences = []
    seqs = list(SeqIO.parse(fasta, "fasta"))
    for seq_rec in seqs:
        seq = _get_dna_from_SeqRecord(seq_rec)
        if not content_only:
            #if seq in sequences:
            #    problems.append(f"{get_name_SeqRecord(seq_rec)} is a repeated sequence.")
            #if get_name_SeqRecord(seq_rec)in names:
            #    problems.append(f"{get_name_SeqRecord(seq_rec)} is a repeated name+description.")
            if seq[:len(seq5)] != seq5:
                problems.append(f"{get_name_SeqRecord(seq_rec)} has incorrect 5' sequence.")
            if seq[-len(seq3):] != seq3:
                problems.append(f"{get_name_SeqRecord(seq_rec)} has incorrect 3' sequence.")
            #sequences.append(seq)
            #names.append(get_name_SeqRecord(seq_rec))
        for n in seq:
            if n not in bases:
                problems.append(f"{get_name_SeqRecord(seq_rec)} has non {bases} base.")
        
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


def get_used_barcodes(fasta, start, end, gu_present = False,num_bp=None,unique=True):
    '''
    from a fasta file return all sequences between start and end inclusive
    '''

    all_seqs = parse_input_sequences(fasta)

    test_barcode = _get_dna_from_SeqRecord(all_seqs[0])[start:end+1]
    print(_get_dna_from_SeqRecord(all_seqs[0]))
    print(f"Confirm {test_barcode} is a correct barcode, otherwise change start and end.")
    barcodes = []
    for record in all_seqs:
        seq = _get_dna_from_SeqRecord(record)
        # end is inclusive
        barcode = seq[start:end+1]
        barcodes.append(str(barcode))
        if gu_present:
            stem5 = barcode[:num_bp]
            stem3 = get_reverse_complement(barcode[-num_bp:])
            if stem5!=stem3:
                other_barcode = stem3+barcode[num_bp:-num_bp]+get_reverse_complement(stem5)
                barcodes.append(other_barcode)
    if unique:
        barcodes = list(set(barcodes))
    return sorted(barcodes)


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
                circularize=False, fraction_use=1, reverse_complement=False,
                viral_prep=False):
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
        reverse_complement (bool): instead of sequence in fasta, use
            reverse complement

    Returns:
        list of SeqRecord with the windows
        if out_fasta specified, also saves these to fasta file
        naming convention is the seqname_start-end with start and
        end inclusive and indexing from 0.
    '''

    print(f'Getting all sliding windows, {window_length}nt every {window_slide}nt.')

    # get sequences and initialize
    seqs = parse_input_sequences(fasta)
    windows = []
    unused_windows = []
    windowsA, windowsB, windowsC, windowsD, windowsE, windowsF = [],[],[],[],[],[]

    for j,seq_rec in enumerate(seqs):
        #print(j)
        # loop through sequence every window_slide nucleotides
        seq = _get_dna_from_SeqRecord(seq_rec)
        if viral_prep:
            circularize = True
            first_third = floor((1/3)*len(seq))
            second_third = floor((2/3)*len(seq))
            third_third = floor((4/3)*len(seq))-window_length+window_slide
        window_limit = floor(fraction_use*len(seq))

        for i in range(0, len(seq), window_slide):

            # when we hit the end of sequence
            if i+window_length > len(seq):

                # add last window
                if not circularize:
                    a, b = len(seq)-window_length, len(seq)
                    new_seq = seq[a:b]
                    namenum = f'{a}-{b-1}'
                    if viral_prep:
                        new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                              f'{get_name_SeqRecord(seq_rec)}_{namenum}',
                                              '', '')
                        new_recrc = SeqIO.SeqRecord(Seq.Seq(get_reverse_complement(new_seq)),
                                              f'{name}_rc', '', '')
                        current_loc = (i + (window_length-1))
                        #if current_loc < first_third:
                        #    windowsA.append(new_rec)
                        #    windowsC.append(new_recrc)
                        if current_loc < second_third:
                            windowsA.append(new_rec)
                            windowsB.append(new_recrc)
                        elif current_loc < third_third:
                            windowsC.append(new_rec)
                            windowsD.append(new_recrc)
                        else: 
                            windowsE.append(new_rec)
                            windowsF.append(new_recrc)
                    else:
                        if reverse_complement:
                            new_rec = SeqIO.SeqRecord(Seq.Seq(get_reverse_complement(new_seq)),
                                                  f'{name}_rc', '', '')
                        else:
                            new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                                  f'{get_name_SeqRecord(seq_rec)}_{namenum}',
                                                  '', '')
                        if fraction_use != 1 and ((i + (window_length-1)) > window_limit):
                            unused_windows.append(new_rec)
                        else:
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
            if len(new_seq)*0.9>new_seq.count('A')+new_seq.count('C')+new_seq.count('T')+new_seq.count('G'):
                print(f'window {new_seq} to uncertain deleting')
                new_seqs = []
            else:
                new_seqs = _fill_in_any_incomplete(new_seq,[new_seq])
            for j,new_seq in enumerate(new_seqs):
                if len(new_seqs) == 1:
                    name = f'{get_name_SeqRecord(seq_rec)}_{namenum}'
                else:
                    name = f'{get_name_SeqRecord(seq_rec)}_amb{j}_{namenum}'
                # save with name inclusive!
                if viral_prep:
                    new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                          f'{get_name_SeqRecord(seq_rec)}_{namenum}',
                                          '', '')
                    new_recrc = SeqIO.SeqRecord(Seq.Seq(get_reverse_complement(new_seq)),
                                          f'{name}_rc', '', '')
                    current_loc = (i + (window_length-1))
                    #if current_loc < first_third:
                    #    windowsA.append(new_rec)
                    #    windowsC.append(new_recrc)
                    if current_loc < second_third:
                        windowsA.append(new_rec)
                        windowsB.append(new_recrc)
                    elif current_loc < third_third:
                        windowsC.append(new_rec)
                        windowsD.append(new_recrc)
                    else:
                        windowsE.append(new_rec)
                        windowsF.append(new_recrc)
                else:
                    if reverse_complement:
                        new_rec = SeqIO.SeqRecord(Seq.Seq(get_reverse_complement(new_seq)),
                                              f'{name}_rc', '', '')
                    else:
                        new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                              f'{get_name_SeqRecord(seq_rec)}_{namenum}',
                                              '', '')
                    if fraction_use != 1 and ((i + (window_length-1)) > window_limit):
                        unused_windows.append(new_rec)
                    else:
                        windows.append(new_rec)

    # remove and save fraction unused
    if fraction_use != 1 and not viral_prep:
        unused_file = f'{out_fasta.rsplit(".",1)[0]}_unused.{out_fasta.rsplit(".",1)[1]}'
        SeqIO.write(unused_windows, unused_file, "fasta")
        print(f'Saved unused windows to {unused_file}.')
    
    # save file
    if out_fasta is not None and not viral_prep:
        SeqIO.write(windows, out_fasta, "fasta")
        print(f'Saved windows to {out_fasta}.')

    if viral_prep:
        SeqIO.write(windowsA, f'{out_fasta.rsplit(".",1)[0]}_A.{out_fasta.rsplit(".",1)[1]}','fasta')
        SeqIO.write(windowsB, f'{out_fasta.rsplit(".",1)[0]}_B.{out_fasta.rsplit(".",1)[1]}','fasta')
        SeqIO.write(windowsC, f'{out_fasta.rsplit(".",1)[0]}_C.{out_fasta.rsplit(".",1)[1]}','fasta')
        SeqIO.write(windowsD, f'{out_fasta.rsplit(".",1)[0]}_D.{out_fasta.rsplit(".",1)[1]}','fasta')
        SeqIO.write(windowsE, f'{out_fasta.rsplit(".",1)[0]}_E.{out_fasta.rsplit(".",1)[1]}','fasta')
        SeqIO.write(windowsF, f'{out_fasta.rsplit(".",1)[0]}_F.{out_fasta.rsplit(".",1)[1]}','fasta')
    return windows

def clean_fasta(fasta,out_fasta):
    seqs = parse_input_sequences(fasta)
    cleaned = []
    for seq in seqs:
        cleaned.append(SeqIO.SeqRecord(Seq.Seq(_get_dna_from_SeqRecord(seq)),get_name_SeqRecord(seq),'',''))
    SeqIO.write(cleaned,out_fasta,'fasta')

def get_name_SeqRecord(seq):
    name = str(seq.name)
    if seq.description != '<unknown description>':
        name += seq.description
    return name.replace('/','__').replace(' ','_')

def get_all_single_mutants(fasta, out_fasta=None, mutational_dict=None, bases=BASES):
    '''
    Get all single mutants from sequences in a fasta file.

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        out_fasta (str): if specified, save mutants to fasta (default None)
        mutational_dict (dict): if specified, nucleotide:list nucleotide to mutate to (default None)
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
    all_WT = parse_input_sequences(fasta)
    all_single_mutants = []

    # TODO cannot handle non ACTG?

    for record in all_WT:
        seq = _get_dna_from_SeqRecord(record)

        # at each position, get single mutants
        for i in range(len(seq)):
            if mutational_dict is None:
                base_list = bases
            else:
                base_list = mutational_dict[seq[i]]
            for mut in base_list:
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
    # print(regionAs,regionBs)

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


def get_bp_mutants(fasta, bp_sets, out_fasta=None,
                    mutate_dict={'AT':['AC','GC'],
                                 'TA':['CA','CG'],
                                 'CG':['CA','TA'],
                                 'GC':['AC','AT'],
                                 'GT':['TT','TA'],
                                 'TG':['TT','AT']},
                    other_mutate=[]):
    '''
    For each base-pair specified get mutants of that base paired as specified by mutate_dict.
    Any base pair not specified in the dictionary will be mutated to the pairs in other_mutate.
    
    # TODO ask RD what he wants as default
    # TODO document
    '''
    print(f"Getting mutants to the base pairs specified according it {mutate_dict}.")

    # get all sequences and initialize
    all_WT = parse_input_sequences(fasta)
    all_mutants = []

    # check user specified set of basepairs for every sequence in the fasta
    if len(all_WT) != len(bp_sets):
        print(f'WARNING: bps must be a list, one for each sequence in fasta, of lists of basepairs to rescue for that sequence. You have {len(all_WT)} inputted sequences and {len(bp_sets)} base-pair sets.')

    for record, bps in zip(all_WT, bp_sets):

        seq = _get_dna_from_SeqRecord(record)

        # at each base pair, get rescue mutants
        for bp in bps:
            current_bp = seq[bp[0]]+seq[bp[1]]
            if current_bp in mutate_dict:
                mutate_list = mutate_dict[current_bp]
            else:
                mutate_list = other_mutate
            for new_bp in mutate_list:
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
                all_mutants.append(new_mut)

    # save file
    if out_fasta is not None:
        SeqIO.write(all_mutants, out_fasta, "fasta")
        print(f'Saved all mutants to {out_fasta}.')

    return all_mutants

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


def delete_used_barcodes(all_barcode,used_barcodes,out_fasta):
    print('barcod')
    all_uids = list(SeqIO.parse(all_barcode, "fasta"))
    all_uids = [(_get_dna_from_SeqRecord(uid.seq),uid.name) for uid in all_uids]
    print('now sorting')
    all_uids = sorted(all_uids)
    used_barcodes = sorted(used_barcodes)
    print('got_all now filtering')
    new_uids = []
    #new_uids_filter = list(filter(lambda i: i[0] not in used_barcodes, all_uids))
    for i,uid in enumerate(all_uids):
        if i % 1000000==0:
            print(i) 
        if len(used_barcodes)>0:       
            while uid[0] > used_barcodes[0]:
                x=used_barcodes.pop(0)
                if len(used_barcodes)==0:
                    break
            if len(used_barcodes)>0:  
                if uid[0] != used_barcodes[0]:
                    new_uids.append(SeqIO.SeqRecord(Seq.Seq(uid[0]),uid[1],'',''))
        else:
            new_uids.append(SeqIO.SeqRecord(Seq.Seq(uid[0]),uid[1],'',''))

    print('now ahuffling and then writing')
    
    #for i,uid in enumerate(new_uids_filter):
    #    if i % 1000000:
    #        print(i)
    #    new_uids.append(SeqIO.SeqRecord(Seq.Seq(uid[0]),uid[1],'',''))
    shuffle(new_uids)
    print(f'deleted {len(all_uids)-len(new_uids)} barcodes')
    SeqIO.write(new_uids, out_fasta, "fasta")


def get_all_barcodes(out_fasta=None, num_bp=8, num5hang=0, num3hang=0,
                     polyA5=0, polyA3=0,
                     loop=TETRALOOP, bases=BASES,
                     used_barcodes=[],
                     shuffle_barcodes=True,
                     used_barcodes_sorted=True):
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
st = time()
get_all_barcodes(num_bp=10)
en = time()
en-st
2.1671266555786133
st = time()
get_all_barcodes(num_bp=12,out_fasta='Temp.txt')
en = time()
en-st
39
    '''

    print("Getting all possible barcodes.")
    all_barcodes = []
    # slight speed up
    # also if ordered gaurenteed can pop!!!!
    if num5hang+num3hang == 0:
        used_barcodes = [barcode[:num_bp] for barcode in used_barcodes]
    # get all possible combinations of bases for random/barcode regions
    i=0

    #if out_fasta is not None:
    #    f = open(out_fasta,'w')
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
        stemB = get_reverse_complement(stemA)

        # put all barcode parts together
        seq = ("A"*polyA5)+hang5+stemA+loop+stemB+hang3+("A"*polyA3)
        # TODO this never actually checks edit distance...
        if used_barcodes_sorted and used_barcodes!=[] and num5hang+num3hang==0:
            if seq[:num_bp] != used_barcodes[0]:
                name = f' stem{stemA}_{hang5}hang{hang3}_{polyA5}polyA{polyA3}'
                seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), name, '', '')
                all_barcodes.append(seq_rec)
            else: 
                while seq[:num_bp] == used_barcodes[0]:
                    used_barcodes.pop(0)
                    if len(used_barcodes)==0:
                        break
        else:
            if seq not in used_barcodes:
                name = f' stem{stemA}_{hang5}hang{hang3}_{polyA5}polyA{polyA3}'
                seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), name, '', '')
                all_barcodes.append(seq_rec)
        i+=1
        #if i%1000000==0: print(i,seq[:num_bp],used_barcodes[0])
    print(len(used_barcodes))
    print("shuffling barcodes")
    if shuffle_barcodes:
        shuffle(all_barcodes)
    # save
    if out_fasta is not None:
        #open(out_fasta,'w').wrte('\n'.join(all_barcodes))
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
            min_length_stem=6, max_length_stem=16,
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

    print('WARNING this code will likely be obsoletely by a code that takes other regions into account when choosing the pad.')

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
                    full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}'
                    struct_results = check_struct_bpp(full_seq,
                                                      regions['unpaired'], [regions['pairedA']],
                                                      [regions['pairedB']], regions['noninteractA'], regions['noninteractB'],
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
                                                      regions['unpaired'], [regions['pairedA']],
                                                      [regions['pairedB']], regions['noninteractA'], regions['noninteractB'],
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
                full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}'
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
                        struct_results = check_struct_bpp(full_seq, regions['unpaired'], [regions['pairedA']],
                                                          [regions['pairedB']], regions['noninteractA'], regions['noninteractB'],
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
                name = f'{get_name_SeqRecord(seq_rec)}_{len(pad5)}pad{len(pad3)}'
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
                                polyA5=num5polyA,shuffle_barcodes=True)
    if len(used_barcodes) != 0:
         if len(all_uids[0]) != len(used_barcodes[0]):
             print('ERROR: used barcodes are not the correct length')
    
    # read sequences, check all same length
    seqs = parse_input_sequences(fasta)
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
        name = get_name_SeqRecord(seq_rec)+'_libraryready'

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
        mutations = _get_mutations_from_name(get_name_SeqRecord(seq_rec))
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
                                              [region_paired_A], [region_paired_B],
                                              [regionA], [regionB],
                                              epsilon_interaction=epsilon_interaction,
                                              epsilon_punpaired=epsilon_punpaired,
                                              epsilon_avg_punpaired=epsilon_avg_punpaired,
                                              epsilon_paired=epsilon_paired,
                                              epsilon_avg_paired=epsilon_avg_paired,
                                              prob_factor=prob_factor,
                                              lines=lines,
                                              save_image=save_image, mutants=mutations)
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


def add_library_elements(fasta, out_fasta=None, 
                              bases=BASES, share_pad='same_length',
                              seq5=SEQ5, seq3=SEQ3,
                              barcode_num_bp=8, barcode_num5hang=0, barcode_num5polyA=4,
                              barcode_loop=TETRALOOP,
                              epsilon_interaction=MAXPROB_NONINTERACT,
                              epsilon_punpaired=MINPROB_UNPAIRED,
                              epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                              epsilon_paired=MINPROB_PAIRED,
                              epsilon_avg_paired=MINAVGPROB_PAIRED,
                              save_image_folder=None, save_bpp_fig=0,
                              punpaired_chunk_size=500, used_barcodes=[],
                              num_barcode_before_reduce=100,
                              percent_reduce_prob=10, min_edit=2,
                              pad_loop=TETRALOOP, pad_hang=0, pad_polyAhang=3, pad_min_num_samples=30,
            pad_max_prop_bad=0.05, pad_side='both',
            min_length_stem=6, max_length_stem=16,
            num_pads_reduce=100, pad_to_length=None,
            barcode_file=None,
            num_replicates=1,
            pad_polyAhang_other_side=0):
    '''
    From a fasta of sequences, add constant regions and barcodes
    Given a fasta of sequence pad all sequence to the same length

    Args:
        fasta (str): fasta file containing sequence to get mutants of
        out_fasta (str): if specified, save library to fasta (default None)

        share_pad (str): either all different pads (none),
            all same for sequences of same length (same_length),
            all same where sequ    Given a fasta of sequence pad all sequence to the same length
ence of different length are just truncated (all)
        pad_loop (str): constant sequence for loop of stem-loop
        pad_hang (int): distance between sequence and the stem will be this or this+1, random sequence (default 3)
        pad_polyAhang (int): distance between sequence and the stem will be this or this+1, polyA (default 3)
        pad_min_num_samples (int): minimum number of samples to test pad structure (default 30)
        pad_max_prop_bad (float): of the samples sequences, proportion that can be bad (default 0.05)
        pad_side (str): whether split pad ('both') or force to be exclusively 5' or 3'
        num_pads_reduce (int): number of pads to try before reducing probability (default 100)
        percent_reduce_prob (float): percent to reduce probabilities each time (default 10)
        pad_to_length (int): minimum length to pad sequences to        

        seq5 (str): the constant sequence to place at 5' end of all sequences
        seq3 (str): the constant sequence to place at 3' end of all sequences

        barcode_num_bp (int): number of base pairs to have in barcode stem (default 8)
        barcode_num5hang (int): number of random nucleotides to put before the stem (default 0)
        barcode_num5polyA (int): number of A to put before the barcode (stem and 5 random hang if applicable) (default 4)
        barcode_loop (str): sequence of loop to have in the hairpin

        min_length_stem (int): minimum number base-pairs to form stem
        max_length_stem (int): maximum number of base-pairs to form stem

        epsilon_interaction (float): Maximum base-pair-probability for 2 regions to be considered non-interacting
        epsilon_punpaired (float): Minimum probability unpaired for region to be unpaired
        epsilon_avg_punpaired (float): Average probability unpaired for region to be unpaired
        epsilon_paired (float): Minimum base-pair-probability for 2 regions to be considered paired
        epsilon_avg_paired (float): Average base-pair-probability for 2 regions to be considered paired

        used_barcodes (list): list of sequences not to use as barcodes (default [])
        num_barcode_before_reduce (int): for each sequence, the number of barcodes to try
            before reducing the probability thresholds by 10% (default 100)
        percent_reduce_prob (float): percent to reduce probability threshold each time (default 10)
        min_edit (int): minimum edit distance of barcodes, base-pairs and num5hang included (default 2)

        save_image_folder (str): folder to save images to
        save_bpp_fig (float): proportion of sequences to save base-pair-probability matrix figure, 0 is none, 1 is all
        punpaired_chunk_size (int): max number of sequences to plot on each p unpaired plot (default 500)

    Returns:
        list of SeqRecord which are library ready, pads, barcodes, and constant regions
        if pad added # pad# where numbers are length of pad
        if out_fasta specified, also saves these to fasta file, _libraryready appended to names        

    TODO 
        fixed pad???
        change prepare_libary over to this function
    '''

    print("IN PROGRESS CODE: combining padding and adding barcode+fixed sequences")

    # STAGE 1 find pad and add pad

    # same as padding code, but need to check if need pad, and then add a random 
    # barcdoe and the 3UT and 5UTR and do a full check

    # if image folder does not exist create it
    if save_image_folder is not None and save_image_folder!="None":
        if not os.path.exists(save_image_folder):
            print(f'{save_image_folder} did not exists, creating.')
            os.makedirs(save_image_folder)

    # get sequences and sort by length
    seqs = parse_input_sequences(fasta)
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
            pad_min_num_samples, len(seq_group)*0.01))
        number_to_select = round(number_to_select)
        selected_sec.append(sample(seq_group, k=number_to_select))
        max_bad_structs[len_seq] = floor(pad_max_prop_bad*number_to_select)
        print(f'Pad for length {len_seq} search using {number_to_select} sequences for structure check for each length.')

    # get max length of pad needed
    if pad_to_length is None:
        desired_len = max(seq_by_length.keys())
    else:
        desired_len = pad_to_length
    porp_reduce = (1-(percent_reduce_prob/100))

    if barcode_file is not None and share_pad not in ['none','same_origin']:
        print("ERROR: if using sbatch, must have share_pad=none when specifying barcode_file.")
    if num_replicates > 1 and share_pad != 'none':
        print("ERROR: multiple pad+barcode replicates currently only implemented for share_pad='none'.")
    # if want all pads to be random
    if share_pad == 'none':
        if min_edit!=2 and barcode_num5hang!=0:
                print("WARNING: the new pad+barcode code for no sharing of barcodes is only currently implemented for stem-only barcodes and edit distance of 2")

        lib = []
        current_uid = 0
        rejected_uids = []
        if barcode_file is None:
            print("Getting random pads for each sequence, not guaranteed to be unique.")

            # get and randomly shuffle all potential barcodes
            all_uids = get_all_barcodes(num_bp=barcode_num_bp, loop=barcode_loop, num5hang=barcode_num5hang,
                                        polyA5=barcode_num5polyA,shuffle_barcodes=True)
            if len(used_barcodes) != 0:
                 if len(all_uids[0]) != len(used_barcodes[0]):
                     print('ERROR: used barcodes are not the correct length')
            
        else:
            all_uids = list(SeqIO.parse(barcode_file, "fasta"))


        for seq_length, seqs in seq_by_length.items():
            # for each length get the structure of the 3' and 5' pad
            pad_length = desired_len-seq_length
            structs, regions = _get_5_3_split(pad_length, pad_hang, pad_polyAhang,
                                              min_length_stem, max_length_stem,
                                              pad_side, pad_loop, seq_length, offset=len(seq5),add_3_3=pad_polyAhang_other_side)
            if structs['3']['N']<pad_polyAhang_other_side:
                pad_polyAhang_other_side=structs['3']['N']
                #print(structs)
                #print(regions)
            # barcode structural regions
            regions_barcode = [len(seq5), len(seq5)+desired_len,
                       len(seq5)+desired_len+barcode_num5hang+barcode_num5polyA+(2*barcode_num_bp)+len(barcode_loop)]
            # barcode should not interact with sequence nor 5' or 3'
            regionA = list(range(regions_barcode[0], regions_barcode[1]))
            regionA2 = list(range(0,regions_barcode[0]))+ list(range(regions_barcode[2],regions_barcode[2]+len(seq3)))
            regionB = list(range(regions_barcode[1], regions_barcode[2]))
            # pad should also not interact with 5' and 3' regions, nor barcode
            for x in regions['noninteractB']:
                x.extend(regionB)
                x.extend(regionA2)
            regionA.extend(regionA2)
            regions['noninteractB'] = [regionA] + regions['noninteractB'] 
            regions['noninteractA'] = [regionB] + regions['noninteractA'] 
            region_unpaired = [list(range(regions_barcode[1],
                                          regions_barcode[1]+barcode_num5hang+barcode_num5polyA)),
                               list(range(regions_barcode[2]-barcode_num_bp-len(barcode_loop),
                                          regions_barcode[2]-barcode_num_bp))]
            region_paired_A = list(range(regions_barcode[1]+barcode_num5hang+barcode_num5polyA,
                                         regions_barcode[2]-barcode_num_bp-len(barcode_loop)))
            region_paired_B = list(range(regions_barcode[2]-barcode_num_bp, regions_barcode[2]))[::-1]
            regions['pairedA'] = [region_paired_A, regions['pairedA']]
            regions['pairedB'] = [region_paired_B, regions['pairedB']]
            regions['unpaired'] = region_unpaired + regions['unpaired']
            print(f"Searching for pad for group len {pad_length}.")
            print(f"Finding a {structs['5']['N']}nt 5' pad with {structs['5']['bp']}bp stem {structs['5']['hang_rand']}nt random hang {structs['5']['hang_polyA']}nt polyA hang.")
            print(f"Finding a {structs['3']['N']}nt 3' pad with {structs['3']['bp']}bp stem {structs['3']['hang_rand']}nt random hang {structs['3']['hang_polyA']}nt polyA hang.")
            #print(regions)
            # for each sequence search for good_pad
            for seq in tqdm(list(seqs)):
                mutate_polyA = False
                for rep in range(num_replicates):
                    good_pad = False
                    # if tried enough relax the probability contstraints
                    prob_factor = {}
                    seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}

                    while not good_pad:
                        # get random pad
                        #print('pad')
                        pad5 = _get_random_barcode(num_bp=structs["5"]['bp'],
                                                   num3hang=structs["5"]['hang_rand'],
                                                   polyA3=structs["5"]['hang_polyA'],
                                                   loop=structs["5"]['loop'])
                        pad3 = _get_random_barcode(num_bp=structs["3"]['bp'],
                                                   num5hang=structs["3"]['hang_rand'],
                                                   polyA5=structs["3"]['hang_polyA'],
                                                   loop=structs["3"]['loop'],
                                                   polyA3=pad_polyAhang_other_side)
                        pad5 = _get_dna_from_SeqRecord(pad5)
                        pad3 = _get_dna_from_SeqRecord(pad3)
                        #print(pad3,pad5)

                        barcode_count = 0
                        uid_good = False
                        num_pads_reduce = 3
                        # TODO random 20 barcode_count < 20 and
                        while not uid_good and barcode_count<20:
                            for type_error, count in seq_count.items():
                                mutiplier = (count // num_pads_reduce)
                                prob_factor[type_error] = porp_reduce**mutiplier
                                if (count-mutiplier) % num_pads_reduce == 0 and count != 0:
                                    seq_count[type_error] += 1
                                    print(f'For {get_name_SeqRecord(seq)}, failed to find pad from {count-mutiplier} pads because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')
                            if mutiplier >= 3 and not mutate_polyA:
                                mutate_polyA = True
                                if barcode_num5polyA != 0:
                                    barcode_num5polyA -= 1
                                    num5hang += 1
                                if structs["3"]['hang_polyA'] != 0:
                                    structs["3"]['hang_rand'] += structs["3"]['hang_polyA']
                                    structs["3"]['hang_polyA'] = 0
                                if structs["5"]['hang_polyA'] != 0:
                                    structs["5"]['hang_rand'] += structs["5"]['hang_polyA']
                                    structs["5"]['hang_polyA'] = 0
                                    seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
                                    #num_barcode = num_barcode_before_reduce*3
                                print(f'{get_name_SeqRecord(seq)}, failed to find barcodes now allowing polyA to be random sequence.')
                                # check barcode folds and barcode does not interact with sequence

                            # TODO need to check uniqueness!!! edit distance and other logic in original code
                            barcode = _get_dna_from_SeqRecord(all_uids[current_uid])

                            full_seq = seq5 + pad5 + _get_dna_from_SeqRecord(seq) + pad3 + barcode + seq3
                            if num_replicates == 1:
                                full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}_libraryready'
                            else:
                                full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}_rep{rep}_libraryready'
                            # check if structure correct, save picture if specified and chance has it
                            if False:#TODO(save_image_folder is not None) and (random() < save_bpp_fig):
                                save_image = f'{save_image_folder}/{full_seq_name}.png'
                                plot_lines = None # TODO lines
                            else:
                                save_image, plot_lines = None, None

                            #print(full_seq)
                            #print(pad3)
                            struct_results = check_struct_bpp(full_seq, regions['unpaired'], 
                                                              regions['pairedA'], regions['pairedB'], 
                                                              regions['noninteractA'], regions['noninteractB'],
                                                              epsilon_interaction=epsilon_interaction,
                                                              epsilon_punpaired=epsilon_punpaired,
                                                              epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                              epsilon_paired=epsilon_paired,
                                                              epsilon_avg_paired=epsilon_avg_paired,
                                                              save_image=save_image,prob_factor=prob_factor)

                            # if barcode fails, get new barcode and start again
                            barcode_fail = False
                            if len(struct_results['paired_fail'])>0:
                                if struct_results['paired_fail'][0][0] == 0:
                                    barcode_fail = True
                            if len(struct_results['unpaired_fail'])>0:
                                if struct_results['unpaired_fail'][0][0] < len(region_unpaired):
                                    barcode_fail = True
                            if len(struct_results['interaction_fail'])>0:
                                if struct_results['interaction_fail'][0][0] == 0:
                                    barcode_fail = True
                                
                            # if barcode is good,
                            if not barcode_fail:
                                uid_good = True
                            else:
                                rejected_uids.append(all_uids[current_uid])

                            
                            good_pad = struct_results["pass"]
                            
                            current_uid += 1
                            barcode_count += 1
                            # if looped through all barcodes, go back to previously rejected barcodes and continue loop
                            if current_uid == len(all_uids):
                                all_uids = rejected_uids
                                current_uid, rejected_uids = 0, []
                        if not good_pad:
                            seq_count = _update_seq_count(seq_count,struct_results)

                            #if not good_pad:
                            #    seq_count = _update_seq_count(seq_count,struct_results)
                            # check pad folds and pad does not interact with sequence or others

                            # if not good structure add to count
                            # stop if reached maximal bad
                            
                        # if its non interact 2 then its barcode needs revamped, otherwise pad                    
                        
                        # get full sequence and check structure


                    # if pass add final sequence
                    padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
                                                 full_seq_name, '', '')
                    lib.append(padded_seq)
        # save
        SeqIO.write(lib, out_fasta, "fasta")
        print(f'Saved all padded sequences to {out_fasta}.')

    elif share_pad == 'same_origin':
        pad_fasta = f"{out_fasta.rsplit('.',1)[0]}_pad.{out_fasta.rsplit('.',1)[1]}"
        lib = []
        pads_by_len = {}
        all_seqs_by_origin = {}
        rejected_uids = []
        if barcode_file is None:
            print("Getting random pads for each sequence, not guaranteed to be unique.")

            # get and randomly shuffle all potential barcodes
            all_uids = get_all_barcodes(num_bp=barcode_num_bp, loop=barcode_loop, num5hang=barcode_num5hang,
                                        polyA5=barcode_num5polyA,shuffle_barcodes=True)
            if len(used_barcodes) != 0:
                 if len(all_uids[0]) != len(used_barcodes[0]):
                     print('ERROR: used barcodes are not the correct length')
            
        else:
            all_uids = list(SeqIO.parse(barcode_file, "fasta"))
        current_uid =0
        for len_group,seqs in seq_by_length.items():
            # get length of pad for this group, if 0 done
            length_to_add = desired_len-len_group
            pads_by_len[len_group] = {}

            # get groups by origin
            seqs_by_origin = {}
            for seq in seqs:
                origin = get_origin_seq_from_name(seq.name)
                if origin in seqs_by_origin:
                    seqs_by_origin[origin].append(seq)
                else:
                    seqs_by_origin[origin] = [seq]
            all_seqs_by_origin[len_group] = seqs_by_origin
            if length_to_add == 0:
                for origin in seqs_by_origin.keys():
                    pads_by_len[len_group][origin] = ['', '']
                continue
            

            # split length of pad to either end of the sequence and get regions
            structs, regions = _get_5_3_split(length_to_add, pad_hang, pad_polyAhang,
                                              min_length_stem, max_length_stem,
                                              pad_side, pad_loop, len_group,offset=len(seq5),add_3_3=pad_polyAhang_other_side)
            #print(structs)
            # barcode structural regions
            regions_barcode = [len(seq5), len(seq5)+desired_len,
                       len(seq5)+desired_len+barcode_num5hang+barcode_num5polyA+(2*barcode_num_bp)+len(barcode_loop)]
            # barcode should not interact with sequence nor 5' or 3'
            regionA = list(range(regions_barcode[0], regions_barcode[1]))
            regionA2 = list(range(0,regions_barcode[0]))+ list(range(regions_barcode[2],regions_barcode[2]+len(seq3)))
            regionB = list(range(regions_barcode[1], regions_barcode[2]))
            # pad should also not interact with 5' and 3' regions, nor barcode
            for x in regions['noninteractB']:
                x.extend(regionB)
                x.extend(regionA2)
            regionA.extend(regionA2)
            regions['noninteractB'] = [regionA] + regions['noninteractB'] 
            regions['noninteractA'] = [regionB] + regions['noninteractA'] 
            region_unpaired = [list(range(regions_barcode[1],
                                          regions_barcode[1]+barcode_num5hang+barcode_num5polyA)),
                               list(range(regions_barcode[2]-barcode_num_bp-len(barcode_loop),
                                          regions_barcode[2]-barcode_num_bp))]
            region_paired_A = list(range(regions_barcode[1]+barcode_num5hang+barcode_num5polyA,
                                         regions_barcode[2]-barcode_num_bp-len(barcode_loop)))
            region_paired_B = list(range(regions_barcode[2]-barcode_num_bp, regions_barcode[2]))[::-1]
            regions['pairedA'] = [region_paired_A, regions['pairedA'], ]
            regions['pairedB'] = [region_paired_B, regions['pairedB']]
            regions['unpaired'] = region_unpaired + regions['unpaired']
            #print(regions)
            print(f"Finding a {structs['5']['N']}nt 5' pad with {structs['5']['bp']}bp stem {structs['5']['hang_rand']}nt random hang {structs['5']['hang_polyA']}nt polyA hang.")
            print(f"Finding a {structs['3']['N']}nt 3' pad with {structs['3']['bp']}bp stem {structs['3']['hang_rand']}nt random hang {structs['3']['hang_polyA']}nt polyA hang.")
            #mutate_polyA = False
            for origin,seqs in seqs_by_origin.items():
                mutate_polyA=False
                print(f"Searching for pad for group len {len_group}, origin {origin}.")
                #for seq in tqdm(seqs):
                #    # TODO TODO

                # loop through to find pad that works
                good_pad = False
                prob_factor = {}
                seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
                if not os.path.exists(f'{save_image_folder}/temp'):
                    os.makedirs(f'{save_image_folder}/temp')
                else:
                    for f in glob(f'{save_image_folder}/temp/*'):
                        os.remove(f)
                max_bad_structs = len(seqs)*0.25 # magic number TODO
                num_pads_reduce =50 # TODO
                # current_pad = 0
                while not good_pad:
                    potential_lib = []
                    # if tried enough relax the probability contstraints
                    # TODO
                    print('pad')
                    # get random pad
                    pad5 = _get_random_barcode(num_bp=structs["5"]['bp'],
                                               num3hang=structs["5"]['hang_rand'],
                                               polyA3=structs["5"]['hang_polyA'],
                                               loop=structs["5"]['loop'])
                    pad3 = _get_random_barcode(num_bp=structs["3"]['bp'],
                                               num5hang=structs["3"]['hang_rand'],
                                               polyA5=structs["3"]['hang_polyA'],
                                               loop=structs["3"]['loop'],
                                               polyA3=pad_polyAhang_other_side)

                    # chek all samples sequences
                    bad_count = 0
                    for i, seq in enumerate(seqs): 
                        if i % 2 == 0 and i != 0:
                            print(i)
                        
                        barcode_count = 0
                        uid_good = False
                        while barcode_count < 20 and not uid_good:
                            # more barcode search... TODO
                            for type_error, count in seq_count.items():
                                mutiplier = (count // num_pads_reduce)
                                prob_factor[type_error] = porp_reduce**mutiplier
                                if (count-mutiplier) % num_pads_reduce == 0 and count != 0:
                                    seq_count[type_error] += 1
                                    print(f'For {len_group}, failed to find pad from {count-mutiplier} pads because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')

                            # when this is 3x, if there is a ployA allow this to mutate
                            if mutiplier >= 3 and not mutate_polyA:
                                mutate_polyA = True
                                if barcode_num5polyA != 0:
                                    barcode_num5polyA -= 1
                                    num5hang += 1
                                if structs["3"]['hang_polyA'] != 0:
                                    structs["3"]['hang_rand'] += structs["3"]['hang_polyA']
                                    structs["3"]['hang_polyA'] = 0
                                if structs["5"]['hang_polyA'] != 0:
                                    structs["5"]['hang_rand'] += structs["5"]['hang_polyA']
                                    structs["5"]['hang_polyA'] = 0
                                seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
                                #num_barcode = num_barcode_before_reduce*3
                                print(f'{get_name_SeqRecord(seq)}, failed to find barcodes now allowing polyA to be random sequence.')

                            # check barcode folds and barcode does not interact with sequence
                            #barcode = _get_random_barcode(num_bp=barcode_num_bp,
                            #                          loop=barcode_loop, 
                            #                          num5hang=barcode_num5hang,
                            #                          polyA5=barcode_num5polyA)
                            barcode = _get_dna_from_SeqRecord(all_uids[current_uid])
                            if mutate_polyA:
                                new_hang = _get_random_barcode(
                                    num_bp=0, num5hang=barcode_num5polyA, loop='').seq
                                # new_loop = _get_random_barcode(num_bp=0, num5hang=len(loop),loop='').seq
                                # uid = new_hang+uid[num5polyA:num5polyA+num_bp]+new_loop+uid[num5polyA+num_bp+len(loop):]
                                barcode = new_hang+barcode[barcode_num5polyA:]
                            full_seq = seq5 + pad5 + _get_dna_from_SeqRecord(seq) + pad3 + barcode + seq3
                            #print(full_seq.seq)
                            full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}_libraryready'
                            #print(pad5,'5')
                            #print(pad3,'3')
                            #print(full_seq.seq,print(regions))
                            # check if structure correct, save picture if specified and chance has it
                            if False: #TDO (save_image_folder is not None) and (random() < save_bpp_fig):
                                save_image = f'{save_image_folder}/temp/{full_seq_name}.png'
                                plot_lines = None # TODO lines

                            else:
                                save_image, plot_lines = None, None
                            # check if structure correct, save picture if specified and chance has it
                            struct_results = check_struct_bpp(full_seq, regions['unpaired'], 
                                                              regions['pairedA'], regions['pairedB'], 
                                                              regions['noninteractA'], regions['noninteractB'],
                                                              epsilon_interaction=epsilon_interaction,
                                                              epsilon_punpaired=epsilon_punpaired,
                                                              epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                              epsilon_paired=epsilon_paired,
                                                              epsilon_avg_paired=epsilon_avg_paired,
                                                              save_image=save_image,prob_factor=prob_factor)

                            # if barcode fails, get new barcode and start again
                            barcode_fail = False
                            if len(struct_results['paired_fail'])>0:
                                if struct_results['paired_fail'][0][0] == 0:
                                    barcode_fail = True
                            if len(struct_results['unpaired_fail'])>0:
                                if struct_results['unpaired_fail'][0][0] < len(region_unpaired):
                                    barcode_fail = True
                            if len(struct_results['interaction_fail'])>0:
                                if struct_results['interaction_fail'][0][0] == 0:
                                    barcode_fail = True
                                
                            # if barcode is good,
                            if not barcode_fail:
                                uid_good = True
                            current_uid += 1
                            good_pad = struct_results["pass"]
                            barcode_count += 1
                            #if not good_pad:
                            #    seq_count = _update_seq_count(seq_count,struct_results)
                            # check pad folds and pad does not interact with sequence or others

                            # if not good structure add to count
                            # stop if reached maximal bad
                        if not good_pad:
                            bad_count += 1
                            seq_count = _update_seq_count(seq_count,struct_results)
                            if bad_count/(i+1) >= max_bad_structs/len(seqs):
                                break
                        #else:
                        #    # good so add to list
                        full_seq = _get_dna_from_SeqRecord(full_seq)
                        padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
                                                 full_seq_name, '', '')
                        potential_lib.append(padded_seq)
                    if len(potential_lib)==len(seqs):
                        good_pad = True

                lib.extend(potential_lib)
                #for f in glob(f'{save_image_folder}/temp/*'):
                os.system(f"mv {save_image_folder}/temp/* ../")#{'/'.join(f.split('/')[:-2]+f.split('/')[-1:])}")
        # save
        SeqIO.write(lib, out_fasta, "fasta")
        print(f'Saved all padded sequences to {out_fasta}.')
        #        pads_by_len[len_group][origin] = [_get_dna_from_SeqRecord(
        #            pad5), _get_dna_from_SeqRecord(pad3)]

        #print('Found pads, now adding pads.')
        #padded_seqs = []

        #pads_by_len = {}
        #all_seqs_by_origin = {}
        #for seq_length, orign_groups in all_seqs_by_origin.items():
        #    for origin, seqs in all_seqs_by_origin[seq_length].items():
        #        pad5, pad3 = pads_by_len[seq_length][origin]
        #        for seq in list(seqs):
        #            full_seq = pad5+_get_dna_from_SeqRecord(seq)+pad3
        #            full_seq_name = f'{seq.name}_{len(pad5)}pad{len(pad3)}'
        #            padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
        #                                         full_seq_name, '', '')
        #            padded_seqs.append(padded_seq)

        # TODO low key just get barcodes as I go? --> faster
        # accept a list of seqs instead of fasta.
        #SeqIO.write(padded_seqs, pad_fasta, "fasta")
        #print(f'Saved all padded sequences to {pad_fasta}.')


        # STAGE 2 with pads now add barcodes
        #lib = add_fixed_seq_and_barcode(fasta=pad_fasta, out_fasta=out_fasta, seq5=seq5, seq3=seq3,
        #                          num_bp=barcode_num_bp, num5hang=barcode_num5hang, num5polyA=barcode_num5polyA,
        #                          loop=barcode_loop,
        #                          epsilon_interaction=epsilon_interaction,
        #                          epsilon_punpaired=epsilon_punpaired,
        #                          epsilon_avg_punpaired=epsilon_avg_punpaired,
        #                          epsilon_paired=epsilon_paired,
        #                          epsilon_avg_paired=epsilon_avg_paired,
        #                          save_image_folder=save_image_folder, save_bpp_fig=save_bpp_fig,
        #                          punpaired_chunk_size=punpaired_chunk_size, used_barcodes=used_barcodes,
        #                          num_barcode_before_reduce=num_barcode_before_reduce,
        #                          percent_reduce_prob=percent_reduce_prob, min_edit=min_edit)

    # if want sequences of same length to share same pad
    elif share_pad == 'same_length':
        pad_fasta = f"{out_fasta.rsplit('.',1)[0]}_pad.{out_fasta.rsplit('.',1)[1]}"

        pads_by_len = {}
        for seqs in selected_sec:
            # get length of pad for this group, if 0 done
            len_group = len(str(seqs[0].seq))
            length_to_add = desired_len-len_group
            if length_to_add == 0:
                pads_by_len[len_group] = ['', '']
                continue

            print(f"Searching for pad for group len {len_group}.")

            # split length of pad to either end of the sequence and get regions
            structs, regions = _get_5_3_split(length_to_add, pad_hang, pad_polyAhang,
                                              min_length_stem, max_length_stem,
                                              pad_side, pad_loop, len_group,offset=len(seq5),add_3_3=pad_polyAhang_other_side)

            # barcode structural regions
            regions_barcode = [len(seq5), len(seq5)+desired_len,
                       len(seq5)+desired_len+barcode_num5hang+barcode_num5polyA+(2*barcode_num_bp)+len(barcode_loop)]
            # barcode should not interact with sequence nor 5' or 3'
            regionA = list(range(regions_barcode[0], regions_barcode[1]))
            regionA2 = list(range(0,regions_barcode[0]))+ list(range(regions_barcode[2],regions_barcode[2]+len(seq3)))
            regionB = list(range(regions_barcode[1], regions_barcode[2]))
            # pad should also not interact with 5' and 3' regions, nor barcode
            for x in regions['noninteractB']:
                x.extend(regionB)
                x.extend(regionA2)
            regionA.extend(regionA2)
            regions['noninteractB'] = [regionA] + regions['noninteractB'] 
            regions['noninteractA'] = [regionB] + regions['noninteractA'] 
            region_unpaired = [list(range(regions_barcode[1],
                                          regions_barcode[1]+barcode_num5hang+barcode_num5polyA)),
                               list(range(regions_barcode[2]-barcode_num_bp-len(barcode_loop),
                                          regions_barcode[2]-barcode_num_bp))]
            region_paired_A = list(range(regions_barcode[1]+barcode_num5hang+barcode_num5polyA,
                                         regions_barcode[2]-barcode_num_bp-len(barcode_loop)))
            region_paired_B = list(range(regions_barcode[2]-barcode_num_bp, regions_barcode[2]))[::-1]
            regions['pairedA'] = [region_paired_A, regions['pairedA'], ]
            regions['pairedB'] = [region_paired_B, regions['pairedB']]
            regions['unpaired'] = region_unpaired + regions['unpaired']

            print(f"Finding a {structs['5']['N']}nt 5' pad with {structs['5']['bp']}bp stem {structs['5']['hang_rand']}nt random hang {structs['5']['hang_polyA']}nt polyA hang.")
            print(f"Finding a {structs['3']['N']}nt 3' pad with {structs['3']['bp']}bp stem {structs['3']['hang_rand']}nt random hang {structs['3']['hang_polyA']}nt polyA hang.")

            # loop through to find pad that works
            good_pad = False
            prob_factor = {}
            seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
            # current_pad = 0
            while not good_pad:

                # if tried enough relax the probability contstraints
                
                
                # get random pad
                pad5 = _get_random_barcode(num_bp=structs["5"]['bp'],
                                           num3hang=structs["5"]['hang_rand'],
                                           polyA3=structs["5"]['hang_polyA'],
                                           loop=structs["5"]['loop'])
                pad3 = _get_random_barcode(num_bp=structs["3"]['bp'],
                                           num5hang=structs["3"]['hang_rand'],
                                           polyA5=structs["3"]['hang_polyA'],
                                           loop=structs["3"]['loop'],
                                           polyA3=pad_polyAhang_other_side)

                # chek all samples sequences
                bad_count = 0
                for i, seq in enumerate(seqs):
                    #if i % 50 == 0 and i != 0:
                    #    print(i)
                    

                    barcode_count = 0
                    uid_good = False
                    while barcode_count < 20 and not uid_good:
                        for type_error, count in seq_count.items():
                            mutiplier = (count // num_pads_reduce)
                            prob_factor[type_error] = porp_reduce**mutiplier
                            if (count-mutiplier) % num_pads_reduce == 0 and count != 0:
                                seq_count[type_error] += 1
                                print(f'For {len_group}, failed to find pad from {count-mutiplier} pads because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')
                        # check barcode folds and barcode does not interact with sequence
                        barcode = _get_random_barcode(num_bp=barcode_num_bp,
                                                  loop=barcode_loop, 
                                                  num5hang=barcode_num5hang,
                                                  polyA5=barcode_num5polyA)
                        full_seq = seq5 + pad5 + _get_dna_from_SeqRecord(seq) + pad3 + barcode + seq3

                        # check if structure correct, save picture if specified and chance has it
                        struct_results = check_struct_bpp(full_seq, regions['unpaired'], 
                                                          regions['pairedA'], regions['pairedB'], 
                                                          regions['noninteractA'], regions['noninteractB'],
                                                          epsilon_interaction=epsilon_interaction,
                                                          epsilon_punpaired=epsilon_punpaired,
                                                          epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                          epsilon_paired=epsilon_paired,
                                                          epsilon_avg_paired=epsilon_avg_paired)

                        # if barcode fails, get new barcode and start again
                        barcode_fail = False
                        if len(struct_results['paired_fail'])>0:
                            if struct_results['paired_fail'][0][0] == 0:
                                barcode_fail = True
                        if len(struct_results['unpaired_fail'])>0:
                            if struct_results['unpaired_fail'][0][0] < len(region_unpaired):
                                barcode_fail = True
                        if len(struct_results['interaction_fail'])>0:
                            if struct_results['interaction_fail'][0][0] == 0:
                                barcode_fail = True
                            
                        # if barcode is good,
                        if not barcode_fail:
                            uid_good = True
                        
                        good_pad = struct_results["pass"]
                        barcode_count += 1
                    
                        # check pad folds and pad does not interact with sequence or others

                        # if not good structure add to count
                        # stop if reached maximal bad
                        if not good_pad:
                            seq_count = _update_seq_count(seq_count,struct_results)
                    if not good_pad:
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
                full_seq_name = f'{get_name_SeqRecord(seq)}_{len(pad5)}pad{len(pad3)}'
                padded_seq = SeqIO.SeqRecord(Seq.Seq(full_seq),
                                             full_seq_name, '', '')
                padded_seqs.append(padded_seq)

        # save TODO adding barcode function should probably
        # accept a list of seqs instead of fasta.
        SeqIO.write(padded_seqs, pad_fasta, "fasta")
        print(f'Saved all padded sequences to {pad_fasta}.')


        # STAGE 2 with pads now add barcodes
        lib = add_fixed_seq_and_barcode(fasta=pad_fasta, out_fasta=out_fasta, seq5=seq5, seq3=seq3,
                                  num_bp=barcode_num_bp, num5hang=barcode_num5hang, num5polyA=barcode_num5polyA,
                                  loop=barcode_loop,
                                  epsilon_interaction=epsilon_interaction,
                                  epsilon_punpaired=epsilon_punpaired,
                                  epsilon_avg_punpaired=epsilon_avg_punpaired,
                                  epsilon_paired=epsilon_paired,
                                  epsilon_avg_paired=epsilon_avg_paired,
                                  save_image_folder=save_image_folder, save_bpp_fig=save_bpp_fig,
                                  punpaired_chunk_size=punpaired_chunk_size, used_barcodes=used_barcodes,
                                  num_barcode_before_reduce=num_barcode_before_reduce,
                                  percent_reduce_prob=percent_reduce_prob, min_edit=min_edit)

    # if want all sequences to share same pad (truncated as needed)
    elif share_pad == 'all':
        print("ERROR not implemented in this function see add_pad where it is implemented")

    else:
        print("ERROR share_pad option not recognized.")

    


    return lib


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
                     lines=[], save_image=None, mutants=[]):
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
        for i,region_unpaired in enumerate(regions_unpaired):
            if region_unpaired != []:
                punpaired_check = p_unpaired[region_unpaired]
                if ((punpaired_check.mean() < epsilon_avg_punpaired*prob_factor['unpaired']) or
                        (punpaired_check.min() < epsilon_punpaired*prob_factor['unpaired'])):
                    results["unpaired_fail"].append([i,region_unpaired])

    if regionA is not None and regionA != [] and regionB != []:
        for i,(regA,regB) in enumerate(zip(regionA,regionB)):
            if regA != [] and regB != []:
                interaction_check = bpp[regA][:, regB]
                if (interaction_check.max() > epsilon_interaction/prob_factor['interaction']):
                    results["interaction_fail"].append([i,regA, regB])

    if region_paired_A is not None and region_paired_A != []:
        for i,(regA, regB) in enumerate(zip(region_paired_A,region_paired_B)):
            if regA != [] and regB != []:
                paired_check = bpp[regA, regB]
                if ((paired_check.mean() < epsilon_avg_paired*prob_factor['paired']) or
                        (paired_check.min() < epsilon_paired*prob_factor['paired'])):
                    results["paired_fail"].append([i,regA, regB])

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
            for l in line:
                plt.vlines(l-0.5, i+0.5, i-0.5,
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
        save_image = f'{save_image_folder}/{get_name_SeqRecord(seq_rec)}.png'
        plot_bpp(bpp, seq, save_image)


###############################################################################
# helpers
###############################################################################


def get_reverse_complement(seq):
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

    names = [get_name_SeqRecord(n) for n in seqsB]
    good_seqs = []
    for seq_rec in seqsA:
        if get_name_SeqRecord(seq_rec) not in names:
            good_seqs.append(seq_rec)
    return good_seqs


def _get_dna_from_SeqRecord(seqrecord):
    '''
    From a sequence, string or SeqRecord, return dna str version
    '''

    if type(seqrecord) == str:
        dna = seqrecord.upper().replace("U", "T")
    elif type(seqrecord) == Seq.Seq:
        dna = str(seqrecord).upper().replace("U", "T")
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
    stemB = get_reverse_complement(stemA)

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
            if '-' in name_part:
                subparts = name_part.split('-')
                if len(subparts)==2:
                    nucnum_A, nuc_B = name_part.split('-')
                    num, nuc_A = nucnum_A[:-1], nucnum_A[-1]
                    if nuc_A in BASES and nuc_B in BASES:
                        nucs.append(int(num))
    return nucs


def _get_pad_length_from_name(name):
    '''
    using the naming convention _#N-N_ where # is nucleotide number
    N are the nucleotides native and mutant, find these and return the number
    '''

    nucs = []
    for name_part in name.split('_'):
        if 'pad' in name_part:
            pad5, pad3 = name_part.split('pad')
            return int(pad5),int(pad3)
    return 0,0


def get_origin_seq_from_name(name):
    '''
    using the naming convention where the original sequence name
    is followed by names added by this code (eg _#N-N_ for mutant,
    _#pad#_ for pad)

    WARNING consider each window (_#-#_) a different origin
    '''

    origin_name = []
    name_parts = '_'.join(name.split()).split('_')
    reach_library_part = False
    while not reach_library_part and len(name_parts)>0:
        next_part = name_parts.pop(0)
        if '-' in next_part:
            subparts = next_part.split('-')
            if len(subparts)==2:
                if subparts[0][-1] in BASES and subparts[1] in BASES and subparts[0][:-1].isnumeric():
                    reach_library_part = True 
                # uncomment if ever want to consider different windows the same origin
                #if subparts[0].isnumeric() and subparts[1].isnumeric():
                #    reach_library_part = True
        if 'pad' in next_part:
            subparts = next_part.split('pad')
            if len(subparts)==2:
                if subparts[0].isnumeric() and subparts[1].isnumeric():
                    reach_library_part = True
        if next_part == 'libraryready':
            reach_library_part = True
        if not reach_library_part:
            origin_name.append(next_part)
    return '_'.join(origin_name)


def XXget_5_3_split(length, hang, polyAhang, min_length_stem, max_length_stem, pad_side, loop, seq_len,offset=0,add_3_3=0):

    # initialize variables
    loop_len = len(loop)
    initial_structure = {"N": 0, 'bp': 0, 'hang_rand': hang,
                         'hang_polyA': polyAhang, 'loop': loop}
    structs = {"5": initial_structure.copy(),
               "3": initial_structure.copy()}
    min_pad_for_stem = min_length_stem*2 + loop_len + hang + polyAhang
    max_pad_for_stem = max_length_stem*2 + loop_len + hang + polyAhang
    unpaired_length = loop_len+hang+polyAhang+add_3_3
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
            structs["3"]['N'] = (length)-add_3_3
            structs["3"]['bp'] = num_bp
        elif pad_side == "5'":
            structs["5"]['N'] = length
            structs["5"]['bp'] = num_bp

    # if pad both, make the split
    elif pad_side == 'both':
        # if too short for 1 stem, split into 2 unstructured regions
        if length < min_pad_for_stem:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = (length//2)-add_3_3
            if length % 2 != 0:
                structs["3"]['N'] += 1

        # if enough for 1 stem but not too much for 1 stem, make just 1 stem
        elif length < max_pad_for_stem:
            structs["3"]['N'] = length -add_3_3
            structs["3"]['bp'] = (length-unpaired_length)//2

        # if too much for 1 stem but not enough for 2 stems, make 1 stem as long as possible
        elif length < 2*min_pad_for_stem:
            structs["5"]['N'] = length-max_pad_for_stem
            structs["3"]['N'] = max_pad_for_stem-add_3_3
            structs["3"]['bp'] = max_length_stem

        # if enough for 2 stems, make 2
        else:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = (length//2)-add_3_3
            if length % 2 != 0:
                structs["3"]['N'] += 1
            structs["5"]['bp'] = (structs["5"]['N']-unpaired_length)//2
            structs["3"]['bp'] = (structs["3"]['N']-unpaired_length)//2
            if (structs["5"]['bp'] > max_length_stem) or (structs["3"]['bp'] > max_length_stem):
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
            #if struct['hang_polyA']#struct['hang_rand'] == 0 and struct['hang_polyA'] != 0:
            #    struct['hang_polyA'] = struct['N'] - \
            #        ((struct['bp']*2)+len(struct['loop']))
            #else:
            struct['hang_rand'] = struct['N'] - \
                    ((struct['bp']*2)+struct['hang_polyA']+len(struct['loop']))

    

    # get parts that should and shouldn't be paired
    regions = {}
    parts = [offset,structs["5"]['bp'], len(structs["5"]['loop']), structs["5"]['bp'],
             structs["5"]['hang_polyA'] + structs["5"]['hang_rand'],
             seq_len, structs["3"]['hang_polyA'] + structs["3"]['hang_rand'],
             structs["3"]['bp'], len(structs["3"]['loop']), structs["3"]['bp']]

    regions['unpaired'] = [list(range(sum(parts[:2]),
                                      sum(parts[:3]))),
                           list(range(sum(parts[:4]),
                                      sum(parts[:5]))),
                           list(range(sum(parts[:6]),
                                      sum(parts[:7]))),
                           list(range(sum(parts[:8]),
                                      sum(parts[:9])))]
    regions['pairedA'] = list(range(sum(parts[:1]),
                                    sum(parts[:2])))
    regions['pairedB'] = list(range(sum(parts[:3]),
                                    sum(parts[:4])))[::-1]
    # 5' pad and 3' pad do not interact with sequence
    if add_3_3 != 0:
        #structs["3"]['N'] += add_3_3
        regions['unpaired'].append(list(range(sum(parts[:10]),sum(parts[:10])+add_3_3)))
        regions['noninteractA'] = [list(range(sum(parts[:1]),
                                             sum(parts[:5]))),
                                   list(range(sum(parts[:6]),
                                                  sum(parts[:10])+add_3_3))]
        regions['noninteractB'] = [list(range(sum(parts[:5]),
                                             sum(parts[:10])+add_3_3)),
                                list(range(sum(parts[:1]),
                                             sum(parts[:6])))]
    else:
        regions['noninteractA'] = [list(range(sum(parts[:1]),
                                             sum(parts[:5]))),
                                   list(range(sum(parts[:6]),
                                                  sum(parts[:10])))]
        regions['noninteractB'] = [list(range(sum(parts[:5]),
                                             sum(parts[:10]))),
                                list(range(sum(parts[:1]),
                                             sum(parts[:6])))]
    regions['pairedA'].extend(list(range(sum(parts[:7]),
                                         sum(parts[:8]))))
    regions['pairedB'].extend(list(range(sum(parts[:9]),
                                         sum(parts[:10])))[::-1])
    

    return structs, regions


def _get_5_3_split(length, hang, polyAhang, min_length_stem, max_length_stem, pad_side, loop, seq_len,offset=0,add_3_3=0):

    # initialize variables
    loop_len = len(loop)
    initial_structure = {"N": 0, 'bp': 0, 'hang_rand': hang,
                         'hang_polyA': polyAhang, 'loop': loop}
    structs = {"5": initial_structure.copy(),
               "3": initial_structure.copy()}
    min_pad_for_stem = min_length_stem*2 + loop_len + hang + polyAhang
    max_pad_for_stem = max_length_stem*2 + loop_len + hang + polyAhang
    unpaired_length = loop_len+hang+polyAhang+add_3_3
    # if only pad one side
    if pad_side == "3'" or pad_side == "5'":
        # make a helix if long enough
        if length -add_3_3>= min_pad_for_stem:
            num_bp = (length-unpaired_length)//2
            if num_bp > max_length_stem:
                print("WARNING: stem too long, fix not implemented.")
        else:
            num_bp = 0

        # assign to correct side
        if pad_side == "3'":
            structs["3"]['N'] = (length)-add_3_3
            structs["3"]['bp'] = num_bp
        elif pad_side == "5'":
            structs["5"]['N'] = length
            structs["5"]['bp'] = num_bp

    # if pad both, make the split
    elif pad_side == 'both':
        # if too short for 1 stem, split into 2 unstructured regions
        if length-add_3_3 < min_pad_for_stem:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = (length//2)-add_3_3
            if length % 2 != 0:
                structs["3"]['N'] += 1

        # if enough for 1 stem but not too much for 2 stema, make just 1 stem
        elif length-add_3_3 < max_pad_for_stem:
            structs["3"]['N'] = length -add_3_3
            structs["3"]['bp'] = (length-unpaired_length)//2

        # if too much for 1 stem but not enough for 2 stems, make 1 stem as long as possible
        elif length -add_3_3< 2*min_pad_for_stem:
            structs["5"]['N'] = length-max_pad_for_stem
            structs["3"]['N'] = max_pad_for_stem-add_3_3
            structs["3"]['bp'] = max_length_stem

        # if enough for 2 stems, make 2
        else:
            structs["5"]['N'] = length//2
            structs["3"]['N'] = (length//2)-add_3_3
            if length % 2 != 0:
                structs["3"]['N'] += 1
            structs["5"]['bp'] = (structs["5"]['N']-unpaired_length)//2
            structs["3"]['bp'] = (structs["3"]['N']-unpaired_length)//2
            if (structs["5"]['bp'] > max_length_stem) or (structs["3"]['bp'] > max_length_stem):
                print("WARNING: stem too long, fix not implemented.")
    else:
        print("ERROR: pad_side not recognized.")

    # get hang and loop lengths
    for side, struct in structs.items():
        if struct['bp'] == 0:
            struct['loop'] = ''
        if struct['N'] <= 0:
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
    parts = [offset,structs["5"]['bp'], len(structs["5"]['loop']), structs["5"]['bp'],
             structs["5"]['hang_polyA'] + structs["5"]['hang_rand'],
             seq_len, structs["3"]['hang_polyA'] + structs["3"]['hang_rand'],
             structs["3"]['bp'], len(structs["3"]['loop']), structs["3"]['bp']]

    regions['unpaired'] = [list(range(sum(parts[:2]),
                                      sum(parts[:3]))),
                           list(range(sum(parts[:4]),
                                      sum(parts[:5]))),
                           list(range(sum(parts[:6]),
                                      sum(parts[:7]))),
                           list(range(sum(parts[:8]),
                                      sum(parts[:9])))]
    regions['pairedA'] = list(range(sum(parts[:1]),
                                    sum(parts[:2])))
    regions['pairedB'] = list(range(sum(parts[:3]),
                                    sum(parts[:4])))[::-1]
    # 5' pad and 3' pad do not interact with sequence
    if add_3_3 != 0:
        print("FD")
        structs["3"]['N'] += add_3_3
        to_Add = min(add_3_3,structs["3"]['N'] )
        regions['unpaired'].append(list(range(sum(parts[:10]),sum(parts[:10])+to_Add)))
        regions['noninteractA'] = [list(range(sum(parts[:1]),
                                             sum(parts[:5]))),
                                   list(range(sum(parts[:6]),
                                                  sum(parts[:10])+to_Add))]
        regions['noninteractB'] = [list(range(sum(parts[:5]),
                                             sum(parts[:10])+to_Add)),
                                list(range(sum(parts[:1]),
                                             sum(parts[:6])))]
    else:
        regions['noninteractA'] = [list(range(sum(parts[:1]),
                                             sum(parts[:5]))),
                                   list(range(sum(parts[:6]),
                                                  sum(parts[:10])))]
        regions['noninteractB'] = [list(range(sum(parts[:5]),
                                             sum(parts[:10]))),
                                list(range(sum(parts[:1]),
                                             sum(parts[:6])))]
    regions['pairedA'].extend(list(range(sum(parts[:7]),
                                         sum(parts[:8]))))
    regions['pairedB'].extend(list(range(sum(parts[:9]),
                                         sum(parts[:10])))[::-1])
    

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
                current_regions[part] = [region,region]
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


def plot_punpaired_from_fasta(fasta, save_image,lenbarcode,len5=len(SEQ5),len3=len(SEQ3),max_per=500,startN=0,just_one=False,plot_lines_from_name=True):
    # NOT well tested
    seqs = list(SeqIO.parse(fasta, "fasta"))
    if just_one:
        endN = startN+1
    else:
        endN = 1+(len(seqs)//max_per)
    for N in range(startN,endN):
        if not os.path.isfile(f'{save_image}_{N}.png'):

            print("Working on image",N)
            p_unpaireds = {}
            seqs_list = []
            muts = []
            pad_lines = []
            for seq_rec in seqs[N*max_per:min(len(seqs),(N+1)*max_per)]:
                seq = str(seq_rec.seq).upper().replace("T", "U")
                bpp = bpps(seq,
                           package='eternafold')
                p_unpaired = 1-bpp.sum(axis=0)
                p_unpaireds[get_name_SeqRecord(seq_rec)] = p_unpaired
                seqs_list.append(seq)
                # TODO bug in reading pads need to check it can be int iffed, otherwise not padd name
                if plot_lines_from_name:
                    pad5,pad3 = _get_pad_length_from_name(get_name_SeqRecord(seq_rec))
                    pad5 += len5
                    pad3 = len(seq) - pad3 - len3 - lenbarcode
                    barcode = len(seq) - len3 - lenbarcode
                    pad_lines.append([len5,pad5,pad3,barcode,len(seq)-len3])
                    mutations = _get_mutations_from_name(get_name_SeqRecord(seq_rec))
                    mutations = [x+pad5 for x in mutations]
                    muts.append(mutations)
            labels = []
            for i in range(len(seqs[0])):
                if i % 10 == 0:
                    labels.append(i)
                else:
                    labels.append('')
            plot_punpaired(p_unpaireds, labels, seqs[N*max_per:min(len(seqs),(N+1)*max_per)], muts, [], pad_lines, f'{save_image}_{N}.png')

#add_library_elements('examples/m2seq_ex_output/example_single_mut.fasta', out_fasta='test.fasta',share_pad='none',save_image_folder='test',save_bpp_fig=1)

#add_library_elements('examples/m2seq_ex_output/example_single_mut.fasta', out_fasta='test.fasta',share_pad='same_length',save_image_folder='test',save_bpp_fig=1)
