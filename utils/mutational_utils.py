import numpy as np
import pandas as pd
from itertools import product,chain
from random import shuffle, sample, choices
from tqdm import tqdm
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
from arnie.bpps import bpps


###############################################################################
# TODO
###############################################################################
# all docstrings
# README with example calls - windows, m2seq, and also a single script for RD
# check all seq high probs stay high (?)
    # all above MINAVGPROB_PAIRED
# check_padding likely remove
# magic numbers in pad
# robustly test various pad lengths
# probably need code to reduce counts eg random select sequences?
# general robustness of pad
# specific mutate rescues
# visuals/coloring could add mutation site(s)
# add bpp options, other packages, linear
# pad, random but share length, all own pad, or directed eg AC*
# do I need to check the 5 and 3 interaction?
# maybe need a min edit distance for the uid?
###############################################################################

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

MAXPROB_NONINTERACT = 0.1
MINPROB_PAIRED = 0.5
MINAVGPROB_PAIRED = 0.8
MINPROB_UNPAIRED = 0.5
MINAVGPROB_UNPAIRED = 0.7

###############################################################################
# utils
###############################################################################


def get_reverse_complement(seq):
    reverse = seq[::-1]
    reverse_complement = reverse.replace("T", "X").replace("A", "T").replace(
        "X", "A").replace("G", "X").replace("C", "G").replace("X", "C")
    return reverse_complement


def combine_fastas(fastas, out_fasta):
    all_seqs = []
    for fasta in fastas:
        all_seqs.extend(list(SeqIO.parse(fasta, "fasta")))
    SeqIO.write(all_seqs, out_fasta, "fasta")


def format_fasta_for_submission(fasta, out_file, file_format='twist'):

    if file_format == 'twist':
        if out_file[-4:] != '.csv':
            print('twist requires .csv files, change out_file')
        names = []
        seqs = []
        for seq in SeqIO.parse(fasta, "fasta"):
            names.append(seq.name)
            seqs.append(str(seq.seq).upper().replace("U", "T"))
        df = pd.DataFrame(np.array([names, seqs]).T,
                          columns=["name", "sequence"])
        df.to_csv(out_file, index=False)

    elif file_format == 'custom_array':
        if out_file[-4:] != '.txt':
            print('custom_array requires .txt files, change out_file')
        with open(out_file, 'w') as f:
            for seq in SeqIO.parse(fasta, "fasta"):
                f.write(str(seq.seq).upper().replace("U", "T"))
                f.write('\n')

    else:
        print('file_format not supported, available: custom_array twist')


###############################################################################
# get desired sequences
###############################################################################


def get_windows(fasta, window_length, window_slide, out_fasta=None,
                circularize=False):
    seqs = list(SeqIO.parse(fasta, "fasta"))
    windows = []
    for seq_rec in seqs:
        seq = seq_rec.seq.upper().replace("U", "T")
        for i in range(0, len(seq), window_slide):
            if i+window_length > len(seq):
                if not circularize:
                    a, b = len(seq)-window_length, len(seq)
                    new_seq = seq[a:b]
                    namenum = f'{a}-{b-1}'
                    break
                else:
                    a, b, c = i, len(seq), window_length-len(seq)+i
                    new_seq = seq[a:b]+seq[:c]
                    namenum = f'{a}-{c-1}'
            else:
                a, b = i, i+window_length
                new_seq = seq[a:b]
                namenum = f'{a}-{b-1}'
            # save with name inclusive!
            new_rec = SeqIO.SeqRecord(Seq.Seq(new_seq),
                                      f'{seq_rec.name}_{namenum}', '', '')
            windows.append(new_rec)
    if out_fasta is not None:
        SeqIO.write(windows, out_fasta, "fasta")
    return windows


def get_all_single_mutants(fasta, out_fasta, bases=BASES):
    all_WT = SeqIO.parse(fasta, "fasta")
    all_single_mutants = []
    for record in all_WT:
        # convert to standard RNA
        seq = record.seq.upper().replace("U", "T")
        # at each position, get single mutants
        for i in range(len(seq)):
            for mut in bases:
                if mut != seq[i]:
                    name = f' {record.id}_{i}{seq[i]}>{mut}'
                    new_seq = seq[:i]+mut+seq[i+1:]
                    new_mut = SeqIO.SeqRecord(Seq.Seq(new_seq), name, '', '')
                    all_single_mutants.append(new_mut)

    SeqIO.write(all_single_mutants, out_fasta, "fasta")


def get_all_double_mutants(fasta, out_fasta, regionAs, regionBs, bases=BASES):
    # 9 per pair * lenregionA * lenregionB

    all_WT = list(SeqIO.parse(fasta, "fasta"))
    all_double_mutants = []  # unlike single mutant, WT not saved
    if len(all_WT) != len(regionAs) or len(all_WT) != len(regionBs):
        print(f'WARNING: you must list regions (as a list of lists) for all sequences in the fasta, number sequences in fasta {len(all_WT)} and number of regions specified {len(regionAs)} {len(regionBs)}')
    for record, regionA, regionB in zip(all_WT, regionAs, regionBs):
        # convert to standard RNA
        seq = record.seq.upper().replace("U", "T")
        # at each position pair, get double mutants
        for i in regionA:
            for mutA in bases:
                if mutA != seq[i]:
                    for j in regionB:
                        for mutB in bases:
                            if mutB != seq[j]:
                                name = f' {record.id}_{i}{seq[i]}>{mutA}_{j}{seq[j]}>{mutB}'
                                if i < j:
                                    new_seq = seq[:i]+mutA + \
                                        seq[i+1:j]+mutB+seq[j+1:]
                                else:
                                    new_seq = seq[:j]+mutB + \
                                        seq[j+1:i]+mutA+seq[i+1:]
                                new_mut = SeqIO.SeqRecord(
                                    Seq.Seq(new_seq), name, '', '')
                                all_double_mutants.append(new_mut)
    SeqIO.write(all_double_mutants, out_fasta, "fasta")


###############################################################################
# add library parts
###############################################################################

def write_all_barcodes(out_fasta=None, num_bp=8, num5hang=0, num3hang=0,
                       loop=TETRALOOP, bases=BASES):
    # probably should add ability to randomly generate but this
    # is fast enough for these small barcode
    all_barcodes = []
    for x in product(bases, repeat=num_bp+num5hang+num3hang):
        uid = ''.join(x)
        
        if num5hang+num3hang!=0:
            hang5,hang3 = '',''
            hang5 = uid[:num5hang]
            if num3hang==0:
                uid = uid[num5hang:]
            else:
                hang3 = uid[-num3hang:]
                uid = uid[num5hang:-num3hang]

            seq = hang5+uid+loop+get_reverse_complement(uid)+hang3
            seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), f' {uid}_{hang5}hang{hang3}', '', '')
        else:
            seq = uid+loop+get_reverse_complement(uid)
            seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), f' {uid}', '', '')
        all_barcodes.append(seq_rec)
    if out_fasta is not None:
        SeqIO.write(all_barcodes, out_fasta, "fasta")
    return all_barcodes


def _get_all_rand_seq(length,bases=BASES):
    all_seq = []
    for x in product(bases, repeat=length):
        seq = ''.join(x)
        seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), f' {seq}', '', '')
        all_seq.append(seq_rec)
    return all_seq


def _get_5_3_split(length):
    # 4+6+6+4 20
    if length < 20:
        pad_5_len, pad_3_len = 0, length
    elif length < 32:
        pad_5_len, pad_3_len = length-20, 20
    elif length % 2 == 0:
        pad_5_len, pad_3_len = length//2, length//2
    else:
        pad_5_len, pad_3_len = length//2, 1+(length//2)
    return pad_5_len, pad_3_len


def _get_stem_pads(pad_length,side="5'",loop=TETRALOOP,min_hang=3,bases=BASES):
    
    if side == "5'":
        # (((....)))....
        # for things that would have <4bp just get random seq
        if pad_length < 15:
            barcodes = _get_all_rand_seq(pad_length,bases)
            return barcodes,0,0,pad_length
        else:
            num_bp = (pad_length-len(loop)-min_hang)//2
            num_hang = pad_length-len(loop)-(2*num_bp)
            barcodes = write_all_barcodes(num_bp=num_bp, num3hang=num_hang,
                       loop=loop, bases=bases)
            return barcodes,num_bp,num_hang,len(loop)
    elif side == "3'":
        # ....(((....)))
        if pad_length < 15:
            barcodes = _get_all_rand_seq(pad_length,bases)
            return barcodes,0,0,pad_length
        else:
            num_bp = (pad_length-len(loop)-min_hang)//2
            num_hang = pad_length-len(loop)-(2*num_bp)
            barcodes = write_all_barcodes(num_bp=num_bp, num5hang=num_hang,
                       loop=loop, bases=bases)
            return barcodes,num_bp,num_hang,len(loop)

    else:
        print("ERROR side must be 5' or 3'")



def add_pad(fasta, out_fasta, bases=BASES, padding_type='SL_same_per_length',
            epsilon_punpaired=MINPROB_UNPAIRED,
            epsilon_interaction=MAXPROB_NONINTERACT,
            epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,loop=TETRALOOP,min_hang=3):

    # rand_same_all SL_same_per_length
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
        number_to_select = min(len(seq_group), max(10, len(seq_group)*0.01))
        selected_sec.append(sample(seq_group, k=number_to_select))

    # get length of pad needed
    desired_len = max(seq_by_length.keys())
    
    pads_by_len = {}
    if padding_type == 'SL_same_per_length':
        for seqs in selected_sec:
            length_to_add = desired_len-len(str(seqs[0].seq))
            if length_to_add==0:
                pads_by_len[len(str(seqs[0].seq))] = ['','']
            else:
                
                pad_5n, pad_3n = _get_5_3_split(length_to_add)

                pad_5s,num_bp5,num_hang5,loop_len5 = _get_stem_pads(pad_5n,"5'",loop,min_hang,bases)
                pad_3s,num_bp3,num_hang3,loop_len3 = _get_stem_pads(pad_3n,"3'",loop,min_hang,bases)

                shuffle(pad_5s)
                shuffle(pad_3s)

                # numbp_loop_numbp_hang_seq_hang_numbp_loop_numbp
                part_lengths = [num_bp5,loop_len5,num_bp5,num_hang5,len(str(seqs[0].seq)),num_hang3,num_bp3,loop_len3,num_bp3]

                region_unpaired = list(range(sum(part_lengths[:1]),sum(part_lengths[:2])))
                region_unpaired.extend(list(range(sum(part_lengths[:3]),sum(part_lengths[:4]))))
                region_unpaired.extend(list(range(sum(part_lengths[:5]),sum(part_lengths[:6]))))
                region_unpaired.extend(list(range(sum(part_lengths[:7]),sum(part_lengths[:8]))))

                region_paired_A = list(range(sum(part_lengths[:0]),sum(part_lengths[:1])))
                region_paired_B = list(range(sum(part_lengths[:2]),sum(part_lengths[:3])))[::-1]
                region_paired_A.extend(list(range(sum(part_lengths[:6]),sum(part_lengths[:7]))))
                region_paired_B.extend(list(range(sum(part_lengths[:8]),sum(part_lengths[:9])))[::-1])

                regionA = list(range(sum(part_lengths[:0]),sum(part_lengths[:4])))
                regionA.extend(list(range(sum(part_lengths[:5]),sum(part_lengths[:9]))))

                regionB = list(range(sum(part_lengths[:4]),sum(part_lengths[:5])))

                # loop through to find one that works
                good_pad = False
                current_pad = 0
                while not good_pad: 
                    
                    for i,seq in enumerate(seqs):
                        if i%2==0 and i!=0:
                            print(i)
                        full_seq = pad_5s[current_pad].seq + str(seq.seq).upper().replace('U','T') + pad_3s[current_pad].seq
                        good_pad, p_unpaired = check_struct_bpp(full_seq,region_unpaired,region_paired_A,region_paired_B,regionA,regionB)
                        if not good_pad:
                            break
                    current_pad += 1
                    if current_pad == len(pad_5s):
                        print('shuffling pads and trying again')
                        shuffle(pad_5s)
                        shuffle(pad_3s)
                        current_pad = 0 

                pads_by_len[len(str(seqs[0].seq))] = [pad_5s[current_pad].seq,pad_3s[current_pad].seq]

        padded_seqs = []
        for seq_length,seqs in seq_by_length.items():
            pad5,pad3 = pads_by_len[seq_length]
            for seq in seqs:
                full_seq = pad5+seq.seq.upper().replace("U", "T")+pad3
                padded_seqs.append(SeqIO.SeqRecord(
                    Seq.Seq(full_seq), f'{seq.name}_{len(pad5)}pad{len(pad3)}', '', ''))

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
                        seq.seq.upper().replace("U", "T") + pad3[:pad_3n]
                    regionA = list(range(pad_5n))
                    regionA += list(range(pad_5n + len(seq.seq),
                                          pad_5n+len(seq.seq)+pad_3n))
                    regionB = list(range(pad_5n, pad_5n+len(seq.seq)))

                    this_good, p_unpaired = check_struct_bpp(full_seq, region_unpaired=regionA,
                                            regionA=regionA,regionB=regionB,
                                            epsilon_punpaired=epsilon_punpaired,
                                                          epsilon_interaction=epsilon_interaction, epsilon_avg_punpaired=epsilon_avg_punpaired)
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
                full_seq = pad5[-pad_5n:] + \
                    seq_rec.seq.upper().replace("U", "T") + pad3[:pad_3n]
                padded_seqs.append(SeqIO.SeqRecord(
                    Seq.Seq(full_seq), f'{seq_rec.name}_{pad_5n}pad{pad_3n}', '', ''))
            else:
                padded_seqs.append(seq_rec)

    if out_fasta is not None:
        SeqIO.write(padded_seqs, out_fasta, "fasta")
    return padded_seqs


def add_fixed_seq_and_barcode(fasta, out_fasta=None, seq5=SEQ5, seq3=SEQ3,
                              num_bp=8, num5hang=2, loop=TETRALOOP,
                              epsilon=MAXPROB_NONINTERACT,
                              epsilon_punpaired=MINPROB_UNPAIRED,
                              epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                              epsilon_paired=MINPROB_PAIRED,
                              epsilon_avg_paired=MINAVGPROB_PAIRED,
                              save_image_folder=None):
    # get and randomly shuffle all potential barcoces
    all_uids = write_all_barcodes(num_bp=num_bp, loop=loop,num5hang=num5hang)
    shuffle(all_uids)
    current_uid = 0
    all_full_seqs = []
    p_unpaireds = {}
    seqs = SeqIO.parse(fasta, "fasta")
    if save_image_folder is not None:
        pad_lines = []
    print("Adding 5', barcode, 3'")
    # TODO do this in reasonable chunks for plotting
    for seq_rec in tqdm(list(seqs)):
        seq = str(seq_rec.seq).upper().replace("U", "T")
        name = seq_rec.name+'_w53barcode'

        # 5_seq_num5hang_numbp_loop_numbp
        regionA = list(range(len(seq5), len(seq5)+len(seq)))
        regionB = list(range(len(seq5)+len(seq), len(seq5) +
                             len(seq)+len(loop)+num5hang+(2*num_bp)))
        region_unpaired = list(range(len(seq5)+len(seq), len(seq5) +
                             len(seq)+num5hang))
        region_unpaired.extend(list(range(len(seq5)+len(seq)+num5hang+num_bp, len(seq5) +
                             len(seq)+num5hang+len(loop))))
        region_paired_A = list(range(len(seq5)+len(seq)+num5hang, len(seq5) +
                             len(seq)+num5hang+num_bp))
        region_paired_B = list(range(len(seq5)+len(seq)+num5hang+num_bp+len(loop), len(seq5) +
                             len(seq)+num5hang+(2*num_bp)+len(loop)))[::-1]

        uid_good = False
        while not uid_good:
            uid = all_uids[current_uid].seq
            full_seq = f'{seq5}{seq}{uid}{seq3}'
            current_uid += 1
            if save_image_folder is not None:

                lines = [len(seq5), len(seq5)+len(seq),
                         len(seq5)+len(seq)+len(loop)+num5hang+(2*num_bp)]
                if 'pad' in name:
                    name.split('pad')
                    pad5_length = int(name.split('pad')[0].split("_")[-1])
                    pad3_length = int(name.split('pad')[1].split("_")[0])
                    new_lines = [len(seq5)+pad5_length,
                                 len(seq5)+len(seq)-pad3_length]
                    lines.extend(new_lines)
                else:
                    new_lines = []
                uid_good, p_unpaired = check_struct_bpp(full_seq,region_unpaired,region_paired_A,region_paired_B,regionA,regionB, epsilon=epsilon,
                                                     epsilon_punpaired=epsilon_punpaired, epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                     epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired,
                                                     lines=lines,
                                                     save_image=f'{save_image_folder}/{name}.png')
            else:
                uid_good, p_unpaired = check_struct_bpp(full_seq,region_unpaired,region_paired_A,region_paired_B,regionA,regionB,
                                                     epsilon=epsilon, epsilon_punpaired=epsilon_punpaired, epsilon_avg_punpaired=epsilon_avg_punpaired,
                                                     epsilon_paired=epsilon_paired, epsilon_avg_paired=epsilon_avg_paired,)
        pad_lines.append(new_lines)
        all_full_seqs.append(SeqIO.SeqRecord(Seq.Seq(full_seq), name, '', ''))
        p_unpaireds[name] = p_unpaired
    if out_fasta is not None:
        SeqIO.write(all_full_seqs, out_fasta, "fasta")
    if save_image_folder is not None:
        plot_punpaired(p_unpaireds, f'{seq5}{" "*len(seq)}{" "*len(uid)}{seq3}',
                       [len(seq5), len(seq5)+len(seq), len(seq5) +
                        len(seq)+len(loop)+num5hang+(2*num_bp)], pad_lines,
                       f'{save_image_folder}/all_p_unpaired.png')
    return all_full_seqs


###############################################################################
# check structures
###############################################################################

def check_struct_bpp(seq,region_unpaired=None,region_paired_A=None,region_paired_B=None,
                  regionA=None,regionB=None,
                  epsilon=MAXPROB_NONINTERACT,
                  epsilon_punpaired=MINPROB_UNPAIRED, epsilon_avg_punpaired=MINAVGPROB_UNPAIRED,
                  epsilon_paired=MINPROB_PAIRED, epsilon_avg_paired=MINAVGPROB_PAIRED,
                  lines=None, save_image=None):

    bpp = bpps(seq.upper().replace("T", "U"),
               package='eternafold')
    p_unpaired = 1-bpp.sum(axis=0)

    bad_struct = False
    if region_unpaired is not None:
        punpaired_check = p_unpaired[region_unpaired]
        bad_struct = bad_struct or (punpaired_check.min() < epsilon_punpaired)
        bad_struct = bad_struct or (punpaired_check.mean() < epsilon_avg_punpaired)

    if regionA is not None:
        interaction_check = bpp[regionA][:, regionB]
        bad_struct = bad_struct or (interaction_check.max() > epsilon) 

    if region_paired_A is not None:
        paired_check = bpp[region_paired_A,region_paired_B]
        bad_struct = bad_struct or (paired_check.mean() < epsilon_avg_paired)
        bad_struct = bad_struct or (paired_check.min() < epsilon_paired)

    if bad_struct:
        return False, p_unpaired

    else:
        if save_image is not None:
            plot_bpp(bpp, seq, lines, save_image)
        return True, p_unpaired


###############################################################################
# visualize structures
###############################################################################

def plot_bpp(bpp, seq, lines, save_image):
    plt.figure(figsize=(len(seq)*0.12, len(seq)*0.12))
    plt.imshow(bpp, origin='upper', cmap='gist_heat_r')
    xlabels = []
    ylabels = []
    for i, s in enumerate(seq):
        if i % 10 == 0:
            ylabels.append(f'{i} {s}')
            xlabels.append(f'{s}\n{i}')
        else:
            xlabels.append(s)
            ylabels.append(s)
    plt.xticks(range(len(seq)), xlabels, size=8)
    plt.yticks(range(len(seq)), ylabels, size=8)
    for line in lines:
        plt.hlines(line, 0, len(seq), color="grey")
        plt.vlines(line, 0, len(seq), color="grey")
    plt.xlim(0, len(seq))
    plt.ylim(len(seq), 0)
    plt.savefig(save_image, bbox_inches='tight', dpi=300)
    plt.close()


def plot_punpaired(p_unpaired, seq, lines, pad_lines, save_image):
    plt.figure(figsize=(len(seq)*0.3, len(p_unpaired)*0.3))
    df = pd.DataFrame(p_unpaired).T
    plt.imshow(df, cmap='gist_heat_r')
    plt.yticks(range(len(df)), df.index, size=8)
    plt.xticks(range(len(seq)), seq, size=8)
    ax = plt.gca()
    y1, y2 = ax.get_ylim()
    for line in lines:
        plt.vlines(line-0.5, y1+1, y2-1, color="cyan",linewidth=8)
    for i, line in enumerate(pad_lines):
        if line != []:
            plt.vlines(line[0]-0.5, i+0.5, i-0.5, color='lime',linewidth=8)
            plt.vlines(line[1]-0.5, i+0.5, i-0.5, color='lime',linewidth=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)
    plt.savefig(save_image, bbox_inches='tight', dpi=300)
    plt.close()


###############################################################################
# example run
###############################################################################

# TODO write_all to get_all
# paralleize idea, run pad search on all single sequences
# then when all done, come together and find one of those that works
# repeat if necessary?
# this may also just be fast enough

#write_all_barcodes('../examples/all_8bp_barcodes.fasta')
#write_all_barcodes('../examples/all_8bp_barcodes_hang.fasta',num_bp=4,num5hang=1, num3hang=2)

'''
get_windows('../examples/example_WT.fasta', 30, 4,
            out_fasta='../examples/example_windows.fasta', circularize=False)
get_all_single_mutants('../examples/example_WT.fasta',
                       '../examples/example_single_mut.fasta')
get_all_double_mutants('../examples/example_WT.fasta',
                       '../examples/example_double_mut.fasta', [[3, 4], [9, 103], [50, 51, 60, 62]], [[10, 11, 12], [30, 31], [120]])
combine_fastas(['../examples/example_WT.fasta', '../examples/example_single_mut.fasta',
                '../examples/example_double_mut.fasta'], '../examples/example_WT_single_double_mut.fasta')


add_pad('../examples/example_WT_single_double_mut.fasta',
        '../examples/example_padded.fasta')

add_fixed_seq_and_barcode('../examples/example_padded.fasta', '../examples/example_finished.fasta',
                          save_image_folder='../examples/bpps')
format_fasta_for_submission('../examples/example_finished.fasta',
                            '../examples/example_finished.csv', file_format='twist')
format_fasta_for_submission('../examples/example_finished.fasta',
                            '../examples/example_finished.txt', file_format='custom_array')
'''

add_fixed_seq_and_barcode('../examples/example_padded_small.fasta', '../examples/example_finished_small.fasta',
                          save_image_folder='../examples/bpps')
format_fasta_for_submission('../examples/example_finished_small.fasta',
                            '../examples/example_finished_small.csv', file_format='twist')
format_fasta_for_submission('../examples/example_finished_small.fasta',
                            '../examples/example_finished_small.txt', file_format='custom_array')
'''
add_pad('../examples/example_WT_single_double_mut.fasta',
        '../examples/example_padded.fasta')
'''
