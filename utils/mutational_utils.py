# function for making mutations from WT sequence

from Bio import SeqIO, Seq
from itertools import product
from random import shuffle, sample, choices
from arnie.bpps import bpps
import matplotlib.pyplot as plt
import pandas as pd

BASES = ['A', 'C', 'G', 'T']
SEQ5 = 'GGGAACGACTCGAGTAGAGTCGAAAA'
SEQ3 = 'AAAAGAAACAACAACAACAAC'  # with the single A spacer added already
TETRALOOP = 'TTCG'
MAXPROB = 0.05
MINPROB = 0.25
MINAVGPROB = 0.7

######### TODO ##########
# specific mutate rescues
# padding
# UID helices
#  26 nt.   5' fixed sequence
# 240 nt.   Region of interest
#  20 nt.   8 bp barcode with UUCG tetraloop (65,536 possible barcodes)
#   1 nt.   Single A spacer
#  20 nt.   3' fixed sequence

#### current procedure ####
# get list of desired sequences: WTs, all single mutants, double mutants of interaction of interest
# if needed generate random pad to fill length (half on 5', half on 3'), select random sampling of desired sequences, iterate finding new random pad until, for every sampled sequence, bpp between sequence and pad is < eps everywhere
# get all possible 8bp barcodes
# for each sequence add 5', selected random barcode (no replacement), and 3', check bpp between sequence and barcode < eps everywhere


def write_all_barcodes(out_fasta=None, num_bp=8, loop=TETRALOOP, bases=BASES):
    # TODO probably should add ability to randomly generate but this
    # is fast enough for now
    all_barcodes = []
    for x in product(bases, repeat=num_bp):
        uid = ''.join(x)
        seq = uid+loop+get_reverse_complement(uid)
        seq_rec = SeqIO.SeqRecord(Seq.Seq(seq), uid, uid, uid)
        all_barcodes.append(seq_rec)
    if out_fasta is not None:
        SeqIO.write(all_barcodes, out_fasta, "fasta")
    return all_barcodes


def get_reverse_complement(seq):
    reverse = seq[::-1]
    reverse_complement = reverse.replace("T", "X").replace("A", "T").replace(
        "X", "A").replace("G", "X").replace("C", "G").replace("X", "C")
    return reverse_complement


def add_pad(fasta, out_fasta, bases=BASES, epsilon_punpaired=MINPROB, epsilon_interaction=MAXPROB, epsilon_avg_punpaired=MINAVGPROB):
    # TODO to add a pad to the sequence that will not interaction
    seqs = list(SeqIO.parse(fasta, "fasta"))
    seq_by_length = {}
    for seq_rec in seqs:
        len_seq = len(seq_rec.seq)
        if len_seq in seq_by_length:
            seq_by_length[len_seq].append(seq_rec)
        else:
            seq_by_length[len_seq] = [seq_rec]
    # for each len randomly select 1% of sequences or 100 number of sequences if < 100
    selected_sec = []
    for len_seq, seq_group in seq_by_length.items():
        number_to_select = min(len(seq_group), max(100, len(seq_group)*0.01))
        selected_sec.extend(sample(seq_group, k=number_to_select))

    # get length of pad needed
    desired_len = max(seq_by_length.keys())
    max_length_to_add = desired_len-min(seq_by_length.keys())
    if max_length_to_add % 2 == 0:
        pad_5_len, pad_3_len = max_length_to_add//2, max_length_to_add//2
    else:
        pad_5_len, pad_3_len = max_length_to_add//2, 1+(max_length_to_add//2)

    good_pad = False
    while not good_pad:
        pad5 = ''.join(choices(bases, k=pad_5_len))
        pad3 = ''.join(choices(bases, k=pad_3_len))
        any_bad = False
        for seq in selected_sec:
            length_to_add = desired_len-len(seq.seq)
            if length_to_add != 0:
                if length_to_add % 2 == 0:
                    pad_5n, pad_3n = length_to_add//2, length_to_add//2
                else:
                    pad_5n, pad_3n = length_to_add//2, 1+(length_to_add//2)
                full_seq = pad5[-pad_5n:] + seq.seq + pad3[:pad_3n]
                regionA = list(range(pad_5n))+list(range(pad_5n +
                                                         len(seq.seq), pad_5n+len(seq.seq)+pad_3n))
                regionB = list(range(pad_5n, pad_5n+len(seq.seq)))
                this_good, p_unpaired = check_punpaired_high(full_seq, regionA, regionB, epsilon_punpaired=epsilon_punpaired,
                                                             epsilon_interaction=epsilon_interaction, epsilon_avg_punpaired=epsilon_avg_punpaired)
                if not this_good:
                    any_bad = True
                    continue
        if not any_bad:
            good_pad = True
    padded_seqs = []
    for seq_rec in seqs:
        length_to_add = desired_len-len(seq_rec.seq)
        if length_to_add != 0:
            if length_to_add % 2 == 0:
                pad_5n, pad_3n = length_to_add//2, length_to_add//2
            else:
                pad_5n, pad_3n = length_to_add//2, 1+(length_to_add//2)
            full_seq = pad5[-pad_5n:] + seq.seq + pad3[:pad_3n]
            padded_seqs.append(SeqIO.SeqRecord(
                Seq.Seq(full_seq), f'{seq_rec.name}_{pad_5n}pad{pad_3n}', '', ''))
        else:
            padded_seqs.append(seq_rec)
    if out_fasta is not None:
        SeqIO.write(padded_seqs, out_fasta, "fasta")
    return padded_seqs


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


def add_fixed_seq_and_barcode(fasta, out_fasta=None, seq5=SEQ5, seq3=SEQ3, num_bp=8, loop=TETRALOOP, epsilon=MAXPROB, save_image_folder=None):
    # TODO maybe need a min edit distance for the uid?
    # get and randomly shuffle all potential barcoces
    all_uids = write_all_barcodes(num_bp=num_bp, loop=loop)
    shuffle(all_uids)
    current_uid = 0
    all_full_seqs = []
    p_unpaireds = {}
    seqs = SeqIO.parse(fasta, "fasta")
    for seq_rec in seqs:
        seq = str(seq_rec.seq)
        name = seq_rec.name+'_w53barcode'

        # TODO do I need to check the 5 and 3 interction?
        # TODO maybe also a p-unpaired check

        regionA = list(range(len(seq5), len(seq5)+len(seq)))
        regionB = list(range(len(seq5)+len(seq), len(seq5) +
                             len(seq)+len(loop)+(2*num_bp)))
        uid_good = False
        while not uid_good:
            uid = all_uids[current_uid].seq
            full_seq = f'{seq5}{seq}{uid}{seq3}'
            current_uid += 1
            if save_image_folder is not None:
                # TODO for lines pad may also go in there
                uid_good, p_unpaired = check_bpp_no_signal(full_seq, regionA, regionB, epsilon=epsilon,
                                                           lines=[len(seq5), len(seq5)+len(seq), len(seq5)+len(seq)+len(loop)+(2*num_bp)], save_image=f'{save_image_folder}/{name}.png')
            else:
                uid_good, p_unpaired = check_bpp_no_signal(
                    full_seq, regionA, regionB, epsilon=epsilon)
        all_full_seqs.append(SeqIO.SeqRecord(Seq.Seq(seq), name, '', ''))
        p_unpaireds[name] = p_unpaired
    if out_fasta is not None:
        SeqIO.write(all_full_seqs, out_fasta, "fasta")
    if save_image_folder is not None:
        plot_punpaired(p_unpaireds, f'{seq5}{" "*len(seq)}{" "*len(uid)}{seq3}', [len(seq5), len(seq5)+len(seq), len(seq5)+len(seq)+len(loop)+(2*num_bp)], f'{save_image_folder}/all_p_unpaired.png')
    return all_full_seqs


def check_bpp_no_signal(seq, regionA, regionB, epsilon=MAXPROB, lines=None, save_image=None):
    # TODO save heatmap of these somewhere, with box of region? -- only good ones
    # TODO likely give optiono of package and linear or not
    bpp = bpps(seq, package='eternafold')  # ,linear=True)
    if save_image is not None:
        plot_bpp(bpp, seq, lines, save_image)
    region_to_check = bpp[regionA][:, regionB]
    p_unpaired = 1-bpp.sum(axis=0)
    if region_to_check.max() > epsilon:
        return False, p_unpaired
    else:
        return True, p_unpaired


def check_punpaired_high(seq, regionA, regionB, epsilon_punpaired=MINPROB, epsilon_interaction=MAXPROB, epsilon_avg_punpaired=MINAVGPROB):
    # TODO save heatmap of these somewhere, with box of region? -- only good ones
    # TODO likely give optiono of package and linear or not
    bpp = bpps(seq, package='eternafold')  # ,linear=True)
    p_unpaired = 1-bpp.sum(axis=0)
    p_unpaired_region_to_check = p_unpaired[regionA]
    region_to_check = bpp[regionA][:, regionB]
    print(p_unpaired_region_to_check.mean())
    if ((p_unpaired_region_to_check.min() < epsilon_punpaired) and
            (region_to_check.max() > epsilon_interaction) and (p_unpaired_region_to_check.mean() < epsilon_avg_punpaired)):
        return False, p_unpaired
    else:
        return True, p_unpaired


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
    # color sequence by 5 region etc, maybe even mut spot
    # or just put hline/vline at these boundaries
    plt.savefig(save_image)


def plot_punpaired(p_unpaired, seq, lines, save_image):
    plt.figure(figsize=(len(p_unpaired)*20, len(p_unpaired)*20))
    # plot desired region not design etc
    df = pd.DataFrame(p_unpaired).T
    plt.imshow(df, cmap='gist_heat_r')
    plt.yticks(range(len(df)), df.index, size=8)
    plt.xticks(range(len(seq)), seq, size=8)
    ax = plt.gca()
    y1,y2=ax.get_ylim()
    for line in lines:
        plt.vlines(line-0.5, y1+1, y2-1, color="blue")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='both', which='both',length=0)
    plt.savefig(save_image)


def format_seqs_for_submission(fasta):
    return None


def get_all_single_mutants(fasta, out_fasta, bases=BASES):
    all_WT = SeqIO.parse(fasta, "fasta")
    all_single_mutants = list(SeqIO.parse(fasta, "fasta"))
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


get_all_single_mutants('../examples/example_WT.fasta',
                       '../examples/example_single_mut.fasta')
write_all_barcodes('../examples/all_8bp_barcodes.fasta')
get_all_double_mutants('../examples/example_WT.fasta',
                       '../examples/example_double_mut.fasta', [[3, 4], [9, 103]], [[10, 11, 12], [30, 31]])
#add_pad('../examples/example_WT.fasta', '../examples/example_padded.fasta')
add_fixed_seq_and_barcode('../examples/example_padded.fasta',
                          save_image_folder='../examples/bpps')
