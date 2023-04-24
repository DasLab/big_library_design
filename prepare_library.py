import argparse
from utils.mutational_utils import *

###############################################################################
# Arguments
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta', type=str, required=True,
                    help='Sequences to make library out of, file in fasta format.')
parser.add_argument('-o', '--output_prefix', type=str,
                    required=True, help='Prefix of the output files.')

programs_ = parser.add_argument_group('programs')
programs = programs_.add_mutually_exclusive_group()
programs.add_argument('--just_library',action='store_true',
                    help='The just_libary program will take all sequences provided and then just prepare the libary using only those sequences.')
programs.add_argument('--window', action='store_true',
                      help='The window program will take each sequence, create sliding windows and then prepare the library sequences.')
programs.add_argument('--m2seq', action='store_true',
                      help='The m2seq program will take each sequence, get all single mutants, create a noninteracting pad to ensure each sequence is the same length if needed, and then prepare the library sequences.')
programs.add_argument('--m2seq_with_double', action='store_true',
                      help='The program runs the m2seq program with the addition of double mutants in user defined regions.')

visuals = parser.add_argument_group('Visuals')
visuals.add_argument('--save_bpp_fig', type=float, default=0,
                     help='Proportion to save images of the base-pair-probability matrices. 0 is none, 1 is all, in between is random seleciton of sequences.')
visuals.add_argument('--save_image_folder', default=None,
                     help="Folder to save images (proabbility unpaired and if specified base-pair-probability matrices), default don't save.")
visuals.add_argument('--max_seq_punpaired_plot', type=int, default=500, help='The maximum number of sequences to put in one probability unpaired plot, will split into seperate images if smaller than total number of sequences.')

struct = parser.add_argument_group('Strctural checks')
struct.add_argument('--Pmax_noninteract', type=float, default=0.05,
                    help='Maximum base-pair-probability for 2 regions to be considered noninteracting.')
struct.add_argument('--Pmin_paired', type=float, default=0.75,
                    help='Minimum base-pair-probability for 2 regions to be considered paired.')
struct.add_argument('--Pavg_paired', type=float, default=0.85,
                    help='Average base-pair-probability for 2 regions to be considered paired.')
struct.add_argument('--Pmin_unpaired', type=float, default=0.75,
                    help='Minimum probability unpaired for region to be unpaired.')
struct.add_argument('--Pavg_unpaired', type=float, default=0.85,
                    help='Average probability unpaired for region to be unpaired.')

library = parser.add_argument_group('Libarary parts')
library.add_argument('--seq5', type=str, default='GGGAACGACTCGAGTAGAGTCGAAAA',
                     help="Constant sequence to place at 5' of every sequence in library.")
library.add_argument('--seq3', type=str, default='AAAAGAAACAACAACAACAAC',
                     help="Constant sequence to place at 3' of every sequence in library.")
library.add_argument('--barcode_numbp', type=int, default=8,
                     help="The length (in bp) of the random barcode stem.")
library.add_argument('--barcode_num5polyA', type=int, default=4,
                     help="The length of polyA stretch placed before the barcode stem, if specified, also put before the random 5'hang..")
library.add_argument('--barcode_num5randomhang', type=int, default=0,
                     help="The legnth of addtional random (single-standed) sequence to place before the barcode stem.")
library.add_argument('--barcode_loop', type=str, default='TTCG',
                     help="Constant loop sequence of the barcode stemloop.")

window = parser.add_argument_group('window')
window.add_argument('--length', type=int, default=100,
                    help='Length of each window to create.')
window.add_argument('--step', type=int, default=10,
                    help='The step size of the sliding window.')
window.add_argument('--circularize', action='store_true',
                    help="Whether to circularize the sequence (at 3' end, don't stop but loop back to 5') or not.")

m2seq = parser.add_argument_group('padding for: m2seq or just_library when length not equal')
m2seq.add_argument('--pad_type', type=str, default='SL_same_per_length', help='If there are sequencees of multiple lengths, to obtain a libarary of equal length some sequences will be padded. This specifies the type of padding with SL_same_per_length (if pad is long enough create a stem-loop, same pad is used for each group of sequences with equal length) or rand_same_all (a single-stranded non-interacting pad is choosen, same for all sequences just using the length of pad required to pad each sequence to same length) as options.')
m2seq.add_argument('--pad_loop',type=str,default='TTCG',help='If padtype is a stem-loop the constant loop to use.')
m2seq.add_argument('--pad_min_hang',type=int,default=3,help='If padtype is a stem-loop the minimum (only +1 possible) to have a random, single-stranded hang between sequence of interest and pad.')
m2seq.add_argument('--pad_num_samples',type=int,default=30,help="Minimum number of sequences to check pad's effect on structure.")

double = parser.add_argument_group('double mutant')
double.add_argument('--doublemut',nargs='+',help='Two region to do double mutagenis of (one mutant in group A one in group B) format: 1-12,15.64-70,72-78 where this would mean one mutant in nucletoides 1to 12 (inclusive) or 15 and one mutant in region 64 to 70 or 72 to 78. 1-12.1-12 would mean all double mutants in 1to 12. If more than one sequence is in the input need to specify the same number of regions seperated by space eg for 2 sequences: 1-12,15.64-70,72-78 34-78.80-85 ')

args = parser.parse_args()


###############################################################################
# just library
###############################################################################

if args.just_library:
    # check pad needed
    if get_same_length(args.input_fasta):
        fasta = args.input_fasta
    else:
        fasta = f'{args.output_prefix}_pad.fasta'
        add_pad(args.input_fasta,
            fasta,
            padding_type=args.pad_type,
            epsilon_interaction=args.Pmax_noninteract,
            epsilon_punpaired=args.Pmin_unpaired,
            epsilon_avg_punpaired=args.Pavg_unpaired,
            epsilon_paired=args.Pmin_paired,
            epsilon_avg_paired=args.Pavg_paired,
            loop=args.pad_loop,
            min_hang=args.pad_min_hang,
            min_num_samples=args.pad_num_samples)
    add_fixed_seq_and_barcode(fasta,
                              f'{args.output_prefix}_library.fasta',
                              seq5=args.seq5,
                              seq3=args.seq3,
                              loop=args.barcode_loop,
                              num_bp=args.barcode_numbp,
                              num5hang=args.barcode_num5randomhang,
                              num5polyA=args.barcode_num5polyA,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

###############################################################################
# Window
###############################################################################

elif args.window:
    get_windows(args.input_fasta, args.length, args.step,
                f'{args.output_prefix}_windowed.fasta',
                circularize=args.circularize)
    add_fixed_seq_and_barcode(f'{args.output_prefix}_windowed.fasta',
                              f'{args.output_prefix}_library.fasta',
                              seq5=args.seq5,
                              seq3=args.seq3,
                              loop=args.barcode_loop,
                              num_bp=args.barcode_numbp,
                              num5hang=args.barcode_num5randomhang,
                              num5polyA=args.barcode_num5polyA,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

###############################################################################
# m2seq
###############################################################################

elif args.m2seq:
    get_all_single_mutants(args.input_fasta, f'{args.output_prefix}_single_mut.fasta')
    combine_fastas([args.input_fasta, f'{args.output_prefix}_single_mut.fasta'], f'{args.output_prefix}_WT_single_mut.fasta')
    add_pad(f'{args.output_prefix}_WT_single_mut.fasta',
            f'{args.output_prefix}_WT_single_mut_pad.fasta',
            padding_type=args.pad_type,
            epsilon_interaction=args.Pmax_noninteract,
            epsilon_punpaired=args.Pmin_unpaired,
            epsilon_avg_punpaired=args.Pavg_unpaired,
            epsilon_paired=args.Pmin_paired,
            epsilon_avg_paired=args.Pavg_paired,
            loop=args.pad_loop,
            min_hang=args.pad_min_hang,
            min_num_samples=args.pad_num_samples)
    add_fixed_seq_and_barcode(f'{args.output_prefix}_WT_single_mut_pad.fasta',
                              f'{args.output_prefix}_library.fasta',
                              seq5=args.seq5,
                              seq3=args.seq3,
                              loop=args.barcode_loop,
                              num_bp=args.barcode_numbp,
                              num5hang=args.barcode_num5randomhang,
                              num5polyA=args.barcode_num5polyA,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

###############################################################################
# m2seq_with_double
###############################################################################

elif args.m2seq_with_double:
    get_all_single_mutants(args.input_fasta, f'{args.output_prefix}_single_mut.fasta')
    regionAs, regionBs = get_regions_for_doublemut(args.doublemut)
    get_all_double_mutants(args.input_fasta, f'{args.output_prefix}_double_mut.fasta',regionAs, regionBs)
    combine_fastas([args.input_fasta, f'{args.output_prefix}_single_mut.fasta',f'{args.output_prefix}_double_mut.fasta'], 
        f'{args.output_prefix}_WT_single_mut.fasta')
    add_pad(f'{args.output_prefix}_WT_single_mut.fasta',
            f'{args.output_prefix}_WT_single_mut_pad.fasta',
            padding_type=args.pad_type,
            epsilon_interaction=args.Pmax_noninteract,
            epsilon_punpaired=args.Pmin_unpaired,
            epsilon_avg_punpaired=args.Pavg_unpaired,
            epsilon_paired=args.Pmin_paired,
            epsilon_avg_paired=args.Pavg_paired,
            loop=args.pad_loop,
            min_hang=args.pad_min_hang,
            min_num_samples=args.pad_num_samples)
    add_fixed_seq_and_barcode(f'{args.output_prefix}_WT_single_mut_pad.fasta',
                              f'{args.output_prefix}_library.fasta',
                              seq5=args.seq5,
                              seq3=args.seq3,
                              loop=args.barcode_loop,
                              num_bp=args.barcode_numbp,
                              num5hang=args.barcode_num5randomhang,
                              num5polyA=args.barcode_num5polyA,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')


else:
    print('ERROR program not yet implemented')
