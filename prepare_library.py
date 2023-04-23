import argparse
from utils.mutational_utils import *

# README with example calls - windows, m2seq, and also a single script for RD
# TODO print in process what is being done

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta', type=str, required=True,
                    help='Sequences to make library out of, file in fasta format.')
parser.add_argument('-o', '--output_prefix', type=str,
                    required=True, help='Prefix of the output files.')

programs_ = parser.add_argument_group('programs')
programs = programs_.add_mutually_exclusive_group()
programs.add_argument('--window', action='store_true',
                      help='The window program will take each sequence, create sliding windows and then prepare the library sequences.')
programs.add_argument('--m2seq', action='store_true',
                      help='The m2seq program will take each sequence, get all single mutants, create a noninteracting pad to ensure each sequence is the same length if needed, and then prepare the library sequences.')

visuals = parser.add_argument_group('Visuals')
visuals.add_argument('--save_bpp_fig', action='store_true',
                     help='Whether to save images of the base-pair-probability matrices.')
visuals.add_argument('--save_image_folder', default=None,
                     help="Folder to save images (proabbility unpaired and if specified base-pair-probability matrices), default don't save.")

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

m2seq = parser.add_argument_group('m2seq')
m2seq.add_argument('--pad_type', type=str, help='')

args = parser.parse_args()

if args.window:
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
                              epsilon=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

elif args.m2seq:
    get_all_single_mutants(args.input_fasta, f'{args.output_prefix}_single_mut.fasta')
    combine_fastas([args.input_fasta, f'{args.output_prefix}_single_mut.fasta'], f'{args.output_prefix}_WT_single_mut.fasta')
    add_pad(f'{args.output_prefix}_WT_single_mut.fasta', f'{args.output_prefix}_WT_single_mut_pad.fasta')
    add_fixed_seq_and_barcode(f'{args.output_prefix}_WT_single_mut_pad.fasta', f'{args.output_prefix}_library.fasta')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

else:
    print('ERROR program not yet implemented')


'''
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

# add_fixed_seq_and_barcode('../examples/example_padded_small.fasta', '../examples/example_finished_small.fasta',
#                          save_image_folder='../examples/bpps')
# format_fasta_for_submission('../examples/example_finished_small.fasta',
#                            '../examples/example_finished_small.csv', file_format='twist')
# format_fasta_for_submission('../examples/example_finished_small.fasta',
#                            '../examples/example_finished_small.txt', file_format='custom_array')
'''
add_pad('../examples/example_WT_single_double_mut.fasta',
        '../examples/example_padded.fasta')
'''
