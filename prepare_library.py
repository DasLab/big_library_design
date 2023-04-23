import argparse
from utils.mutational_utils import *

parser = argparse.ArgumentParser()

programs = parser.add_mutually_exclusive_group()

programs.add_argument('--window',action='store_true',help='')
programs.add_argument('--m2seq',action='store_true',help='')

parser.add_argument('-i','--input_fasta',type=str)
parser.add_argument('-o','--output_prefix',type=str)

window = parser.add_argument_group('window')
window.add_argument('--length',type=int,default=100,help='')
window.add_argument('--step',type=int,default=10,help='')
window.add_argument('--circularize',action='store_true',help='')

m2seq = parser.add_argument_group('m2seq')
m2seq.add_argument('--pad_type',type=str,help='')

args = parser.parse_args()

if args.window:
    get_windows(args.input_fasta,args.length,args.step,f'{args.output_prefix}_windowed.fasta',circularize=args.circularize)
    add_fixed_seq_and_barcode(f'{args.output_prefix}_windowed.fasta',f'{args.output_prefix}_library.fasta')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.csv',file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.txt',file_format='custom_array')

elif args.m2seq:
    get_all_single_mutants(args.input_fasta,f'{args.output_prefix}_single_mut.fasta')
    combine_fastas([args.input_fasta,f'{args.output_prefix}_single_mut.fasta'], f'{args.output_prefix}_WT_single_mut.fasta')
    add_pad(f'{args.output_prefix}_WT_single_mut.fasta',f'{args.output_prefix}_WT_single_mut_pad.fasta')
    add_fixed_seq_and_barcode(f'{args.output_prefix}_WT_single_mut_pad.fasta',f'{args.output_prefix}_library.fasta')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.csv',file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.txt',file_format='custom_array')

else:
    print('ERROR program not yet implemented')


###############################################################################
# example run
###############################################################################



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

#add_fixed_seq_and_barcode('../examples/example_padded_small.fasta', '../examples/example_finished_small.fasta',
#                          save_image_folder='../examples/bpps')
#format_fasta_for_submission('../examples/example_finished_small.fasta',
#                            '../examples/example_finished_small.csv', file_format='twist')
#format_fasta_for_submission('../examples/example_finished_small.fasta',
#                            '../examples/example_finished_small.txt', file_format='custom_array')
'''
add_pad('../examples/example_WT_single_double_mut.fasta',
        '../examples/example_padded.fasta')
'''

