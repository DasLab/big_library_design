import argparse
from utils.mutational_utils import get_windows,add_fixed_seq_and_barcode
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

if args.m2seq:
    get_windows(args.input_fasta,args.length,args.step,f'{args.output_prefix}_windowed.fasta',circularize=args.circularize)
    add_fixed_seq_and_barcode(f'{args.output_prefix}_windowed.fasta',f'{args.output_prefix}_library.fasta')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.csv',file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta',f'{args.output_prefix}_library.txt',file_format='custom_array')

else:
    print('ERROR program not yet implemented')

