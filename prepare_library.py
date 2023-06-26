import argparse
from utils.mutational_utils import *

###############################################################################
# Arguments
###############################################################################

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

required = parser.add_argument_group('Required')
required.add_argument('-i', '--input_fasta', type=str, required=True,
                    help='Sequences to make library out of, file in fasta format.')
required.add_argument('-o', '--output_prefix', type=str,
                    required=True, help='Prefix of the output files.')

programs_ = parser.add_argument_group('programs')
programs = programs_.add_mutually_exclusive_group()
programs.add_argument('--check_library', action='store_true',
                      help='The check_libary program will take one fasta and make sure it is ready for submission, reporting any problems.')
programs.add_argument('--just_library', action='store_true',
                      help='The just_libary program will take all sequences provided and then just prepare the library using only those sequences.')
programs.add_argument('--window', action='store_true',
                      help='The window program will take each sequence, create sliding windows and then prepare the library sequences.')
programs.add_argument('--m2seq', action='store_true',
                      help='The m2seq program will take each sequence, get all single mutants, create a non-interacting pad to ensure each sequence is the same length if needed, and then prepare the library sequences.')
programs.add_argument('--m2seq_with_double', action='store_true',
                      help='The program runs the m2seq program with the addition of double mutants in user defined regions.')

visuals = parser.add_argument_group('Visuals')
visuals.add_argument('--save_bpp_fig', type=float, default=0,
                     help='Proportion to save images of the base-pair-probability matrices. 0 is none, 1 is all, in between is random selection of sequences.')
visuals.add_argument('--save_image_folder', default=None,
                     help="Folder to save images (probability unpaired and if specified base-pair-probability matrices), default don't save.")
visuals.add_argument('--max_seq_punpaired_plot', type=int, default=500,
                     help='The maximum number of sequences to put in one probability unpaired plot, will split into separate images if smaller than total number of sequences.')

struct = parser.add_argument_group('Structural checks')
struct.add_argument('--Pmax_noninteract', type=float, default=0.05,
                    help='Maximum base-pair-probability for 2 regions to be considered non-interacting.')
struct.add_argument('--Pmin_paired', type=float, default=0.75,
                    help='Minimum base-pair-probability for 2 regions to be considered paired.')
struct.add_argument('--Pavg_paired', type=float, default=0.85,
                    help='Average base-pair-probability for 2 regions to be considered paired.')
struct.add_argument('--Pmin_unpaired', type=float, default=0.75,
                    help='Minimum probability unpaired for region to be unpaired.')
struct.add_argument('--Pavg_unpaired', type=float, default=0.85,
                    help='Average probability unpaired for region to be unpaired.')

library = parser.add_argument_group('Library parts')
library.add_argument('--seq5', type=str, default='GGGAACGACTCGAGTAGAGTCGAAAA',
                     help="Constant sequence to place at 5' of every sequence in library.")
library.add_argument('--seq3', type=str, default='AAAAGAAACAACAACAACAAC',
                     help="Constant sequence to place at 3' of every sequence in library.")
library.add_argument('--barcode_numbp', type=int, default=8,
                     help="The length (in bp) of the random barcode stem.")
library.add_argument('--barcode_num5polyA', type=int, default=4,
                     help="The length of polyA stretch placed before the barcode stem, if specified, also put before the random 5'hang..")
library.add_argument('--barcode_num5randomhang', type=int, default=0,
                     help="The length of additional random (single-standed) sequence to place before the barcode stem.")
library.add_argument('--barcode_loop', type=str, default='TTCG',
                     help="Constant loop sequence of the barcode stem-loop.")
library.add_argument('--num_barcodes_reduce_prob', type=float, default=100,
                     help='Number of barcodes to try before reducing probability thresholds by 10 percent.')
library.add_argument('--avoid_barcodes_files', nargs='+',
                     help='Fasta files which contains sequences with barcodes to remove.')
library.add_argument('--avoid_barcodes_start', type=int,
                     help='First nucleotide in barcode.')
library.add_argument('--avoid_barcodes_end', type=int,
                     help='Last nucleotide in barcode.')
library.add_argument('--min_edit', type=int, default=2,
                     help='Minimum edit distance for barcodes (2 or 3 should suffice.')

window = parser.add_argument_group('window')
window.add_argument('--length', type=int, default=100,
                    help='Length of each window to create.')
window.add_argument('--step', type=int, default=10,
                    help='The step size of the sliding window.')
window.add_argument('--circularize', action='store_true',
                    help="Whether to circularize the sequence (at 3' end, don't stop but loop back to 5') or not.")
window.add_argument('--reverse_complement', action='store_true',
                    help="Whether to use the negative sense strand instead.")
window.add_argument('--prop_windows_keep', type=float, default=1.0,
                    help='The proportion of windows to keep (from the start), others not added to library.')
window.add_argument('--prop_windows_keep_random', type=float, default=1.0,
                    help='The proportion of windows to keep (randomly selected), others not added to library.')

m2seq = parser.add_argument_group(
    'padding for: m2seq or just_library when length not equal')
m2seq.add_argument('--share_pad', type=str, default='same_length',
                   help='If there are sequences of multiple lengths, to obtain a library of equal length some sequences will be padded. This specifies which sequence share the same pad, all, none or same_length sequences.')
m2seq.add_argument('--pad_loop', type=str, default='TTCG',
                   help='If pad has stem the constant loop to use.')
m2seq.add_argument('--pad_polyAhang', type=int, default=3,
                   help='Number nucleotides (only +1 possible) to have a polyA, single-stranded hang between sequence of interest and pad.')
m2seq.add_argument('--pad_hang', type=int, default=0,
                   help='Number nucleotides (only +1 possible) to have a random, single-stranded hang between sequence of interest and pad.')
m2seq.add_argument('--pad_num_samples', type=int, default=30,
                   help="Minimum number of sequences to check pad's effect on structure.")
m2seq.add_argument('--pad_to_length', type=int, default=None,
                     help='Length of sequence to pad to, if nothing, longest length in fasta file.')

double = parser.add_argument_group('double mutant')
double.add_argument('--doublemut', nargs='+',
                    help='Two region to do double mutants of (one mutant in group A one in group B) format: 1-12,15.64-70,72-78 where this would mean one mutant in nucleotides 1 to 12 (inclusive) or 15 and one mutant in region 64 to 70 or 72 to 78. 1-12.1-12 would mean all double mutants in 1to 12. If more than one sequence is in the input need to specify the same number of regions separated by space eg for 2 sequences: 1-12,15.64-70,72-78 34-78.80-85 ')

args = parser.parse_args()

###############################################################################
# check library
###############################################################################

if args.check_library:
    barcodes = get_used_barcodes(args.input_fasta, args.avoid_barcodes_start, args.avoid_barcodes_end)
    print('Confirm these are barcodes, otherwise change --avoid_barcodes_start and --avoid_barcodes_end',sample(barcodes,5))
    names = [s.name for s in list(SeqIO.parse(args.input_fasta, "fasta"))]
    if get_same_length(args.input_fasta)[0]:
        print("All sequences are of the same length.")
    else:
        print('ERROR: seuqences are not all the same lengh')
    seq_good, problems = check_sequences_contents(args.input_fasta,seq5=args.seq5, seq3=args.seq3,bases=['A','C','G','T'])
    if seq_good:
        print(f"All sequences have correct 5' {args.seq5} and 3' {args.seq3} regions and no unkown bases.")
    else:
        for p in problems:
            print('ERROR', p)
    all_good = True
    for i,x in enumerate(barcodes):
        for j,y in enumerate(barcodes[i+1:]):
            if edit_distance(x,y) < args.min_edit:
                all_good = False
                print(f'ERROR {names[i]} has close barcode {x}, with {names[i+j+1]}, {y}, distance is {edit_distance(x,y)} which is less than minimum specified {args.min_edit}.')
    if all_good:
        print(f'All barcodes are unqiue to the specfied edit distance of {args.min_edit}.')


###############################################################################
# just library
###############################################################################

elif args.just_library:
    '''
    # check pad needed
    same_length, length = get_same_length(args.input_fasta)
    if same_length and length==args.pad_to_length:
        fasta = args.input_fasta
    else:
        fasta = f'{args.output_prefix}_pad.fasta'
        add_pad(args.input_fasta,
                fasta,
                share_pad=args.share_pad,
                epsilon_interaction=args.Pmax_noninteract,
                epsilon_punpaired=args.Pmin_unpaired,
                epsilon_avg_punpaired=args.Pavg_unpaired,
                epsilon_paired=args.Pmin_paired,
                epsilon_avg_paired=args.Pavg_paired,
                loop=args.pad_loop,
                hang=args.pad_hang,
                min_num_samples=args.pad_num_samples,
                pad_to_length=args.pad_to_length)
    '''
    used_barcodes = []
    if args.avoid_barcodes_files is not None:
        for file in args.avoid_barcodes_files:
            used_barcodes.extend(get_used_barcodes(file, args.avoid_barcodes_start, args.avoid_barcodes_end))
    '''
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
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit)
    '''

    add_library_elements(args.input_fasta, out_fasta=f'{args.output_prefix}_library.fasta', 
                              share_pad=args.share_pad,
                              seq5=args.seq5, seq3=args.seq3,
                              barcode_num_bp=args.barcode_numbp, 
                              barcode_num5hang=args.barcode_num5randomhang, 
                              barcode_num5polyA=args.barcode_num5polyA,
                              barcode_loop=args.barcode_loop,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit,
                              pad_loop=args.pad_loop,
                              pad_hang=args.pad_hang,
                              # pad_polyAhang, num_barcode_before_reduce, percent_reduce_prob
                              #, ,pad_max_prop_bad, pad_side
                              #min_length_stem,  num_pads_reduce  max_length_stem
                              pad_polyAhang = args.pad_polyAhang,
                              pad_min_num_samples=args.pad_num_samples,
                              pad_to_length=args.pad_to_length)

    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

###############################################################################
# Window
###############################################################################

elif args.window:
    get_windows(args.input_fasta, args.length, args.step,
                f'{args.output_prefix}_windowed.fasta',
                circularize=args.circularize,
                reverse_complement=args.reverse_complement,
                fraction_use=args.prop_windows_keep)
    if args.prop_windows_keep_random != 1:
        randomly_select_seqs(f'{args.output_prefix}_windowed.fasta', f'{args.output_prefix}_windowed_rejected.fasta', args.prop_windows_keep_random)

    used_barcodes = []
    if args.avoid_barcodes_files is not None:
        for file in args.avoid_barcodes_files:
            used_barcodes.extend(get_used_barcodes(file, args.avoid_barcodes_start, args.avoid_barcodes_end))

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
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit)

    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')

###############################################################################
# m2seq
###############################################################################

elif args.m2seq:
    get_all_single_mutants(args.input_fasta, f'{args.output_prefix}_single_mut.fasta')
    combine_fastas([args.input_fasta, f'{args.output_prefix}_single_mut.fasta'], f'{args.output_prefix}_WT_single_mut.fasta')

    '''
    # check pad needed
    same_length, length = get_same_length(args.input_fasta)
    if same_length and length==args.pad_to_length:
        fasta = f'{args.output_prefix}_WT_single_mut.fasta'
    else:
        fasta = f'{args.output_prefix}_WT_single_mut_pad.fasta'
        add_pad(f'{args.output_prefix}_WT_single_mut.fasta',
                fasta,
                share_pad=args.share_pad,
                epsilon_interaction=args.Pmax_noninteract,
                epsilon_punpaired=args.Pmin_unpaired,
                epsilon_avg_punpaired=args.Pavg_unpaired,
                epsilon_paired=args.Pmin_paired,
                epsilon_avg_paired=args.Pavg_paired,
                loop=args.pad_loop,
                hang=args.pad_hang,
                min_num_samples=args.pad_num_samples,
                pad_to_length=args.pad_to_length)
    '''
    used_barcodes = []
    if args.avoid_barcodes_files is not None:
        for file in args.avoid_barcodes_files:
            used_barcodes.extend(get_used_barcodes(file, args.avoid_barcodes_start, args.avoid_barcodes_end))
    '''
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
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit)
    '''
    add_library_elements(f'{args.output_prefix}_WT_single_mut.fasta', out_fasta=f'{args.output_prefix}_library.fasta', 
                              share_pad=args.share_pad,
                              seq5=args.seq5, seq3=args.seq3,
                              barcode_num_bp=args.barcode_numbp, 
                              barcode_num5hang=args.barcode_num5randomhang, 
                              barcode_num5polyA=args.barcode_num5polyA,
                              barcode_loop=args.barcode_loop,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit,
                              pad_loop=args.pad_loop,
                              pad_hang=args.pad_hang,
                              # pad_polyAhang, num_barcode_before_reduce, percent_reduce_prob
                              #, ,pad_max_prop_bad, pad_side
                              #min_length_stem,  num_pads_reduce  max_length_stem
                              pad_polyAhang = args.pad_polyAhang,
                              pad_min_num_samples=args.pad_num_samples,
                              pad_to_length=args.pad_to_length)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')


###############################################################################
# m2seq_with_double
###############################################################################


elif args.m2seq_with_double:
    get_all_single_mutants(args.input_fasta, f'{args.output_prefix}_single_mut.fasta')
    regionAs, regionBs = get_regions_for_doublemut(args.doublemut)
    get_all_double_mutants(args.input_fasta, regionAs, regionBs,out_fasta=f'{args.output_prefix}_double_mut.fasta')
    combine_fastas([args.input_fasta, f'{args.output_prefix}_single_mut.fasta', f'{args.output_prefix}_double_mut.fasta'],
                   f'{args.output_prefix}_WT_single_double_mut.fasta')
    '''
    # check pad needed
    same_length, length = get_same_length(args.input_fasta)
    if same_length and length==args.pad_to_length:
        fasta = f'{args.output_prefix}_WT_single_double_mut.fasta'
    else:
        fasta = f'{args.output_prefix}_WT_single_double_mut_pad.fasta'
        add_pad(f'{args.output_prefix}_WT_single_double_mut.fasta',
                fasta,
                share_pad=args.share_pad,
                epsilon_interaction=args.Pmax_noninteract,
                epsilon_punpaired=args.Pmin_unpaired,
                epsilon_avg_punpaired=args.Pavg_unpaired,
                epsilon_paired=args.Pmin_paired,
                epsilon_avg_paired=args.Pavg_paired,
                loop=args.pad_loop,
                hang=args.pad_hang,
                min_num_samples=args.pad_num_samples,
                pad_to_length=args.pad_to_length)
    '''
    used_barcodes = []
    if args.avoid_barcodes_files is not None:
        for file in args.avoid_barcodes_files:
            used_barcodes.extend(get_used_barcodes(file, args.avoid_barcodes_start, args.avoid_barcodes_end))
    '''
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
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit)
    '''
    add_library_elements(f'{args.output_prefix}_WT_single_double_mut.fasta', out_fasta=f'{args.output_prefix}_library.fasta', 
                              share_pad=args.share_pad,
                              seq5=args.seq5, seq3=args.seq3,
                              barcode_num_bp=args.barcode_numbp, 
                              barcode_num5hang=args.barcode_num5randomhang, 
                              barcode_num5polyA=args.barcode_num5polyA,
                              barcode_loop=args.barcode_loop,
                              epsilon_interaction=args.Pmax_noninteract,
                              epsilon_punpaired=args.Pmin_unpaired,
                              epsilon_avg_punpaired=args.Pavg_unpaired,
                              epsilon_paired=args.Pmin_paired,
                              epsilon_avg_paired=args.Pavg_paired,
                              save_image_folder=args.save_image_folder,
                              save_bpp_fig=args.save_bpp_fig,
                              punpaired_chunk_size=args.max_seq_punpaired_plot,
                              num_barcode_before_reduce=args.num_barcodes_reduce_prob,
                              used_barcodes=used_barcodes,
                              min_edit=args.min_edit,
                              pad_loop=args.pad_loop,
                              pad_hang=args.pad_hang,
                              # pad_polyAhang, num_barcode_before_reduce, percent_reduce_prob
                              #, ,pad_max_prop_bad, pad_side
                              #min_length_stem,  num_pads_reduce  max_length_stem
                              pad_polyAhang = args.pad_polyAhang,
                              pad_min_num_samples=args.pad_num_samples,
                              pad_to_length=args.pad_to_length)
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.csv', file_format='twist')
    format_fasta_for_submission(f'{args.output_prefix}_library.fasta', f'{args.output_prefix}_library.txt', file_format='custom_array')


else:
    print('ERROR program not yet implemented')
