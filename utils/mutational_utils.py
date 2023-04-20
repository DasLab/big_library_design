# function for making mutations from WT sequence

from Bio import SeqIO,Seq
from itertools import product
from random import shuffle
from arnie.bpps import bpps
import matplotlib.pyplot as plt

BASES = ['A','C','G','T']
SEQ5 = 'GGGAACGACTCGAGTAGAGTCGAAAA'
SEQ3 = 'AAAAGAAACAACAACAACAAC'# with the single A spacer added already
TETRALOOP = 'TTCG'
MAXPROB = 0.05

######### TODO ##########
# specific mutate rescues
# padding
# UID helices
#  26 nt.   5' fixed sequence 
# 240 nt.   Region of interest
#  20 nt.   8 bp barcode with UUCG tetraloop (65,536 possible barcodes)
#   1 nt.   Single A spacer
#  20 nt.   3' fixed sequence

def write_all_barcodes(out_fasta=None,num_bp=8,loop=TETRALOOP,bases=BASES):
	# TODO probably should add ability to randomly generate but this
	# is fast enough for now
	all_barcodes = []
	for x in product(bases,repeat=num_bp):
		uid = ''.join(x)
		seq = uid+loop+get_reverse_complement(uid)
		seq_rec = SeqIO.SeqRecord(Seq.Seq(seq),uid,uid,uid)
		all_barcodes.append(seq_rec)
	if out_fasta is not None:
		SeqIO.write(all_barcodes,out_fasta,"fasta")
	return all_barcodes

def get_reverse_complement(seq):
	reverse = seq[::-1]
	reverse_complement = reverse.replace("T","X").replace("A","T").replace("X","A").replace("G","X").replace("C","G").replace("X","C")
	return reverse_complement

def add_pad(seqs):
	# TODO to add a pad to the sequence that will not interaction
	return None

def get_all_double_mutants(fasta,out_fasta,regionAs,regionBs,bases=BASES):
	# 9 per pair * lenregionA * lenregionB

	all_WT = list(SeqIO.parse(fasta,"fasta"))
	all_double_mutants = [] # unlike single mutant, WT not saved
	if len(all_WT) != len(regionAs) or len(all_WT) != len(regionBs):
		print(f'WARNING: you must list regions (as a list of lists) for all sequences in the fasta, number sequences in fasta {len(all_WT)} and number of regions specified {len(regionAs)} {len(regionBs)}')
	for record,regionA,regionB in zip(all_WT,regionAs,regionBs):
		# convert to standard RNA
		seq = record.seq.upper().replace("U","T")
		# at each position pair, get double mutants
		for i in regionA:
			for mutA in bases:
					if mutA != seq[i]:
						for j in regionB:
							for mutB in bases:
								if mutB != seq[j]:
									name = f' {record.id}_{i}{seq[i]}>{mutA}_{j}{seq[j]}>{mutB}'
									if i < j:
										new_seq = seq[:i]+mutA+seq[i+1:j]+mutB+seq[j+1:]
									else:
										new_seq = seq[:j]+mutB+seq[j+1:i]+mutA+seq[i+1:]
									new_mut = SeqIO.SeqRecord(Seq.Seq(new_seq),name,'','')
									all_double_mutants.append(new_mut)
	SeqIO.write(all_double_mutants,out_fasta,"fasta")

def add_fixed_seq_and_barcode(fasta,out_fasta=None,seq5=SEQ5,seq3=SEQ3,num_bp=8,loop=TETRALOOP,epsilon=MAXPROB,save_image_folder=None):
	# TODO maybe need a min edit distance for the uid?
	# get and randomly shuffle all potential barcoces
	all_uids = write_all_barcodes(num_bp=num_bp,loop=loop)
	shuffle(all_uids)
	current_uid = 0
	all_full_seqs = []
	seqs = SeqIO.parse(fasta,"fasta")
	for seq_rec in seqs:
		seq = str(seq_rec.seq)
		name = seq_rec.name+'_w53barcode'

		# TODO do I need to check the 5 and 3 interction?

		regionA = list(range(len(seq5),len(seq5)+len(seq)))
		regionB = list(range(len(seq5)+len(seq),len(seq5)+len(seq)+len(loop)+(2*num_bp)))
		uid_good = False
		while not uid_good:
			uid = all_uids[current_uid].seq
			full_seq = f'{seq5}{seq}{uid}{seq3}'
			current_uid += 1
			if save_image_folder is not None:
				# TODO for lines pad may also go in there
				uid_good = check_bpp_no_signal(full_seq,regionA,regionB,epsilon=epsilon,
					lines=[len(seq5),len(seq5)+len(seq),len(seq5)+len(seq)+len(loop)+(2*num_bp)],save_image=f'{save_image_folder}/{name}.png')
			else:
				uid_good = check_bpp_no_signal(full_seq,regionA,regionB,epsilon=epsilon)
		all_full_seqs.append(SeqIO.SeqRecord(Seq.Seq(seq),name,name,name))
	if out_fasta is not None:
		SeqIO.write(all_full_seqs,out_fasta,"fasta")
	return all_full_seqs

def check_bpp_no_signal(seq,regionA,regionB,epsilon=MAXPROB,lines=None,save_image=None):
	# TODO save heatmap of these somewhere, with box of region? -- only good ones
	# TODO likely give optiono of package and linear or not
	bpp = bpps(seq,package='eternafold')#,linear=True)
	print(bpp.max())
	if save_image is not None:
		plot_bpp(bpp,seq,lines,save_image)
	region_to_check = bpp[regionA][:,regionB]
	if region_to_check.max()>epsilon:
		return False
	else:
		return True

def plot_bpp(bpp,seq,lines,save_image):
	plt.figure(figsize=(len(seq)*0.12,len(seq)*0.12))
	plt.imshow(bpp, origin='upper', cmap='gist_heat_r')
	xlabels = []
	ylabels = []
	for i,s in enumerate(seq):
		if i%10==0:
			ylabels.append(f'{i} {s}')
			xlabels.append(f'{s}\n{i}')
		else:
			xlabels.append(s)
			ylabels.append(s)
	plt.xticks(range(len(seq)),xlabels,size=8)
	plt.yticks(range(len(seq)),ylabels,size=8)
	for line in lines:
		plt.hlines(line,0,len(seq),color="grey")
		plt.vlines(line,0,len(seq),color="grey")
	plt.xlim(0,len(seq))
	plt.ylim(len(seq),0)
	# color sequence by 5 region etc, maybe even mut spot
	# or just put hline/vline at these boundaries
	plt.savefig(save_image)

def format_seqs_for_submission(fasta):
	return None

def get_all_single_mutants(fasta,out_fasta,bases=BASES):
	all_WT = SeqIO.parse(fasta,"fasta")
	all_single_mutants = list(SeqIO.parse(fasta,"fasta"))
	for record in all_WT:
		# convert to standard RNA
		seq = record.seq.upper().replace("U","T")
		# at each position, get single mutants
		for i in range(len(seq)):
			for mut in bases:
				if mut != seq[i]:
					name = f' {record.id}_{i}{seq[i]}>{mut}'
					new_seq = seq[:i]+mut+seq[i+1:]
					new_mut = SeqIO.SeqRecord(Seq.Seq(new_seq),name,'','')
					all_single_mutants.append(new_mut)

	SeqIO.write(all_single_mutants,out_fasta,"fasta")

get_all_single_mutants('../examples/example_WT.fasta','../examples/example_single_mut.fasta')
write_all_barcodes('../examples/all_8bp_barcodes.fasta')
get_all_double_mutants('../examples/example_WT.fasta','../examples/example_double_mut.fasta',[[3,4],[9,103]],[[10,11,12],[30,31]])
add_fixed_seq_and_barcode('../examples/example_WT.fasta',save_image_folder='../examples/bpps')

#### current procedure ####
# get list of desired sequences: WTs, all single mutants, double mutants of interaction of interest
# if needed generate random pad (half on 5' half on 3'), select random sampling of desired sequences, iterate finding random pad until bpp between sequences and pad is < 0.1 everywhere
# get all possible 8bp barcodes
# for each sequence add 5', selected random barcode (no replacement), and 3', check bpp between sequence and barcode < 0.1 everywhere

