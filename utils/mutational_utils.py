# function for making mutations from WT sequence

from Bio import SeqIO,Seq
from itertools import product

BASES = ['A','C','G','T']
SEQ5 = 'GGGAACGACTCGAGTAGAGTCGAAAA'
SEQ3 = 'AAAAGAAACAACAACAACAAC'# with the single A spacer added already
TETRALOOP = 'TTCG'
MAXPROB = 0.1

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
	return None

def add_fixed_seq_and_barcode(seqs):
	return None

def check_bpp(seq,region,epsilon=MAXPROB):
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
					new_mut = SeqIO.SeqRecord(Seq.Seq(new_seq),name,name,name)
					all_single_mutants.append(new_mut)

	SeqIO.write(all_single_mutants,out_fasta,"fasta")

get_all_single_mutants('../examples/example_WT.fasta','../examples/example_single_mut.fasta')
write_all_barcodes('../examples/all_8bp_barcodes.fasta')