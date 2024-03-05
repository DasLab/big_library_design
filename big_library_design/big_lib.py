import os
from itertools import product
from random import shuffle

from Bio import SeqIO
from arnie.bpps import bpps

from big_library_design.mutational_utils import get_reverse_complement, plot_punpaired, _update_seq_count

BASES = ["A", "C", "G", "U"]
INCOMPLETE_SEQ = {'N': ['A', 'C', 'T', 'G'],
                  'R': ['A', 'G'],
                  'Y': ['C', 'T'],
                  'S': ['C', 'G'],
                  'W': ['A', 'T'],
                  'K': ['T', 'G'],
                  'M': ['A', 'C'],
                  'B': ['C', 'T', 'G'],
                  'D': ['A', 'T', 'G'],
                  'H': ['A', 'C', 'T'],
                  'V': ['A', 'C', 'G'], }


class LibSeq:
    def __init__(self, seq, name=None, description=None):

        # save sequence and name
        if type(seq) == str:
            assert name is not None, "ERROR when creating a LibSeq by string, need to include name"
            self.soi = seq
            self.name = name
            self.description = description
        elif isinstance(seq, SeqIO.SeqRecord):
            self.soi = str(seq.seq)
            self.name = seq.id
            if seq.id.strip() != seq.description.strip():
                self.description = seq.description
            else:
                self.description = ""
        else:
            assert True, "ERROR did recongize seq type"

        # save padding and barcode if they exist
        # TODO already has barcode or pad
        self.pad3 = ""
        self.pad5 = ""
        self.barcode = ""

        # standardize and check sequence
        self.standardize_seq()
        incorrect_bases = self.check_seq_validity()
        if len(incorrect_bases) > 0:
            print(f"ERROR non {incorrect_bases} base found.")

        # initialize variables not yet
        self.bpp = None
        self.complete = False
        self.problems = {}
        self.struct_str = ""

    def standardize_seq(self):
        # standardize all sequence --> upper case, RNA not DNA
        self.soi = self.soi.upper().replace("T", "U")
        self.pad3 = self.pad3.upper().replace("T", "U")
        self.pad5 = self.pad5.upper().replace("T", "U")
        self.barcode = self.barcode.upper().replace("T", "U")

    def check_seq_validity(self):
        # check if non standard letters, return them
        nts = set(self.get_seq())
        for b in BASES:
            if b in nts:
                nts.remove(b)
        return nts

    def get_seq(self, region="all"):
        # TODO other options
        return self.pad5 + self.soi + self.pad3 + self.barcode

    def parse_struct_str(self):
        regions = {'unpaired': [], 'pairedA': [], 'pairedB': [],
                   'noninteractA': [], 'noninteractB': []}
        for i, s in enumerate(self.struct_str):
            if s in ['l', "L"]:
                regions['unpaired'].append(i)
            elif s in ['a', '(', '[', '{', '<']:
                regions['pairedA'].append(i)
            elif s in ['b', ')', ']', '}', '>']:
                regions['pairedB'].append(i)
            
            # TODO warning if str weird?
            # TODO this only catches soi - other interactions, ok?
            if s in ['0', '1']:
                regions['noninteractA'].append(i)
            elif s != 'C':
                regions['noninteractB'].append(i)
        # TODO need to fix this properly
        regions['pairedB'] = [regions['pairedB'][::-1]]
        regions['pairedA'] = [regions['pairedA']]
        return regions

    def get_bpp(self, const5, const3, save=True, package="eternafold"):
        # get base pair probability matrix and save
        full_seq = const5 + self.get_seq() + const3
        bpp = bpps(full_seq, package=package)
        if save:
            self.bpp = bpp
        return bpp

    # TODO clean-up with a probability dictionary, like factor dict
    def check_struct(self, regions, epsilon_prob={'unpaired': 0.75, 'avg_unpaired': 0.85, 'paired': 0.75, 'avg_paired': 0.85, 'interaction': 0.09},
                     prob_factor={'unpaired': 1, 'paired': 1, 'interaction': 1}):

        
        if self.bpp is None:
            self.get_bpp('', '')  # TODO figure out const5,const3

        p_unpaired = 1-self.bpp.sum(axis=0)

        # for regions specified, check if structure is ok
        results = {"unpaired_fail": [], "paired_fail": [],
                   "interaction_fail": [], "p_unpaired": p_unpaired}
        if regions['unpaired'] is not None:
            for i, region_unpaired in enumerate(regions['unpaired']):
                if region_unpaired != []:
                    punpaired_check = p_unpaired[region_unpaired]
                    if ((punpaired_check.mean() < epsilon_prob['avg_unpaired']*prob_factor['unpaired']) or
                            (punpaired_check.min() < epsilon_prob['unpaired']*prob_factor['unpaired'])):
                        results["unpaired_fail"].append([i, region_unpaired])

        if regions['noninteractA'] is not None and regions['noninteractA'] != [] and regions['noninteractB'] != []:
            # TODO only check soi and other interaction not 
            #for i, (regA, regB) in enumerate(zip(regions['noninteractA'], regions['noninteractB'])):
            #    if regA != [] and regB != []:
            interaction_check = self.bpp[regions['noninteractA']][:, regions['noninteractB']] #self.bpp[regA][:, regB]
            if (interaction_check.max() > epsilon_prob['interaction']/prob_factor['interaction']):
                results["interaction_fail"].append([regions['noninteractA'],regions['noninteractB']])#[i, regA, regB])

        if regions['pairedA'] is not None and regions['pairedA'] != []:
            for i, (regA, regB) in enumerate(zip(regions['pairedA'], regions['pairedB'])):
                if regA != [] and regB != []:
                    paired_check = self.bpp[regA, regB]
                    if ((paired_check.mean() < epsilon_prob['avg_paired']*prob_factor['paired']) or
                            (paired_check.min() < epsilon_prob['paired']*prob_factor['paired'])):
                        results["paired_fail"].append([i, regA, regB])

        if results["unpaired_fail"] + results["interaction_fail"] + results["paired_fail"] == []:
            results['pass'] = True
        else:
            results['pass'] = False

        return results

    def prepare_seq(self, length, const5, const3, unused_barcodes, bc_str, epsilon_prob,prob_factor={'unpaired': 1, 'paired': 1, 'interaction': 1}):
        # TODO padding
        self.soi = self.soi[:length]  # TODO real padding...
        good_seq = False
        self.struct_str = ("C"*len(const5)) + \
            ("0"*len(self.soi)) + bc_str + ("C"*len(const3))
        # TODO padding, (,[,{,<>}]) and L for unpaired
        # TODO these parameters should be changeabl
        numbc_reduce_prob = 10
        percent_reduce_prob = 10
        regions = self.parse_struct_str()

        seq_count = {'unpaired': 0, 'paired': 0, 'interaction': 0}
        while not self.complete:
            self.barcode = unused_barcodes.pop(0)
            self.get_bpp(const5, const3)
            struct = self.check_struct(regions,epsilon_prob,prob_factor)
            if struct['pass']:
                self.complete = True
                # TODO register any problems, implement after educe probs stuff?
            else:
                unused_barcodes.append(self.barcode)
                seq_count = _update_seq_count(seq_count,struct)
                for type_error, count in seq_count.items():
                    mutiplier = (count // numbc_reduce_prob)
                    prob_factor[type_error] = (1-(percent_reduce_prob/100))**mutiplier
                    if (count-mutiplier) % numbc_reduce_prob == 0 and count != 0:
                        seq_count[type_error] += 1
                        print(f'For {self.name}, failed to find pad from {count-mutiplier} barcodes because of {type_error}, reducing probabilities needed by a factor of {prob_factor[type_error]}.')

        return unused_barcodes
        # TODO when to give up
        # TODO need to prepend if not good
        # TODO images to save?

class Library:

    def __init__(self, input_file, lib_length, barcode_added=False, const5="GGGAACGACTCGAGTAGAGTCGAAAA", const3="AAAAGAAACAACAACAACAAC", barcode_numbp=8, barcode_loop="UUCG", barcode_num5randomhang=0, barcode_num5polyA=2, min_edit=2, barcode_stem_gu_present=False, num_replicates=1, Pmax_noninteract=0.05, Pmin_paired=0.75, Pavg_paired=0.85, Pmin_unpaired=0.75, Pavg_unpaired=0.85):

        # initialize
        self.libseqs = []
        self.unused_barcodes = None
        self.length = lib_length

        # save information about library
        self.const5 = const5.upper().replace("T", "U")
        self.const3 = const3.upper().replace("T", "U")
        # TODO check validity

        # save data about structure of the barcode
        self.bc_added = barcode_added
        self.bc_numbp = barcode_numbp
        self.bc_loop = barcode_loop
        self.bc_num5randomghang = barcode_num5randomhang
        self.bc_num5polyA = barcode_num5polyA
        self.bc_minedit = min_edit  # TODO only really works with 2 currently?
        self.bc_gu_present = barcode_stem_gu_present

        # set barcode string
        self.bc_str = ("l"*(self.bc_num5randomghang+self.bc_num5polyA)) + \
            ("a"*self.bc_numbp) + ("l"*len(self.bc_loop)) + ("b"*self.bc_numbp)

        # save information about sequence of interest
        self.num_replicates = num_replicates

        # save information about structure checks
        self.struct_pmax_noninteract = Pmax_noninteract
        self.struct_pmin_paired = Pmin_paired
        self.struct_pavg_paired = Pavg_paired
        self.struct_pmin_unpaired = Pmin_unpaired
        self.struct_pavg_unpaired = Pavg_unpaired

        if ".fa" in input_file or ".fasta" in input_file:
            self.seqs_from_fasta(input_file)
        else:
            print("ERROR not yet implemented input type")

    def seqs_from_fasta(self, fasta):
        assert os.path.isfile(fasta), "ERROR file does not extist."
        new_seqs = list(SeqIO.parse(fasta, "fasta"))
        new_seqs = [LibSeq(s) for s in new_seqs]
        self.libseqs.extend(new_seqs)

    # TODO this whole function??
    def seqs_from_spreasheet(self, sheet):
        assert os.path.isfile(fasta), "ERROR file does not extist."
        if ".tsv" in sheet:

            df = pd.read_table(seqs)
            cols = df.columns
            name_seq_columns = [['id', 'sequence']]
            for name_col, seq_col in name_seq_columns:
                if name_col in cols and seq_col in cols:
                    df[seq_col] = df[seq_col].str.strip()
                    return df.apply(lambda row: SeqIO.SeqRecord(Seq.Seq(_get_dna_from_SeqRecord(row[seq_col].strip())),
                                                                name=str(row[name_col])), axis=1).to_list()
            print("ERROR recognized tsv but did not recognize column names, please submmit and issue to request new input format.")
        else:
            print("ERROR unrecognized input file extension, please submit an issue to have the input format recognized.")

    def get_all_bpps(self):
        for seq in self.libseqs:
            seq.get_bpp(self.const5, self.const3)

    def get_barcodes(self):
        # TODO need way to ignore barcodes
        print("Getting all possible barcodes.")
        all_barcodes = []
        for x in product(BASES, repeat=self.bc_num5randomghang+self.bc_numbp):
            uid = ''.join(x)

            # split barcode in stem and hang regions
            hang5 = uid[:self.bc_num5randomghang]
            stemA = uid[self.bc_num5randomghang:]
            stemB = get_reverse_complement(stemA).replace("T","U")

            # put all barcode parts together
            seq = ("A"*self.bc_num5polyA)+hang5+stemA+self.bc_loop+stemB
            all_barcodes.append(seq)
        shuffle(all_barcodes)
        self.unused_barcodes = all_barcodes

    def prepare(self):
        assert not self.bc_added, "ERROR nothing to prepare, user specified library already preparred"
        if self.unused_barcodes is None:
            unused_barcodes = self.get_barcodes()
            # TODO save to file instead???
        for seq in self.libseqs:

            self.unused_barcodes = seq.prepare_seq(self.length, self.const5, self.const3, self.unused_barcodes,
                                                   self.bc_str,{'unpaired': self.struct_pmin_unpaired, 'avg_unpaired': self.struct_pavg_unpaired, 'paired': self.struct_pmin_paired, 'avg_paired': self.struct_pavg_paired, 
                'interaction': self.struct_pmax_noninteract})

    def plot_punpaired(self, save_image):
        p_unpaired = {}
        seqs = []
        muts = []  # location of 1 in struct_str
        pad_lines = []  # if empty don't draw anything
        # TODO save mutation name etc in common format?
        for seq in self.libseqs:
            seqs.append(self.const5 + seq.get_seq() + self.const3)
            p_unpaired[seq.name+seq.description] = 1-seq.bpp.sum(axis=0)
            # TODO only if type is M2seq...???
            muts.append([i for i, x in enumerate(seq.struct_str) if x == "1"])
            if "L" in seq.struct_str:  # has pads
                # TODO or 1 ...
                pads_lines.append(
                    [seq.struct_str.index("0"), seq.struct_str.rindex("0")])
            else:
                pad_lines.append([])
        lines = [len(self.const5), len(self.const5)+self.length, len(self.const5) +
                 self.length+len(self.bc_str)]  # locatino of barcode, const
        xlabels = [x if x%10==0 else '' for x in range(1,1+len(seqs[0]))]
        plot_punpaired(p_unpaired, xlabels, seqs, muts,
                       lines, pad_lines, save_image)

    # TODO output need to change U to T again, probably an Experiment level property?


class WindowLibrary(Library):
    def __init__(self, const5, const3):
        Library.__init__(self, const5, const3)
        # TODO need to windowize


class Expirement:
    def __init__(self):
        self.libs = []


'''
        # TODO check validite --> A,T,C,G       
    def _fill_in_any_incomplete(seq,seqs):

                    if len(new_seq)*0.9>new_seq.count('A')+new_seq.count('C')+new_seq.count('T')+new_seq.count('G'):
                print(f'window {new_seq} to uncertain deleting')
                new_seqs = []
            else:
                new_seqs = _fill_in_any_incomplete(new_seq,[new_seq])
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
'''
