# Das Lab Library Design

## This repository is currently under construction, please expect large changes to occur frequently while this message is here. Please submit issues or pull requests for any bugs, thanks!

A collection of functionalities for preparing library for RNA structure probing.

## Example functionalities

Major functionalities have been scripted in `prepare_library.py`. Please see `python prepare_library.py -h` for details on the options and example runs below. Note these scripted functionalities are a subset of what is possible with the code base. Please submit an issue if you want an additional library preperation method to be added to this main script.

```
usage: prepare_library.py [-h] -i INPUT_FASTA -o OUTPUT_PREFIX
                          [--just_library | --window | --m2seq | --m2seq_with_double]
                          [--save_bpp_fig SAVE_BPP_FIG] [--save_image_folder SAVE_IMAGE_FOLDER]
                          [--max_seq_punpaired_plot MAX_SEQ_PUNPAIRED_PLOT]
                          [--Pmax_noninteract PMAX_NONINTERACT] [--Pmin_paired PMIN_PAIRED]
                          [--Pavg_paired PAVG_PAIRED] [--Pmin_unpaired PMIN_UNPAIRED]
                          [--Pavg_unpaired PAVG_UNPAIRED] [--seq5 SEQ5] [--seq3 SEQ3]
                          [--barcode_numbp BARCODE_NUMBP]
                          [--barcode_num5polyA BARCODE_NUM5POLYA]
                          [--barcode_num5randomhang BARCODE_NUM5RANDOMHANG]
                          [--barcode_loop BARCODE_LOOP] [--length LENGTH] [--step STEP]
                          [--circularize] [--pad_type PAD_TYPE] [--pad_loop PAD_LOOP]
                          [--pad_min_hang PAD_MIN_HANG] [--pad_num_samples PAD_NUM_SAMPLES]
                          [--doublemut DOUBLEMUT [DOUBLEMUT ...]]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Sequences to make library out of, file in fasta format.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix of the output files.

programs:
  --just_library        The just_libary program will take all sequences provided and then just
                        prepare the libary using only those sequences.
  --window              The window program will take each sequence, create sliding windows and
                        then prepare the library sequences.
  --m2seq               The m2seq program will take each sequence, get all single mutants,
                        create a noninteracting pad to ensure each sequence is the same length
                        if needed, and then prepare the library sequences.
  --m2seq_with_double   The program runs the m2seq program with the addition of double mutants
                        in user defined regions.

Visuals:
  --save_bpp_fig SAVE_BPP_FIG
                        Proportion to save images of the base-pair-probability matrices. 0 is
                        none, 1 is all, in between is random seleciton of sequences.
  --save_image_folder SAVE_IMAGE_FOLDER
                        Folder to save images (proabbility unpaired and if specified base-pair-
                        probability matrices), default don't save.
  --max_seq_punpaired_plot MAX_SEQ_PUNPAIRED_PLOT
                        The maximum number of sequences to put in one probability unpaired plot,
                        will split into seperate images if smaller than total number of
                        sequences.

Strctural checks:
  --Pmax_noninteract PMAX_NONINTERACT
                        Maximum base-pair-probability for 2 regions to be considered
                        noninteracting.
  --Pmin_paired PMIN_PAIRED
                        Minimum base-pair-probability for 2 regions to be considered paired.
  --Pavg_paired PAVG_PAIRED
                        Average base-pair-probability for 2 regions to be considered paired.
  --Pmin_unpaired PMIN_UNPAIRED
                        Minimum probability unpaired for region to be unpaired.
  --Pavg_unpaired PAVG_UNPAIRED
                        Average probability unpaired for region to be unpaired.

Libarary parts:
  --seq5 SEQ5           Constant sequence to place at 5' of every sequence in library.
  --seq3 SEQ3           Constant sequence to place at 3' of every sequence in library.
  --barcode_numbp BARCODE_NUMBP
                        The length (in bp) of the random barcode stem.
  --barcode_num5polyA BARCODE_NUM5POLYA
                        The length of polyA stretch placed before the barcode stem, if
                        specified, also put before the random 5'hang..
  --barcode_num5randomhang BARCODE_NUM5RANDOMHANG
                        The legnth of addtional random (single-standed) sequence to place before
                        the barcode stem.
  --barcode_loop BARCODE_LOOP
                        Constant loop sequence of the barcode stemloop.

window:
  --length LENGTH       Length of each window to create.
  --step STEP           The step size of the sliding window.
  --circularize         Whether to circularize the sequence (at 3' end, don't stop but loop back
                        to 5') or not.

padding for: m2seq or just_library when length not equal:
  --pad_type PAD_TYPE   If there are sequencees of multiple lengths, to obtain a libarary of
                        equal length some sequences will be padded. This specifies the type of
                        padding with SL_same_per_length (if pad is long enough create a stem-
                        loop, same pad is used for each group of sequences with equal length) or
                        rand_same_all (a single-stranded non-interacting pad is choosen, same
                        for all sequences just using the length of pad required to pad each
                        sequence to same length) as options.
  --pad_loop PAD_LOOP   If padtype is a stem-loop the constant loop to use.
  --pad_min_hang PAD_MIN_HANG
                        If padtype is a stem-loop the minimum (only +1 possible) to have a
                        random, single-stranded hang between sequence of interest and pad.
  --pad_num_samples PAD_NUM_SAMPLES
                        Minimum number of sequences to check pad's effect on structure.

double mutant:
  --doublemut DOUBLEMUT [DOUBLEMUT ...]
                        Two region to do double mutagenis of (one mutant in group A one in group
                        B) format: 1-12,15.64-70,72-78 where this would mean one mutant in
                        nucletoides 1to 12 (inclusive) or 15 and one mutant in region 64 to 70
                        or 72 to 78. 1-12.1-12 would mean all double mutants in 1to 12. If more
                        than one sequence is in the input need to specify the same number of
                        regions seperated by space eg for 2 sequences: 1-12,15.64-70,72-78
                        34-78.80-85
```

### Creating library from already prepared sequence list

```
python prepare_library.py --just_library -i examples/example_WT.fasta -o examples/just_library_ex_output/example --save_bpp_fig 1 --save_image_folder examples/just_library_ex_output/figs/
```

Example output is in `examples/just_library_ex_output`, note you should not expect to see same padding or barcodes as these are randomly generated. 

### Creating libary of sliding windows

```
python prepare_library.py --window -i examples/example_WT.fasta -o examples/window_ex_output/example --length 100 --step 10 --Pmax_noninteract 0.15 --save_bpp_fig 1 --save_image_folder examples/window_ex_output/figs/ --max_seq_punpaired_plot 20
```

Example output is in `examples/window_ex_output`, note you should not expect exact same sequences as barcodes are randomly generated. Only 1 example base-pair-probaility matrix figure uploaded to this github, but this call will generate one for all sequences.

### Creating library for M2Seq (all single-mutants)

```
python prepare_library.py --m2seq -i examples/example_WT.fasta -o examples/m2seq_ex_output/example --Pmax_noninteract 0.15 --Pmin_paired 0.7 --Pavg_paired 0.8 --Pmin_unpaired 0.6 --Pmin_unpaired 0.8 --save_bpp_fig 1 --save_image_folder examples/m2seq_ex_output/figs/
```

Example output is in `examples/m2seq_ex_output`, note you should expect to see the same mutants in the region of interest, but you should not expect to see same pads or barcodes. Only 1 example base-pair-probability matrix figure is uploaded tothis github, but this call will generate one for every sequence in the libaray.

### Creating library with all single-mutants and select double

```
python prepare_library.py --m2seq_with_double -i examples/example_WT.fasta -o examples/double_ex_output/example --Pmax_noninteract 0.15 --Pmin_paired 0.7 --Pavg_paired 0.8 --Pmin_unpaired 0.6 --Pmin_unpaired 0.8 --doublemut 20-23.24-28 56.82,90 130,140.150-151 --save_bpp_fig 1 --save_image_folder examples/double_ex_output/figs/
```

Example output is in `examples/double_ex_output`, note you should expect to see the same mutants in the region of interest, but you should not expect to see same pads or barcodes. Only 1 example base-pair-probability matrix figure is uploaded to this github, but this call will generate one for every sequence in the libaray.
