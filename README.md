# Das Lab Library Design

A collection of functionalities for preparing library for RNA structure probing.

## Example functionalities

### Creating libary of sliding windows

```
python prepare_library.py --window -i examples/example_WT.fasta -o examples/window_ex_output/example --length 100 --step 10 --Pmax_noninteract 0.15 --save_bpp_fig --save_image_folder examples/window_ex_output/figs/
```

Example output is in `examples/window_ex_output`, note only 1 example base-pair-probaility matrix uploaded to this github.

### Creating library for M2Seq (all single-mutants)


