# Current Steps to Run Uniwig

### Input Bed File

Currently, Uniwig accepts a single `.bed` file. It should be sorted by chromosome.

The below script can be used to create a sorted bed file from a directory of bed files:

```shell
#!/bin/sh
# directory for the raw data (bed files)
RAWDATA_DIR="./data/raw/"
# directory for combined data
COMBDATA_DIR="./data/combined/"
# raw data filename
raw="*.bed"
# unsorted combined data filename
unsorted="combined_unsort.bed"
# chrsorted combined data filename
chrsorted="combined_chrsort.bed"
cat $RAWDATA_DIR$raw > $COMBDATA_DIR$unsorted
sort -k1,1V $COMBDATA_DIR$unsorted > $COMBDATA_DIR$chrsorted
```
### Running uniwig

Once you have your single, sorted bedfile, you can run uniwig with the following command:

```
cargo run uniwig -b /home/drc/Downloads/uniwig_testing_19apr2024/sourcefiles/test_30_lines_sorted.bed -c /home/drc/Downloads/uniwig_testing_19apr2024/sourcefiles/hg38.chrom.sizes -m 5 -t 1 -l /home/drc/Downloads/uniwig_testing_19apr2024/wiggles_created_with_rust/final_wiggles/ -y wig

```

Note that we provide a chrom.sizes reference file (hg38) in the testing folder -> `genimtools/tests/hg38.chrom.sizes`

### Usage
```
Usage: genimtools uniwig --bed <bed> --chromref <chromref> --smoothsize <smoothsize> --stepsize <stepsize> --fileheader <fileheader> --outputtype <outputtype>

Options:
  -b, --bed <bed>                Path to the combined bed file we want to tranforms
  -c, --chromref <chromref>      Path to chromreference
  -m, --smoothsize <smoothsize>  Integer value for smoothing
  -t, --stepsize <stepsize>      Integer value for stepsize
  -l, --fileheader <fileheader>  Name of the file
  -y, --outputtype <outputtype>  Output as wiggle or CSV
  -h, --help                     Print help

```

### Create bigwig files from wiggle files

Once you have created wiggle files, you can convert them to bigWig files using `wigToBigWig` (see: https://genome.ucsc.edu/goldenPath/help/bigWig.html, https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils/wigToBigWig):

```
./wigToBigWig ./test_rust_wig/_end.wig ./sourcefiles/hg38.chrom.sizes ./end_rust.bw
```

### Export types

Currently only `.wig` is supported as an output type. 