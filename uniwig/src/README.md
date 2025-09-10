# Current Steps to Run Uniwig

### Input Bed File

Currently, Uniwig accepts a single `.bed` file, `.narrowPeak` file, `.bam` file. It should be sorted by chromosome. This single file will be used to create 3 output files:
`_start` -> accumulations of start coordinates
`_end` -> accumulations of end coordinates
`_core` -> accumulations of peaks (between starts and ends)

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
awk 'NF {print} END {print ""}' $RAWDATA_DIR$raw > $COMBDATA_DIR$unsorted
sort -k1,1V $COMBDATA_DIR$unsorted | grep '.' > $COMBDATA_DIR$chrsorted
```
### Running uniwig

Once you have your single, sorted bedfile, you can run uniwig with the following command:

```
cargo run uniwig -f /home/drc/Downloads/uniwig_testing_19apr2024/sourcefiles/test_30_lines_sorted.bed -c /home/drc/Downloads/uniwig_testing_19apr2024/sourcefiles/hg38.chrom.sizes -m 5 -s 1 -l /home/drc/Downloads/uniwig_testing_19apr2024/wiggles_created_with_rust/final_wiggles/ -y wig

```

Note that we provide a chrom.sizes reference file (hg38) in the testing folder -> `genimtools/tests/hg38.chrom.sizes`
The chrom.sizes reference is an optional argument. Uniwig will default to using the combined bed file's last chromosome end position to determine chrom size by default.


### Usage
```
Create accumulation files from a BED or BAM file

Usage: gtars uniwig [OPTIONS] --file <file> --smoothsize <smoothsize> --stepsize <stepsize> --fileheader <fileheader> --outputtype <outputtype>

Options:
  -f, --file <file>              Path to the combined bed file we want to transform or a sorted bam file
  -t, --filetype <filetype>      Input file type, 'bed' 'bam' or 'narrowpeak' [default: bed]
  -c, --chromref <chromref>      Path to chromreference
  -m, --smoothsize <smoothsize>  Integer value for smoothing
  -s, --stepsize <stepsize>      Integer value for stepsize
  -l, --fileheader <fileheader>  Name of the file
  -y, --outputtype <outputtype>  Output as wiggle or npy
  -u, --counttype <counttype>    Select to only output start, end, or core. Select `shift` for bam workflows. Defaults to all. [default: all]
  -p, --threads <threads>        Number of rayon threads to use for parallel processing [default: 6]
  -o, --score                    Count via score (narrowPeak only!)
  -a, --no-bamshift              Set bam shift to False, i.e. uniwig will count raw reads without considering read direction.
  -z, --zoom <zoom>              Number of zoom levels (for bw file output only [default: 1]
  -d, --debug                    Print more verbose debug messages?
  -h, --help                     Print help

```

### Processing bam files to bw

Example command
```
gtars uniwig -f "test1_chr1_chr2.bam" -m 5 -s 1 -l /myoutput/directory/test_file_name -y bw -t bam -p 6 -c /genome/alias/hg38/fasta/default/hg38.chrom.sizes -u all

```


### Export types

For Input types: `.bed` and `.narrowPeak`
Output types include `.wig`, `.npy`, `.bedGraph`, and `.bw`

For Input Types: `.bam`
Output types include `.bw` and `.bed`
