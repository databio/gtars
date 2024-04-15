# Current Manual testing

Full command example:
```
cargo run uniwig -s -b /home/drc/GITHUB/genimtools/genimtools/tests/data/test_sorted_small.bed -c /home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes -m 5 -t 1 -l /home/drc/Downloads/test -y wig

```

# Uniwig

Given a set of bed files, we want to produce 2 [BigWig](http://genome.ucsc.edu/goldenPath/help/bigWig.html) files: one track of the start coordinates, one track of the end coordinates, and one track for core coordinates.

# Usage

CLI or Python Bindings