Attempting to replicate IGD in Rust from C:
https://github.com/databio/IGD
https://academic.oup.com/bioinformatics/article/37/1/118/6050710

Current manual test:

Input: /home/drc/IGD_TEST/bedfiles/
Output: /home/drc/IGD_TEST/output/

Full command:

Create
```
cargo run igd create --output /home/drc/IGD_TEST/output/ --filelist /home/drc/IGD_TEST/bedfiles/
```

Search
```
cargo run igd search -d /home/drc/IGD_TEST/output/igd_database.igd -q /home/drc/IGD_TEST/bedfiles/test_small_bed_file.bed

```