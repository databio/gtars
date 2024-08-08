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

temp comparison
```
cargo run igd create --output /home/drc/IGD_TEST_2/igd_rust_output/ --filelist /home/drc/IGD_TEST_2/source_bedfiles/
```

Search
```
cargo run igd search --database /home/drc/IGD_TEST_2/igd_rust_output/igd_database.igd --query /home/drc/IGD_TEST_2/query_bed_file/igd_query_test.bed
```