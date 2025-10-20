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


USING BLOOM FILTERS

CREATE COMMAND EXAMPLE

```
./target/release/gtars igd bloom --action create 
--universe "/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed" 
--bedfilesuniverse "/home/drc/Downloads/bloom_testing/test1/two_real_bed_files/" 
--bloomdirectory "/home/drc/Downloads/bloom_testing/test1/" 
--bloomname "test" 
--numitems 10000 
--falsepositive 0.001
```

SEARCH COMMAND EXAMPLE

```
./target/release/gtars igd bloom --action search 
--universe "/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed" 
--bloomdirectory "/home/drc/Downloads/bloom_testing/test1/" 
--bloomname "test" 
--querybed "/home/drc/Downloads/bloom_testing/test1/query2.bed"

```

