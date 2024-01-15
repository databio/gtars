Attempting to replicate IGD in Rust from C:
https://github.com/databio/IGD
https://academic.oup.com/bioinformatics/article/37/1/118/6050710

Current manual test:

Input: /home/drc/IGD_TEST/bedfiles/
Output: /home/drc/IGD_TEST/output/

Full command:
```
cargo run igd --output /home/drc/IGD_TEST/output/ --filelist /home/drc/IGD_TEST/bedfiles/
```