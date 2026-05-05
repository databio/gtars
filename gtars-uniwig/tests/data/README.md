# gtars-uniwig test fixtures

## `vn_less_pg.bam` / `.bai`

A small chr22 BAM derived from `tests/data/test_chr22_small.bam` (workspace
root) with a hand-crafted `@PG\tID:custom_tool\tPN:custom_tool` record (no
`VN:` tag) injected after the `@SQ` lines. This reproduces the noodles
strict-parse panic that motivated the `bam_header::read_header_lenient`
helper.

Regenerate with:

```bash
samtools view -h tests/data/test_chr22_small.bam > /tmp/test.sam
sed -i '/^@SQ/a @PG\tID:custom_tool\tPN:custom_tool' /tmp/test.sam
samtools view -bS /tmp/test.sam > tests/data/vn_less_pg.bam
samtools index tests/data/vn_less_pg.bam
samtools view -H tests/data/vn_less_pg.bam | awk '/^@SQ/ {
    for (i=1;i<=NF;i++) {
        if ($i ~ /^SN:/) sn=substr($i,4);
        if ($i ~ /^LN:/) ln=substr($i,4);
    }
    print sn"\t"ln
}' > tests/data/vn_less_pg.chrom.sizes
```
