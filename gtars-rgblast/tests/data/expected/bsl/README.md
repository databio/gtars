`reference.bsl` was written by BLAST+ 2.17.0 from `ids.txt`:

    blastdb_aliastool -seqid_file_in ids.txt -seqid_title reference -seqid_file_out reference.bsl

It covers `SQ.` ids, `gnl|db|tag` and `lcl|x` spellings (normalized to `db:tag` and `x`), a
numeric tag, a digest starting with `-`, a repeated id, and a 300-byte id (`0xFF` + u32 length).
`tests/taxonomy.rs` (`bsl_matches_blastdb_aliastool`) compares `write_bsl` output with it, the
date taken from the reference.
