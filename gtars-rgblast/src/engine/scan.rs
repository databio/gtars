//! Subject scanning: find lookup-table words in a packed (ncbi2na) subject.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_nascan.c` (`s_MBScanSubject_Any`,
//! `s_BlastMBLookupRetrieve`, `s_BlastSmallNaScanSubject_Any`, `s_BlastSmallNaRetrieveHits`,
//! `s_BlastNaScanSubject_Any`, `s_BlastNaScanSubject_8_4`, `s_BlastLookupGetNumHits`,
//! `s_BlastLookupRetrieve`, `s_NaChooseScanSubject`).
//! The width-specialised kernels (`s_MBScanSubject_9_1` .. `_11_3Mod4`,
//! `s_BlastSmallNaScanSubject_4_1` .. `_8_3Mod4`, `_8_4`) are `mb_scan_w`, `small_na_scan_w`
//! and `small_na_scan_8_4`, chosen by `mb_scan` / `small_na_scan` as `s_MBChooseScanSubject` /
//! `s_SmallNaChooseScanSubject` do; they emit the same seeds in the same order as the generic
//! kernels and stop at the same offset. `subject` must have at least 4 readable bytes past the
//! last word (the C over-reads too; the arena pads every subject).

use super::lookup::{MbTable, NaTable, SmallNaTable};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OffsetPair {
    pub q_off: i32,
    pub s_off: i32,
}

/// `s_MBScanSubject_Any`. `scan_range` = [start, end] (end inclusive); on return
/// `scan_range[0]` is where the next call resumes. Returns the number of pairs written.
pub fn mb_scan_any(
    t: &MbTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mask = t.hashsize - 1;
    let mut total_hits: i32 = 0;
    let w = t.lut_word_length;
    let mut scan_step = t.scan_step;
    let max_hits = max_hits - t.longest_chain as i32;
    let retrieve = |index: u64, s_off: i32, total: i32, out: &mut [OffsetPair]| -> i32 {
        let mut i = 0i32;
        let mut q_off = t.hashtable[index as usize];
        while q_off != 0 {
            out[(total + i) as usize] = OffsetPair {
                q_off: q_off - 1,
                s_off,
            };
            i += 1;
            q_off = t.next_pos[q_off as usize];
        }
        i
    };
    if scan_step % 4 == 0 {
        let s_end = (scan_range[1] / 4) as usize;
        let shift = 2 * (12 - w);
        let mut s = (scan_range[0] / 4) as usize;
        scan_step /= 4;
        while s <= s_end {
            let index =
                ((subject[s] as u64) << 16 | (subject[s + 1] as u64) << 8 | subject[s + 2] as u64)
                    >> shift;
            if t.has_hits(index) {
                if total_hits >= max_hits {
                    break;
                }
                total_hits += retrieve(index, (s * 4) as i32, total_hits, offset_pairs);
            }
            s += scan_step as usize;
        }
        scan_range[0] = (s * 4) as i32;
    } else if w == 16 {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (4 - (scan_range[0] % 4));
            let s = (scan_range[0] / 4) as usize;
            let wv = (subject[s] as u64) << 32
                | (subject[s + 1] as u64) << 24
                | (subject[s + 2] as u64) << 16
                | (subject[s + 3] as u64) << 8
                | subject[s + 4] as u64;
            let index = (wv >> shift) & mask;
            if t.has_hits(index) {
                if total_hits >= max_hits {
                    break;
                }
                total_hits += retrieve(index, scan_range[0], total_hits, offset_pairs);
            }
            scan_range[0] += scan_step;
        }
    } else if w > 9 {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (16 - (scan_range[0] % 4 + w));
            let s = (scan_range[0] / 4) as usize;
            let v = (subject[s] as u64) << 24
                | (subject[s + 1] as u64) << 16
                | (subject[s + 2] as u64) << 8
                | subject[s + 3] as u64;
            let index = (v >> shift) & mask;
            if t.has_hits(index) {
                if total_hits >= max_hits {
                    break;
                }
                total_hits += retrieve(index, scan_range[0], total_hits, offset_pairs);
            }
            scan_range[0] += scan_step;
        }
    } else {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (12 - (scan_range[0] % 4 + w));
            let s = (scan_range[0] / 4) as usize;
            let v =
                (subject[s] as u64) << 16 | (subject[s + 1] as u64) << 8 | subject[s + 2] as u64;
            let index = (v >> shift) & mask;
            if t.has_hits(index) {
                if total_hits >= max_hits {
                    break;
                }
                total_hits += retrieve(index, scan_range[0], total_hits, offset_pairs);
            }
            scan_range[0] += scan_step;
        }
    }
    total_hits
}

/// `s_BlastSmallNaScanSubject_Any`.
pub fn small_na_scan_any(
    t: &SmallNaTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mut total_hits: i32 = 0;
    let mut scan_step = t.scan_step;
    let mask = t.mask;
    let w = t.lut_word_length;
    let max_hits = max_hits - t.longest_chain as i32;
    let retrieve = |index: i32, s_off: i32, total: i32, out: &mut [OffsetPair]| -> i32 {
        if index >= 0 {
            out[total as usize] = OffsetPair {
                q_off: index,
                s_off,
            };
            return 1;
        }
        let mut n = 0i32;
        let mut src = (-index) as usize;
        let mut idx = t.overflow[src] as i32;
        src += 1;
        loop {
            out[(total + n) as usize] = OffsetPair { q_off: idx, s_off };
            n += 1;
            idx = t.overflow[src] as i32;
            src += 1;
            if idx < 0 {
                break;
            }
        }
        n
    };
    if w > 5 {
        if scan_step % 4 == 0 {
            let s_end = (scan_range[1] / 4) as usize;
            let shift = 2 * (8 - w);
            let mut s = (scan_range[0] / 4) as usize;
            scan_step /= 4;
            while s <= s_end {
                let v = (subject[s] as i32) << 8 | subject[s + 1] as i32;
                let index = t.final_backbone[(v >> shift) as usize] as i32;
                if index != -1 {
                    if total_hits > max_hits {
                        break;
                    }
                    total_hits += retrieve(index, (s * 4) as i32, total_hits, offset_pairs);
                }
                s += scan_step as usize;
            }
            scan_range[0] = (s * 4) as i32;
        } else {
            while scan_range[0] <= scan_range[1] {
                let shift = 2 * (12 - (scan_range[0] % 4 + w));
                let s = (scan_range[0] / 4) as usize;
                let v = (subject[s] as i32) << 16
                    | (subject[s + 1] as i32) << 8
                    | subject[s + 2] as i32;
                let index = t.final_backbone[((v >> shift) & mask) as usize] as i32;
                if index != -1 {
                    if total_hits > max_hits {
                        break;
                    }
                    total_hits += retrieve(index, scan_range[0], total_hits, offset_pairs);
                }
                scan_range[0] += scan_step;
            }
        }
    } else {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (8 - (scan_range[0] % 4 + w));
            let s = (scan_range[0] / 4) as usize;
            let v = (subject[s] as i32) << 8 | subject[s + 1] as i32;
            let index = t.final_backbone[((v >> shift) & mask) as usize] as i32;
            if index != -1 {
                if total_hits > max_hits {
                    break;
                }
                total_hits += retrieve(index, scan_range[0], total_hits, offset_pairs);
            }
            scan_range[0] += scan_step;
        }
    }
    total_hits
}

/// `s_BlastLookupRetrieve` into `out[total..]`; returns the number written.
#[inline]
fn na_retrieve(t: &NaTable, index: i32, s_off: i32, total: i32, out: &mut [OffsetPair]) -> i32 {
    let pos = t.positions(index);
    for (i, &q) in pos.iter().enumerate() {
        out[total as usize + i] = OffsetPair { q_off: q, s_off };
    }
    pos.len() as i32
}

/// `s_BlastNaScanSubject_8_4` (`s_NaChooseScanSubject` picks it for width 8, stride 4).
pub fn na_scan_8_4(
    t: &NaTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mut total_hits = 0i32;
    let mut s = (scan_range[0] / 4) as usize;
    let s_end = (scan_range[1] / 4) as usize;
    while s <= s_end {
        let index = (subject[s] as i32) << 8 | subject[s + 1] as i32;
        let num_hits = t.num_hits(index);
        if num_hits == 0 {
            s += 1;
            continue;
        }
        if num_hits > max_hits - total_hits {
            break;
        }
        total_hits += na_retrieve(t, index, (s * 4) as i32, total_hits, offset_pairs);
        s += 1;
    }
    scan_range[0] = (s * 4) as i32;
    total_hits
}

/// `s_BlastNaScanSubject_Any`.
pub fn na_scan_any(
    t: &NaTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mut total_hits = 0i32;
    let mask = t.mask;
    let mut scan_step = t.scan_step;
    let w = t.lut_word_length;
    if w > 5 {
        if scan_step % 4 == 0 {
            let s_end = (scan_range[1] / 4) as usize;
            let shift = 2 * (8 - w);
            let mut s = (scan_range[0] / 4) as usize;
            scan_step /= 4;
            while s <= s_end {
                let index = ((subject[s] as i32) << 8 | subject[s + 1] as i32) >> shift;
                let num_hits = t.num_hits(index);
                if num_hits == 0 {
                    s += scan_step as usize;
                    continue;
                }
                if num_hits > max_hits - total_hits {
                    break;
                }
                total_hits += na_retrieve(t, index, (s * 4) as i32, total_hits, offset_pairs);
                s += scan_step as usize;
            }
            scan_range[0] = (s * 4) as i32;
        } else {
            while scan_range[0] <= scan_range[1] {
                let shift = 2 * (12 - (scan_range[0] % 4 + w));
                let s = (scan_range[0] / 4) as usize;
                let v = (subject[s] as i32) << 16
                    | (subject[s + 1] as i32) << 8
                    | subject[s + 2] as i32;
                let index = (v >> shift) & mask;
                let num_hits = t.num_hits(index);
                if num_hits == 0 {
                    scan_range[0] += scan_step;
                    continue;
                }
                if num_hits > max_hits - total_hits {
                    break;
                }
                total_hits += na_retrieve(t, index, scan_range[0], total_hits, offset_pairs);
                scan_range[0] += scan_step;
            }
        }
    } else {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (8 - (scan_range[0] % 4 + w));
            let s = (scan_range[0] / 4) as usize;
            let v = (subject[s] as i32) << 8 | subject[s + 1] as i32;
            let index = (v >> shift) & mask;
            let num_hits = t.num_hits(index);
            if num_hits == 0 {
                scan_range[0] += scan_step;
                continue;
            }
            if num_hits > max_hits - total_hits {
                break;
            }
            total_hits += na_retrieve(t, index, scan_range[0], total_hits, offset_pairs);
            scan_range[0] += scan_step;
        }
    }
    total_hits
}
