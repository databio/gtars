//! Symmetric DUST, ported from BLAST+ 2.17.0 `algo/dustmask/symdust.cpp` / `symdust.hpp`
//! (`CSymDustMasker`, its `triplets` window, `find_perfect`, `save_masked_regions`,
//! `CIupac2Ncbi2na_converter`) and applied as `algo/blast/api/dust_filter.cpp` does.
//!
//! The converter draws from a default-seeded NCBI `CRandom` for every `N`, one generator per
//! query, exactly as the C++ does (`CSymDustMasker` is constructed per query).

use crate::convert::NcbiRandom;
use std::collections::VecDeque;

const TRIPLET_MASK: u8 = 0x3F;

#[derive(Debug, Clone, Copy)]
struct Perfect {
    first: usize,
    second: usize,
    score: u32,
    len: usize,
}

struct Triplets {
    start: usize,
    stop: usize,
    max_size: usize,
    low_k: u8,
    l: usize,
    r_w: u32,
    r_v: u32,
    c_w: [u8; 64],
    c_v: [u8; 64],
    num_diff: i32,
    /// front = newest triplet.
    list: VecDeque<u8>,
}

#[inline]
fn add_triplet_info(r: &mut u32, c: &mut [u8; 64], t: u8) {
    *r += c[t as usize] as u32;
    c[t as usize] += 1;
}

#[inline]
fn rem_triplet_info(r: &mut u32, c: &mut [u8; 64], t: u8) {
    c[t as usize] -= 1;
    *r -= c[t as usize] as u32;
}

impl Triplets {
    fn new(window: usize, low_k: u8) -> Self {
        Triplets {
            start: 0,
            stop: 0,
            max_size: window - 2,
            low_k,
            l: 0,
            r_w: 0,
            r_v: 0,
            c_w: [0; 64],
            c_v: [0; 64],
            num_diff: 0,
            list: VecDeque::new(),
        }
    }

    fn shift_high(&mut self, t: u8, p: &mut Vec<Perfect>) -> bool {
        let s = self.list.pop_back().unwrap();
        rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
        if self.c_w[s as usize] == 0 {
            self.num_diff -= 1;
        }
        self.start += 1;
        self.list.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        self.stop += 1;
        if self.num_diff <= 1 {
            p.insert(
                0,
                Perfect {
                    first: self.start,
                    second: self.stop + 1,
                    score: 0,
                    len: 0,
                },
            );
            false
        } else {
            true
        }
    }

    fn shift_window(&mut self, t: u8, p: &mut Vec<Perfect>) -> bool {
        if self.list.len() >= self.max_size {
            if self.num_diff <= 1 {
                return self.shift_high(t, p);
            }
            let s = self.list.pop_back().unwrap();
            rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
            if self.c_w[s as usize] == 0 {
                self.num_diff -= 1;
            }
            if self.l == self.start {
                self.l += 1;
                rem_triplet_info(&mut self.r_v, &mut self.c_v, s);
            }
            self.start += 1;
        }
        self.list.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        add_triplet_info(&mut self.r_v, &mut self.c_v, t);
        if self.c_v[t as usize] > self.low_k {
            let mut off = self.list.len() - (self.l - self.start) - 1;
            loop {
                let x = self.list[off];
                rem_triplet_info(&mut self.r_v, &mut self.c_v, x);
                self.l += 1;
                let same = self.list[off] == t;
                off = off.wrapping_sub(1);
                if same {
                    break;
                }
            }
        }
        self.stop += 1;
        if self.list.len() >= self.max_size && self.num_diff <= 1 {
            p.clear();
            p.insert(
                0,
                Perfect {
                    first: self.start,
                    second: self.stop + 1,
                    score: 0,
                    len: 0,
                },
            );
            false
        } else {
            true
        }
    }

    fn needs_processing(&self, thresholds: &[u32]) -> bool {
        let count = self.stop - self.l;
        count < self.list.len() && 10 * self.r_w > thresholds[count]
    }

    fn find_perfect(&self, p: &mut Vec<Perfect>, thresholds: &[u32]) {
        let mut counts = self.c_v;
        let mut count = self.stop - self.l;
        let mut score = self.r_v;
        let mut perfect_iter = 0usize;
        let mut max_perfect_score: u32 = 0;
        let mut max_len: usize = 0;
        // pos = L - 1 (may be "negative" in C only if L == 0 and processing happens, which
        // needs_processing() rules out: count < size means L > start >= 0 ... keep wrapping).
        let mut pos = self.l.wrapping_sub(1);
        for idx in count..self.list.len() {
            let t = self.list[idx];
            let cnt = counts[t as usize];
            add_triplet_info(&mut score, &mut counts, t);
            if cnt > 0 && score * 10 > thresholds[count] {
                while perfect_iter < p.len() && pos <= p[perfect_iter].first {
                    let pi = &p[perfect_iter];
                    if max_perfect_score == 0
                        || (max_len as u64) * (pi.score as u64)
                            > (max_perfect_score as u64) * (pi.len as u64)
                    {
                        max_perfect_score = pi.score;
                        max_len = pi.len;
                    }
                    perfect_iter += 1;
                }
                if max_perfect_score == 0
                    || (score as u64) * (max_len as u64)
                        >= (max_perfect_score as u64) * (count as u64)
                {
                    max_perfect_score = score;
                    max_len = count;
                    p.insert(
                        perfect_iter,
                        Perfect {
                            first: pos,
                            second: self.stop + 1,
                            score: max_perfect_score,
                            len: count,
                        },
                    );
                }
            }
            count += 1;
            pos = pos.wrapping_sub(1);
        }
    }
}

/// `CSymDustMasker`.
pub struct SymDustMasker {
    level: u32,
    window: usize,
    linker: usize,
    low_k: u8,
    thresholds: Vec<u32>,
    rng: NcbiRandom,
}

impl SymDustMasker {
    pub fn new(level: u32, window: usize, linker: usize) -> Self {
        let level = if (2..=64).contains(&level) { level } else { 20 };
        let window = if (8..=64).contains(&window) {
            window
        } else {
            64
        };
        let linker = if (1..=32).contains(&linker) {
            linker
        } else {
            1
        };
        let mut thresholds = Vec::with_capacity(window - 2);
        thresholds.push(1);
        for i in 1..(window - 2) {
            thresholds.push(i as u32 * level);
        }
        SymDustMasker {
            level,
            window,
            linker,
            low_k: (level / 5) as u8,
            thresholds,
            rng: NcbiRandom::default_seeded(),
        }
    }

    /// `CIupac2Ncbi2na_converter::operator()` on an IUPAC letter (uppercase as CSeqVector gives it).
    #[inline]
    fn convert(&mut self, r: u8) -> u8 {
        match r {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => (self.rng.next() & 0x3) as u8,
            _ => 0,
        }
    }

    fn save_masked_regions(
        &self,
        res: &mut Vec<(usize, usize)>,
        p: &mut Vec<Perfect>,
        wstart: usize,
        start: usize,
    ) {
        if let Some(b) = p.last().copied() {
            if b.first < wstart {
                let b1 = (b.first + start, b.second + start);
                if let Some(last) = res.last_mut() {
                    let s = last.1;
                    if s + self.linker >= b1.0 {
                        last.1 = s.max(b1.1);
                    } else {
                        res.push(b1);
                    }
                } else {
                    res.push(b1);
                }
                while let Some(l) = p.last() {
                    if l.first < wstart {
                        p.pop();
                    } else {
                        break;
                    }
                }
            }
        }
    }

    /// `operator()(seq, start, stop)`: masked intervals, inclusive, in `seq` coordinates.
    /// `seq` holds IUPAC letters (uppercase).
    pub fn mask(&mut self, seq: &[u8]) -> Vec<(usize, usize)> {
        let mut res: Vec<(usize, usize)> = Vec::new();
        if seq.is_empty() {
            return res;
        }
        let mut start = 0usize;
        let stop = seq.len() - 1;
        let mut p: Vec<Perfect> = Vec::new();
        while stop > 2 + start {
            p.clear();
            let mut w = Triplets::new(self.window, self.low_k);
            let c1 = seq[start];
            let c2 = seq[start + 1];
            let mut t = (self.convert(c1) << 2) + self.convert(c2);
            let mut it = start + w.stop + 2;
            let mut done = false;
            while !done && it <= stop {
                self.save_masked_regions(&mut res, &mut p, w.start, start);
                t = ((t << 2) & TRIPLET_MASK) + (self.convert(seq[it]) & 0x3);
                it += 1;
                if w.shift_window(t, &mut p) {
                    if w.needs_processing(&self.thresholds) {
                        w.find_perfect(&mut p, &self.thresholds);
                    }
                } else {
                    while it <= stop {
                        self.save_masked_regions(&mut res, &mut p, w.start, start);
                        t = ((t << 2) & TRIPLET_MASK) + (self.convert(seq[it]) & 0x3);
                        if w.shift_window(t, &mut p) {
                            done = true;
                            break;
                        }
                        it += 1;
                    }
                }
            }
            {
                let mut wstart = w.start;
                while !p.is_empty() {
                    self.save_masked_regions(&mut res, &mut p, wstart, start);
                    wstart += 1;
                }
            }
            if w.start > 0 {
                start += w.start;
            } else {
                break;
            }
        }
        let _ = self.level;
        res
    }
}

/// `Blast_FindDustFilterLoc` for one query: dust intervals on the plus strand as
/// `(left, right)` inclusive, merged (`CSeq_loc::Merge` with `fMerge_All|fSort`: overlapping
/// or abutting intervals become one).
pub fn dust_query(iupac_upper: &[u8], level: u32, window: usize, linker: usize) -> Vec<(i32, i32)> {
    let mut m = SymDustMasker::new(level, window, linker);
    let raw = m.mask(iupac_upper);
    let mut out: Vec<(i32, i32)> = Vec::with_capacity(raw.len());
    for (a, b) in raw {
        let (a, b) = (a as i32, b as i32);
        if let Some(last) = out.last_mut()
            && a <= last.1 + 1
        {
            last.1 = last.1.max(b);
            continue;
        }
        out.push((a, b));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_a_is_masked_and_random_is_not() {
        let s = b"ACGTTGCAAGCTTAGCGATCGATTACGAGCTTAGGCATCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTTAGCGGATCCGATAGCTAGC";
        let m = dust_query(s, 20, 64, 1);
        assert_eq!(m.len(), 1, "{m:?}");
        // the A run spans 42..=81
        assert!(m[0].0 <= 44 && m[0].1 >= 79, "{m:?}");
        let r = b"ACGTTGCAAGCTTAGCGATCGATTACGAGCTTAGGCATCGACGTTAGCGGATCCGATAGCTAGCTTGACCAGT";
        assert!(dust_query(r, 20, 64, 1).is_empty());
    }

    #[test]
    fn short_sequences() {
        assert!(dust_query(b"", 20, 64, 1).is_empty());
        assert!(dust_query(b"ACG", 20, 64, 1).is_empty());
        assert!(dust_query(b"AAAAAAAAAA", 20, 64, 1).len() <= 1);
    }
}
