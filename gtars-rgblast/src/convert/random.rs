//! makeblastdb's random fill for ambiguous bases: NCBI `CRandom` (util/random_gen.cpp),
//! a lagged Fibonacci generator with lags 33 and 13. WriteDB seeds it with the
//! sequence length and draws once per ambiguous base, in order.

pub struct NcbiRandom {
    s: [u32; 33],
    rj: i32,
    rk: i32,
}

impl NcbiRandom {
    pub fn new(seed: u32) -> Self {
        let mut s = [0u32; 33];
        s[0] = seed;
        for i in 1..33 {
            s[i] = 1103515245u32.wrapping_mul(s[i - 1]).wrapping_add(12345);
        }
        let mut r = NcbiRandom { s, rj: 12, rk: 32 };
        for _ in 0..330 {
            r.next();
        }
        r
    }

    /// `CRandom::CRandom()` (default-constructed, `Reset()`): the fixed state table
    /// `sm_State` with `m_RJ = 12`, `m_RK = 32`. `CSymDustMasker` draws from one of these for
    /// every `N` it meets.
    pub fn default_seeded() -> Self {
        NcbiRandom {
            s: [
                0xd53f1852, 0xdfc78b83, 0x4f256096, 0xe643df7, 0x82c359bf, 0xc7794dfa, 0xd5e9ffaa,
                0x2c8cb64a, 0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589, 0xa5853646, 0xe71576e2,
                0xdae30df, 0xb09ce711, 0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8, 0xd788d765,
                0x67498087, 0x9d7aba26, 0x261351d4, 0x411ee7ea, 0x393a263, 0x2c5a5835, 0xc115fcd8,
                0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98, 0x6e37bd55,
            ],
            rj: 12,
            rk: 32,
        }
    }

    /// `CRandom::GetRand()`: the 32-bit LFG output with the low bit dropped.
    #[inline(always)]
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> u32 {
        let r = self.s[self.rk as usize].wrapping_add(self.s[self.rj as usize]);
        self.s[self.rk as usize] = r;
        self.rj -= 1;
        self.rk -= 1;
        if self.rk < 0 {
            self.rk = 32;
        } else if self.rj < 0 {
            self.rj = 32;
        }
        r >> 1
    }
}
