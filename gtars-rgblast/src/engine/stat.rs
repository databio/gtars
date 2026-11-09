//! Karlin-Altschul statistics for nucleotide searches.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_stat.c` (`BlastScoreBlkNuclMatrixCreate`,
//! `BlastScoreBlkMaxScoreSet`, `Blast_ResFreqStdComp`, `Blast_ResFreqString`,
//! `Blast_ScoreFreqNew`, `BlastScoreFreqCalc`, `NlmKarlinLambdaNR`, `Blast_KarlinLambdaNR`,
//! `BlastKarlinLtoH`, `BlastKarlinLHtoK`, `Blast_KarlinBlkUngappedCalc`,
//! `s_GetNuclValuesArray`, `s_SplitArrayOf8`, `s_AdjustGapParametersByGcd`,
//! `Blast_KarlinBlkNuclGappedCalc`, `BLAST_CheckRewardPenaltyScores`, `Blast_GetNuclAlphaBeta`,
//! `BLAST_ComputeLengthAdjustment`, `BLAST_KarlinStoE_simple`, `BLAST_KarlinEtoS_simple`,
//! `BLAST_GapDecayDivisor`), `ncbi_math.c` (`BLAST_Expm1`, `BLAST_Gcd`, `BLAST_Nint`,
//! `BLAST_Powi`) and `blast_setup.c` (`Blast_ScoreBlkMatrixInit`, `BlastSetup_ScoreBlkInit`,
//! `Blast_ScoreBlkKbpGappedCalc`, `BLAST_CalcEffLengths`).
//!
//! Every formula keeps the C operation order so the doubles come out bit-identical.

use crate::error::Error;

/// blastna -> ncbi4na bit set (blast_encoding.c).
pub const BLASTNA_TO_NCBI4NA: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;
const BLAST_KARLIN_K_ITER_MAX: i32 = 100;

// ---- ncbi_math.c ----

#[rustfmt::skip] // deeply nested Horner form makes rustfmt run for hours
pub fn expm1(x: f64) -> f64 {
    let absx = x.abs();
    if absx > 0.33 {
        return x.exp() - 1.;
    }
    if absx < 1.0e-16 {
        return x;
    }
    x * (1.
        + x * (1. / 2.
            + x * (1. / 6.
                + x * (1. / 24.
                    + x * (1. / 120.
                        + x * (1. / 720.
                            + x * (1. / 5040.
                                + x * (1. / 40320.
                                    + x * (1. / 362880.
                                        + x * (1. / 3628800.
                                            + x * (1. / 39916800.
                                                + x * (1. / 479001600. + x / 6227020800.))))))))))))
}

pub fn gcd(mut a: i32, mut b: i32) -> i32 {
    b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

pub fn nint(mut x: f64) -> i64 {
    x += if x >= 0. { 0.5 } else { -0.5 };
    x as i64
}

pub fn powi(mut x: f64, mut n: i32) -> f64 {
    if n == 0 {
        return 1.;
    }
    if x == 0. {
        if n < 0 {
            return f64::INFINITY;
        }
        return 0.;
    }
    if n < 0 {
        x = 1. / x;
        n = -n;
    }
    let mut y = 1.;
    while n > 0 {
        if n & 1 != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}

// ---- score matrix ----

/// `BlastScoreBlkNuclMatrixCreate`: the 16x16 blastna matrix. Row 15 (gap) is `INT4_MIN/2`.
pub fn nucl_matrix(reward: i32, penalty: i32) -> [[i32; 16]; 16] {
    let mut m = [[0i32; 16]; 16];
    let mut degeneracy = [0i32; 16];
    for d in degeneracy.iter_mut().take(4) {
        *d = 1;
    }
    for i in 4..16 {
        let mut degen = 0;
        for j in 0..4 {
            if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                degen += 1;
            }
        }
        degeneracy[i] = degen;
    }
    for i in 0..16 {
        for j in i..16 {
            if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                m[i][j] = nint(
                    (((degeneracy[j] - 1) * penalty + reward) as f64) / (degeneracy[j] as f64),
                ) as i32;
                if i != j {
                    m[j][i] = m[i][j];
                }
            } else {
                m[i][j] = penalty;
                m[j][i] = penalty;
            }
        }
    }
    // "The value of 15 is a gap, which is a sentinel between strands": both the row and the
    // column (the C sets matrix[15][*] and matrix[*][15])
    for v in m[15].iter_mut() {
        *v = i32::MIN / 2;
    }
    for row in m.iter_mut() {
        row[15] = i32::MIN / 2;
    }
    m
}

/// `BlastScoreBlkMaxScoreSet`: (loscore, hiscore) over the matrix, ignoring out-of-range cells.
pub fn max_score_set(m: &[[i32; 16]; 16]) -> (i32, i32) {
    let mut lo = BLAST_SCORE_MAX;
    let mut hi = BLAST_SCORE_MIN;
    for row in m {
        for &s in row {
            if s <= BLAST_SCORE_MIN || s >= BLAST_SCORE_MAX {
                continue;
            }
            lo = lo.min(s);
            hi = hi.max(s);
        }
    }
    (lo.max(BLAST_SCORE_MIN), hi.min(BLAST_SCORE_MAX))
}

/// `Blast_ResFreqStdComp` for nucleotides: 0.25 for A, C, G, T.
pub fn std_res_freq() -> [f64; 16] {
    let mut p = [0.0; 16];
    // nt_prob holds 25.00 each; Blast_ResFreqNormalize divides by the sum (100) and multiplies by 1.
    let sum = 25.00 * 4.0;
    for v in p.iter_mut().take(4) {
        *v = 25.00 / sum * 1.0;
    }
    p
}

/// `Blast_ResFreqString`: composition of one query context in blastna, with the ambiguous
/// residues 'N' (14) and '-' (15) zeroed (`BLAST_ScoreSetAmbigRes` in `Blast_ScoreBlkMatrixInit`).
pub fn res_freq_string(seq: &[u8]) -> [f64; 16] {
    let mut comp = [0i64; 16];
    for &b in seq {
        comp[(b & 0x0f) as usize] += 1;
    }
    comp[14] = 0;
    comp[15] = 0;
    let mut sum = 0.0f64;
    for &c in &comp {
        sum += c as f64;
    }
    let mut p = [0.0; 16];
    if sum == 0. {
        return p;
    }
    for i in 0..16 {
        p[i] = comp[i] as f64 / sum;
    }
    p
}

/// `Blast_ScoreFreq`: score probabilities, `sprob[score - score_min]`.
#[derive(Debug, Clone)]
pub struct ScoreFreq {
    pub score_min: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    sprob: Vec<f64>,
}

impl ScoreFreq {
    #[inline]
    pub fn prob(&self, score: i32) -> f64 {
        self.sprob[(score - self.score_min) as usize]
    }
}

fn score_chk(lo: i32, hi: i32) -> bool {
    !(lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
        && (hi - lo) <= (BLAST_SCORE_MAX - BLAST_SCORE_MIN)
}

/// `Blast_ScoreFreqNew` + `BlastScoreFreqCalc`.
pub fn score_freq_calc(
    m: &[[i32; 16]; 16],
    loscore: i32,
    hiscore: i32,
    rfp1: &[f64; 16],
    rfp2: &[f64; 16],
) -> Option<ScoreFreq> {
    if !score_chk(loscore, hiscore) {
        return None;
    }
    let range = (hiscore - loscore + 1) as usize;
    let mut sfp = ScoreFreq {
        score_min: loscore,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        sprob: vec![0.0; range],
    };
    for i in 0..16 {
        for j in 0..16 {
            let score = m[i][j];
            if score >= loscore {
                sfp.sprob[(score - loscore) as usize] += rfp1[i] * rfp2[j];
            }
        }
    }
    let mut score_sum = 0.;
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;
    for score in loscore..=hiscore {
        let p = sfp.prob(score);
        if p > 0. {
            score_sum += p;
            obs_max = score;
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }
    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;
    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in obs_min..=obs_max {
            let idx = (score - loscore) as usize;
            sfp.sprob[idx] /= score_sum;
            score_avg += score as f64 * sfp.sprob[idx];
        }
    }
    sfp.score_avg = score_avg;
    Some(sfp)
}

// ---- Karlin blocks ----

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct KarlinBlk {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub h: f64,
}

/// `NlmKarlinLambdaNR`.
#[allow(clippy::too_many_arguments)]
fn nlm_karlin_lambda_nr(
    probs: &dyn Fn(i32) -> f64,
    d: i32,
    low: i32,
    high: i32,
    lambda0: f64,
    tolx: f64,
    itmax: i32,
    max_newton: i32,
) -> f64 {
    let (mut a, mut b) = (0.0f64, 1.0f64);
    let mut f = 4.0f64;
    let mut is_newton = false;
    let x0 = (-lambda0).exp();
    let mut x = if 0. < x0 && x0 < 1. { x0 } else { 0.5 };
    for k in 0..itmax {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;
        let mut g = 0.0;
        f = probs(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + probs(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + probs(0) - 1.;
        let mut i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + probs(i);
            i += d;
        }
        if f > 0. {
            a = x;
        } else if f < 0. {
            b = x;
        } else {
            break;
        }
        if b - a < 2. * a * (1. - b) * tolx {
            x = (a + b) / 2.;
            break;
        }
        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0. {
            x = (a + b) / 2.;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1. - x) {
                    break;
                }
            }
        }
    }
    -x.ln() / d as f64
}

/// `Blast_KarlinLambdaNR`.
pub fn karlin_lambda_nr(sfp: &ScoreFreq, initial_guess: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    if sfp.score_avg >= 0. {
        return -1.0;
    }
    if !score_chk(low, high) {
        return -1.;
    }
    let mut d = -low;
    let mut i = 1;
    while i <= high - low && d > 1 {
        if sfp.prob(i + low) != 0.0 {
            d = gcd(d, i);
        }
        i += 1;
    }
    nlm_karlin_lambda_nr(
        &|s| sfp.prob(s),
        d,
        low,
        high,
        initial_guess,
        BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
        20,
        20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT,
    )
}

/// `BlastKarlinLtoH`.
pub fn karlin_lto_h(sfp: &ScoreFreq, lambda: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    if lambda < 0. {
        return -1.;
    }
    if !score_chk(low, high) {
        return -1.;
    }
    let etonlam = (-lambda).exp();
    let mut sum = low as f64 * sfp.prob(low);
    for score in (low + 1)..=high {
        sum = score as f64 * sfp.prob(score) + etonlam * sum;
    }
    let scale = powi(etonlam, high);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * high as f64 + sum.ln()).exp()
    }
}

/// `BlastKarlinLHtoK`.
pub fn karlin_lh_to_k(sfp: &ScoreFreq, mut lambda: f64, h: f64) -> f64 {
    if lambda <= 0. || h <= 0. {
        return -1.;
    }
    if sfp.score_avg >= 0.0 {
        return -1.;
    }
    let mut low = sfp.obs_min;
    let mut high = sfp.obs_max;
    let mut range = high - low;
    // probArrayStartLow[i] == sprob[low + i]
    let pl = |i: i32| sfp.prob(low + i);
    let mut divisor = -low;
    let mut i = 1;
    while i <= range && divisor > 1 {
        if pl(i) != 0.0 {
            divisor = gcd(divisor, i);
        }
        i += 1;
    }
    let orig_low = low;
    high /= divisor;
    low /= divisor;
    lambda *= divisor as f64;
    range = high - low;
    let mut first_term_closed_form = h / lambda;
    let exp_minus_lambda = (-lambda).exp();
    if low == -1 && high == 1 {
        let pl0 = sfp.prob(low * divisor);
        let ph0 = sfp.prob(high * divisor);
        return (pl0 - ph0) * (pl0 - ph0) / pl0;
    }
    if low == -1 || high == 1 {
        if high != 1 {
            let score_avg = sfp.score_avg / divisor as f64;
            first_term_closed_form = (score_avg * score_avg) / first_term_closed_form;
        }
        return first_term_closed_form * (1.0 - exp_minus_lambda);
    }
    let sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    let iterlimit = BLAST_KARLIN_K_ITER_MAX;
    let mut asp = vec![0.0f64; (iterlimit * range + 1) as usize];
    let mut outer_sum = 0.;
    let mut low_alignment_score = 0i32;
    let mut high_alignment_score = 0i32;
    asp[0] = 1.;
    let mut inner_sum = 1.;
    let mut iter_counter = 0;
    // probArrayStartLow with the *original* low: the C pointer is set before the division.
    let pl_orig = |i: i32| sfp.prob(orig_low + i);
    while iter_counter < iterlimit && inner_sum > sumlimit {
        let mut first = range;
        let mut last = range;
        low_alignment_score += low;
        high_alignment_score += high;
        let mut p = high_alignment_score - low_alignment_score; // index of ptrP
        while p >= 0 {
            let mut ptr1 = p - first;
            let ptr1e = p - last;
            let mut ptr2 = first;
            inner_sum = 0.;
            while ptr1 >= ptr1e {
                inner_sum += asp[ptr1 as usize] * pl_orig(ptr2);
                ptr1 -= 1;
                ptr2 += 1;
            }
            asp[p as usize] = inner_sum;
            if first != 0 {
                first -= 1;
            }
            if p <= range {
                last -= 1;
            }
            p -= 1;
        }
        // ptrP is now at index -1; "*++ptrP" reads index 0.
        let mut p = 0usize;
        inner_sum = asp[p];
        let mut i = low_alignment_score + 1;
        while i < 0 {
            p += 1;
            inner_sum = asp[p] + inner_sum * exp_minus_lambda;
            i += 1;
        }
        inner_sum *= exp_minus_lambda;
        while i <= high_alignment_score {
            p += 1;
            inner_sum += asp[p];
            i += 1;
        }
        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }
    -(-2.0 * outer_sum).exp() / (first_term_closed_form * expm1(-lambda))
}

/// `Blast_KarlinBlkUngappedCalc`; `None` when any parameter fails.
pub fn karlin_blk_ungapped_calc(sfp: &ScoreFreq) -> Option<KarlinBlk> {
    let lambda = karlin_lambda_nr(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
    if lambda < 0. {
        return None;
    }
    let h = karlin_lto_h(sfp, lambda);
    if h < 0. {
        return None;
    }
    let k = karlin_lh_to_k(sfp, lambda, h);
    if k < 0. {
        return None;
    }
    Some(KarlinBlk {
        lambda,
        k,
        log_k: k.ln(),
        h,
    })
}

// ---- gapped nucleotide tables ----

type ArrayOf8 = [f64; 8];

// Only the tables of the task defaults: 1/-3 (blastn-short), 1/-2 (megablast), 2/-3 (blastn).
const BLASTN_VALUES_1_3: &[ArrayOf8] = &[
    [0., 0., 1.374, 0.711, 1.31, 1.05, 0., 100.],
    [2., 2., 1.37, 0.70, 1.2, 1.1, 0., 99.],
    [1., 2., 1.35, 0.64, 1.1, 1.2, -1., 98.],
    [0., 2., 1.25, 0.42, 0.83, 1.5, -2., 91.],
    [2., 1., 1.34, 0.60, 1.1, 1.2, -1., 97.],
    [1., 1., 1.21, 0.34, 0.71, 1.7, -2., 88.],
];
const BLASTN_VALUES_1_2: &[ArrayOf8] = &[
    [0., 0., 1.28, 0.46, 0.85, 1.5, -2., 96.],
    [2., 2., 1.33, 0.62, 1.1, 1.2, 0., 99.],
    [1., 2., 1.30, 0.52, 0.93, 1.4, -2., 97.],
    [0., 2., 1.19, 0.34, 0.66, 1.8, -3., 89.],
    [3., 1., 1.32, 0.57, 1.0, 1.3, -1., 99.],
    [2., 1., 1.29, 0.49, 0.92, 1.4, -1., 96.],
    [1., 1., 1.14, 0.26, 0.52, 2.2, -5., 85.],
];
const BLASTN_VALUES_2_3: &[ArrayOf8] = &[
    [0., 0., 0.55, 0.21, 0.46, 1.2, -5., 87.],
    [4., 4., 0.63, 0.42, 0.84, 0.75, -2., 99.],
    [2., 4., 0.615, 0.37, 0.72, 0.85, -3., 97.],
    [0., 4., 0.55, 0.21, 0.46, 1.2, -5., 87.],
    [3., 3., 0.615, 0.37, 0.68, 0.9, -3., 97.],
    [6., 2., 0.63, 0.42, 0.84, 0.75, -2., 99.],
    [5., 2., 0.625, 0.41, 0.78, 0.8, -2., 99.],
    [4., 2., 0.61, 0.35, 0.68, 0.9, -3., 96.],
    [2., 2., 0.515, 0.14, 0.33, 1.55, -9., 81.],
];

/// `s_GetNuclValuesArray`'s result: the affine rows (`normal`), the linear-gap row when the
/// table has one, the "anything more stringent" limits and the `round_down` flag.
#[derive(Debug, Clone)]
pub struct NuclValues {
    pub normal: Vec<ArrayOf8>,
    pub linear: Option<ArrayOf8>,
    pub gap_open_max: i32,
    pub gap_extend_max: i32,
    pub round_down: bool,
}

/// `s_GetNuclValuesArray` for the task defaults' reward/penalty pairs (all coprime, so
/// `s_AdjustGapParametersByGcd` never applies); `Err` carries blastn's "Substitution scores
/// R and P are not supported".
pub fn nucl_values_array(reward: i32, penalty: i32) -> Result<NuclValues, Error> {
    let (table, gap_open_max, gap_extend_max, round_down): (&[ArrayOf8], i32, i32, bool) =
        match (reward, penalty) {
            (1, -3) => (BLASTN_VALUES_1_3, 2, 2, false),
            (1, -2) => (BLASTN_VALUES_1_2, 2, 2, false),
            (2, -3) => (BLASTN_VALUES_2_3, 6, 4, true),
            _ => {
                return Err(Error::InvalidOptions(format!(
                    "Substitution scores {reward} and {penalty} are not supported"
                )));
            }
        };
    // s_SplitArrayOf8
    let (normal, linear) = if table[0][0] == 0. && table[0][1] == 0. {
        (table[1..].to_vec(), Some(table[0]))
    } else {
        (table.to_vec(), None)
    };
    Ok(NuclValues {
        normal,
        linear,
        gap_open_max,
        gap_extend_max,
        round_down,
    })
}

/// `Blast_KarlinBlkNuclGappedCalc`: the gapped block for the given gap costs, plus the
/// table's `round_down` flag. `Err` carries blastn's message text.
pub fn karlin_blk_nucl_gapped_calc(
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    kbp_ungap: &KarlinBlk,
) -> Result<(KarlinBlk, bool), Error> {
    let v = nucl_values_array(reward, penalty)?;
    if gap_open == 0 && gap_extend == 0 {
        if let Some(l) = v.linear {
            return Ok((
                KarlinBlk {
                    lambda: l[2],
                    k: l[3],
                    log_k: l[3].ln(),
                    h: l[4],
                },
                v.round_down,
            ));
        }
    }
    for row in &v.normal {
        if row[0] == gap_open as f64 && row[1] == gap_extend as f64 {
            return Ok((
                KarlinBlk {
                    lambda: row[2],
                    k: row[3],
                    log_k: row[3].ln(),
                    h: row[4],
                },
                v.round_down,
            ));
        }
    }
    if gap_open >= v.gap_open_max && gap_extend >= v.gap_extend_max {
        return Ok((*kbp_ungap, v.round_down));
    }
    let mut msg = format!(
        "Gap existence and extension values {gap_open} and {gap_extend} are not supported for substitution scores {reward} and {penalty}\n"
    );
    for row in &v.normal {
        msg.push_str(&format!(
            "{} and {} are supported existence and extension values\n",
            row[0] as i64, row[1] as i64
        ));
    }
    msg.push_str(&format!(
        "{} and {} are supported existence and extension values\n",
        v.gap_open_max, v.gap_extend_max
    ));
    msg.push_str(&format!(
        "Any values more stringent than {} and {} are supported\n",
        v.gap_open_max, v.gap_extend_max
    ));
    Err(Error::InvalidOptions(msg))
}

fn ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        -2.
    } else {
        0.
    }
}

/// `Blast_GetNuclAlphaBeta`.
pub fn nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    kbp_ungap: &KarlinBlk,
    gapped: bool,
) -> Result<(f64, f64), Error> {
    let v = nucl_values_array(reward, penalty)?;
    if gapped && !v.normal.is_empty() {
        if gap_open == 0 && gap_extend == 0 {
            if let Some(l) = v.linear {
                return Ok((l[5], l[6]));
            }
        } else {
            for row in &v.normal {
                if row[0] == gap_open as f64 && row[1] == gap_extend as f64 {
                    return Ok((row[5], row[6]));
                }
            }
        }
    }
    Ok((
        kbp_ungap.lambda / kbp_ungap.h,
        ungapped_beta(reward, penalty),
    ))
}

/// `BLAST_ComputeLengthAdjustment`; returns (length_adjustment, converged).
pub fn compute_length_adjustment(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
) -> (i32, bool) {
    let max_iterations = 20;
    let m = query_length as f64;
    let n = db_length as f64;
    let nn = db_num_seqs as f64;
    let mut ell_min = 0.0f64;
    let ell_max;
    let mut converged = false;
    let mut ell_next = 0.0f64;
    {
        let a = nn;
        let mb = m * nn + n;
        let c = n * m - m.max(n) / k;
        if c < 0. {
            return (0, false);
        } else {
            ell_max = 2. * c / (mb + (mb * mb - 4. * a * c).sqrt());
        }
    }
    let mut ell_max = ell_max;
    for i in 1..=max_iterations {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;
        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }
        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.
            };
        }
    }
    let mut length_adjustment;
    if converged {
        length_adjustment = ell_min as i32;
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                length_adjustment = ell as i32;
            }
        }
    } else {
        length_adjustment = ell_min as i32;
    }
    (length_adjustment, converged)
}

/// `BLAST_KarlinStoE_simple`.
#[inline]
pub fn karlin_stoe_simple(s: i32, kbp: &KarlinBlk, searchsp: i64) -> f64 {
    if kbp.lambda < 0. || kbp.k < 0. || kbp.h < 0. {
        return -1.;
    }
    (searchsp as f64) * ((-kbp.lambda * s as f64) + kbp.log_k).exp()
}

/// `BLAST_KarlinEtoS_simple`.
pub fn karlin_etos_simple(e: f64, kbp: &KarlinBlk, searchsp: i64) -> i32 {
    if kbp.lambda < 0. || kbp.k < 0. || kbp.h < 0. {
        return BLAST_SCORE_MIN;
    }
    if e == 0.0 {
        return BLAST_SCORE_MAX;
    }
    let s = ((kbp.k * searchsp as f64 / e).ln() / kbp.lambda).ceil();
    (BLAST_SCORE_MIN as f64).max(s) as i32
}

/// `BLAST_GapDecayDivisor`.
pub fn gap_decay_divisor(decayrate: f64, nsegs: u32) -> f64 {
    (1. - decayrate) * powi(decayrate, nsegs as i32 - 1)
}

// ---- the score block ----

/// `BlastScoreBlk` for blastn: the matrix, score range and per-context Karlin blocks.
#[derive(Debug, Clone)]
pub struct ScoreBlk {
    pub reward: i32,
    pub penalty: i32,
    pub matrix: [[i32; 16]; 16],
    pub loscore: i32,
    pub hiscore: i32,
    /// `kbp_std[context]`: None for an invalid context.
    pub kbp_std: Vec<Option<KarlinBlk>>,
    /// `kbp_gap_std[context]`.
    pub kbp_gap_std: Vec<Option<KarlinBlk>>,
    pub round_down: bool,
}

impl ScoreBlk {
    /// `BlastSetup_ScoreBlkInit` for blastn: `Blast_ScoreBlkMatrixInit`,
    /// `Blast_ScoreBlkKbpUngappedCalc` over every valid context (composition from the query
    /// letters), then `Blast_ScoreBlkKbpGappedCalc`. `contexts` gives each context's blastna
    /// letters (None = invalid). Returns the block and the contexts that failed (marked
    /// invalid, as the C does).
    pub fn new(
        reward: i32,
        penalty: i32,
        gap_open: i32,
        gap_extend: i32,
        gapped: bool,
        contexts: &[Option<&[u8]>],
    ) -> Result<(ScoreBlk, Vec<bool>), Error> {
        let matrix = nucl_matrix(reward, penalty);
        let (loscore, hiscore) = max_score_set(&matrix);
        let stdrfp = std_res_freq();
        let mut kbp_std = Vec::with_capacity(contexts.len());
        let mut valid = Vec::with_capacity(contexts.len());
        let mut any_valid = false;
        for c in contexts {
            match c {
                Some(seq) => {
                    let rfp = res_freq_string(seq);
                    let kbp = score_freq_calc(&matrix, loscore, hiscore, &rfp, &stdrfp)
                        .and_then(|sfp| karlin_blk_ungapped_calc(&sfp));
                    valid.push(kbp.is_some());
                    any_valid |= kbp.is_some();
                    kbp_std.push(kbp);
                }
                None => {
                    valid.push(false);
                    kbp_std.push(None);
                }
            }
        }
        if !any_valid && !contexts.is_empty() {
            // Blast_ScoreBlkKbpUngappedCalc returns 1; BLAST_MainSetUp reports it.
            return Err(Error::InvalidOptions(
                "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options".into(),
            ));
        }
        let mut kbp_gap_std = vec![None; contexts.len()];
        let mut round_down = false;
        if gapped {
            for (i, kbp) in kbp_std.iter().enumerate() {
                if let Some(k) = kbp {
                    let (g, rd) =
                        karlin_blk_nucl_gapped_calc(gap_open, gap_extend, reward, penalty, k)?;
                    kbp_gap_std[i] = Some(g);
                    round_down = rd;
                }
            }
        }
        Ok((
            ScoreBlk {
                reward,
                penalty,
                matrix,
                loscore,
                hiscore,
                kbp_std,
                kbp_gap_std,
                round_down,
            },
            valid,
        ))
    }

    /// `BLAST_CalcEffLengths` for one context: (length_adjustment, eff_searchsp).
    /// `db_length` is already `-dbsize` or the real length; `db_num_seqs` the real count.
    pub fn eff_lengths(
        &self,
        context: usize,
        gap_open: i32,
        gap_extend: i32,
        gapped: bool,
        query_length: i32,
        db_length: i64,
        db_num_seqs: i32,
    ) -> Result<(i32, i64), Error> {
        let kbp_ungap = match &self.kbp_std[context] {
            Some(k) => k,
            None => return Ok((0, 0)),
        };
        let kbp = if gapped {
            self.kbp_gap_std[context].as_ref().unwrap_or(kbp_ungap)
        } else {
            kbp_ungap
        };
        if query_length <= 0 {
            return Ok((0, 0));
        }
        let (alpha, beta) = nucl_alpha_beta(
            self.reward,
            self.penalty,
            gap_open,
            gap_extend,
            kbp_ungap,
            gapped,
        )?;
        let (length_adjustment, _) = compute_length_adjustment(
            kbp.k,
            kbp.log_k,
            alpha / kbp.lambda,
            beta,
            query_length,
            db_length,
            db_num_seqs,
        );
        let mut effective_db_length = db_length - (db_num_seqs as i64 * length_adjustment as i64);
        if effective_db_length <= 0 {
            effective_db_length = 1;
        }
        let eff_searchsp = effective_db_length * (query_length - length_adjustment) as i64;
        Ok((length_adjustment, eff_searchsp))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_1_2() {
        let m = nucl_matrix(1, -2);
        assert_eq!(m[0][0], 1);
        assert_eq!(m[0][1], -2);
        // A vs R (A|G): degeneracy 2 -> nint((1*-2 + 1)/2) = nint(-0.5) = -1
        assert_eq!(m[0][4], -1);
        assert_eq!(m[4][0], -1);
        // A vs N: degeneracy 4 -> nint((3*-2+1)/4) = nint(-1.25) = -1
        assert_eq!(m[0][14], -1);
        assert_eq!(m[15][3], i32::MIN / 2);
        assert_eq!(max_score_set(&m), (-2, 1));
    }

    #[test]
    fn ungapped_karlin_1_2_matches_blastn_footer() {
        // blastn 2.17 -reward 1 -penalty -2 prints "Lambda K H: 1.28 0.460 0.850" gapped from
        // the table and ungapped "1.33 0.621 1.12" for an ACGT-balanced query.
        let m = nucl_matrix(1, -2);
        let (lo, hi) = max_score_set(&m);
        let q: Vec<u8> = (0..400u32).map(|i| (i % 4) as u8).collect();
        let sfp = score_freq_calc(&m, lo, hi, &res_freq_string(&q), &std_res_freq()).unwrap();
        let k = karlin_blk_ungapped_calc(&sfp).unwrap();
        assert!((k.lambda - 1.33).abs() < 0.005, "{k:?}");
        assert!((k.k - 0.621).abs() < 0.0005, "{k:?}");
        assert!((k.h - 1.12).abs() < 0.005, "{k:?}");
        let (g, rd) = karlin_blk_nucl_gapped_calc(0, 0, 1, -2, &k).unwrap();
        assert_eq!((g.lambda, g.k, g.h, rd), (1.28, 0.46, 0.85, false));
    }

    #[test]
    fn ungapped_karlin_2_3_general_path() {
        // -task blastn defaults: reward 2 penalty -3 gap 5/2; footer "Lambda K H 0.625 0.410 0.780"
        // gapped and "0.634 0.408 0.912" ungapped.
        let m = nucl_matrix(2, -3);
        let (lo, hi) = max_score_set(&m);
        let q: Vec<u8> = (0..400u32).map(|i| (i % 4) as u8).collect();
        let sfp = score_freq_calc(&m, lo, hi, &res_freq_string(&q), &std_res_freq()).unwrap();
        let k = karlin_blk_ungapped_calc(&sfp).unwrap();
        assert!((k.lambda - 0.634).abs() < 0.0015, "{k:?}");
        assert!((k.k - 0.408).abs() < 0.0015, "{k:?}");
        assert!((k.h - 0.912).abs() < 0.0015, "{k:?}");
        let (g, rd) = karlin_blk_nucl_gapped_calc(5, 2, 2, -3, &k).unwrap();
        assert_eq!((g.lambda, g.k, g.h, rd), (0.625, 0.41, 0.78, true));
    }

    #[test]
    fn unsupported_scores() {
        assert!(nucl_values_array(1, -7).is_err());
        let kbp = KarlinBlk {
            lambda: 1.,
            k: 1.,
            log_k: 0.,
            h: 1.,
        };
        let e = karlin_blk_nucl_gapped_calc(1, 0, 1, -2, &kbp).unwrap_err();
        assert!(
            e.to_string()
                .starts_with("Gap existence and extension values 1 and 0")
        );
        // more stringent than the table's maximum: the ungapped block is used
        let (g, _) = karlin_blk_nucl_gapped_calc(7, 3, 1, -2, &kbp).unwrap();
        assert_eq!(g, kbp);
    }

    #[test]
    fn gcd_and_math() {
        assert_eq!(gcd(6, -4), 2);
        assert_eq!(gcd(1, -2), 1);
        assert_eq!(nint(-0.5), -1);
        assert_eq!(nint(0.5), 1);
        assert_eq!(powi(2.0, 10), 1024.0);
        assert!((expm1(0.1) - 0.1f64.exp_m1()).abs() < 1e-15);
    }
}
