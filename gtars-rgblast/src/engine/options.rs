//! Search options (task defaults) and the derived parameters (cutoffs).
//!
//! Source: BLAST+ 2.17.0 `algo/blast/api/blast_nucl_options.cpp`
//! (`SetTraditionalMegablastDefaults`, `SetTraditionalBlastnDefaults` and the `Set*Defaults`
//! they call), `algo/blast/core/blast_options.c` (the `*OptionsNew` constructors for the fields
//! the handle leaves alone), `blast_options.h` constants, and `blast_parameters.c`
//! (`BlastInitialWordParametersNew/Update`, `BlastExtensionParametersNew`,
//! `BlastHitSavingParametersNew/Update`, `BlastScoringParametersNew`,
//! `BlastEffectiveLengthsParametersNew`) with `blast_stat.c::BLAST_Cutoffs` and
//! `blast_setup.c::BLAST_GapAlignSetUp` (the order things are derived in).

use super::stat::{self, KarlinBlk, ScoreBlk};
use crate::error::Error;
use serde::{Deserialize, Serialize};

pub const NCBIMATH_LN2: f64 = 0.693_147_180_559_945_309_417_232_121_458_18;
pub const MAX_DBSEQ_LEN: i32 = 5_000_000;
pub const DBSEQ_CHUNK_OVERLAP: i32 = 100;
const CUTOFF_E_BLASTN: f64 = 0.05;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Task {
    #[default]
    Megablast,
    Blastn,
    BlastnShort,
}

impl Task {
    pub fn parse(s: &str) -> Option<Task> {
        Some(match s {
            "megablast" => Task::Megablast,
            "blastn" => Task::Blastn,
            "blastn-short" => Task::BlastnShort,
            _ => return None,
        })
    }
    pub fn name(&self) -> &'static str {
        match self {
            Task::Megablast => "megablast",
            Task::Blastn => "blastn",
            Task::BlastnShort => "blastn-short",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Strand {
    #[default]
    Both,
    Plus,
    Minus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LookupType {
    /// `eMBLookupTable`
    Mb,
    /// `eNaLookupTable` (the SmallNa variant is chosen at build time)
    Na,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PrelimGapExt {
    GreedyScoreOnly,
    DynProgScoreOnly,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TbackExt {
    GreedyTbck,
    DynProgTbck,
}

/// `-dust`: `no`, `yes` (20 64 1) or explicit level/window/linker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct DustSpec {
    pub level: i32,
    pub window: i32,
    pub linker: i32,
}

impl Default for DustSpec {
    fn default() -> Self {
        DustSpec {
            level: 20,
            window: 64,
            linker: 1,
        }
    }
}

/// The nucleotide fields of `BlastLookupTableOptions`, `QuerySetUpOptions`,
/// `BlastInitialWordOptions`, `BlastExtensionOptions`, `BlastScoringOptions`,
/// `BlastHitSavingOptions` and `BlastEffectiveLengthsOptions`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlastOptions {
    pub task: Task,
    // lookup table
    pub word_size: i32,
    pub lut_type: LookupType,
    // query setup
    pub dust: Option<DustSpec>,
    pub mask_at_hash: bool,
    pub strand: Strand,
    // initial word
    pub x_dropoff: f64,
    pub gap_trigger: f64,
    // extension
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub prelim_gap_ext: PrelimGapExt,
    pub tback_ext: TbackExt,
    // scoring
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gapped_calculation: bool,
    // hit saving
    pub hitlist_size: i32,
    pub expect_value: f64,
    pub percent_identity: f64,
    /// `-max_hsps` (`SetMaxHspsPerSubject`); 0 = unlimited.
    pub max_hsps_per_subject: i32,
    pub min_diag_separation: i32,
    pub mask_level: i32,
    pub cutoff_score: i32,
}

impl BlastOptions {
    /// `CBlastOptionsFactory::CreateTask("megablast")`.
    pub fn megablast() -> Self {
        BlastOptions {
            task: Task::Megablast,
            word_size: 28,
            lut_type: LookupType::Mb,
            dust: Some(DustSpec::default()),
            mask_at_hash: true,
            strand: Strand::Both,
            // `CBlastNucleotideOptionsHandle::SetMBInitialWordOptionsDefaults` sets only the
            // window size, so the megablast task runs with `x_dropoff == 0`;
            // `BlastInitialWordParametersNew` then gets `x_dropoff_init == 0` and
            // `BlastInitialWordParametersUpdate` uses the ungapped cutoff score as the x-drop
            // (`blast_parameters.c` 380-381). The blastn handle sets 20.
            x_dropoff: 0.0,
            gap_trigger: 27.0,
            gap_x_dropoff: 25.0,
            gap_x_dropoff_final: 100.0,
            prelim_gap_ext: PrelimGapExt::GreedyScoreOnly,
            tback_ext: TbackExt::GreedyTbck,
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
            gapped_calculation: true,
            hitlist_size: 500,
            expect_value: 10.0,
            percent_identity: 0.0,
            max_hsps_per_subject: 0,
            min_diag_separation: 6,
            mask_level: 101,
            cutoff_score: 0,
        }
    }

    /// `CBlastOptionsFactory::CreateTask("blastn")`.
    pub fn blastn() -> Self {
        BlastOptions {
            task: Task::Blastn,
            word_size: 11,
            lut_type: LookupType::Na,
            x_dropoff: 20.0,
            gap_x_dropoff: 30.0,
            gap_x_dropoff_final: 100.0,
            prelim_gap_ext: PrelimGapExt::DynProgScoreOnly,
            tback_ext: TbackExt::DynProgTbck,
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            min_diag_separation: 50,
            ..Self::megablast()
        }
    }

    /// `CBlastOptionsFactory::CreateTask("blastn-short")`: blastn with word 7, no dust,
    /// evalue 1000, reward/penalty 1/-3, gap 5/2.
    pub fn blastn_short() -> Self {
        BlastOptions {
            task: Task::BlastnShort,
            word_size: 7,
            dust: None,
            expect_value: 1000.0,
            reward: 1,
            penalty: -3,
            ..Self::blastn()
        }
    }

    pub fn for_task(task: Task) -> Self {
        match task {
            Task::Megablast => Self::megablast(),
            Task::Blastn => Self::blastn(),
            Task::BlastnShort => Self::blastn_short(),
        }
    }

    /// The checks the remaining user flags can fail (scores and gap costs are the task's).
    pub fn validate(&self) -> Result<(), Error> {
        if self.word_size < 4 {
            return Err(Error::InvalidOptions(
                "Word-size must be 4 or greater for nucleotide comparison".into(),
            ));
        }
        if self.hitlist_size < 1 {
            return Err(Error::InvalidOptions("No hits are being saved".into()));
        }
        if self.expect_value <= 0.0 {
            return Err(Error::InvalidOptions(
                "expect value must be greater than zero".into(),
            ));
        }
        Ok(())
    }
}

/// `BlastUngappedCutoffs` (per context).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct UngappedCutoffs {
    pub x_dropoff_init: i32,
    pub x_dropoff: i32,
    pub cutoff_score: i32,
    pub reduced_nucl_cutoff_score: i32,
}

/// `BlastGappedCutoffs` (per context).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct GappedCutoffs {
    pub cutoff_score: i32,
    pub cutoff_score_max: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiagContainer {
    DiagArray,
    DiagHash,
}

/// Per-context statistics the search uses (`BlastContextInfo`'s statistical fields).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct ContextStats {
    pub query_index: u32,
    pub frame: i32,
    pub query_length: i32,
    pub is_valid: bool,
    pub length_adjustment: i32,
    pub eff_searchsp: i64,
    pub kbp_ungapped: Option<KarlinBlk>,
    pub kbp_gapped: Option<KarlinBlk>,
}

/// Everything `BLAST_GapAlignSetUp` and `BlastInitialWordParametersNew` derive.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Parameters {
    // BlastEffectiveLengthsParameters
    /// What `BLAST_CalcEffLengths` used for `db_length` / `db_num_seqs`.
    pub db_length: i64,
    pub db_num_seqs: i32,
    // BlastScoringParameters
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    // BlastExtensionParameters
    pub gap_x_dropoff: i32,
    pub gap_x_dropoff_final: i32,
    // BlastHitSavingParameters
    pub hit_cutoffs: Vec<GappedCutoffs>,
    pub hit_cutoff_score_min: i32,
    pub prelim_evalue: f64,
    pub mask_level: i32,
    // BlastInitialWordParameters
    pub word_cutoffs: Vec<UngappedCutoffs>,
    pub word_cutoff_score_min: i32,
    pub x_dropoff_max: i32,
    pub container_type: DiagContainer,
    pub nucl_score_table: Vec<i32>,
    pub ungapped_extension: bool,
}

/// `BLAST_Cutoffs`: returns (S, E) after adjustment.
pub fn blast_cutoffs(
    s_in: i32,
    e_in: f64,
    kbp: &KarlinBlk,
    searchsp: i64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> (i32, f64) {
    let mut s = s_in;
    let mut e = e_in;
    let esave = e;
    let mut s_changed = false;
    let mut es = 1;
    if e > 0. {
        if dodecay && gap_decay_rate > 0. && gap_decay_rate < 1. {
            e *= stat::gap_decay_divisor(gap_decay_rate, 1);
        }
        es = stat::karlin_etos_simple(e, kbp, searchsp);
    }
    if es > s {
        s_changed = true;
        s = es;
    }
    if esave <= 0. || !s_changed {
        e = stat::karlin_stoe_simple(s, kbp, searchsp);
        if dodecay && gap_decay_rate > 0. && gap_decay_rate < 1. {
            e /= stat::gap_decay_divisor(gap_decay_rate, 1);
        }
    } else {
        e = esave;
    }
    (s, e)
}

/// `BlastHitSavingParametersNew` + `Update` for gapped blastn without sum statistics.
fn hit_saving_params(
    o: &BlastOptions,
    sbp: &ScoreBlk,
    ctx: &[ContextStats],
) -> (Vec<GappedCutoffs>, i32) {
    let mut cutoffs = vec![GappedCutoffs::default(); ctx.len()];
    let mut cutoff_min = i32::MAX;
    if o.cutoff_score > 0 {
        for c in cutoffs.iter_mut() {
            c.cutoff_score = o.cutoff_score;
            c.cutoff_score_max = o.cutoff_score;
        }
        return (cutoffs, o.cutoff_score);
    }
    for (i, c) in ctx.iter().enumerate() {
        if !c.is_valid {
            cutoffs[i].cutoff_score = i32::MAX;
            continue;
        }
        let kbp = if o.gapped_calculation {
            sbp.kbp_gap_std[i].as_ref()
        } else {
            sbp.kbp_std[i].as_ref()
        };
        let Some(kbp) = kbp else {
            cutoffs[i].cutoff_score = i32::MAX;
            continue;
        };
        let (new_cutoff, _) = blast_cutoffs(1, o.expect_value, kbp, c.eff_searchsp, false, 0.0);
        cutoffs[i].cutoff_score = new_cutoff;
        cutoffs[i].cutoff_score_max = new_cutoff;
    }
    for (i, c) in ctx.iter().enumerate() {
        if c.is_valid {
            cutoff_min = cutoff_min.min(cutoffs[i].cutoff_score);
        }
    }
    (cutoffs, cutoff_min)
}

impl Parameters {
    /// `BLAST_GapAlignSetUp` then `BlastInitialWordParametersNew`, for blastn.
    /// `ctx` holds `query_index`, `frame`, `query_length` and `is_valid` per context; this
    /// fills `length_adjustment`, `eff_searchsp` and the Karlin blocks.
    /// `db_length`/`num_seqs` are the `BLAST_GetSubjectTotals` pair (STATS already applied);
    /// `avg_subject_length` is `BlastSeqSrcGetAvgSeqLen`.
    pub fn new(
        o: &BlastOptions,
        sbp: &ScoreBlk,
        ctx: &mut [ContextStats],
        total_length: i64,
        num_seqs: i32,
        avg_subject_length: i32,
    ) -> Result<Parameters, Error> {
        // BlastEffectiveLengthsParametersNew + BLAST_CalcEffLengths (no -dbsize, -searchsp
        // or dbseq_num override: the totals are used as given)
        let db_length = total_length;
        let db_num_seqs = num_seqs;
        for (i, c) in ctx.iter_mut().enumerate() {
            c.kbp_ungapped = sbp.kbp_std[i];
            c.kbp_gapped = sbp.kbp_gap_std[i];
            if !c.is_valid || c.query_length <= 0 {
                c.length_adjustment = 0;
                c.eff_searchsp = 0;
                continue;
            }
            let (la, ss) = sbp.eff_lengths(
                i,
                o.gap_open,
                o.gap_extend,
                o.gapped_calculation,
                c.query_length,
                db_length,
                db_num_seqs,
            )?;
            c.length_adjustment = la;
            c.eff_searchsp = ss;
        }
        // BlastExtensionParametersNew
        let mut gap_x_dropoff = 0;
        let mut gap_x_dropoff_final = 0;
        if o.gapped_calculation {
            let mut min_lambda = i32::MAX as f64;
            for k in sbp.kbp_gap_std.iter().flatten() {
                if k.lambda > 0. && k.k > 0. && k.h > 0. && min_lambda > k.lambda {
                    min_lambda = k.lambda;
                }
            }
            gap_x_dropoff = (o.gap_x_dropoff * NCBIMATH_LN2 / min_lambda) as i32;
            gap_x_dropoff_final = (o.gap_x_dropoff_final * NCBIMATH_LN2 / min_lambda)
                .max(gap_x_dropoff as f64) as i32;
        }
        // BlastHitSavingParametersNew (min_subject_length = total_length / num_seqs)
        let _ = avg_subject_length;
        let (hit_cutoffs, hit_cutoff_score_min) = hit_saving_params(o, sbp, ctx);
        // BlastInitialWordParametersNew + Update
        let k_query_len_for_hash_table = 8000;
        let concat_len: i64 = ctx.iter().map(|c| c.query_length as i64 + 1).sum();
        let container_type = if concat_len > k_query_len_for_hash_table {
            DiagContainer::DiagHash
        } else {
            DiagContainer::DiagArray
        };
        let mut word_cutoffs = vec![UngappedCutoffs::default(); ctx.len()];
        let mut cutoff_min = i32::MAX;
        let mut xdrop_max = 0;
        for (i, c) in ctx.iter().enumerate() {
            if !c.is_valid {
                word_cutoffs[i].cutoff_score = i32::MAX;
                continue;
            }
            let kbp = sbp.kbp_std[i].expect("valid context has an ungapped block");
            let wc = &mut word_cutoffs[i];
            wc.x_dropoff_init = (o.x_dropoff * NCBIMATH_LN2 / kbp.lambda).ceil() as i32;
            let gap_trigger = ((o.gap_trigger * NCBIMATH_LN2 + kbp.log_k) / kbp.lambda) as i32;
            let mut new_cutoff = if o.gapped_calculation {
                gap_trigger
            } else {
                let mut query_length = c.query_length as u64 * 2;
                if query_length == 0 {
                    query_length = 0;
                }
                let kb = sbp.kbp_std[i].unwrap();
                let searchsp =
                    (avg_subject_length as u64).min(query_length) * avg_subject_length as u64;
                let (s, _) = blast_cutoffs(1, CUTOFF_E_BLASTN, &kb, searchsp as i64, true, 0.0);
                s
            };
            new_cutoff = new_cutoff.min(hit_cutoffs[i].cutoff_score_max);
            wc.cutoff_score = new_cutoff;
            wc.x_dropoff = if wc.x_dropoff_init == 0 {
                new_cutoff
            } else {
                wc.x_dropoff_init
            };
            cutoff_min = cutoff_min.min(new_cutoff);
            xdrop_max = xdrop_max.max(wc.x_dropoff);
            wc.reduced_nucl_cutoff_score = (0.8 * new_cutoff as f64) as i32;
        }
        let mut nucl_score_table = vec![0i32; 256];
        for (i, t) in nucl_score_table.iter_mut().enumerate() {
            let mut score = 0;
            score += if i & 3 != 0 { o.penalty } else { o.reward };
            score += if (i >> 2) & 3 != 0 {
                o.penalty
            } else {
                o.reward
            };
            score += if (i >> 4) & 3 != 0 {
                o.penalty
            } else {
                o.reward
            };
            score += if i >> 6 != 0 { o.penalty } else { o.reward };
            *t = score;
        }
        Ok(Parameters {
            db_length,
            db_num_seqs,
            reward: o.reward,
            penalty: o.penalty,
            gap_open: o.gap_open,
            gap_extend: o.gap_extend,
            gap_x_dropoff,
            gap_x_dropoff_final,
            hit_cutoffs,
            hit_cutoff_score_min,
            prelim_evalue: o.expect_value,
            mask_level: if o.mask_level >= 0 { o.mask_level } else { 101 },
            word_cutoffs,
            word_cutoff_score_min: cutoff_min,
            x_dropoff_max: xdrop_max,
            container_type,
            nucl_score_table,
            ungapped_extension: true,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn task_defaults() {
        let m = BlastOptions::megablast();
        assert_eq!(
            (m.word_size, m.reward, m.penalty, m.gap_open, m.gap_extend),
            (28, 1, -2, 0, 0)
        );
        assert_eq!(m.min_diag_separation, 6);
        let b = BlastOptions::blastn();
        assert_eq!(
            (b.word_size, b.reward, b.penalty, b.gap_open, b.gap_extend),
            (11, 2, -3, 5, 2)
        );
        assert_eq!((b.gap_x_dropoff, b.min_diag_separation), (30.0, 50));
        assert_eq!(Task::parse("dc-megablast"), None);
        m.validate().unwrap();
    }

    #[test]
    fn megablast_parameters_hg38_like() {
        // A 1000-base balanced query against 3.1 Gbp in 25 sequences.
        let o = BlastOptions::megablast();
        let q: Vec<u8> = (0..1000u32).map(|i| (i % 4) as u8).collect();
        let (sbp, _) = ScoreBlk::new(1, -2, 0, 0, true, &[Some(&q), Some(&q)]).unwrap();
        let mut ctx = vec![
            ContextStats {
                query_index: 0,
                frame: 1,
                query_length: 1000,
                is_valid: true,
                ..Default::default()
            },
            ContextStats {
                query_index: 0,
                frame: -1,
                query_length: 1000,
                is_valid: true,
                ..Default::default()
            },
        ];
        let p = Parameters::new(&o, &sbp, &mut ctx, 3_100_000_000, 25, 124_000_000).unwrap();
        // gap x-drop: 25*ln2/1.28 = 13.5 -> 13; final 100*ln2/1.28 = 54.1 -> 54
        assert_eq!((p.gap_x_dropoff, p.gap_x_dropoff_final), (13, 54));
        // ungapped cutoff: gap_trigger (27 bits) in raw score, (27*ln2 + logK)/lambda = 13
        assert_eq!(p.word_cutoffs[0].cutoff_score, 13);
        // megablast leaves the ungapped x-drop option at 0, so the x-drop is the cutoff score
        // (blastn 2.17.0 on hg38 extends through a run of 6 mismatches, -12, with cutoff 13)
        assert_eq!(p.word_cutoffs[0].x_dropoff_init, 0);
        assert_eq!(p.word_cutoffs[0].x_dropoff, 13);
        assert_eq!(p.container_type, DiagContainer::DiagArray);
        // the blastn task sets 20 bits: ceil(20*ln2/lambda) with the ungapped lambda of 2/-3
        let b = BlastOptions::blastn();
        let (sbp_b, _) = ScoreBlk::new(2, -3, 5, 2, true, &[Some(&q), Some(&q)]).unwrap();
        let pb = Parameters::new(&b, &sbp_b, &mut ctx, 3_100_000_000, 25, 124_000_000).unwrap();
        assert_eq!(
            pb.word_cutoffs[0].x_dropoff_init,
            (20.0 * NCBIMATH_LN2 / sbp_b.kbp_std[0].unwrap().lambda).ceil() as i32
        );
        assert_eq!(
            pb.word_cutoffs[0].x_dropoff,
            pb.word_cutoffs[0].x_dropoff_init
        );
        assert!(pb.word_cutoffs[0].x_dropoff > 0);
        assert!(ctx[0].length_adjustment > 0 && ctx[0].eff_searchsp > 0);
        assert_eq!(ctx[0].eff_searchsp, ctx[1].eff_searchsp);
        assert_eq!(p.nucl_score_table[0], 4);
        assert_eq!(p.nucl_score_table[255], -8);
    }
}
