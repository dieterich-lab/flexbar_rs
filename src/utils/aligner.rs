use bio::alignment::pairwise::*;
use crate::types::TrimEnd;

pub struct AlignerFactory {
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
}

impl AlignerFactory {
    pub fn new(
        trim_end: TrimEnd,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix) = match trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (MIN_SCORE, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, MIN_SCORE, 0, 0),
        };

        Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
        }
    }

    pub fn get(&self, m: usize, n: usize) -> Aligner<impl Fn(u8, u8) -> i32> {
        let match_score = self.match_score;
        let mismatch_score = self.mismatch_score;

        let scoring = Scoring {
            match_fn: move |a, b| {
                if a == b || a == b'N' || b == b'N' {
                    match_score
                } else {
                    mismatch_score
                }
            },
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        };

        Aligner::with_capacity_and_scoring(m, n, scoring)
    }
}
