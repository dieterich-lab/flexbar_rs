
mod generic;
mod bio_rust;
mod result;

use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::{BioRustAlignerConfig, BioRustPairwiseAligner};
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::AlignmentResult;

enum SpecificAligner<'r> {
    BioRustPairwiseAligner(BioRustPairwiseAligner<'r>),
}

impl<'r> GenericAligner<'r> for SpecificAligner<'r> {
    fn align<'s>(&mut self, sequence: &'s [u8], reference_index: usize) 
        -> AlignmentResult<'r, 's> 
    {
        match self {
            SpecificAligner::BioRustPairwiseAligner(aligner) => 
                aligner.align(sequence, reference_index),
        }
    }
}

pub struct MultiAligner<'r> {
    references: &'r Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    cached_aligners: thread_local::ThreadLocal<std::cell::RefCell<Vec<SpecificAligner<'r>>>>,
}


impl<'r> MultiAligner<'r> {
    pub fn new(
        references: &'r Vec<SequenceWithId>,
        trim_end: TrimEnd,
        error_rate: f32,
        min_overlap: usize,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        Self {
            references,
            trim_end,
            error_rate,
            min_overlap,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            cached_aligners: thread_local::ThreadLocal::new(),
        }
    }

    pub fn align<'s>(&mut self, sequence: &'s [u8]) -> Option<AlignmentResult<'r, 's>> {
        let aligners_cell = self.cached_aligners.get_or(
            || std::cell::RefCell::new(MultiAligner::<'r>::create_aligners(self))
        );
        let mut aligners = aligners_cell.borrow_mut();
        
        let mut best: Option<AlignmentResult<'r, 's>> = None;
        for index in 0..self.references.len() {
            for aligner in aligners.iter_mut() {
                let result = aligner.align(sequence, index);
                if !result.validate(self.error_rate, self.min_overlap) { continue; }
                let old_score = match best.as_ref() {
                    Some(b) => b.score(),
                    None => i32::MAX,
                };
                if old_score >= result.score() { continue; }
                best = Some(result);
                break;
            }
        }
        best
    }

    fn create_aligners(&self) -> Vec<SpecificAligner<'r>> {
        let mut aligners = Vec::new();
        let required_modes = match self.trim_end {
            TrimEnd::Any => vec![
                AlignerMode::Local, AlignerMode::LeftOverlap, AlignerMode::RightOverlap
            ],
            TrimEnd::Left => vec![AlignerMode::Local, AlignerMode::LeftOverlap],
            TrimEnd::Right => vec![AlignerMode::Local, AlignerMode::RightOverlap],
            TrimEnd::LTail => vec![AlignerMode::LeftOverlap],
            TrimEnd::RTail => vec![AlignerMode::RightOverlap],
        };
        for mode in required_modes {
            let config = self.get_aligner_config(mode);
            aligners.push(SpecificAligner::BioRustPairwiseAligner(
                BioRustPairwiseAligner::new(config, self.references)
            ))
        }
        aligners
    }

    fn get_aligner_config(&self, mode: AlignerMode) -> BioRustAlignerConfig {
        BioRustAlignerConfig {
            mode,
            match_score: self.match_score,
            mismatch_score: self.mismatch_score,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
        }
    }
}