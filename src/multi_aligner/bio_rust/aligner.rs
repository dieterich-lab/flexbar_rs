use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::{AlignmentResult, RawAlignment};

pub struct BioRustAlignerConfig {
    pub mode: AlignerMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

pub type BoxedAligner = Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>;

pub struct BioRustPairwiseAligner<'a> {
    aligner: Option<BoxedAligner>,
    references: &'a [SequenceWithId],
    config: BioRustAlignerConfig,
    x_capacity: usize,
    y_capacity: usize,
}

impl<'a> BioRustPairwiseAligner<'a> {
    pub fn new(config: BioRustAlignerConfig, references: &'a [SequenceWithId]) -> Self {
        let y_capacity = references.iter().map(|r| r.sequence.len()).max().unwrap();
        BioRustPairwiseAligner {
            aligner: None,
            references,
            config,
            x_capacity: 0,
            y_capacity
        }
    }
    fn setup_aligner(&mut self, x_capacity: usize) {
        match self.aligner {
            Some(_) => {
                if self.x_capacity >= x_capacity {
                    return;
                } else {
                    self.x_capacity = x_capacity + x_capacity / 10;
                }
            }
            None => {
                self.x_capacity = x_capacity;
            }
        }
        let new_aligner = self.build_aligner();
        self.aligner = Some(new_aligner);
    }

    fn build_aligner(&self) -> BoxedAligner {
        let match_score = self.config.match_score;
        let mismatch_score = self.config.mismatch_score;

        let match_fn: Box<dyn Fn(u8, u8) -> i32 + Send + Sync> = Box::new(move |a, b| {
            if a == b || a == b'N' || b == b'N' {
                match_score
            } else {
                mismatch_score
            }
        });
        let xclip_prefix = match self.config.mode {
            AlignerMode::LeftOverlap => MIN_SCORE,
            _ => 0,
        };
        let xclip_suffix = match self.config.mode {
            AlignerMode::RightOverlap => MIN_SCORE,
            _ => 0,
        };
        let yclip_suffix = match self.config.mode {
            AlignerMode::RightOverlap => 0,
            _ => MIN_SCORE,
        };
        let yclip_prefix = match self.config.mode {
            AlignerMode::LeftOverlap => MIN_SCORE,
            _ => 0,
        };

        let scoring = Scoring {
            match_fn,
            match_scores: None,
            gap_open: self.config.gap_open,
            gap_extend: self.config.gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
        };

        Aligner::with_capacity_and_scoring(self.x_capacity, self.y_capacity, scoring)
    }
}

impl<'a> GenericAligner<'a> for BioRustPairwiseAligner<'a> {
    fn align<'s>(&mut self, sequence: &'s [u8], reference_index: usize)
         -> AlignmentResult<'a, 's>
    {
        let reference = &self.references[reference_index];
        self.setup_aligner(sequence.len());
        let aligner = self.aligner.as_mut().unwrap();

        let seq_len = sequence.len();
        let ref_len = reference.sequence.len();
        let mode = self.config.mode;
        let offset = if mode == AlignerMode::RightOverlap
            && seq_len > ref_len { seq_len - ref_len } else { 0 };
        let length = if self.config.mode != AlignerMode::Local
            && seq_len > ref_len { ref_len } else { seq_len };
        let relevant_seq = &sequence[offset..offset + length];

        let raw = aligner.custom(relevant_seq, reference.sequence.as_bytes());
        AlignmentResult { raw: RawAlignment::BioRust(raw), reference, sequence, mode, offset }
    }
}