use bio::alignment::pairwise::{Scoring, MIN_SCORE};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::{MultiAlignerResult, RawAlignment};

pub struct BioRustAlignerConfigBase<ExtraConfig> {
    pub mode: AlignerMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub extra: ExtraConfig,
}

pub trait BioRustAlignerSetup<ExtraConfig> {
    type AlignerType;
    type ExtraAlignerData;

    fn build_extra_data(config: &BioRustAlignerConfigBase<ExtraConfig>, references: &[SequenceWithId]) -> Self::ExtraAlignerData;
    
    fn build_aligner(
        config: &BioRustAlignerConfigBase<ExtraConfig>,
        x_capacity: usize,
        y_capacity: usize,
        extra_data: &Self::ExtraAlignerData,
    ) -> Self::AlignerType;

    fn perform_alignment(
        aligner: &mut Self::AlignerType,
        relevant_seq: &[u8],
        reference: &SequenceWithId,
        extra_data: &Self::ExtraAlignerData,
        reference_index: usize,
    ) -> bio::alignment::Alignment;
}

pub struct BioRustAligner<ExtraConfig, Setup: BioRustAlignerSetup<ExtraConfig>> {
    pub(crate) aligner: Option<Setup::AlignerType>,
    pub(crate) references: Vec<SequenceWithId>,
    pub(crate) extra_data: Setup::ExtraAlignerData,
    pub(crate) config: BioRustAlignerConfigBase<ExtraConfig>,
    pub(crate) x_capacity: usize,
    pub(crate) y_capacity: usize,
}

impl<ExtraConfig, Setup: BioRustAlignerSetup<ExtraConfig>> BioRustAligner<ExtraConfig, Setup> {
    pub fn new(config: BioRustAlignerConfigBase<ExtraConfig>, references: Vec<SequenceWithId>) -> Self {
        let y_capacity = references.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        let extra_data = Setup::build_extra_data(&config, &references);
        
        BioRustAligner {
            aligner: None,
            references,
            extra_data,
            config,
            x_capacity: 0,
            y_capacity,
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
        let new_aligner = Setup::build_aligner(&self.config, self.x_capacity, self.y_capacity, &self.extra_data);
        self.aligner = Some(new_aligner);
    }
}

pub fn get_scoring<ExtraConfig>(config: &BioRustAlignerConfigBase<ExtraConfig>) -> Scoring<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>> {
    let match_score = config.match_score;
    let mismatch_score = config.mismatch_score;

    let match_fn: Box<dyn Fn(u8, u8) -> i32 + Send + Sync> = Box::new(move |a, b| {
        if a == b || a == b'N' || b == b'N' {
            match_score
        } else {
            mismatch_score
        }
    });

    let xclip_prefix = match config.mode {
        AlignerMode::LeftOverlap => MIN_SCORE,
        _ => 0,
    };
    let xclip_suffix = match config.mode {
        AlignerMode::RightOverlap => MIN_SCORE,
        _ => 0,
    };
    let yclip_suffix = match config.mode {
        AlignerMode::RightOverlap => 0,
        _ => MIN_SCORE,
    };
    let yclip_prefix = match config.mode {
        AlignerMode::LeftOverlap => 0,
        _ => MIN_SCORE,
    };

    Scoring {
        match_fn,
        match_scores: None,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        xclip_prefix,
        xclip_suffix,
        yclip_prefix,
        yclip_suffix,
    }
}

impl<ExtraConfig, Setup: BioRustAlignerSetup<ExtraConfig>> GenericAligner for BioRustAligner<ExtraConfig, Setup> {
    fn align(&mut self, sequence: &[u8], reference_index: usize) -> MultiAlignerResult {
        self.setup_aligner(sequence.len());
        let reference = &self.references[reference_index];
        let aligner = self.aligner.as_mut().unwrap();

        let seq_len = sequence.len();
        let ref_len = reference.sequence.len();
        let mode = self.config.mode;
        
        let offset = if mode == AlignerMode::RightOverlap && seq_len > ref_len {
            seq_len - ref_len
        } else {
            0
        };
        
        let length = if mode != AlignerMode::Local && seq_len > ref_len {
            ref_len
        } else {
            seq_len
        };
        
        let relevant_seq = &sequence[offset..offset + length];

        let raw = Setup::perform_alignment(aligner, relevant_seq, reference, &self.extra_data, reference_index);
        
        MultiAlignerResult {
            raw: RawAlignment::BioRust(raw),
            offset,
            reference_index,
        }
    }
}
