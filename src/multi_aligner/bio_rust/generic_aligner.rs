use bio::alignment::pairwise::{Scoring, MIN_SCORE};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::AlignerMode;
use crate::multi_aligner::result::{MultiAlignerResult, RawAlignment};

pub type BoxedScoreFunction = Box<dyn Fn(u8, u8) -> i32 + Send + Sync>;

#[derive(Clone)]
pub struct BioRustAlignerParameters {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_extend: i32,
    pub gap_open: i32,
    pub error_rate: f32,
    pub min_overlap: Option<usize>,
}

#[derive(Clone)]
pub struct BioRustAlignerConfig<ExtraConfig: Clone> {
    pub mode: AlignerMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub extra: ExtraConfig,
}

pub struct BioRustAligner<Raw: Send, ExtraConfig: Clone> {
    pub config: BioRustAlignerConfig<ExtraConfig>,
    pub x_capacity: usize,
    pub y_capacity: usize,
    pub reference_lengths: Vec<usize>,
    pub raw: Option<Raw>,
}

pub trait BioRustRawAligner<ExtraConfig: Clone> {
    type AuxData<'s>;
    fn create_raw(
        config: &BioRustAlignerConfig<ExtraConfig>,
        x_capacity: usize, y_capacity: usize
    ) -> Self;

    fn create_aux_data<'s>(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &'s [SequenceWithId]
    ) -> Self::AuxData<'s>;

    fn perform_align<'s>(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        aux_data: &Self::AuxData<'s>,
    ) -> bio::alignment::Alignment;
}

impl<'s, ExtraConfig: Clone, Raw: BioRustRawAligner<ExtraConfig>+Send>
    BioRustAligner<Raw, ExtraConfig>
{
    pub fn new(config: &BioRustAlignerConfig<ExtraConfig>, references: &'s [SequenceWithId]) -> Self {
        let reference_lengths: Vec<usize> = references.iter().map(|r| r.sequence.len()).collect();
        let y_capacity = *reference_lengths.iter().max().unwrap_or(&0);
        BioRustAligner {
            config: config.clone(),
            x_capacity: 0,
            y_capacity,
            reference_lengths,
            raw: None,
        }
    }

    pub fn setup_aligner(&mut self, x_capacity: usize) {
        if self.raw.is_some() && x_capacity <= self.x_capacity { return; }
        
        let new_x_capacity = if self.x_capacity == 0 {
            x_capacity
        } else {
            x_capacity + x_capacity / 10
        };

        self.x_capacity = new_x_capacity;
        self.raw = Some(Raw::create_raw(&self.config, self.x_capacity, self.y_capacity));
    }

    pub fn align<'a>(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        aux_data: &'a Raw::AuxData<'a>
    ) -> MultiAlignerResult {
        let aligner = self.raw.as_mut().expect("Aligner not setup. Call setup_aligner first.");
        let seq_len = sequence.len();
        let ref_len = self.reference_lengths[reference_index];
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
        let raw = aligner.perform_align(relevant_seq, reference_index, &aux_data);
        MultiAlignerResult {
            raw: RawAlignment::BioRust(raw),
            offset,
            reference_index,
        }
    }
}

pub fn get_scoring<ExtraConfig: Clone>(config: &BioRustAlignerConfig<ExtraConfig>) -> Scoring<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>> {
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
