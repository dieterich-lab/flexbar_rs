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
    pub min_overlap: Option<usize>,
    pub extra: ExtraConfig,
}

pub struct BioRustAligner<Raw: Send, ExtraConfig: Clone> {
    pub config: BioRustAlignerConfig<ExtraConfig>,
    pub y_capacity: usize,
    pub allowed_errors: usize,
    pub reference_lengths: Vec<usize>,
    _phantom: std::marker::PhantomData<Raw>,
}

pub struct BioRustAlignerState<Raw: Send, Extra: Send> {
    pub x_capacity: usize,
    pub raw: Option<Raw>,
    pub extra: Extra,
}

pub trait BioRustRawAligner<ExtraConfig: Clone> {
    type AuxData;
    type StateExtra: Send;

    fn create_raw(
        config: &BioRustAlignerConfig<ExtraConfig>,
        x_capacity: usize, y_capacity: usize, allowed_errors: usize
    ) -> Self;

    fn create_state_extra(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra;

    fn create_aux_data(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &[SequenceWithId]
    ) -> Self::AuxData;

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        aux_data: &Self::AuxData,
        state_extra: &mut Self::StateExtra,
    ) -> (bio::alignment::Alignment, std::ops::Range<usize>);
}

impl<ExtraConfig: Clone, Raw: BioRustRawAligner<ExtraConfig>+Send>
    BioRustAligner<Raw, ExtraConfig>
{
    pub fn new(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &[SequenceWithId],
        error_rate: f32
    ) -> Self {
        let reference_lengths: Vec<usize> = references.iter().map(|r| r.sequence.len()).collect();
        let y_capacity = *reference_lengths.iter().max().unwrap_or(&0);
        let allowed_errors = (error_rate * y_capacity as f32) as usize;

        BioRustAligner {
            config: config.clone(),
            y_capacity,
            allowed_errors,
            reference_lengths,
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn create_state(&self, references: &[SequenceWithId]) -> BioRustAlignerState<Raw, Raw::StateExtra> {
        BioRustAlignerState {
            x_capacity: 0,
            raw: None,
            extra: Raw::create_state_extra(&self.config, references),
        }
    }

    pub fn setup_aligner(
        &self,
        state: &mut BioRustAlignerState<Raw, Raw::StateExtra>,
        sequence_len: usize,
    ) {
        if state.raw.is_some() && sequence_len <= state.x_capacity { return; }

        let new_x_capacity = if self.config.mode != AlignerMode::Local {
            self.y_capacity
        } else if state.x_capacity == 0 {
            sequence_len
        } else {
            sequence_len + sequence_len / 10
        };

        state.x_capacity = new_x_capacity;
        state.raw = Some(Raw::create_raw(
            &self.config, state.x_capacity, self.y_capacity, self.allowed_errors
        ));
    }

    pub fn align(
        &self,
        state: &mut BioRustAlignerState<Raw, Raw::StateExtra>,
        sequence: &[u8],
        reference_index: usize,
        aux_data: &Raw::AuxData
    ) -> MultiAlignerResult {
        let aligner = state.raw.as_mut().expect("Aligner not setup. Call setup_aligner first.");
        let seq_len = sequence.len();
        let ref_len = self.reference_lengths[reference_index];
        let mode = self.config.mode;

        let offset = if mode == AlignerMode::RightOverlap && seq_len > ref_len {
            seq_len - ref_len
        } else if mode == AlignerMode::RightMin && seq_len > ref_len {
            seq_len - ref_len
        } else {
            0
        };

        let length = if (mode == AlignerMode::LeftOverlap || mode == AlignerMode::RightOverlap || mode == AlignerMode::LeftMin || mode == AlignerMode::RightMin) && seq_len > ref_len {
            ref_len
        } else {
            seq_len
        };

        let relevant_seq = &sequence[offset..offset + length];
        let (raw, reference_range) = aligner.perform_align(relevant_seq, reference_index, &aux_data, &mut state.extra);
        MultiAlignerResult {
            raw: RawAlignment::BioRust(raw),
            sequence_offset: offset,
            reference_index,
            reference_range,
        }
    }
}

pub fn get_scoring<ExtraConfig: Clone>(config: &BioRustAlignerConfig<ExtraConfig>) -> Scoring<BoxedScoreFunction> {
    let match_score = config.match_score;
    let mismatch_score = config.mismatch_score;

    let match_fn: BoxedScoreFunction = Box::new(move |a, b| {
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
