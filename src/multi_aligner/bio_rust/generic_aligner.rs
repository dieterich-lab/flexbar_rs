use bio::alignment::Alignment;
use bio::alignment::pairwise::{Scoring, MIN_SCORE};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode, ReferenceWithRange};
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
    _phantom: std::marker::PhantomData<Raw>,
}

pub struct BioRustAlignerState<Raw: Send, Extra: Send> {
    pub x_capacity: usize,
    pub raw: Option<Raw>,
    pub extra: Extra,
}

pub trait BioRustRawAligner<ExtraConfig: Clone> {
    type StateExtra: Send;

    fn create_raw(
        config: &BioRustAlignerConfig<ExtraConfig>,
        x_capacity: usize, y_capacity: usize, allowed_errors: usize
    ) -> Self;

    fn create_state_extra(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra;

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        state_extra: &mut Self::StateExtra,
    ) -> Alignment;

    // The get_reference signature is awkward. Unfortunately, the references
    // are kept in aligner's "state_extra", because they might be tightly linked with
    // other structures in the state, such as the kmer hashes for the banded aligner.
    fn get_reference<'r, 's: 'r, 'e: 'r>(&'s self, state_extra: &'e Self::StateExtra, reference_index: usize)
        -> &'r ReferenceWithRange;
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
    ) -> MultiAlignerResult {
        let aligner = state.raw.as_mut().expect("Aligner not setup. Call setup_aligner first.");
        let mode = self.config.mode;
        let r = aligner.get_reference(&state.extra, reference_index);
        let seq_len = sequence.len();
        let full_ref_len = r.reference.sequence.len();
        let seq_start = match mode {
            AlignerMode::RightOverlap | AlignerMode::RightMin => if seq_len > full_ref_len
                { seq_len - full_ref_len + r.get_aligner_start_offset() } else { 0 },
            _ => if seq_len > r.get_aligner_start_offset()
                { r.get_aligner_start_offset() } else { 0 },
        };
        let seq_end = match mode {
            AlignerMode::LeftOverlap | AlignerMode::LeftMin => if seq_len > full_ref_len
                { full_ref_len - r.get_aligner_end_offset() } else { seq_len },
            _ => if seq_len > seq_start + r.get_aligner_end_offset()
                { seq_len - r.get_aligner_end_offset() } else { seq_len },
        };
        let relevant_seq = &sequence[seq_start..seq_end];
        let reference_range = r.aligner_range.clone();
        let stripped_leading_n_count = r.get_aligner_start_offset();
        let stripped_trailing_n_count = r.get_aligner_end_offset();
        let raw = aligner.perform_align(relevant_seq, reference_index, &mut state.extra);

        MultiAlignerResult {
            raw: RawAlignment::BioRust(raw),
            sequence_offset: seq_start,
            reference_index,
            reference_range,
            stripped_leading_n_count,
            stripped_trailing_n_count,
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

pub fn create_references_with_range<Extra: Clone>(
    config: &BioRustAlignerConfig<Extra>,
    references: &[SequenceWithId],
) -> Vec<ReferenceWithRange> {
    let mut refs_with_range = Vec::new();
    for r in references {
        let range = match config.mode {
            AlignerMode::LeftMin => {
                let min_overlap = config.min_overlap.unwrap_or(r.sequence.len());
                (r.sequence.len() - min_overlap)..r.sequence.len()
            },
            AlignerMode::RightMin => {
                let min_overlap = config.min_overlap.unwrap_or(r.sequence.len());
                0..min_overlap
            },
            _ => 0..r.sequence.len(),
        };

        let ref_bytes = r.sequence.as_bytes();
        let mut core_start = range.start;
        while core_start < range.end && ref_bytes[core_start] == b'N' {
            core_start += 1;
        }
        let mut core_end = range.end;
        while core_end > core_start && ref_bytes[core_end - 1] == b'N' {
            core_end -= 1;
        }
        let aligner_range = core_start..core_end;

        refs_with_range.push(
            ReferenceWithRange { reference: r.clone(), range, aligner_range }
        );
    }
    refs_with_range
}