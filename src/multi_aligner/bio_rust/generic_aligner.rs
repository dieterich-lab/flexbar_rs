use bio::alignment::pairwise::{Scoring, MIN_SCORE};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode};
use crate::multi_aligner::result::{MultiAlignerResult, RawAlignment};

pub type BoxedScoreFunction = Box<dyn Fn(u8, u8) -> i32 + Send + Sync>;

pub struct BioRustAlignerConfig<ExtraConfig> {
    pub mode: AlignerMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub extra: ExtraConfig,
}

pub struct BioRustAligner<Raw: Send, ExtraConfig, AuxData> {
    pub config: BioRustAlignerConfig<ExtraConfig>,
    pub references: Vec<SequenceWithId>,
    pub aux_data: AuxData,
    pub y_capacity: usize,
    pub _phantom: std::marker::PhantomData<Raw>,
}

pub struct BioRustAlignerState<Raw> {
    pub raw: Option<Raw>,
    pub x_capacity: usize,
}

impl<Raw> Default for BioRustAlignerState<Raw> {
    fn default() -> Self {
        Self { raw: None, x_capacity: 0 }
    }
}

pub trait BioRustRawAligner<ExtraConfig, AuxData> {
    fn create_raw(
        config: &BioRustAlignerConfig<ExtraConfig>,
        x_capacity: usize, y_capacity: usize
    ) -> Self;

    fn create_aux_data(
        config: &BioRustAlignerConfig<ExtraConfig>,
        references: &Vec<SequenceWithId>
    ) -> AuxData;

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        references: &[SequenceWithId],
        extra_data: &AuxData,
    ) -> bio::alignment::Alignment;
}


impl<ExtraConfig, AuxData, Raw: BioRustRawAligner<ExtraConfig, AuxData>+Send>
    BioRustAligner<Raw, ExtraConfig, AuxData>
{
    pub fn new(config: BioRustAlignerConfig<ExtraConfig>, references: Vec<SequenceWithId>) -> Self {
        let aux_data = Raw::create_aux_data(&config, &references);
        let y_capacity = references.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        BioRustAligner {
            references,
            aux_data,
            config,
            y_capacity,
            _phantom: std::marker::PhantomData,
        }
    }

    fn setup_aligner(&self, state: &mut BioRustAlignerState<Raw>, x_capacity: usize) {
        if state.raw.is_some() && x_capacity <= state.x_capacity { return; }
        
        let new_x_capacity = if state.x_capacity == 0 {
            x_capacity
        } else {
            x_capacity + x_capacity / 10
        };

        state.x_capacity = new_x_capacity;
        state.raw = Some(Raw::create_raw(&self.config, state.x_capacity, self.y_capacity));
    }

    pub fn align(&self, state: &mut BioRustAlignerState<Raw>, sequence: &[u8], reference_index: usize) -> MultiAlignerResult {
        self.setup_aligner(state, sequence.len());
        let aligner = state.raw.as_mut().unwrap();
        let reference = &self.references[reference_index];
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
        let raw = aligner.perform_align(relevant_seq, reference_index, &self.references, &self.aux_data);
        MultiAlignerResult {
            raw: RawAlignment::BioRust(raw),
            offset,
            reference_index,
        }
    }
}

pub fn get_scoring<ExtraConfig>(config: &BioRustAlignerConfig<ExtraConfig>) -> Scoring<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>> {
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
