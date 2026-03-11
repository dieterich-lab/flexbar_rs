use std::ops::Range;
use bio::alignment::Alignment;
use bio::alignment::pairwise::Aligner;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

use crate::multi_aligner::generic::ReferenceWithRange;
use crate::multi_aligner::bio_rust::generic_aligner::create_references_with_range;

pub type BioRustPairwiseRawAligner = Aligner<BoxedScoreFunction>;

#[derive(Clone)]
pub struct AlignerState {
    pub references: Vec<ReferenceWithRange>,
}

impl BioRustRawAligner<()> for BioRustPairwiseRawAligner {
    type StateExtra = AlignerState;

    fn create_raw(
        config: &BioRustAlignerConfig<()>,
        x_capacity: usize, y_capacity: usize, _allowed_errors: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Aligner::with_capacity_and_scoring(x_capacity, y_capacity, scoring)
    }

    fn create_state_extra(
        config: &BioRustAlignerConfig<()>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra {
        AlignerState { references: create_references_with_range(config, references) }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        state_extra: &mut Self::StateExtra,
    ) -> (Alignment, Range<usize>) {
        let entry = &state_extra.references[reference_index];
        let reference = &entry.reference.sequence.as_bytes()[entry.range.clone()];
        (self.custom(sequence, reference), entry.range.clone())
    }
}

pub type BioRustPairwiseAligner = BioRustAligner<Aligner<BoxedScoreFunction>,()>;