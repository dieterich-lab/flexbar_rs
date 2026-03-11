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
    ) -> Alignment {
        let reference = &state_extra.references[reference_index].get_aligner_sequence();
        self.custom(sequence, reference)
    }

    fn get_reference<'r, 's: 'r, 'e: 'r>(
        &'s self, state_extra: &'e Self::StateExtra, reference_index: usize
    ) -> &'r ReferenceWithRange {
        &state_extra.references[reference_index]
    }
}

pub type BioRustPairwiseAligner = BioRustAligner<Aligner<BoxedScoreFunction>,()>;