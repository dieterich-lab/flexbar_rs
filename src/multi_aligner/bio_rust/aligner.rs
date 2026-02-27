use bio::alignment::Alignment;
use bio::alignment::pairwise::Aligner;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

pub type BioRustPairwiseRawAligner = Aligner<BoxedScoreFunction>;

#[derive(Clone)]
pub struct AlignerData<'s> {
    references: &'s [SequenceWithId],
}

impl BioRustRawAligner<()> for BioRustPairwiseRawAligner {
    type AuxData<'s> = AlignerData<'s>;

    fn create_raw(
        config: &BioRustAlignerConfig<()>,
        x_capacity: usize, y_capacity: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Aligner::with_capacity_and_scoring(x_capacity, y_capacity, scoring)
    }

    fn create_aux_data<'s>(
        _config: &BioRustAlignerConfig<()>, references: &'s [SequenceWithId]
    ) -> Self::AuxData<'s> {
        AlignerData { references }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        extra_data: &AlignerData
    ) -> Alignment {
        let reference = &extra_data.references[reference_index].sequence.as_bytes();
        self.custom(sequence, reference)
    }
}

pub type BioRustPairwiseAligner = BioRustAligner<Aligner<BoxedScoreFunction>,()>;