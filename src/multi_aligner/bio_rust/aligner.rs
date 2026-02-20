use bio::alignment::pairwise::Aligner;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAligner, BioRustAlignerConfigBase, BioRustAlignerSetup, get_scoring};

pub type BioRustAlignerConfig = BioRustAlignerConfigBase<()>;

pub type BoxedAligner = Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>;

pub struct StandardSetup;

impl BioRustAlignerSetup<()> for StandardSetup {
    type AlignerType = BoxedAligner;
    type ExtraAlignerData = ();

    fn build_extra_data(_config: &BioRustAlignerConfig, _references: &[SequenceWithId]) -> Self::ExtraAlignerData {
        ()
    }

    fn build_aligner(
        config: &BioRustAlignerConfig,
        x_capacity: usize,
        y_capacity: usize,
        _extra_data: &Self::ExtraAlignerData,
    ) -> Self::AlignerType {
        let scoring = get_scoring(config);
        Aligner::with_capacity_and_scoring(x_capacity, y_capacity, scoring)
    }

    fn perform_alignment(
        aligner: &mut Self::AlignerType,
        relevant_seq: &[u8],
        reference: &SequenceWithId,
        _extra_data: &Self::ExtraAlignerData,
        _reference_index: usize,
    ) -> bio::alignment::Alignment {
        aligner.custom(relevant_seq, reference.sequence.as_bytes())
    }
}

pub type BioRustPairwiseAligner = BioRustAligner<(), StandardSetup>;