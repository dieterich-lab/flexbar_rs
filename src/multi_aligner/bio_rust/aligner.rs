use std::ops::Range;
use bio::alignment::Alignment;
use bio::alignment::pairwise::Aligner;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

pub type BioRustPairwiseRawAligner = Aligner<BoxedScoreFunction>;

#[derive(Clone)]
pub struct ReferenceWithRange {
    pub reference: SequenceWithId,
    pub range: Range<usize>,
}

#[derive(Clone)]
pub struct AlignerData {
    pub references: Vec<ReferenceWithRange>,
}

impl BioRustRawAligner<()> for BioRustPairwiseRawAligner {
    type AuxData = AlignerData;
    type StateExtra = ();

    fn create_raw(
        config: &BioRustAlignerConfig<()>,
        x_capacity: usize, y_capacity: usize, _allowed_errors: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Aligner::with_capacity_and_scoring(x_capacity, y_capacity, scoring)
    }

    fn create_state_extra(
        _config: &BioRustAlignerConfig<()>,
        _references: &[SequenceWithId],
    ) -> Self::StateExtra {
        ()
    }

    fn create_aux_data(
        config: &BioRustAlignerConfig<()>, references: &[SequenceWithId]
    ) -> Self::AuxData {
        let mut refs_with_range = Vec::new();
        for r in references {
            let range = match config.mode {
                crate::multi_aligner::generic::AlignerMode::LeftMin => {
                    let min_overlap = config.min_overlap.unwrap_or(r.sequence.len());
                    (r.sequence.len() - min_overlap)..r.sequence.len()
                },
                crate::multi_aligner::generic::AlignerMode::RightMin => {
                    let min_overlap = config.min_overlap.unwrap_or(r.sequence.len());
                    0..min_overlap
                },
                _ => 0..r.sequence.len(),
            };
            refs_with_range.push(ReferenceWithRange { reference: r.clone(), range });
        }
        AlignerData { references: refs_with_range }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        extra_data: &AlignerData,
        _state_extra: &mut Self::StateExtra,
    ) -> (Alignment, Range<usize>) {
        let entry = &extra_data.references[reference_index];
        let reference = &entry.reference.sequence.as_bytes()[entry.range.clone()];
        (self.custom(sequence, reference), entry.range.clone())
    }
}

pub type BioRustPairwiseAligner = BioRustAligner<Aligner<BoxedScoreFunction>,()>;