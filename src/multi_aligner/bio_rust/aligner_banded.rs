use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAligner, BioRustAlignerConfigBase, BioRustAlignerSetup, get_scoring};

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
}

pub type BioRustBandedAlignerConfig = BioRustAlignerConfigBase<BandedExtraConfig>;

pub type BoxedBandedAligner = Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>;

pub struct BandedSetup;

impl BioRustAlignerSetup<BandedExtraConfig> for BandedSetup {
    type AlignerType = BoxedBandedAligner;
    type ExtraAlignerData = Vec<HashMapFx<Vec<u8>, Vec<u32>>>;

    fn build_extra_data(config: &BioRustBandedAlignerConfig, references: &[SequenceWithId]) -> Self::ExtraAlignerData {
        references.iter()
            .map(|r| {
                hash_kmers(r.sequence.as_bytes(), config.extra.k)
                    .into_iter()
                    .map(|(k, v)| (k.to_vec(), v))
                    .collect()
            })
            .collect()
    }

    fn build_aligner(
        config: &BioRustBandedAlignerConfig,
        x_capacity: usize,
        y_capacity: usize,
        _extra_data: &Self::ExtraAlignerData,
    ) -> Self::AlignerType {
        let scoring = get_scoring(config);
        // Using w = k * 2 as a heuristic for window size
        Aligner::with_capacity_and_scoring(x_capacity, y_capacity, scoring, config.extra.k, config.extra.k * 2)
    }

    fn perform_alignment(
        aligner: &mut Self::AlignerType,
        relevant_seq: &[u8],
        reference: &SequenceWithId,
        extra_data: &Self::ExtraAlignerData,
        reference_index: usize,
    ) -> bio::alignment::Alignment {
        let ref_hash = &extra_data[reference_index];
        let temp_hash: HashMapFx<&[u8], Vec<u32>> = ref_hash.iter()
            .map(|(k, v)| (k.as_slice(), v.clone()))
            .collect();
        
        aligner.custom_with_prehash(relevant_seq, reference.sequence.as_bytes(), &temp_hash)
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<BandedExtraConfig, BandedSetup>;
