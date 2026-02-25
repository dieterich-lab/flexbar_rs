use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
}

pub struct BandedAlignerData {
    kmer_hash: Vec<HashMapFx<Vec<u8>, Vec<u32>>>,
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;
impl BioRustRawAligner<BandedExtraConfig> for BioRustBandedRawAligner {
    type AuxData = BandedAlignerData;

    fn create_raw(
        config: &BandedBioRustAlignerConfig,
        x_capacity: usize, y_capacity: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Self::with_capacity_and_scoring(
            x_capacity, y_capacity, scoring, config.extra.k, config.extra.k * 2
        )
    }

    fn create_aux_data(
        config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &Vec<SequenceWithId>
    ) -> Self::AuxData {
        BandedAlignerData {
            kmer_hash: references.iter()
                .map( | r| {
                    hash_kmers(r.sequence.as_bytes(), config.extra.k)
                    .into_iter()
                    .map( | (k, v) | (k.to_vec(), v))
                    .collect()
                })
                .collect()
        }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        references: &[SequenceWithId],
        extra_data: &BandedAlignerData,
    ) -> bio::alignment::Alignment {
        let ref_hash = &extra_data.kmer_hash[reference_index];
        let temp_hash: HashMapFx<&[u8], Vec<u32>> = ref_hash.iter()
            .map(|(k, v)| (k.as_slice(), v.clone()))
            .collect();
        let reference = references[reference_index].sequence.as_bytes();
        self.custom_with_prehash(sequence, reference, &temp_hash)
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
