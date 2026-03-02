use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
}

#[derive(Clone)]
pub struct BandedAlignerData {
    references: Vec<SequenceWithId>,
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;

pub struct BandedStateExtra {
    pub kmer_hashes: Vec<HashMapFx<Vec<u8>, Vec<u32>>>,
}

impl BioRustRawAligner<BandedExtraConfig> for BioRustBandedRawAligner {
    type AuxData = BandedAlignerData;
    type StateExtra = BandedStateExtra;

    fn create_raw(
        config: &BandedBioRustAlignerConfig,
        x_capacity: usize, y_capacity: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Self::with_capacity_and_scoring(
            x_capacity, y_capacity, scoring, config.extra.k, config.extra.k * 2
        )
    }

    fn create_state_extra(
        config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra {
        let mut kmer_hash: Vec<HashMapFx<Vec<u8>, Vec<u32>>> = Vec::new();
        for r in references.iter() {
            let hash = hash_kmers(r.sequence.as_bytes(), config.extra.k);
            let mut owned_hash = HashMapFx::default();
            for (k, v) in hash {
                owned_hash.insert(k.to_vec(), v);
            }
            kmer_hash.push(owned_hash);
        }
        BandedStateExtra {
            kmer_hashes: kmer_hash,
        }
    }

    fn create_aux_data(
        _config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &[SequenceWithId]
    ) -> Self::AuxData {
        BandedAlignerData {
            references: references.to_vec(),
        }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        extra_data: &BandedAlignerData,
        state_extra: &mut Self::StateExtra,
    ) -> bio::alignment::Alignment {
        let reference = &extra_data.references[reference_index].sequence.as_bytes();
        let ref_hash = &state_extra.kmer_hashes[reference_index];
        
        let mut temp_hash = HashMapFx::default();
        for (k, v) in ref_hash {
            temp_hash.insert(k.as_slice(), v.clone());
        }
        self.custom_with_prehash(sequence, reference, &temp_hash)
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
