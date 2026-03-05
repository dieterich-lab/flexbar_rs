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
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;

pub struct BandedStateExtra {
    pub owned_references: Vec<Vec<u8>>,
    pub kmer_hashes: Vec<HashMapFx<&'static [u8], Vec<u32>>>,
}

unsafe impl Send for BandedStateExtra {}

impl BioRustRawAligner<BandedExtraConfig> for BioRustBandedRawAligner {
    type AuxData = BandedAlignerData;
    type StateExtra = BandedStateExtra;

    fn create_raw(
        config: &BandedBioRustAlignerConfig,
        x_capacity: usize, y_capacity: usize, allowed_errors: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Self::with_capacity_and_scoring(
            x_capacity, y_capacity, scoring, config.extra.k, allowed_errors
        )
    }

    fn create_state_extra(
        config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra {
        let mut owned_references = Vec::new();
        for r in references.iter() {
            owned_references.push(r.sequence.as_bytes().to_vec());
        }

        let mut kmer_hashes: Vec<HashMapFx<&'static [u8], Vec<u32>>> = Vec::new();
        for r in owned_references.iter() {
            let hash = hash_kmers(r, config.extra.k);
            let mut pre_hash = HashMapFx::default();
            for (k, v) in hash {
                // SAFETY: We are casting the reference to the owned sequence to 'static.
                // This is safe because BandedStateExtra owns both 'owned_references' and 'kmer_hashes',
                // and 'owned_references' is never modified or moved individually after creation.
                let static_k: &'static [u8] = unsafe { std::mem::transmute(k) };
                pre_hash.insert(static_k, v);
            }
            kmer_hashes.push(pre_hash);
        }
        BandedStateExtra {
            owned_references,
            kmer_hashes,
        }
    }

    fn create_aux_data(
        _config: &BioRustAlignerConfig<BandedExtraConfig>,
        _references: &[SequenceWithId]
    ) -> Self::AuxData {
        BandedAlignerData {}
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        _extra_data: &BandedAlignerData,
        state_extra: &mut Self::StateExtra,
    ) -> bio::alignment::Alignment {
        let reference = &state_extra.owned_references[reference_index];
        let kmer_hash = &state_extra.kmer_hashes[reference_index];
        self.custom_with_prehash(sequence, reference, kmer_hash)
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
