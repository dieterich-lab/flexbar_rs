use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
}

#[derive(Clone)]
pub struct BandedAlignerData<'s> {
    references: &'s [SequenceWithId],
    kmer_hashes: Vec<HashMapFx<&'s [u8], Vec<u32>>>,
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;
impl BioRustRawAligner<BandedExtraConfig> for BioRustBandedRawAligner {
    type AuxData<'a> = BandedAlignerData<'a>;

    fn create_raw(
        config: &BandedBioRustAlignerConfig,
        x_capacity: usize, y_capacity: usize
    ) -> Self {
        let scoring = get_scoring(config);
        Self::with_capacity_and_scoring(
            x_capacity, y_capacity, scoring, config.extra.k, config.extra.k * 2
        )
    }

    fn create_aux_data<'s>(
        config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &'s [SequenceWithId]
    ) -> Self::AuxData<'s> {
        let mut kmer_hash: Vec<HashMapFx<&'s [u8], Vec<u32>>> = Vec::new();
        for r in references.iter() {
            kmer_hash.push(hash_kmers(r.sequence.as_bytes(), config.extra.k));
        }
        BandedAlignerData {
            kmer_hashes: kmer_hash,
            references,
        }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        extra_data: &BandedAlignerData,
    ) -> bio::alignment::Alignment {
        let reference = &extra_data.references[reference_index].sequence.as_bytes();
        let ref_hash = &extra_data.kmer_hashes[reference_index];
        self.custom_with_prehash(sequence, reference, ref_hash)
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
