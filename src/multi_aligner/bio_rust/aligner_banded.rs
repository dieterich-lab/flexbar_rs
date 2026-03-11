use std::ops::Range;
use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};
use crate::multi_aligner::bio_rust::generic_aligner::create_references_with_range;
use crate::multi_aligner::generic::{AlignerMode, ReferenceWithRange};

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
    pub min_overlap: Option<usize>,
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;

pub struct BandedStateExtra {
    pub references: Vec<ReferenceWithRange>,
    pub kmer_hashes: Vec<HashMapFx<&'static [u8], Vec<u32>>>,
}

unsafe impl Send for BandedStateExtra {}

impl BioRustRawAligner<BandedExtraConfig> for BioRustBandedRawAligner {
    type StateExtra = BandedStateExtra;

    fn create_raw(
        config: &BandedBioRustAlignerConfig,
        x_capacity: usize, y_capacity: usize, allowed_errors: usize
    ) -> Self {
        let scoring = get_scoring(config);
        let window_size = if config.mode == AlignerMode::Local {
            allowed_errors
        } else {
            // I do not really understand why this value has to be so huge
            // in the case of semi-global alignments. Bug in the library?
            y_capacity - config.extra.min_overlap.unwrap_or(y_capacity)
        };
        Self::with_capacity_and_scoring(
            x_capacity, y_capacity, scoring, config.extra.k, window_size
        )
    }

    fn create_state_extra(
        config: &BioRustAlignerConfig<BandedExtraConfig>,
        references: &[SequenceWithId],
    ) -> Self::StateExtra {
        let references_with_range = create_references_with_range(config, references);

        let mut kmer_hashes: Vec<HashMapFx<&'static [u8], Vec<u32>>> = Vec::new();
        for entry in references_with_range.iter() {
            let relevant_ref = &entry.reference.sequence.as_bytes()[entry.range.clone()];
            let hash = hash_kmers(relevant_ref, config.extra.k);
            let mut pre_hash = HashMapFx::default();
            for (k, v) in hash {
                // SAFETY: We are casting the reference to the owned sequence to 'static.
                // This is safe because BandedStateExtra owns both 'references' and 'kmer_hashes',
                // and 'references' is never modified or moved individually after creation.
                let static_k: &'static [u8] = unsafe { std::mem::transmute(k) };
                pre_hash.insert(static_k, v);
            }
            kmer_hashes.push(pre_hash);
        }
        BandedStateExtra {
            references: references_with_range,
            kmer_hashes,
        }
    }

    fn perform_align(
        &mut self,
        sequence: &[u8],
        reference_index: usize,
        state_extra: &mut Self::StateExtra,
    ) -> (bio::alignment::Alignment, Range<usize>) {
        let entry = &state_extra.references[reference_index];
        let reference = &entry.reference.sequence.as_bytes()[entry.range.clone()];
        let kmer_hash = &state_extra.kmer_hashes[reference_index];
        (self.custom_with_prehash(sequence, reference, kmer_hash), entry.range.clone())
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
