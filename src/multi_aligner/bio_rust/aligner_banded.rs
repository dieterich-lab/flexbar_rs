use std::ops::Range;
use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner, get_scoring, BoxedScoreFunction, BioRustAligner};
use crate::multi_aligner::generic::AlignerMode;

#[derive(Clone, Copy)]
pub struct BandedExtraConfig {
    pub k: usize,
    pub min_overlap: Option<usize>,
}

#[derive(Clone)]
pub struct BandedAlignerData {
}

pub type BandedBioRustAlignerConfig = BioRustAlignerConfig<BandedExtraConfig>;

pub type BioRustBandedRawAligner = Aligner<BoxedScoreFunction>;

pub struct OwnedReferenceWithRange {
    pub reference: Vec<u8>,
    pub range: Range<usize>,
}

pub struct BandedStateExtra {
    pub owned_references: Vec<OwnedReferenceWithRange>,
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
        let mut owned_references = Vec::new();
        for r in references.iter() {
            let range = match config.mode {
                AlignerMode::LeftMin => {
                    let min_overlap = config.extra.min_overlap.unwrap_or(r.sequence.len());
                    (r.sequence.len() - min_overlap)..r.sequence.len()
                },
                AlignerMode::RightMin => {
                    let min_overlap = config.extra.min_overlap.unwrap_or(r.sequence.len());
                    0..min_overlap
                },
                _ => 0..r.sequence.len(),
            };
            owned_references.push(OwnedReferenceWithRange { reference: r.sequence.as_bytes().to_vec(), range });
        }

        let mut kmer_hashes: Vec<HashMapFx<&'static [u8], Vec<u32>>> = Vec::new();
        for entry in owned_references.iter() {
            let relevant_ref = &entry.reference[entry.range.clone()];
            let hash = hash_kmers(relevant_ref, config.extra.k);
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
    ) -> (bio::alignment::Alignment, Range<usize>) {
        let entry = &state_extra.owned_references[reference_index];
        let reference = &entry.reference[entry.range.clone()];
        let kmer_hash = &state_extra.kmer_hashes[reference_index];
        (self.custom_with_prehash(sequence, reference, kmer_hash), entry.range.clone())
    }
}

pub type BioRustBandedPairwiseAligner = BioRustAligner<
    Aligner<BoxedScoreFunction>,
    BandedExtraConfig
>;
