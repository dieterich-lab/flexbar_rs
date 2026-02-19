use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::pairwise::{Scoring, MIN_SCORE};
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::{MultiAlignerResult, RawAlignment};

pub struct BioRustBandedAlignerConfig {
    pub mode: AlignerMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub k: usize,
}

pub type BoxedBandedAligner = Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>;

pub struct BioRustBandedPairwiseAligner {
    aligner: Option<BoxedBandedAligner>,
    references: Vec<SequenceWithId>,
    reference_hashes: Vec<HashMapFx<Vec<u8>, Vec<u32>>>,
    config: BioRustBandedAlignerConfig,
    x_capacity: usize,
    y_capacity: usize,
}

impl BioRustBandedPairwiseAligner {
    pub fn new(config: BioRustBandedAlignerConfig, references: Vec<SequenceWithId>) -> Self {
        let y_capacity = references.iter().map(|r| r.sequence.len()).max().unwrap();
        let reference_hashes = references.iter()
            .map(|r| {
                hash_kmers(r.sequence.as_bytes(), config.k)
                    .into_iter()
                    .map(|(k, v)| (k.to_vec(), v))
                    .collect()
            })
            .collect();

        BioRustBandedPairwiseAligner {
            aligner: None,
            references,
            reference_hashes,
            config,
            x_capacity: 0,
            y_capacity
        }
    }

    fn setup_aligner(&mut self, x_capacity: usize) {
        match self.aligner {
            Some(_) => {
                if self.x_capacity >= x_capacity {
                    return;
                } else {
                    self.x_capacity = x_capacity + x_capacity / 10;
                }
            }
            None => {
                self.x_capacity = x_capacity;
            }
        }
        let new_aligner = self.build_aligner();
        self.aligner = Some(new_aligner);
    }

    fn build_aligner(&self) -> BoxedBandedAligner {
        let match_score = self.config.match_score;
        let mismatch_score = self.config.mismatch_score;

        let match_fn: Box<dyn Fn(u8, u8) -> i32 + Send + Sync> = Box::new(move |a, b| {
            if a == b || a == b'N' || b == b'N' {
                match_score
            } else {
                mismatch_score
            }
        });
        let xclip_prefix = match self.config.mode {
            AlignerMode::LeftOverlap => MIN_SCORE,
            _ => 0,
        };
        let xclip_suffix = match self.config.mode {
            AlignerMode::RightOverlap => MIN_SCORE,
            _ => 0,
        };
        let yclip_suffix = match self.config.mode {
            AlignerMode::RightOverlap => 0,
            _ => MIN_SCORE,
        };
        let yclip_prefix = match self.config.mode {
            AlignerMode::LeftOverlap => 0,
            _ => MIN_SCORE,
        };

        let scoring = Scoring {
            match_fn,
            match_scores: None,
            gap_open: self.config.gap_open,
            gap_extend: self.config.gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
        };

        // Using w = k * 2 as a heuristic for window size
        Aligner::with_capacity_and_scoring(self.x_capacity, self.y_capacity, scoring, self.config.k, self.config.k * 2)
    }

    pub fn estimate_k(references: &[SequenceWithId], error_rate: f32) -> usize {
        let mut min_k = usize::MAX;
        for r in references.iter() {
            let seq = r.sequence.as_bytes();
            let mut max_no_n = 0;
            let mut current_no_n = 0;
            for &b in seq {
                if b != b'N' {
                    current_no_n += 1;
                } else {
                    if current_no_n > max_no_n {
                        max_no_n = current_no_n;
                    }
                    current_no_n = 0;
                }
            }
            if current_no_n > max_no_n {
                max_no_n = current_no_n;
            }

            let errors = (seq.len() as f32 * error_rate) as usize;
            let k = if errors == 0 { max_no_n } else { max_no_n / (errors + 1) };
            if k < min_k {
                min_k = k;
            }
        }
        
        if min_k > 8 { 8 } else { min_k }
    }
}

impl GenericAligner for BioRustBandedPairwiseAligner {
    fn align(&mut self, sequence: &[u8], reference_index: usize)
         -> MultiAlignerResult
    {
        self.setup_aligner(sequence.len());
        let reference = &self.references[reference_index];
        let ref_hash = &self.reference_hashes[reference_index];
        
        let temp_hash: HashMapFx<&[u8], Vec<u32>> = ref_hash.iter()
            .map(|(k, v)| (k.as_slice(), v.clone()))
            .collect();

        let aligner = self.aligner.as_mut().unwrap();

        let seq_len = sequence.len();
        let ref_len = reference.sequence.len();
        let mode = self.config.mode;
        let offset = if mode == AlignerMode::RightOverlap
            && seq_len > ref_len { seq_len - ref_len } else { 0 };
        let length = if mode != AlignerMode::Local
            && seq_len > ref_len { ref_len } else { seq_len };
        let relevant_seq = &sequence[offset..offset + length];

        let raw = aligner.custom_with_prehash(relevant_seq, reference.sequence.as_bytes(), &temp_hash);
        MultiAlignerResult { raw: RawAlignment::BioRust(raw), mode, offset, reference_index }
    }
}
