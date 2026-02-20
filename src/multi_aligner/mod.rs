
mod generic;
mod bio_rust;
pub mod result;

use std::cmp::min;
use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::{BioRustAlignerConfig, BioRustPairwiseAligner};
use crate::multi_aligner::bio_rust::aligner_banded::{BandedExtraConfig, BioRustBandedAlignerConfig, BioRustBandedPairwiseAligner};
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::MultiAlignerResult;
use log::{debug, info, warn};

fn get_minimum_kmer_k() -> usize {
    // If k becomes too low, two problems may arise:
    // - The hash may become huge.
    // - The kmer approach may become inefficient.
    4
}

enum SpecificAligner {
    BioRustPairwiseAligner(BioRustPairwiseAligner),
    BioRustBandedPairwiseAligner(BioRustBandedPairwiseAligner),
}

impl GenericAligner for SpecificAligner {
    fn align(&mut self, sequence: &[u8], reference_index: usize) 
        -> MultiAlignerResult
    {
        match self {
            SpecificAligner::BioRustPairwiseAligner(aligner) => 
                aligner.align(sequence, reference_index),
            SpecificAligner::BioRustBandedPairwiseAligner(aligner) =>
                aligner.align(sequence, reference_index),
        }
    }
}

pub struct MultiAligner {
    name: String,
    references: Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: Option<usize>,
    error_rate: f32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    kmer_k: usize,
    cached_aligners: thread_local::ThreadLocal<std::cell::RefCell<Vec<SpecificAligner>>>,
}


impl MultiAligner {
    pub fn new(
        name: String,
        references: Vec<SequenceWithId>,
        trim_end: TrimEnd,
        error_rate: f32,
        min_overlap: Option<usize>,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        k_override: Option<usize>,
    ) -> Self {
        let kmer_k = determine_kmer_k(&name, &references, error_rate, min_overlap, k_override);
        if kmer_k > 0 {
            debug!("{}: Using banded aligner with k={} with trim end {:?}.", name, kmer_k, trim_end);
        } else {
            debug!("{}: Using non-banded aligner with trim end {:?}.", name, trim_end);
        }
        Self {
            name,
            references,
            trim_end,
            error_rate,
            min_overlap,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            kmer_k,
            cached_aligners: thread_local::ThreadLocal::new(),
        }
    }

    pub fn align<'s>(&self, sequence: &[u8]) -> Option<MultiAlignerResult> {
        let aligners_cell = self.cached_aligners.get_or(
            || std::cell::RefCell::new(self.create_aligners())
        );
        let mut aligners = aligners_cell.borrow_mut();
        
        let mut best: Option<MultiAlignerResult> = None;
        for index in 0..self.references.len() {
            for aligner in aligners.iter_mut() {
                let result = aligner.align(sequence, index);
                let reference = &self.references[index].sequence.as_bytes();
                if !result.validate(
                    self.error_rate, sequence, reference, self.min_overlap, result.offset
                ) {
                    continue;
                }
                let old_score = match best.as_ref() {
                    Some(b) => b.score(),
                    None => i32::MIN,
                };
                if old_score >= result.score() { continue; }
                best = Some(result);
                break;
            }
        }
        best
    }

    fn create_aligners(&self) -> Vec<SpecificAligner> {
        let mut aligners = Vec::new();
        let required_modes = match self.trim_end {
            TrimEnd::Any => vec![
                AlignerMode::Local, AlignerMode::LeftOverlap, AlignerMode::RightOverlap
            ],
            TrimEnd::Left => vec![AlignerMode::Local, AlignerMode::LeftOverlap],
            TrimEnd::Right => vec![AlignerMode::Local, AlignerMode::RightOverlap],
            TrimEnd::LTail => vec![AlignerMode::LeftOverlap],
            TrimEnd::RTail => vec![AlignerMode::RightOverlap],
        };

        for mode in required_modes {
            let refs = self.references.clone();
            if self.kmer_k > 0 {
                let config = BioRustBandedAlignerConfig {
                    mode,
                    match_score: self.match_score,
                    mismatch_score: self.mismatch_score,
                    gap_open: self.gap_open,
                    gap_extend: self.gap_extend,
                    extra: BandedExtraConfig { k: self.kmer_k },
                };
                aligners.push(SpecificAligner::BioRustBandedPairwiseAligner(
                    BioRustBandedPairwiseAligner::new(config, refs)
                ));
            } else {
                debug!("{}: Using standard aligner for mode {:?}.", self.name, mode);
                let config = self.get_aligner_config(mode);
                aligners.push(SpecificAligner::BioRustPairwiseAligner(
                    BioRustPairwiseAligner::new(config, refs)
                ));
            }
        }
        aligners
    }

    fn get_aligner_config(&self, mode: AlignerMode) -> BioRustAlignerConfig {
        BioRustAlignerConfig {
            mode,
            match_score: self.match_score,
            mismatch_score: self.mismatch_score,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            extra: (),
        }
    }
}

fn determine_kmer_k(
    name: &str,
    references: &Vec<SequenceWithId>,
    error_rate: f32,
    min_overlap: Option<usize>,
    k_override: Option<usize>
) -> usize {
    let mut highest_possible_k = usize::MAX;
    for r in references.iter() {
        let possible_k = determine_kmer_k_for_sequence(&r.sequence, error_rate, min_overlap);
        if possible_k < highest_possible_k {
            highest_possible_k = possible_k;
        }
    }
    match k_override {
        Some(0) => {
            info!("{}: Using non-banded aligner on user request.", name);
            0
        },
        Some(val) => {
            if val > highest_possible_k {
                warn!(
                        "Warning: The k-mer size ({}, <{}) for {} might be too big.",
                        val, highest_possible_k, name
                    );
                warn!("Warning: Valid alignments for the given error rate/overlap may be missed.");
            }
            val
        },
        None => {
            if highest_possible_k > get_minimum_kmer_k() {
                highest_possible_k
            } else {
                warn!(
                    "Warning: The k-mer size ({}) for {} is too small - using non-banded aligner.",
                    highest_possible_k, name
                );
                warn!("Warning: Performance might be improved using a lower error rate and/or a larger overlap");
                0
            }
        },
    }
}

fn determine_kmer_k_for_sequence(
    sequence: &str,
    error_rate: f32,
    min_overlap: Option<usize>,
) -> usize {
    let mut max_clean_length: usize = 0;
    let mut present_clean_length: usize = 0;
    let mut n_count = 0;
    for c in sequence.chars() {
        if c == 'N' {
            if present_clean_length > max_clean_length {
                max_clean_length = present_clean_length;
            }
            present_clean_length = 0;
            n_count += 1;
        } else {
            present_clean_length += 1;
        }
    }
    if present_clean_length > max_clean_length {
        max_clean_length = present_clean_length;
    }
    let allowed_errors = (max_clean_length as f32 * error_rate) as usize;
    let min_match_length = max_clean_length / (allowed_errors + 1);
    match min_overlap {
        Some(min_overlap) => {
            if min_overlap < n_count {
                // That is a pathological case, which will give a warning anyway.
                0
            } else {
                min(min_match_length, min_overlap - n_count)
            }
        },
        None => min_match_length,
    }
}

#[cfg(test)]
mod tests {
    use crate::context::SequenceWithId;
    use crate::types::TrimEnd;
    use super::MultiAligner;
    use super::result::MultiAlignerResult;

    fn run_alignment(
        trim_end: TrimEnd,
        error_rate: f32,
        min_overlap: Option<usize>,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        seq_x: &[u8],
        references: Vec<SequenceWithId>
    ) -> Option<MultiAlignerResult> {
        let aligner = MultiAligner::new(
            "Test Aligner".to_string(),
            references,
            trim_end,
            error_rate,
            min_overlap,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            None,
        );
        aligner.align(seq_x)
    }

    fn get_as_references(raw: &[u8]) -> Vec<SequenceWithId> {
       vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: String::from_utf8_lossy(raw).to_string(),
       }]
    }

    #[test]
    fn test_aligner_trivial() {
        let references = get_as_references(b"ACGT");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(3), 1, -1, -1, -1,
            b"ACGT", references
        );
        assert!(result.is_some());
        let alignment = result.unwrap();
        assert_eq!(alignment.score(), 4);
        assert_eq!(alignment.get_start_index(), 0);
        assert_eq!(alignment.get_end_index(), 4);
    }

    #[test]
    fn test_aligner_right_overhang() {
        let references = get_as_references(b"NNNTTTTCCC");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(3), 1, -1, -1, -1,
            b"AAAAATTT", references
        );
        let alignment = result.unwrap();
        // Matching prefix NNNTTT at position 2
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.get_start_index(), 2);
        assert_eq!(alignment.get_end_index(), 8);
    }

    #[test]
    fn test_aligner_right_overhang_ambiguous() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(3), 1, -1, -1, -1,
            b"AAAAATTT", references
        );
        let alignment = result.unwrap();
        // Matching prefix TTT at position 5
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 8);
    }

    #[test]
    fn test_aligner_right_overhang_with_error() {
        let references = get_as_references(b"TTTTTTT");
        let result = run_alignment(
            TrimEnd::Right, 0.2,Some(5), 1, -1, -1, -1,
            b"AAAAATTATT", references
        );
        let alignment = result.unwrap();
        // Matching prefix TTxTT with an error starting at x[5]
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded() {
        let references = get_as_references(b"TTTTT");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(5), 1, -1, -1, -1,
            b"AAAAATTTTTAAA", references
        );
        let alignment = result.unwrap();
        // Prefect match att position 5
        assert_eq!(alignment.score(), 5);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded_with_error() {
        let references = get_as_references(b"TTTTT");
        let result = run_alignment(
            TrimEnd::Right,  0.2, Some(5), 1, -1, -1, -1,
            b"AAAAATTATTAAA", references
        );
        let alignment = result.unwrap();
        // Match with error TTxTT at position 5
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded_with_wildcard() {
        let references = get_as_references(b"TTNTT");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(5),1, -1, -1, -1,
            b"AAAAATTATTAAA", references
        );
        let alignment = result.unwrap();
        // Matching TTNTT at position 5
        assert_eq!(alignment.score(), 5);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_wrong_side_overhang() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Right, 0.2, Some(3), 1, -1, -1, -1,
            b"TTTAAAA", references
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_aligner_right_partial_embedded() {
        // This a strange case: The multi_aligner will score an invalid barcode match
        // in the middle of the string higher than an exact prefix match at the end.
        // So the same read may be mapped to a barcode with TrimEnd::RTail but not
        // with TrimEnd::Right.
        let references = get_as_references(b"NNNTTTTTCCCC");
        let result = run_alignment(
            TrimEnd::Right, 0.0, Some(6), 1, -1, -1, -1,
            b"AAAAATTTTTAAATTT", references
        );
        let alignment = result.unwrap();
        // Matching NNNTTTTTT (prefix) at position 2
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.get_start_index(), 10);
        assert_eq!(alignment.get_end_index(), 16);
    }

    #[test]
    fn test_aligner_left_overhang() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Left, 0.0, Some(1), 1, -1, -1, -1,
            b"TTTAAAAA", references
        );
        let alignment = result.unwrap();
        // Match TTT (suffix) starting x[0]
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 0);
        assert_eq!(alignment.get_end_index(), 3);
    }

    #[test]
    fn test_aligner_rtail_overhang() {
        // See also: test_aligner_right_partial_embedded()!!!
        let references = get_as_references(b"NNNTTTTTCCCC");
        let result = run_alignment(
            TrimEnd::RTail, 0.0, Some(1), 1, -1, -1, -1,
            b"AAAAATTTTTAAATTT", references
        );
        let alignment = result.unwrap();
        // Match NNNTTT starting x[10]
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.get_start_index(), 10);
        assert_eq!(alignment.get_end_index(), 16);
    }

    #[test]
    fn test_aligner_k_override_zero() {
        let references = get_as_references(b"ACGTACGTACGT"); // Length 12
        // Default estimate_k for error_rate 0.0 is max_no_n = 12
        // 12 > 8, so default k would be 8 (banded aligner used since 8 >= 4)
        
        let aligner = MultiAligner::new(
            "Test Aligner".to_string(),
            references,
            TrimEnd::Right,
            0.0,
            Some(4),
            1, -1, -1, -1,
            Some(0), // Override k to 0 should force old aligner
        );
        
        let aligners = aligner.create_aligners();
        for a in aligners {
            match a {
                super::SpecificAligner::BioRustPairwiseAligner(_) => (),
                _ => panic!("Expected BioRustPairwiseAligner when k_override is 0"),
            }
        }
    }

    #[test]
    fn test_aligner_k_override_banded() {
        let references = get_as_references(b"ACGT");
        // Default estimate_k for length 4 is 4.
        
        let aligner = MultiAligner::new(
            "Test Aligner".to_string(),
            references,
            TrimEnd::Right,
            0.0,
            Some(4),
            1, -1, -1, -1,
            Some(5), // Override k to 5 should use banded aligner even if k=4 was estimated
        );
        
        let aligners = aligner.create_aligners();
        for a in aligners {
            match a {
                super::SpecificAligner::BioRustBandedPairwiseAligner(_) => (),
                _ => panic!("Expected BioRustBandedPairwiseAligner when k_override is 5"),
            }
        }
    }

    #[test]
    fn test_determine_k() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, None, None
        );
        // We need to detect matches of length 22 (28-6)
        assert_eq!(k, 22);
    }

    #[test]
    fn test_determine_k_small_overlap() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, Some(11), None
        );
        // With an overlap of 11 and 6 'N's we must be able to detect
        // matches of length 5. So K can be at most 5.
        assert_eq!(k, 5);
    }

    #[test]
    fn test_determine_k_with_min_overlap_low() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTGC".to_string(), // length 9, 6 Ns, max_no_n = 3
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, Some(8), None
        );
        // With only 3 bases other than 'N', we need a k below 4 - so we disable
        // the banded aligner.
        assert_eq!(k, 0);
    }

    #[test]
    fn test_determine_k_with_min_overlap_high() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, None, Some(0)
        );
        // Normale 22 - but the user disabled the banded aligner.
        assert_eq!(k, 0);
    }

    #[test]
    fn test_determine_k_with_error_rate() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.1, None, None
        );
        // With error rate 0.0 we would need to detect matches of length 22 (28-6).
        // With an error rate of 0.1, we may see 2 errors in 22 bases, splitting
        // the alignment worst case into 3 parts of about 7 bases.
        assert_eq!(k, 7);
    }

    #[test]
    fn test_format_alignment_wildcard() {
        let x = b"ACGTATCG";
        let y = b"ACGTNTCG";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.0, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);
        assert!(formatted.contains("ACGTATCG"));
        assert!(formatted.contains("||||.|||"));
        assert!(formatted.contains("ACGTNTCG"));
    }

    #[test]
    fn test_format_alignment_embedded() {
        let x = b"AAAAACCCCTTTTT";
        let y = b"CCCC";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.0, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);
        // Expected alignment:
        // AAAAACCCCTTTTT
        //      ||||
        //      CCCC
        assert!(formatted.contains("AAAAACCCCTTTTT"));
        assert!(formatted.contains("||||"));
        assert!(formatted.contains("CCCC"));
    }

    #[test]
    fn test_format_alignment_prefix_suffix() {
        let x = b"AAACCC";
        let y = b"CCCGGG";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.0, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        // AAACCC---
        //    |||
        // ---CCCGGG

        assert!(formatted.contains("AAACCC"));
        assert!(formatted.contains("|||"));
        assert!(formatted.contains("CCCGGG"));
    }

    #[test]
    fn test_format_alignment_suffix_prefix() {
        let x = b"CCCGGG";
        let y = b"AAACCC";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.0, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        // ---CCCGGG
        //    |||
        // AAACCC---

        assert!(formatted.contains("CCCGGG"));
        assert!(formatted.contains("|||"));
        assert!(formatted.contains("AAACCC"));
    }

    #[test]
    fn test_format_alignment_partial_with_error() {
        let x = b"CCCCCGGG";
        let y = b"AAACCACC";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.2, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CCCCCGGG"));
        assert!(formatted.contains("|\\|"));
        assert!(formatted.contains("AAACCACC"));
    }

    #[test]
    fn test_format_alignment_partial_with_deletion() {
        let x = b"CAATCCTTTT";
        let y = b"CCCCAACCT";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.2, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CAATCCTTT"));
        assert!(formatted.contains("CCCCAA CCT"));
    }

    #[test]
    fn test_format_alignment_partial_with_insertion() {
        let x = b"CAACCTTTT";
        let y = b"CCCCAATCCT";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Any, 0.2, Some(1), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CAA CCTTTT"));
        assert!(formatted.contains("CCCCAATCCT"));
    }

    #[test]
    fn test_umi_tag() {
        let x = b"CAATCCTTT";
        let y = b"NNNTTT";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Right, 0.0, Some(6), 1, -1, -1, -1, x, references);
        let alignment = result.unwrap();
        assert_eq!(alignment.get_umi_tags(x, y), vec![String::from("TCC")])
    }

    #[test]
    fn test_offset_handling() {
        let x = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let y =b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);
        let result = run_alignment(TrimEnd::Right, 0.0, Some(11), 1, -1, -5, -1, x, references);
        assert!(result.is_some());
        let alignment = result.unwrap();
        let umi_tags = alignment.get_umi_tags(x, y);
        
        assert_eq!(umi_tags.len(), 1);
        assert_eq!(umi_tags[0], "GTGCCC");

        // Check format
        let formatted = alignment.format(x, y);
        assert_eq!(formatted, concat!(
            "Alignment with score 21:\n",
            "  CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA \n",
            "  ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´......||||||||||||||||||||| \n",
            "                                              NNNNNNTAGCTCTGTAGGCACCATCAAT\n"
        ));
    }

    #[test]
    fn invalid_alignment_at_end() {
        let x = b"TATTCGCTTTGGGCGCGCGCCGCGGCTGGACGAGGCGCCGCCGGCATCATAGCTCTGTAGGACCATCAATA";
        let y = b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);

        let result = run_alignment(TrimEnd::Right, 0.0, Some(11), 1, -1, -5, -1, x, references);
        assert!(result.is_none());
    }
}