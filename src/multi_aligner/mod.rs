
mod generic;
mod bio_rust;
pub mod result;

use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::{BioRustAlignerConfig, BioRustPairwiseAligner};
use crate::multi_aligner::bio_rust::aligner_banded::{BioRustBandedAlignerConfig, BioRustBandedPairwiseAligner};
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::MultiAlignerResult;

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
    references: Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: Option<usize>,
    error_rate: f32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    k_override: Option<usize>,
    cached_aligners: thread_local::ThreadLocal<std::cell::RefCell<Vec<SpecificAligner>>>,
}


impl MultiAligner {
    pub fn new(
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
        Self {
            references,
            trim_end,
            error_rate,
            min_overlap,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            k_override,
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

        let k = match self.k_override {
            Some(0) => 0,
            Some(val) => val,
            None => BioRustBandedPairwiseAligner::estimate_k(&self.references, self.error_rate),
        };

        for mode in required_modes {
            let refs = self.references.clone();
            if k >= 4 {
                let config = BioRustBandedAlignerConfig {
                    mode,
                    match_score: self.match_score,
                    mismatch_score: self.mismatch_score,
                    gap_open: self.gap_open,
                    gap_extend: self.gap_extend,
                    k,
                };
                aligners.push(SpecificAligner::BioRustBandedPairwiseAligner(
                    BioRustBandedPairwiseAligner::new(config, refs)
                ));
            } else {
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
        }
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
    fn test_aliagner_rtail_overhang() {
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