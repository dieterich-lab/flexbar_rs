
mod generic;
mod bio_rust;
mod result;

use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::{BioRustAlignerConfig, BioRustPairwiseAligner};
use crate::multi_aligner::generic::{AlignerMode, GenericAligner};
use crate::multi_aligner::result::AlignmentResult;

enum SpecificAligner<'r> {
    BioRustPairwiseAligner(BioRustPairwiseAligner<'r>),
}

impl<'r> GenericAligner<'r> for SpecificAligner<'r> {
    fn align<'s>(&mut self, sequence: &'s [u8], reference_index: usize) 
        -> AlignmentResult<'r, 's> 
    {
        match self {
            SpecificAligner::BioRustPairwiseAligner(aligner) => 
                aligner.align(sequence, reference_index),
        }
    }
}

pub struct MultiAligner<'r> {
    references: &'r Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    cached_aligners: thread_local::ThreadLocal<std::cell::RefCell<Vec<SpecificAligner<'r>>>>,
}


impl<'r> MultiAligner<'r> {
    pub fn new(
        references: &'r Vec<SequenceWithId>,
        trim_end: TrimEnd,
        error_rate: f32,
        min_overlap: usize,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
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
            cached_aligners: thread_local::ThreadLocal::new(),
        }
    }

    pub fn align<'s>(&mut self, sequence: &'s [u8]) -> Option<AlignmentResult<'r, 's>> {
        let aligners_cell = self.cached_aligners.get_or(
            || std::cell::RefCell::new(self.create_aligners())
        );
        let mut aligners = aligners_cell.borrow_mut();
        
        let mut best: Option<AlignmentResult<'r, 's>> = None;
        for index in 0..self.references.len() {
            for aligner in aligners.iter_mut() {
                let result = aligner.align(sequence, index);
                if !result.validate(self.error_rate, self.min_overlap) { continue; }
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

    fn create_aligners(&self) -> Vec<SpecificAligner<'r>> {
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
            let config = self.get_aligner_config(mode);
            aligners.push(SpecificAligner::BioRustPairwiseAligner(
                BioRustPairwiseAligner::new(config, self.references)
            ))
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
    use std::rc::Rc;
    use bio::alignment::Alignment;
    use crate::types::TrimEnd;
    use crate::context::SequenceWithId;
    use crate::multi_aligner::generic::AlignerMode;
    use super::MultiAligner;
    use super::result::{AlignmentResult, RawAlignment};

    fn run_alignment<'r, 's>(
        trim_end: TrimEnd,
        error_rate: f32,
        min_overlap: usize,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        seq_x: &'s [u8],
        references: &'r Vec<SequenceWithId>
    ) -> Option<AlignmentResult<'r, 's>> {
        let mut aligner = MultiAligner::new(
            references,
            trim_end,
            error_rate,
            min_overlap,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
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
            TrimEnd::Right, 0.0, 3, 1, -1, -1, -1,
            b"ACGT", &references
        );
        assert!(result.is_some());
        let alignment = result.unwrap();
        assert_eq!(alignment.score(), 4);
        assert_eq!(alignment.get_start_index(), 0);
        assert_eq!(alignment.get_end_index(), 4);
    }

    // #[test]
    // fn test_aligner_right_overhang() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.0,3, 1, -1, -1, -1,
    //         b"AAAAATTT", b"NNNTTTTCCC"
    //     );
    //     // Matching prefix NNNTTT at position 2
    //     assert_eq!(alignment.score, 6);
    //     assert_eq!(alignment.xstart, 2);
    //     assert_eq!(alignment.xend, 8);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 6);
    // }
    //
    // #[test]
    // fn test_aligner_right_overhang_ambiguous() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.0, 3, 1, -1, -1, -1,
    //         b"AAAAATTT", b"TTTT"
    //     );
    //     // Matching prefix TTT at position 5
    //     assert_eq!(alignment.score, 3);
    //     assert_eq!(alignment.xstart, 5);
    //     assert_eq!(alignment.xend, 8);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 3);
    // }
    //
    // #[test]
    // fn test_aligner_right_overhang_with_error() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.2,5, 1, -1, -1, -1,
    //         b"AAAAATTATT", b"TTTTTTT"
    //     );
    //     // Matching prefix TTxTT with an error starting at x[5]
    //     assert_eq!(alignment.score, 3);
    //     assert_eq!(alignment.xstart, 5);
    //     assert_eq!(alignment.xend, 10);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 5);
    // }
    //
    // #[test]
    // fn test_aligner_right_embedded() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.0., 5, 1, -1, -1, -1,
    //         b"AAAAATTTTTAAA", b"TTTTT"
    //     );
    //     // Prefect match att position 5
    //     assert_eq!(alignment.score, 5);
    //     assert_eq!(alignment.xstart, 5);
    //     assert_eq!(alignment.xend, 10);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 5);
    // }
    //
    // #[test]
    // fn test_aligner_right_embedded_with_error() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right,  0.2, 5, 1, -1, -1, -1,
    //         b"AAAAATTATTAAA", b"TTTTT"
    //     );
    //     // Match with error TTxTT at position 5
    //     assert_eq!(alignment.score, 3);
    //     assert_eq!(alignment.xstart, 5);
    //     assert_eq!(alignment.xend, 10);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 5);
    // }
    //
    // #[test]
    // fn test_aligner_right_embedded_with_wildcard() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.0, 5,1, -1, -1, -1,
    //         b"AAAAATTATTAAA", b"TTNTT"
    //     );
    //     // Matching TTNTT at position 5
    //     assert_eq!(alignment.score, 5);
    //     assert_eq!(alignment.xstart, 5);
    //     assert_eq!(alignment.xend, 10);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 5);
    // }
    //
    // #[test]
    // fn test_aligner_right_wrong_side_overhang() {
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 0.2, 3, 1, -1, -1, -1,
    //         b"TTTAAAA", b"TTTT"
    //     );
    //     // Match TTT (prefix) starting x[0]
    //     assert_eq!(alignment.score, 3);
    //     assert_eq!(alignment.xstart, 0);
    //     assert_eq!(alignment.xend, 3);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 3);
    // }
    //
    // #[test]
    // fn test_aligner_right_partial_embedded() {
    //     // This a strange case: The multi_aligner will score an invalid barcode match
    //     // in the middle of the string higher than an exact prefix match at the end.
    //     // So the same read may be mapped to a barcode with TrimEnd::RTail but not
    //     // with TrimEnd::Right.
    //     let alignment = run_alignment(
    //         TrimEnd::Right, 1, -1, -1, -1,
    //         b"AAAAATTTTTAAATTT", b"NNNTTTTTCCCC"
    //     );
    //     // Matching NNNTTTTTT (prefix) at position 2
    //     assert_eq!(alignment.score, 8);
    //     assert_eq!(alignment.xstart, 2);
    //     assert_eq!(alignment.xend, 10);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 8);
    // }
    //
    // #[test]
    // fn test_aligner_left_overhang() {
    //     let alignment = run_alignment(
    //         TrimEnd::Left, 1, -1, -1, -1,
    //         b"TTTAAAAA", b"TTTT"
    //     );
    //     // Match TTT (suffix) starting x[0]
    //     assert_eq!(alignment.score, 3);
    //     assert_eq!(alignment.xstart, 0);
    //     assert_eq!(alignment.xend, 3);
    //     assert_eq!(alignment.ystart, 1);
    //     assert_eq!(alignment.yend, 4);
    // }
    //
    // #[test]
    // fn test_aliagner_rtail_overhang() {
    //     // See also: test_aligner_right_partial_embedded()!!!
    //     let alignment = run_alignment(
    //         TrimEnd::RTail, 1, -1, -1, -1,
    //         b"AAAAATTTTTAAATTT", b"NNNTTTTTCCCC"
    //     );
    //     // Match NNNTTT starting x[10]
    //     assert_eq!(alignment.score, 6);
    //     assert_eq!(alignment.xstart, 10);
    //     assert_eq!(alignment.xend, 16);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 6);
    // }
    //
    // #[test]
    // fn test_format_alignment_wildcard() {
    //     let x = b"ACGTNTCG";
    //     let y = b"ACGTATCG";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //     assert_eq!(l1, "ACGTNTCG");
    //     assert_eq!(l2, "||||.|||");
    //     assert_eq!(l3, "ACGTATCG");
    // }
    //
    // #[test]
    // fn test_format_alignment_embedded() {
    //     let x = b"AAAAACCCCTTTTT";
    //     let y = b"CCCC";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //     // Expected alignment:
    //     // AAAAACCCCTTTTT
    //     //      ||||
    //     //      CCCC
    //     // Note: Global alignment with Any might include clips.
    //     // Actually run_alignment with TrimEnd::Any uses 0 for all clips.
    //
    //     // Let's check what actually happens.
    //     assert_eq!(l1, "AAAAACCCCTTTTT");
    //     assert_eq!(l2, "     ||||     ");
    //     assert_eq!(l3, "     CCCC     ");
    // }
    //
    // #[test]
    // fn test_format_alignment_prefix_suffix() {
    //     let x = b"AAACCC";
    //     let y = b"CCCGGG";
    //     // To get a prefix of y aligned to suffix of x, we might need TrimEnd::Any
    //     // or a specific mode. TrimEnd::Any sets all clips to 0.
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //
    //     // AAACCC---
    //     //    |||
    //     // ---CCCGGG
    //
    //     assert_eq!(l1, "AAACCC   ");
    //     assert_eq!(l2, "   |||   ");
    //     assert_eq!(l3, "   CCCGGG");
    // }
    //
    // #[test]
    // fn test_format_alignment_suffix_prefix() {
    //     let x = b"CCCGGG";
    //     let y = b"AAACCC";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //
    //     // ---CCCGGG
    //     //    |||
    //     // AAACCC---
    //
    //     assert_eq!(l1, "   CCCGGG");
    //     assert_eq!(l2, "   |||   ");
    //     assert_eq!(l3, "AAACCC   ");
    // }
    //
    // #[test]
    // fn test_format_alignment_partial_with_error() {
    //     let x = b"CCCCCGGG";
    //     let y = b"AAACCACC";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //
    //     // ---CCCGGG
    //     //    |||
    //     // AAACCC---
    //
    //     assert_eq!(l1, "   CCCCCGGG");
    //     assert_eq!(l2, "   ||\\||   ");
    //     assert_eq!(l3, "AAACCACC   ");
    // }
    //
    // #[test]
    // fn test_format_alignment_partial_with_deletion() {
    //     let x = b"CAATCCTTTT";
    //     let y = b"CCCCAACCT";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //
    //     // ---CCCGGG
    //     //    |||
    //     // AAACCC---
    //
    //     assert_eq!(l1, "   CAATCCTTTT");
    //     assert_eq!(l2, "   ||| |||   ");
    //     assert_eq!(l3, "CCCCAA CCT   ");
    // }
    //
    // #[test]
    // fn test_format_alignment_partial_with_insertion() {
    //     let x = b"CAACCTTTT";
    //     let y = b"CCCCAATCCT";
    //     let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
    //     let (l1, l2, l3) = format_alignment(&alignment, x, y);
    //
    //     // ---CCCGGG
    //     //    |||
    //     // AAACCC---
    //
    //     assert_eq!(l1, "   CAA CCTTTT");
    //     assert_eq!(l2, "   ||| |||   ");
    //     assert_eq!(l3, "CCCCAATCCT   ");
    // }
}