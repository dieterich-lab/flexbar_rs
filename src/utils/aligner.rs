use bio::alignment::{Alignment, AlignmentOperation};
use bio::alignment::pairwise::*;
use log::trace;
use crate::types::TrimEnd;


struct Capacity {
    n: usize,
    m: usize,
}

struct CachedAligner {
    aligner: Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>,
    capacity: Capacity,
}

pub struct AlignerFactory {
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    cache: thread_local::ThreadLocal<std::cell::RefCell<Option<CachedAligner>>>,
}

impl AlignerFactory {
    pub fn new(
        trim_end: TrimEnd,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        let xclip_prefix = match trim_end {
            TrimEnd::LTail => MIN_SCORE,
            _ => 0,
        };
        let xclip_suffix = match trim_end {
            TrimEnd::RTail => MIN_SCORE,
            _ => 0,
        };
        let yclip_prefix = match trim_end {
            TrimEnd::Right => MIN_SCORE,
            TrimEnd::RTail => MIN_SCORE,
            _ => 0,
        };
        let yclip_suffix= match trim_end {
            TrimEnd::Left => MIN_SCORE,
            TrimEnd::LTail => MIN_SCORE,
            _ => 0,
        };

        Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            cache: thread_local::ThreadLocal::new(),
        }
    }

    pub fn get(&self, n: usize, m: usize) -> std::cell::RefMut<'_, Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>> {
        let cache_cell = self.cache.get_or(|| std::cell::RefCell::new(None));
        let mut cache = cache_cell.borrow_mut();
        let new_cap = self.get_capacity_requirement(&cache, n, m);

        match new_cap {
            Some(cap) => {
                let cached_aligner = self.build_aligner(cap);
                *cache = Some(cached_aligner);
            },
            None => ()
        }

        std::cell::RefMut::map(cache, |opt| {
            &mut opt.as_mut().unwrap().aligner
        })
    }

    fn get_capacity_requirement(&self, cache: &Option<CachedAligner>, n: usize, m: usize)
                                -> Option<Capacity> {
        match cache {
            Some(c) => {
                if n <= c.capacity.n && m <= c.capacity.m {
                    None
                } else {
                    let new_n = if n > c.capacity.n { n + 16 } else { c.capacity.n };
                    let new_m = if m > c.capacity.m { m + 16 } else { c.capacity.m };
                    Some(Capacity { n: new_n, m: new_m }) 
                }
            },
            None => Some(Capacity { n, m })
        }
    }

    fn build_aligner(&self, capacity: Capacity) -> CachedAligner {
        let match_score = self.match_score;
        let mismatch_score = self.mismatch_score;

        let match_fn: Box<dyn Fn(u8, u8) -> i32 + Send + Sync> = Box::new(move |a, b| {
            if a == b || a == b'N' || b == b'N' {
                match_score
            } else {
                mismatch_score
            }
        });

        let scoring = Scoring {
            match_fn,
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        };

        let aligner = Aligner::with_capacity_and_scoring(capacity.m, capacity.n, scoring);
        CachedAligner { aligner, capacity }
    }
}

pub fn format_alignment(alignment: &Alignment, x: &[u8], y: &[u8]) -> (String, String, String) {
    let mut line1 = String::new();
    let mut line2 = String::new();
    let mut line3 = String::new();

    let mut x_pos = 0;
    let mut y_pos = 0;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {
                let x_char = x[x_pos];
                let y_char = y[y_pos];
                line1.push(x_char as char);
                line3.push(y_char as char);
                if x_char == b'N' || y_char == b'N' {
                    line2.push('.');
                } else {
                    line2.push('|');
                }
                x_pos += 1;
                y_pos += 1;
            }
            AlignmentOperation::Subst => {
                let x_char = x[x_pos];
                let y_char = y[y_pos];
                line1.push(x_char as char);
                line3.push(y_char as char);
                if x_char == b'N' || y_char == b'N' {
                    line2.push('.');
                } else {
                    line2.push('\\');
                }
                x_pos += 1;
                y_pos += 1;
            }
            AlignmentOperation::Ins => {
                line1.push(x[x_pos] as char);
                line2.push(' ');
                line3.push(' ');
                x_pos += 1;
            }
            AlignmentOperation::Del => {
                line1.push(' ');
                line2.push(' ');
                line3.push(y[y_pos] as char);
                y_pos += 1;
            }
            AlignmentOperation::Xclip(len) => {
                for _ in 0..*len {
                    line1.push(x[x_pos] as char);
                    line2.push(' ');
                    line3.push(' ');
                    x_pos += 1;
                }
            }
            AlignmentOperation::Yclip(len) => {
                for _ in 0..*len {
                    line1.push(' ');
                    line2.push(' ');
                    line3.push(y[y_pos] as char);
                    y_pos += 1;
                }
            }
        }
    }
    (line1, line2, line3)
}

pub fn trace_alignment(alignment: &Alignment, x: &[u8], y: &[u8]) {
    let (line1, line2, line3) = format_alignment(alignment, x, y);

    trace!(
        "Alignment with score {}:\n  {}\n  {}\n  {}",
        alignment.score, line1, line2, line3
    );
}

mod tests {
    use bio::alignment::Alignment;
    use crate::types::TrimEnd;
    use super::format_alignment;

    fn run_alignment(
        trim_end: TrimEnd,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        seq_x: &[u8],
        seq_y: &[u8]
    ) -> Alignment {
        let aligner_factory = super::AlignerFactory::new(
            trim_end, match_score, mismatch_score, gap_open, gap_extend
        );
        let mut aligner = aligner_factory.get(seq_x.len(), seq_y.len());
        aligner.custom(seq_x, seq_y)
    }

    #[test]
    fn test_aligner_trivial() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"ACGT", b"ACGT"
        );
        assert_eq!(alignment.score, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 4);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 4);
    }

    #[test]
    fn test_aligner_right_overhang() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTT", b"NNNTTTTCCC"
        );
        // Matching prefix NNNTTT at position 2
        assert_eq!(alignment.score, 6);
        assert_eq!(alignment.xstart, 2);
        assert_eq!(alignment.xend, 8);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 6);
    }

    #[test]
    fn test_aligner_right_overhang_ambiguous() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTT", b"TTTT"
        );
        // Matching prefix TTT at position 5
        assert_eq!(alignment.score, 3);
        assert_eq!(alignment.xstart, 5);
        assert_eq!(alignment.xend, 8);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 3);
    }

    #[test]
    fn test_aligner_right_overhang_with_error() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTATT", b"TTTTTTT"
        );
        // Matching prefix TTxTT with an error starting at x[5]
        assert_eq!(alignment.score, 3);
        assert_eq!(alignment.xstart, 5);
        assert_eq!(alignment.xend, 10);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 5);
    }

    #[test]
    fn test_aligner_right_embedded() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTTTTAAA", b"TTTTT"
        );
        // Prefect match att position 5
        assert_eq!(alignment.score, 5);
        assert_eq!(alignment.xstart, 5);
        assert_eq!(alignment.xend, 10);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 5);
    }

    #[test]
    fn test_aligner_right_embedded_with_error() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTATTAAA", b"TTTTT"
        );
        // Match with error TTxTT at position 5
        assert_eq!(alignment.score, 3);
        assert_eq!(alignment.xstart, 5);
        assert_eq!(alignment.xend, 10);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 5);
    }

    #[test]
    fn test_aligner_right_embedded_with_wildcard() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTATTAAA", b"TTNTT"
        );
        // Matching TTNTT at position 5
        assert_eq!(alignment.score, 5);
        assert_eq!(alignment.xstart, 5);
        assert_eq!(alignment.xend, 10);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 5);
    }

    #[test]
    fn test_aligner_right_wrong_side_overhang() {
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"TTTAAAA", b"TTTT"
        );
        // Match TTT (prefix) starting x[0]
        assert_eq!(alignment.score, 3);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 3);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 3);
    }

    #[test]
    fn test_aligner_right_partial_embedded() {
        // This a strange case: The aligner will score an invalid barcode match
        // in the middle of the string higher than an exact prefix match at the end.
        // So the same read may be mapped to a barcode with TrimEnd::RTail but not
        // with TrimEnd::Right.
        let alignment = run_alignment(
            TrimEnd::Right, 1, -1, -1, -1,
            b"AAAAATTTTTAAATTT", b"NNNTTTTTCCCC"
        );
        // Matching NNNTTTTTT (prefix) at position 2
        assert_eq!(alignment.score, 8);
        assert_eq!(alignment.xstart, 2);
        assert_eq!(alignment.xend, 10);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 8);
    }

    #[test]
    fn test_aligner_left_overhang() {
        let alignment = run_alignment(
            TrimEnd::Left, 1, -1, -1, -1,
            b"TTTAAAAA", b"TTTT"
        );
        // Match TTT (suffix) starting x[0]
        assert_eq!(alignment.score, 3);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 3);
        assert_eq!(alignment.ystart, 1);
        assert_eq!(alignment.yend, 4);
    }

    #[test]
    fn test_aliagner_rtail_overhang() {
        // See also: test_aligner_right_partial_embedded()!!!
        let alignment = run_alignment(
            TrimEnd::RTail, 1, -1, -1, -1,
            b"AAAAATTTTTAAATTT", b"NNNTTTTTCCCC"
        );
        // Match NNNTTT starting x[10]
        assert_eq!(alignment.score, 6);
        assert_eq!(alignment.xstart, 10);
        assert_eq!(alignment.xend, 16);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 6);
    }

    #[test]
    fn test_format_alignment_wildcard() {
        let x = b"ACGTNTCG";
        let y = b"ACGTATCG";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);
        assert_eq!(l1, "ACGTNTCG");
        assert_eq!(l2, "||||.|||");
        assert_eq!(l3, "ACGTATCG");
    }

    #[test]
    fn test_format_alignment_embedded() {
        let x = b"AAAAACCCCTTTTT";
        let y = b"CCCC";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);
        // Expected alignment:
        // AAAAACCCCTTTTT
        //      ||||
        //      CCCC
        // Note: Global alignment with Any might include clips.
        // Actually run_alignment with TrimEnd::Any uses 0 for all clips.
        
        // Let's check what actually happens.
        assert_eq!(l1, "AAAAACCCCTTTTT");
        assert_eq!(l2, "     ||||     ");
        assert_eq!(l3, "     CCCC     ");
    }

    #[test]
    fn test_format_alignment_prefix_suffix() {
        let x = b"AAACCC";
        let y = b"CCCGGG";
        // To get a prefix of y aligned to suffix of x, we might need TrimEnd::Any
        // or a specific mode. TrimEnd::Any sets all clips to 0.
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);
        
        // AAACCC---
        //    |||
        // ---CCCGGG
        
        assert_eq!(l1, "AAACCC   ");
        assert_eq!(l2, "   |||   ");
        assert_eq!(l3, "   CCCGGG");
    }

    #[test]
    fn test_format_alignment_suffix_prefix() {
        let x = b"CCCGGG";
        let y = b"AAACCC";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);
        
        // ---CCCGGG
        //    |||
        // AAACCC---

        assert_eq!(l1, "   CCCGGG");
        assert_eq!(l2, "   |||   ");
        assert_eq!(l3, "AAACCC   ");
    }

    #[test]
    fn test_format_alignment_partial_with_error() {
        let x = b"CCCCCGGG";
        let y = b"AAACCACC";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);

        // ---CCCGGG
        //    |||
        // AAACCC---

        assert_eq!(l1, "   CCCCCGGG");
        assert_eq!(l2, "   ||\\||   ");
        assert_eq!(l3, "AAACCACC   ");
    }

    #[test]
    fn test_format_alignment_partial_with_deletion() {
        let x = b"CAATCCTTTT";
        let y = b"CCCCAACCT";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);

        // ---CCCGGG
        //    |||
        // AAACCC---

        assert_eq!(l1, "   CAATCCTTTT");
        assert_eq!(l2, "   ||| |||   ");
        assert_eq!(l3, "CCCCAA CCT   ");
    }

    #[test]
    fn test_format_alignment_partial_with_insertion() {
        let x = b"CAACCTTTT";
        let y = b"CCCCAATCCT";
        let alignment = run_alignment(TrimEnd::Any, 1, -1, -1, -1, x, y);
        let (l1, l2, l3) = format_alignment(&alignment, x, y);

        // ---CCCGGG
        //    |||
        // AAACCC---

        assert_eq!(l1, "   CAA CCTTTT");
        assert_eq!(l2, "   ||| |||   ");
        assert_eq!(l3, "CCCCAATCCT   ");
    }
}