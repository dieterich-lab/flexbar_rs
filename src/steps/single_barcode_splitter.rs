use std::cell::RefMut;
use std::cmp::min;
use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::{Alignment, AlignmentOperation};
use bio::alignment::pairwise::Aligner;
use crate::context::SequenceWithId;
use crate::utils::aligner::{trace_alignment, AlignerFactory};

pub struct SingleBarcodeSplitter {
    barcodes: Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    aligner_factory: AlignerFactory,
    tail_aligner_factory: Option<AlignerFactory>,
    barcode_keep: bool,
    barcode_unassigned: bool,
    removal_tags: bool,
    umi_tags: bool,
    max_barcode_len: usize,
    trace: bool,
}

struct BarCodeAlignment<'a> {
    alignment: Alignment,
    barcode: &'a SequenceWithId,
}

type BoxedAligner<'a> = RefMut<'a, Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>>;

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitter {
    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        let mut best_alignment: Option<BarCodeAlignment> = None;
        let r_seq = rec.get_sequence();
        let seq_len = r_seq.len();
        let mut aligners = self.build_aligners(seq_len);

        for barcode in &self.barcodes {
            match self.try_alignment(&mut aligners, r_seq, barcode) {
                Some(candidate) => {
                    match best_alignment {
                        Some(ref best) => {
                            if candidate.alignment.score > best.alignment.score {
                                best_alignment = Some(candidate);
                            }
                        },
                        None => { best_alignment = Some(candidate); }
                    }
                }
                None => ()
            }
        }

        self.apply_best_alignment(&mut rec, &mut best_alignment);
        rec
    }
}

impl SingleBarcodeSplitter {
    pub fn new(args: &Args, barcodes: &Vec<SequenceWithId>) -> Self {
        let max_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .max().unwrap_or(0);
        let barcode_min_overlap = args.barcode_min_overlap.unwrap_or(max_barcode_len);

        let aligner_factory = AlignerFactory::new(
            args.barcode_trim_end,
            args.barcode_match,
            args.barcode_mismatch,
            args.barcode_gap_open,
            args.barcode_gap_extend,
        );

        // In the case of a trim end of RIGHT or LEFT, we might find an alignment
        // in the middle of the sequence, which has a good score but is invalid.
        // That may prevent the aligner from finding a valid alignment at the end
        // with the same or a lower score. To prevent that, we check again with
        // the 'TAIL' version of the aligner specifically for alignments at the end.
        let tail_aligner_factory = match args.barcode_trim_end {
            TrimEnd::Right => Some(AlignerFactory::new(
                TrimEnd::RTail,
                args.barcode_match,
                args.barcode_mismatch,
                args.barcode_gap_open,
                args.barcode_gap_extend,
            )),
            TrimEnd::Left => Some(AlignerFactory::new(
                TrimEnd::LTail,
                args.barcode_match,
                args.barcode_mismatch,
                args.barcode_gap_open,
                args.barcode_gap_extend,
            )),
            _ => None,
        };

        Self {
            barcodes: (*barcodes).clone(),
            trim_end: args.barcode_trim_end,
            min_overlap: barcode_min_overlap,
            error_rate: args.barcode_error_rate,
            aligner_factory,
            tail_aligner_factory,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            max_barcode_len,
            trace: args.trace,
        }
    }

    fn build_aligners(&self, sequence_len: usize)
        -> Vec<RefMut<Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>>>
    {
        let aligner = self.aligner_factory.get(sequence_len, self.max_barcode_len);
        let tail_aligner = self.tail_aligner_factory.as_ref().map(
            |f| f.get(sequence_len, self.max_barcode_len)
        );
        match tail_aligner {
            None => vec![aligner],
            Some(extra) => vec![aligner, extra],
        }
    }

    fn try_alignment<'a>(
        &self,
        aligners: &mut [BoxedAligner],
        r_seq: &[u8],
        barcode: &'a SequenceWithId
    )
        -> Option<BarCodeAlignment<'a>> 
    {
        let b_seq = barcode.sequence.as_bytes();
        for aligner in aligners.iter_mut() {
            let alignment = aligner.custom(r_seq, b_seq);
            if self.validate_alignment(&alignment, barcode, r_seq.len()) {
                return Some(BarCodeAlignment { alignment, barcode });
            }
        }
        None
    }

    fn validate_alignment(
        &self,
        alignment: &Alignment,
        barcode: &SequenceWithId,
        sequence_len: usize,
    ) -> bool {
        let overlap = alignment.xend - alignment.xstart;
        let barcode_overlap = alignment.yend - alignment.ystart;
        if overlap < self.min_overlap || barcode_overlap < self.min_overlap {
            if self.trace { trace!("Skipping alignment: Overlap too short.") }
            return false;
        }
        let allowed_errors = (overlap as f32 * self.error_rate) as usize;
        let required_overlap = min(
            sequence_len - alignment.xstart,
            barcode.sequence.len()
        );
        if !self.check_match(
            alignment,
            barcode.sequence.as_ref(),
            required_overlap,
            allowed_errors
        ) {
            if self.trace { trace!("Skipping alignment: Too many errors.") }
            return false;
        }
        true
    }

    fn check_match(
        &self, alignment: &Alignment,
        barcode_seq: &[u8],
        overlap: usize,
        allowed_errors: usize
    ) -> bool {
        let mut barcode_pos = alignment.ystart;
        let mut bases_to_check = overlap;
        let mut free_errors: i32 = allowed_errors as i32;

        for op in &alignment.operations {
            match op {
                AlignmentOperation::Match => {
                    barcode_pos += 1;
                    bases_to_check -= 1;
                }
                AlignmentOperation::Subst => {
                    if barcode_seq[barcode_pos] != b'N' {
                        free_errors -= 1;
                    }
                    barcode_pos += 1;
                    bases_to_check -= 1;
                }
                AlignmentOperation::Ins => {
                    bases_to_check -= 1;
                    free_errors -= 1;
                }
                AlignmentOperation::Del => {
                    barcode_pos += 1;
                    free_errors -= 1;
                }
                _ => {}
            }
            if free_errors < 0 {
                return false
            }
            if bases_to_check <= 0 {
                return true
            }
        }
        bases_to_check <= free_errors as usize
    }

    fn apply_best_alignment(
        &self,
        rec:  &mut RecordInTransformation<SingeRecordDto>,
        best_alignment: &mut Option<BarCodeAlignment>
    ) {
        match best_alignment {
            Some(best) => {
                if self.trace {
                    trace!(
                            "Barcode splitter on {} selected the following alignment:", 
                            rec.record.id()
                        );
                    trace_alignment(
                        &best.alignment, rec.get_sequence(),best.barcode.sequence.as_bytes()
                    );
                }
                self.assign_barcode_tag(rec, format!("barcode_{}", best.barcode.id));
                self.handle_umi_tags(rec, &best);
                self.remove_barcode(rec, &best);
            },
            None => {
                if self.barcode_unassigned {
                    self.assign_barcode_tag(rec, "barcode_unassigned_".to_string());
                } else {
                    Self::drop_record(rec);
                }
            }

        }
    }

    fn assign_barcode_tag(&self, rec: &mut RecordInTransformation<SingeRecordDto>, tag: String) {
        if self.trace {
            trace!("Barcode splitter on {}: Assigned tag {}", rec.record.id(), tag);
        }
        rec.alignment_results.barcode_tag = Some(tag);
    }

    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        best: &BarCodeAlignment
    ) {
        if !self.umi_tags { return; }
        let new_umi_tags = extract_umi_tags(
            &best.alignment, best.barcode.sequence.as_bytes(), rec.get_sequence()
        );
        if self.trace {
            trace!(
                "Barcode splitter on {}: Extracted UMI tags {:?}",
                rec.record.id(), new_umi_tags.clone()
            );
        }
        rec.alignment_results.umi_tags.extend(new_umi_tags);
    }

    fn remove_barcode(&self, rec: &mut RecordInTransformation<SingeRecordDto>, best: &BarCodeAlignment) {
        if self.barcode_keep { return; }
        let (rel_start, rel_end) = get_trim_range(
            rec.get_sequence().len(), &best.alignment, self.trim_end
        );
        let old_start = rec.alignment_results.seq_start;
        rec.alignment_results.seq_start = old_start + rel_start;
        rec.alignment_results.seq_end = old_start + rel_end;
        if self.removal_tags {
            rec.alignment_results.removal_tags.push("Flexbar_removal".to_string());
        }
        if self.trace {
            trace!(
                "Barcode splitter on {} shortened sequence to:\n{}",
                rec.record.id(), rec.get_sequence_as_str()
            );
            if self.removal_tags {
                trace!("Barcode splitter on {}: Added removal tag", rec.record.id());
            }
        }
    }

    fn drop_record(rec: &mut RecordInTransformation<SingeRecordDto>) {
        rec.stats.total_discarded_barcode += 1;
        rec.alignment_results.record_dropped = true;
        trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record.id());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::steps::generic::Transformable;

    fn process_record(record: SingeRecordDto, args: Args) -> RecordInTransformation<SingeRecordDto>
    {
        let barcodes = vec![
            SequenceWithId {
                id: "bc1".to_string(),
                sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string()
            }
        ];
        let mut splitter = SingleBarcodeSplitter::new(&args, &barcodes);
        splitter.apply(
            record.get_record_in_transformation()
        )
    }

    fn get_record_not_matching_barcode() -> SingeRecordDto {
        SingeRecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII")
    }

    #[test]
    fn test_barcode_unassigned_false() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4"])
        );
        
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4", "--barcode-unassigned"])
        );
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned_".to_string()));
        assert_eq!(res.alignment_results.seq_end - res.alignment_results.seq_start, 8);
    }

    fn get_record_with_barcode_at_end() -> SingeRecordDto {
        SingeRecordDto::with_attrs(
            "read1",
            None,
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAAT",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIIII"
        )
    }

    #[test]
    fn test_barcode_trimming_right() {
        let res = process_record(
            get_record_with_barcode_at_end(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 44);
        assert_eq!(res.get_sequence_as_str(), "CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    fn get_record_with_barcode_at_start() -> SingeRecordDto {
        SingeRecordDto::with_attrs(
            "read1",
            None,
            b"AAAAAATAGCTCTGTAGGCACCATCAATCATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIIII"
        )
    }

    #[test]
    fn test_barcode_trimming_left() {
        let res = process_record(
            get_record_with_barcode_at_start(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "LEFT",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 28);
        assert_eq!(res.alignment_results.seq_end, 72);
        assert_eq!(res.get_sequence_as_str(), "CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    #[test]
    fn test_barcode_trimming_any_right() {
        let res = process_record(
            get_record_with_barcode_at_end(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "ANY",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 44);
    }

    #[test]
    fn test_barcode_trimming_any_left() {
        let res = process_record(
            get_record_with_barcode_at_start(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "ANY",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 28);
        assert_eq!(res.alignment_results.seq_end, 72);
    }

    #[test]
    fn test_barcode_trimming_ltail() {
        let res = process_record(
            get_record_with_barcode_at_start(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "LTAIL",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 72);
    }

    #[test]
    fn test_barcode_trimming_rtail() {
        let res = process_record(
            get_record_with_barcode_at_end(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "RTAIL",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 72);
    }

    fn get_record_with_last_barcode_base_missing_in_sequence() -> SingeRecordDto {
        // Scenario: Barcode matches with one deletion in the sequence.
        // Sequence: CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA
        // Barcode:                                              NNNNNNTAGCTCTGTAGGCACCATCAAT
        SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:28608:1038",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII"
        )
    }

    #[test]
    fn test_with_last_barcode_base_missing_in_sequence() {
        let res = process_record(
            get_record_with_last_barcode_base_missing_in_sequence(),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    fn get_record_with_with_13_bases_barcode_overlap() -> SingeRecordDto {
        // Scenario: Barcode matches with an overlap of 11 bases.
        // Sequence: CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCT
        // Barcode:                                              NNNNNNTAGCTCTGTAGGCACCATCAAT
        SingeRecordDto::with_attrs(
            "generic",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCT",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIII"
        )
    }

    #[test]
    fn test_record_with_with_13_bases_barcode_overlap() {
        let res = process_record(
            get_record_with_with_13_bases_barcode_overlap(),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags",
            ])
        );
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    fn get_record_with_a_deletion_in_the_barcode() -> SingeRecordDto {
        // Scenario: Barcode matches with one deletion in the sequence.
        // Sequence: GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCTCCCTGGTAGCTCTG AGGCACCATCAATAGATCGGAA
        // Barcode:                                     NNNNNNTAGCTCTGTAGGCACCATCAAT
        SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:37924:1000",
            None,
            b"GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCTCCCTGGTAGCTCTGAGGCACCATCAATAGATCGGAA",
            b"99I999-99-9III--I9IIIIIIIIIIIIIII9IIIIIIIIIIII9IIIIIIIII9III-IIIII9II9I"
        )
    }

    #[test]
    fn test_record_error_no_tolerance() {
        let res = process_record(
            get_record_with_a_deletion_in_the_barcode(),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(res.alignment_results.record_dropped);
    }

    #[test]
    fn test_record_error_with_tolerance() {
        let res = process_record(
            get_record_with_a_deletion_in_the_barcode(),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--barcode-error-rate", "0.1",
                "--umi-tags"
            ])
        );
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(), b"GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCT");
    }

    fn get_partial_overlap_record() -> SingeRecordDto {
        // Scenario: Barcode matches with more than "barcode-min-overlap" bases (11),
        // but has a deletion later in the sequence.
        // Seq:     GATAGCAACTCTAGCATCACCAGTATTAGAGGCTTTCTGTAGCTCTGTAGGCAC-ATCAATAGATCGGAAGA
        // W/o bar: GATAGCAACTCTAGCATCACCAGTATTAGAGGC
        // Matching:                                 |||||||||||||||||||||-||||||
        // Barcode:                                  NNNNNNTAGCTCTGTAGGCACCATCAAT
        SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:39591:1000",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"GATAGCAACTCTAGCATCACCAGTATTAGAGGCTTTCTGTAGCTCTGTAGGCACATCAATAGATCGGAAGA",
            b"I9I9III9IIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        )
    }

    #[test]
    fn test_partial_overlap_record_no_tolerance() {
        let res = process_record(
            get_partial_overlap_record(),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(res.alignment_results.record_dropped);
    }

    #[test]
    fn test_two_pass_alignment_fallback() {
        // This test specifically targets the case where RIGHT fails check_match
        // but RTAIL might succeed.
        
        // Sequence: GATAGCAACTCTAGCATCACCAGTATTAGAGGCTTTCTGTAGCTCTGTAGGCACATCAATAGATCGGAAGA
        // Barcode:                                  NNNNNNTAGCTCTGTAGGCACCATCAAT
        
        // With RIGHT, it might find a high scoring alignment that is NOT at the end,
        // but if it's invalid (e.g. because of check_match), it should fall back to RTAIL.
        
        // Let's use a case where there is a good-ish match in the middle but it fails check_match,
        // and an exact match at the end that RTAIL will find.
        
        // Actually, the existing test_partial_overlap_record_no_tolerance fails because
        // the alignment found by RIGHT (which includes the mismatch/deletion) is invalid.
        // If our fix works, it should now PASS if RTAIL can find a valid alignment.
        // Wait, in get_partial_overlap_record():
        // Seq:     GATAGCAACTCTAGCATCACCAGTATTAGAGGCTTTCTGTAGCTCTGTAGGCAC-ATCAATAGATCGGAAGA
        // Barcode:                                  NNNNNNTAGCTCTGTAGGCACCATCAAT
        // There is a '-' (deletion) in the sequence relative to barcode.
        // If we use RIGHT, it finds this alignment. check_match fails.
        // Then we try RTAIL. RTAIL also finds the same alignment because it's at the end anyway.
        // So it still fails.
        
        // To test the fallback, we need a case where RIGHT finds a DIFFERENT (invalid) alignment
        // than RTAIL (which finds a valid one).
        
        let barcodes = vec![
            SequenceWithId {
                id: "bc1".to_string(),
                sequence: "NNNTTTTTCCCC".to_string()
            }
        ];
        let mut splitter = SingleBarcodeSplitter::new(&Args::parse_from(&[
            "flexbar_rs", "--reads", "test.fq", 
            "--barcode-trim-end", "RIGHT", 
            "--barcode-min-overlap", "6"
        ]), &barcodes);
        
        let rec = SingeRecordDto::with_attrs("read1", None, b"AAAAATTTTTAAATTT", b"IIIIIIIIIIIIIIII");
        let res = splitter.apply(rec.get_record_in_transformation());
        
        // In this case:
        // Primary aligner (RIGHT) finds:
        // AAAAATTTTTAAATTT
        //   ||||||||
        //   NNNTTTTTCCCC
        // overlap = 8 (TTTTTAAA vs TTTTTCCC) -> 3 mismatches.
        // allowed_errors (error_rate=0.1) = 0.
        // check_match fails.
        
        // Tail aligner (RTAIL) finds:
        // AAAAATTTTTAAATTT
        //           ||||||
        //           NNNTTTTTCCCC
        // overlap = 6 (TTT vs NNN) -> 0 errors (N matches anything).
        // check_match succeeds!
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
    }
}
