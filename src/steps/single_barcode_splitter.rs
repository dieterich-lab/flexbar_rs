use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::pairwise::*;
use bio::alignment::{Alignment, AlignmentOperation};
use crate::context::SequenceWithId;

pub struct SingleBarcodeSplitter {
    barcodes: Vec<SequenceWithId>,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    gap_open: i32,
    gap_extend: i32,
    barcode_match: i32,
    barcode_mismatch: i32,
    barcode_keep: bool,
    barcode_unassigned: bool,
    removal_tags: bool,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    umi_tags: bool,
    max_barcode_len: usize,
    trace: bool,
}

struct BarCodeAlignment<'a> {
    alignment: Alignment,
    barcode: &'a SequenceWithId,
}

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitter {
    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        let mut aligner = Aligner::with_capacity_and_scoring(
            rec.get_sequence().len(),
            self.max_barcode_len,
            self.get_scoring(),
        );
        let mut best_alignment: Option<BarCodeAlignment> = None;
        for barcode in &self.barcodes {
            let alignment = aligner.semiglobal(
                barcode.sequence.as_bytes(), rec.get_sequence()
            );
            let candidate = BarCodeAlignment { alignment, barcode: &barcode };
            self.update_best_alignment(&mut best_alignment, candidate);
        }
        self.apply_best_alignment(&mut rec, &mut best_alignment);
        rec
    }
}

impl SingleBarcodeSplitter {
    pub fn new(args: &Args, barcodes: &Vec<SequenceWithId>) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match args.barcode_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (MIN_SCORE, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, MIN_SCORE, 0, 0),
        };

        let max_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .max().unwrap_or(0);
        let barcode_min_overlap = args.barcode_min_overlap.unwrap_or(max_barcode_len);

        Self {
            barcodes: (*barcodes).clone(),
            trim_end: args.barcode_trim_end,
            min_overlap: barcode_min_overlap,
            error_rate: args.barcode_error_rate,
            gap_open: args.barcode_gap_open,
            gap_extend: args.barcode_gap_extend,
            barcode_match: args.barcode_match,
            barcode_mismatch: args.barcode_mismatch,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            umi_tags: args.umi_tags,
            max_barcode_len,
            trace: args.trace
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.barcode_match;
        let mm = self.barcode_mismatch;
        Scoring {
            match_fn: move |a, b|
                if a == b'N' || b == b'N' { 0 }
                else if a == b { m }
                else { mm },
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    fn update_best_alignment<'a>(
        &self,
        best_alignment: &mut Option<BarCodeAlignment<'a>>,
        candidate: BarCodeAlignment<'a>,
    ) {
        let overlap = candidate.alignment.yend - candidate.alignment.ystart;
        let barcode_overlap = candidate.alignment.xend - candidate.alignment.xstart;
        if overlap < self.min_overlap || barcode_overlap < self.min_overlap { return; }
        let allowed_errors = (overlap as f32 * self.error_rate) as usize;
        if !self.check_match(
            &candidate.alignment,
            candidate.barcode.sequence.as_ref(),
            overlap,
            allowed_errors
        ) { return; }

        match best_alignment {
            Some(best) => {
                if candidate.alignment.score > best.alignment.score {
                    *best_alignment = Some(candidate);
                }
            },
            None => { *best_alignment = Some(candidate); }
        }
    }

    fn check_match(
        &self, alignment: &Alignment,
        barcode_seq: &[u8],
        overlap: usize,
        allowed_errors: usize
    ) -> bool {
        let mut barcode_pos = alignment.xstart;
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
                    barcode_pos += 1;
                    free_errors -= 1;
                }
                AlignmentOperation::Del => {
                    bases_to_check -= 1;
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
        false
    }

    fn apply_best_alignment(
        &self,
        rec:  &mut RecordInTransformation<SingeRecordDto>,
        best_alignment: &mut Option<BarCodeAlignment>
    ) {
        match best_alignment {
            Some(best) => {
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
        let splitter = SingleBarcodeSplitter::new(&args, &barcodes);
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

    fn get_record_with_with_11_bases_barcode_overlap() -> SingeRecordDto {
        // Scenario: Barcode matches with an overlap of 11 bases.
        // Sequence: CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCT
        // Barcode:                                              NNNNNNTAGCTCTGTAGGCACCATCAAT
        SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:28608:1038",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCT",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII"
        )
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
}
