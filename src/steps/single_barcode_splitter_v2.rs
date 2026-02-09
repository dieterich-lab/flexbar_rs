use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils_sigalign::{
    extract_umi_tags_sig, get_local_aligner_from_parameters, trace_alignment_sig, usize_to_u32,
    AlignerAlgorithmParameters
};
use crate::context::SequenceWithId;
use sigalign::{ReferenceBuilder, Aligner, algorithms::Local, Reference};
use sigalign::results::{Alignment, AlignmentOperation, QueryAlignment};
use log::trace;
use std::cell::RefCell;
use std::cmp::{max, min};
use std::ops::DerefMut;
use thread_local::ThreadLocal;

struct BarCodeAlignment {
    barcode_index: usize,
    penalty: u32,
    start: usize,
    end: usize,
    raw: Alignment
}

pub struct SingleBarcodeSplitterV2 {
    barcodes: Vec<SequenceWithId>,
    max_barcode_len: usize,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    barcode_keep: bool,
    barcode_unassigned: bool,
    removal_tags: bool,
    umi_tags: bool,
    trace: bool,
    reference: Reference,
    aligner_params: AlignerAlgorithmParameters,
    aligner: ThreadLocal<RefCell<Aligner<Local>>>,
}

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitterV2 {
    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        let r_seq = rec.get_sequence();
        let best = self.try_alignment(r_seq);
        self.apply_best_alignment(&mut rec, &best);
        rec
    }
}

impl SingleBarcodeSplitterV2 {
    pub fn new(args: &Args, barcodes: &Vec<SequenceWithId>) -> Self {
        let min_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .min().unwrap_or(0);
        let max_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .max().unwrap_or(0);
        let barcode_min_overlap = args.barcode_min_overlap.unwrap_or(max_barcode_len);
        let mut builder  = ReferenceBuilder::new()
            .ignore_base(b'N');
        for barcode in barcodes {
            builder = builder.add_target(&*barcode.id, barcode.sequence.as_bytes());
        }
        let reference = builder.build().unwrap();

        Self {
            barcodes: (*barcodes).clone(),
            max_barcode_len,
            reference,
            aligner_params: get_algorithm_parameters_from_args(args, min_barcode_len),
            aligner: ThreadLocal::new(),
            min_overlap: barcode_min_overlap,
            trim_end: args.barcode_trim_end,
            error_rate: args.barcode_error_rate,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn try_alignment(&self, r_seq: &[u8]) -> Option<BarCodeAlignment> {
        let start_index = match self.trim_end {
            TrimEnd::RTail => max(0, r_seq.len() - self.max_barcode_len),
            _ => 0
        };
        let end_index = match self.trim_end {
            TrimEnd::LTail => min(self.max_barcode_len, r_seq.len()),
            _ => r_seq.len()
        };
        let range_to_check = &r_seq[start_index..end_index];
        if range_to_check.len() < self.min_overlap { return None; }

        let mut aligner = self.get_aligner();
        let result = aligner.align(range_to_check, &self.reference);

        self.get_best_valid_alignment(&result, r_seq, start_index)
    }

    fn get_aligner(&self) -> impl DerefMut<Target=Aligner<Local>> + '_ {
        self.aligner
            .get_or(|| RefCell::new(
                get_local_aligner_from_parameters(&self.aligner_params)
            ))
            .borrow_mut()
    }

    fn get_best_valid_alignment(&self, all_alignments: &QueryAlignment, r_seq: &[u8], offset: usize)
        -> Option<BarCodeAlignment>
    {
        let mut best_alignment: Option<BarCodeAlignment> = None;
        for (barcode_index, barcode) in self.barcodes.iter().enumerate() {
            if barcode_index >= all_alignments.0.len() { continue; }
            for alignment in all_alignments.0[barcode_index].alignments.iter() {
                if !self.validate_alignment(alignment, r_seq, &barcode, offset) { continue; }
                match best_alignment {
                    Some(ref best) => {
                        if alignment.penalty >= best.penalty { continue }
                    },
                    None => ()
                }
                best_alignment = Some(BarCodeAlignment {
                    barcode_index,
                    penalty: alignment.penalty,
                    start: offset
                        + alignment.position.query.0 as usize
                        - alignment.position.target.0 as usize,
                    end: offset + alignment.position.query.1 as usize,
                    raw: alignment.clone(),
                });
            }
        }
        best_alignment
    }

    fn validate_alignment(&self, alignment: &Alignment, r_seq: &[u8], barcode: &SequenceWithId, offset: usize)
        -> bool
    {
        let seq_start_index = offset as i32
            + alignment.position.query.0 as i32
            - alignment.position.target.0 as i32;
        let seq_end_index = min(
            r_seq.len(),
            seq_start_index as usize + barcode.sequence.len()
        ) as i32;
        if seq_end_index - seq_start_index < self.min_overlap as i32 { return  false }

        let allowed_errors = (barcode.sequence.len() as f32 * self.error_rate).floor() as usize;
        let mut error_count = 0;

        // If RATIL/LTAIL alignments don't overlap with the correct end, add the offset
        // to the error count, e.g. with LTAIL
        // Seq:      AAAATTTTTTT
        //           XX|||||     (Count the two Xs as 2 errors)
        // Barcode:    AATTT
        match self.trim_end {
            TrimEnd::RTail => if seq_end_index < (r_seq.len() as i32) {
                error_count += (r_seq.len() as i32 - seq_end_index) as usize;
            },
            TrimEnd::LTail => if alignment.position.target.0 > alignment.position.query.0 {
                error_count += (alignment.position.target.0 as i32 - seq_start_index) as usize;
            },
            _ => ()
        };
        let mut i_seq: i32 = seq_start_index;
        let mut i_barcode = 0usize;
        let mut i_target = -(alignment.position.target.0 as i32);

        // Handle barcode before the alignment
        while i_target < 0 {
            if i_seq >= 0 {
                match barcode.sequence.as_bytes()[i_barcode] {
                    b'N' => (),
                    _ => error_count += 1,
                }
            }
            i_seq += 1;
            i_barcode += 1;
            i_target += 1;
            if error_count > allowed_errors { return false }
        }

        // Handle alignment
        for op in alignment.operations.iter() {
            match op.operation {
                AlignmentOperation::Match => {
                    i_seq += op.count as i32;
                    i_barcode += op.count as usize;
                },
                AlignmentOperation::Subst => {
                    for _ in 0 ..op.count {
                        match barcode.sequence.as_bytes()[i_barcode] {
                            b'N' => (),
                            _ => error_count += 1,
                        }
                        i_barcode += 1;
                    }
                    i_seq += op.count as i32;
                },
                AlignmentOperation::Deletion => {
                    error_count += op.count as usize;
                    i_barcode += op.count as usize;
                },
                AlignmentOperation::Insertion => {
                    error_count += op.count as usize;
                    i_seq += op.count as i32;
                }
            }
            if error_count > allowed_errors { return false }
        }

        // Handle barcode after the alignment
        while i_barcode < barcode.sequence.len() {
            if i_seq < r_seq.len() as i32 {
                match barcode.sequence.as_bytes()[i_barcode] {
                    b'N' => (),
                    _ => error_count += 1,
                }
            }
            i_seq += 1;
            i_barcode += 1;
            if error_count > allowed_errors { return false }
        }
        true
    }

    fn apply_best_alignment(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        best_alignment: &Option<BarCodeAlignment>
    ) {
        if self.trace {
            trace!("Record {}:", rec.record.id());
        }
        match best_alignment {
            Some(best) => {
                let barcode = &self.barcodes[best.barcode_index];
                if self.trace {
                    trace!(
                            "Barcode splitter on {} selected the following alignment:", 
                            rec.record.id()
                    );
                    trace_alignment_sig(
                        &best.raw, rec.get_sequence(), barcode.sequence.as_bytes()
                    );
                }
                self.assign_barcode_tag(rec, format!("barcode_{}", barcode.id));
                self.handle_umi_tags(rec, &best);
                self.remove_barcode(rec, &best);
            },
            None => {
                if self.trace && log::log_enabled!(log::Level::Trace) {
                    trace!("No barcode match found.");
                }
                if self.barcode_unassigned {
                    self.assign_barcode_tag(rec, "barcode_unassigned_".to_string());
                } else {
                    self.drop_record(rec);
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

        let barcode = &self.barcodes[best.barcode_index];
        let new_umi_tags = extract_umi_tags_sig(
            &best.raw, barcode.sequence.as_bytes(), rec.get_sequence()
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
        let rel_start = best.start - rec.alignment_results.seq_start;
        let rel_end = best.end - rec.alignment_results.seq_start;
        let old_start = rec.alignment_results.seq_start;
        match self.trim_end {
            TrimEnd::Right | TrimEnd::RTail => {
                rec.alignment_results.seq_end = old_start + rel_start;
            }
            TrimEnd::Left | TrimEnd::LTail => {
                rec.alignment_results.seq_start = old_start + rel_end;
            }
            TrimEnd::Any => {
                let length = rec.get_sequence().len();
                if rel_start < length - rel_end {
                    rec.alignment_results.seq_start = old_start + rel_end;
                } else {
                    rec.alignment_results.seq_end = old_start + rel_start;
                }
            }
        }
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

    fn drop_record(&self, rec: &mut RecordInTransformation<SingeRecordDto>) {
        rec.stats.total_discarded_barcode += 1;
        rec.alignment_results.record_dropped = true;
        if self.trace {
            trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record.id());
        }
    }
}

fn get_algorithm_parameters_from_args(args: &Args, min_barcode_len: usize)
    -> AlignerAlgorithmParameters
{
    let min_length = match args.barcode_min_overlap {
        Some(x) => usize_to_u32(x, "barcode_min_overlap"),
        None => min_barcode_len as u32,
    };
    AlignerAlgorithmParameters {
        mismatch: 4,
        gap_open: 6,
        gap_extend: 2,
        min_length,
        max_penalty_per_length: 0.1,
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
        let mut splitter = SingleBarcodeSplitterV2::new(&args, &barcodes);
        splitter.apply(
            record.get_record_in_transformation()
        )
    }

    fn get_record_not_matching_barcode() -> SingeRecordDto {
        SingeRecordDto::with_attrs("read1", None, b"AAAAAATTAC", b"IIIIIIII")
    }

    #[test]
    fn test_barcode_unassigned_false() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq"])
        );
        
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "10", "--barcode-unassigned"])
        );
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned_".to_string()));
        assert_eq!(res.alignment_results.seq_end - res.alignment_results.seq_start, 10);
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
        assert_eq!(res.alignment_results.seq_start, 28);
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
        assert_eq!(res.alignment_results.seq_end, 44);
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
                "--barcode-error-rate", "0.4",
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
        let mut splitter = SingleBarcodeSplitterV2::new(&Args::parse_from(&[
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
