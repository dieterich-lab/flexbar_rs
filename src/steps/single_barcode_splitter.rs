use crate::args::Args;
use crate::types::{SingeRecordDto, TrimEnd};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup};
use crate::multi_aligner::result::MultiAlignerResult;

use std::any::Any;

pub struct SingleBarcodeSplitter {
    multi_aligner: MultiAligner,
}

#[derive(Clone)]
pub struct SingleBarcodeSplitterSetup {
    pub barcodes: Vec<SequenceWithId>,
    pub multi_aligner_setup: MultiAlignerSetup,
    pub trim_end: TrimEnd,
    pub barcode_keep: bool,
    pub barcode_unassigned: bool,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitter {
    type Setup = SingleBarcodeSplitterSetup;

    fn create_setup(args: &Args, context: &Context) -> SingleBarcodeSplitterSetup {
        let barcodes = context.barcodes.clone();
        let multi_aligner_setup = MultiAlignerSetup::new(
            "Barcode Splitter".to_string(),
            args.barcode_trim_end,
            barcodes.clone(),
            crate::multi_aligner::AlignerConfig {
                match_score: args.barcode_match,
                mismatch_score: args.barcode_mismatch,
                gap_open: args.barcode_gap_open,
                gap_extend: args.barcode_gap_extend,
                error_rate: args.barcode_error_rate,
                min_overlap: args.barcode_min_overlap,
            },
            args.barcode_kmer
        );

        SingleBarcodeSplitterSetup {
            barcodes,
            multi_aligner_setup,
            trim_end: args.barcode_trim_end,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn new(setup: SingleBarcodeSplitterSetup) -> Self {
        Self {
            multi_aligner: setup.multi_aligner_setup.create_aligner(),
        }
    }

    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &SingleBarcodeSplitterSetup) -> RecordInTransformation<SingeRecordDto> {
        let sequence = rec.get_sequence();
        let best_alignment = self.multi_aligner.align(&setup.multi_aligner_setup, &sequence);
        
        match best_alignment {
            Some(best) => {
                let barcode = &setup.barcodes[best.reference_index];
                if setup.trace {
                    let bar_seq = barcode.sequence.as_bytes();
                    trace!(
                            "Barcode splitter on {} selected the following alignment:", 
                            rec.record.id()
                        );
                    trace!("{}", best.format(sequence, bar_seq));
                }
                self.assign_barcode_tag(&mut rec, format!("barcode_{}", barcode.id), setup);
                self.handle_umi_tags(&mut rec, &best, barcode, setup);
                self.remove_barcode(&mut rec, &best, setup);
            },
            None => {
                if setup.barcode_unassigned {
                    self.assign_barcode_tag(&mut rec, "barcode_unassigned_".to_string(), setup);
                } else {
                    Self::drop_record(&mut rec);
                }
            }
        }
        rec
    }
}

impl<'a> SingleBarcodeSplitter {
    pub fn new(setup: SingleBarcodeSplitterSetup) -> Self {
        Self {
            multi_aligner: setup.multi_aligner_setup.create_aligner(),
        }
    }

    fn assign_barcode_tag(&self, rec: &mut RecordInTransformation<SingeRecordDto>, tag: String, setup: &SingleBarcodeSplitterSetup) {
        if setup.trace {
            trace!("Barcode splitter on {}: Assigned tag {}", rec.record.id(), tag);
        }
        rec.alignment_results.barcode_tag = Some(tag);
    }

    fn handle_umi_tags(
        &self, 
        rec: &mut RecordInTransformation<SingeRecordDto>, 
        best: &MultiAlignerResult,
        barcode: &SequenceWithId,
        setup: &SingleBarcodeSplitterSetup,
    ) {
        if !setup.umi_tags { return; }
        let new_umi_tags = best.get_umi_tags(
            rec.get_sequence(), barcode.sequence.as_bytes()
        );
        if setup.trace {
            trace!(
                "Barcode splitter on {}: Extracted UMI tags {:?}",
                rec.record.id(), new_umi_tags.clone()
            );
        }
        rec.alignment_results.umi_tags.extend(new_umi_tags);
    }

    fn remove_barcode(
        &self, 
        rec: &mut RecordInTransformation<SingeRecordDto>, 
        best: &MultiAlignerResult,
        setup: &SingleBarcodeSplitterSetup,
    ) {
        if setup.barcode_keep { return; }
        rec.trim(setup.trim_end, best.get_start_index(), best.get_end_index());
        if setup.removal_tags {
            rec.alignment_results.removal_tags.push("Flexbar_removal".to_string());
        }
        if setup.trace {
            trace!(
                "Barcode splitter on {} shortened sequence to:\n{}",
                rec.record.id(), rec.get_sequence_as_str()
            );
            if setup.removal_tags {
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
        let context = Context {
            barcodes: barcodes.clone(),
            ..Context::empty()
        };
        let setup = SingleBarcodeSplitter::create_setup(&args, &context);
        let mut splitter = SingleBarcodeSplitter::new(setup.clone());
        splitter.apply(
            record.get_record_in_transformation(),
            &setup
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
        
        let args = Args::parse_from(&[
            "flexbar_rs", "--reads", "test.fq", 
            "--barcode-trim-end", "RIGHT", 
            "--barcode-min-overlap", "6"
        ]);
        let barcodes = vec![
            SequenceWithId {
                id: "bc1".to_string(),
                sequence: "NNNTTTTTCCCC".to_string()
            }
        ];
        let context = Context {
            barcodes: barcodes.clone(),
            ..Context::empty()
        };
        let setup = SingleBarcodeSplitter::create_setup(&args, &context);
        let mut splitter = SingleBarcodeSplitter::new(setup.clone());
        
        let rec = SingeRecordDto::with_attrs("read1", None, b"AAAAATTTTTAAATTT", b"IIIIIIIIIIIIIIII");
        let res = splitter.apply(rec.get_record_in_transformation(), &setup);
        
        // In this case:
        // Primary multi_aligner (RIGHT) finds:
        // AAAAATTTTTAAATTT
        //   ||||||||
        //   NNNTTTTTCCCC
        // overlap = 8 (TTTTTAAA vs TTTTTCCC) -> 3 mismatches.
        // allowed_errors (error_rate=0.1) = 0.
        // check_match fails.
        
        // Tail multi_aligner (RTAIL) finds:
        // AAAAATTTTTAAATTT
        //           ||||||
        //           NNNTTTTTCCCC
        // overlap = 6 (TTT vs NNN) -> 0 errors (N matches anything).
        // check_match succeeds!
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
    }
}
