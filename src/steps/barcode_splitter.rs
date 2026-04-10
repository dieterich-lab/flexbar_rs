use crate::args::Args;
use crate::types::{TrimEnd, RunType};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{AlignerConfig, MultiAligner, MultiAlignerSetup, MultiAlignerState};
use crate::multi_aligner::result::MultiAlignerResult;

pub struct BarcodeSplitter {
    pub barcodes: [Vec<SequenceWithId>; 2],
    pub multi_aligners: [Option<MultiAligner>; 2],
    pub multi_aligner_setups: [Option<MultiAlignerSetup>; 2],
    pub trim_end: TrimEnd,
    pub barcode_keep: bool,
    pub barcode_unassigned: bool,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

pub struct BarcodeSplitterState {
    pub multi_aligner_states: [Option<MultiAlignerState>; 2],
}

impl BarcodeSplitter {
    fn create_state_by_index(&self, read_idx: usize) -> Option<MultiAlignerState> {
        match (&self.multi_aligner_setups[read_idx], &self.multi_aligners[read_idx]) {
            (Some(setup), Some(multi_aligner))
                => Some(setup.create_state(multi_aligner)),
            _ => None,
        }
    }
}

fn get_aligner_setup(name: &str, args: &Args, barcodes: &[SequenceWithId]) -> MultiAlignerSetup {
    MultiAlignerSetup::new(
        name.to_string(),
        args.barcode_trim_end,
        barcodes,
        get_barcode_aligner_config(args),
        args.barcode_kmer
    )
}

fn get_barcode_aligner_config(args: &Args) -> AlignerConfig {
    AlignerConfig {
        match_score: args.barcode_match,
        mismatch_score: args.barcode_mismatch,
        gap_open: args.barcode_gap_open,
        gap_extend: args.barcode_gap_extend,
        error_rate: args.barcode_error_rate,
        min_overlap: args.barcode_min_overlap,
        cycles: if args.barcode_keep { 1 } else { args.barcode_cycles },
    }
}

pub fn assign_barcode_tag(
    rec: &mut RecordInTransformation,
    tag: String,
    trace: bool
) {
    if trace {
        trace!("Barcode splitter on {}: Assigned tag {}", rec.record_id(), tag);
    }
    rec.alignment_results.set_barcode_tag(Some(tag));
}

impl GenericStep for BarcodeSplitter {
    type State = BarcodeSplitterState;

    fn new(args: &Args, context: &Context) -> Self {
        let mut barcodes: [Vec<SequenceWithId>; 2] = [Vec::new(), Vec::new()];
        let mut multi_aligner_setups: [Option<MultiAlignerSetup>; 2] = [None, None];
        match context.run_type {
            RunType::Single | RunType::SingleBarcoded => {
                barcodes[0] = context.barcodes.clone();
                if !barcodes[0].is_empty() {
                    multi_aligner_setups[0] = Some(
                        get_aligner_setup("Barcode Splitter", args, &barcodes[0])
                    )
                }
            },
            RunType::Paired | RunType::PairedBarcoded => {
                barcodes[0] = context.barcodes.clone();
                barcodes[1] = context.barcodes2.clone();
                for i in 0..2 {
                    if !barcodes[i].is_empty() {
                        multi_aligner_setups[i] = Some(
                            get_aligner_setup(&format!("Barcode Splitter R{}", i + 1), args, &barcodes[i])
                        )
                    }
                }
            },
        }

        let mut multi_aligners: [Option<MultiAligner>; 2] = [None, None];
        for i in 0..2 {
            if let Some(setup) = &multi_aligner_setups[i] {
                multi_aligners[i] = Some(setup.create_aligner());
            }
        }

        Self {
            barcodes,
            multi_aligners,
            multi_aligner_setups,
            trim_end: args.barcode_trim_end,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn create_state(&self) -> Self::State {
        Self::State {
            multi_aligner_states: [
                self.create_state_by_index(0),
                self.create_state_by_index(1),
            ]
        }
    }

    fn apply(&self, mut rec: RecordInTransformation, state: &mut Self::State) -> RecordInTransformation {
        let mut best_overall: Option<(usize, MultiAlignerResult)> = None;

        if let Some(barcode_read) = &rec.record.barcode_read {
            // Use separate barcode read for detection
            if let Some(setup) = &self.multi_aligner_setups[0] {
                let aligner = self.multi_aligners[0].as_ref().unwrap();
                let (tag, bar_seq) = {
                    let seq = barcode_read.seq();
                    let inner_state = state.multi_aligner_states[0].as_mut().unwrap();
                    if let Some(best) = aligner.align(inner_state, &setup, seq) {
                        let barcode = &self.barcodes[0][best.reference_index];
                        let tag = format!("barcode_{}", barcode.id);
                        let bar_seq = if self.trace { Some((barcode.sequence.as_bytes().to_vec(), best.format(seq, barcode.sequence.as_bytes()))) } else { None };
                        (Some(tag), bar_seq)
                    } else {
                        (None, None)
                    }
                };

                if let Some(tag) = tag {
                    assign_barcode_tag(&mut rec, tag, self.trace);
                    if let Some((_, formatted)) = bar_seq {
                        trace!(
                            "Barcode splitter on {} selected the following alignment from Barcode Read:",
                            rec.record_id()
                        );
                        trace!("{}", formatted);
                    }
                    // No barcode removal for separate barcode reads
                    return rec;
                }
            }
        } else {
            let num_reads = if rec.record.reads[1].is_some() { 2 } else { 1 };

            for read_idx in 0..num_reads {
                if self.multi_aligner_setups[read_idx].is_none() { continue; }
                let setup = self.multi_aligner_setups[read_idx].as_ref().unwrap();
                let aligner = self.multi_aligners[read_idx].as_ref().unwrap();
                let mut state = state.multi_aligner_states[read_idx].as_mut().unwrap();
                let seq = rec.get_sequence(read_idx);
                if let Some(best) = aligner.align(&mut state, &setup, &seq) {
                    match &best_overall {
                        Some((_, current_best)) if best.score() > current_best.score() => {
                            best_overall = Some((read_idx, best));
                        }
                        None => {
                            best_overall = Some((read_idx, best));
                        }
                        _ => {}
                    }
                }
            }
        }

        match best_overall {
            Some((read_idx, best)) => {
                let barcode = &self.barcodes[read_idx][best.reference_index];
                let tag = format!("barcode_{}", barcode.id);
                assign_barcode_tag(&mut rec, tag, self.trace);
                if self.trace {
                    let bar_seq = barcode.sequence.as_bytes();
                    let sequence: &[u8] = rec.get_sequence(read_idx);
                    trace!(
                        "Barcode splitter on {} selected the following alignment from Read {}:",
                        rec.record_id(),
                        read_idx + 1
                    );
                    trace!("{}", best.format(sequence, bar_seq));
                }
                self.handle_umi_tags(&mut rec, &best, barcode, read_idx);
                self.remove_barcode(&mut rec, &best, read_idx);
            }
            None => {
                self.handle_unassigned(&mut rec);
            }
        }

        rec
    }
}

impl BarcodeSplitter {
    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation,
        best: &MultiAlignerResult,
        barcode: &SequenceWithId,
        read_idx: usize,
    ) {
        if !self.umi_tags { return; }
        
        let sequence: &[u8] = rec.get_sequence(read_idx);
        let bar_seq = barcode.sequence.as_bytes();
        let new_umi_tags = best.get_umi_tags(sequence, bar_seq);
        
        if self.trace {
            trace!(
                "Barcode splitter on {}: Extracted UMI tags {:?}",
                rec.record_id(), new_umi_tags.clone()
            );
        }
        rec.alignment_results.per_read[read_idx].as_mut().unwrap().umi_tags.extend(new_umi_tags);
    }

    fn remove_barcode(
        &self,
        rec: &mut RecordInTransformation,
        best: &MultiAlignerResult,
        read_idx: usize,
    ) {
        if self.barcode_keep { return; }
        
        let old_len = rec.get_read_length(read_idx);
        rec.trim(read_idx, self.trim_end, best.get_start_index(), best.get_end_index());
        let new_len = rec.get_read_length(read_idx);
        rec.alignment_results.per_read[read_idx].as_mut().unwrap().bases_barcode_removal += old_len - new_len;

        if self.removal_tags {
            rec.alignment_results.per_read[read_idx].as_mut().unwrap().removal_tags.push("Flexbar_removal".to_string());
        }

        if self.trace {
            trace!(
                "Barcode splitter on {} shortened sequence of Read {} to length {}",
                rec.record_id(), read_idx + 1, new_len
            );
        }
    }

    fn handle_unassigned(&self, rec: &mut RecordInTransformation) {
        if self.barcode_unassigned {
            assign_barcode_tag(rec, "barcode_unassigned".to_string(), self.trace);
        } else {
            let num_reads = if rec.record.reads[1].is_some() { 2 } else { 1 };
            for read_idx in 0..num_reads {
                if !rec.alignment_results.per_read[read_idx].as_ref().unwrap().dropped {
                    rec.alignment_results.per_read[read_idx].as_mut().unwrap().dropped = true;
                    rec.stats.total_discarded_barcode += 1;
                }
            }
            if self.trace {
                trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record_id());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::types::{RecordDto, SingeRecordDto};

    fn process_record(record: RecordDto, args: Args) -> RecordInTransformation
    {
        let barcodes = vec![
            SequenceWithId {
                id: "bc1".to_string(),
                sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string()
            }
        ];
        let context = Context {
            barcodes: barcodes.clone(),
            run_type: if record.reads[1].is_none() { RunType::Single } else { RunType::Paired },
            ..Context::empty()
        };
        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();
        step.apply(
            RecordInTransformation::from_dto(record),
            &mut state
        )
    }

    fn get_record_not_matching_barcode_single() -> RecordDto {
        RecordDto::single(SingeRecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII"))
    }

    #[test]
    fn test_barcode_unassigned_false_single() {
        let res = process_record(
            get_record_not_matching_barcode_single(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4"])
        );
        
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true_single() {
        let res = process_record(
            get_record_not_matching_barcode_single(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4", "--barcode-unassigned"])
        );
        
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end - res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 8);
    }

    fn get_record_with_barcode_at_end_single() -> RecordDto {
        RecordDto::single(SingeRecordDto::with_attrs(
            "read1",
            None,
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAAT",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIIII"
        ))
    }

    #[test]
    fn test_barcode_trimming_right_single() {
        let res = process_record(
            get_record_with_barcode_at_end_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 0);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 44);
        assert_eq!(res.get_sequence(0), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    fn get_record_with_barcode_at_start_single() -> RecordDto {
        RecordDto::single(SingeRecordDto::with_attrs(
            "read1",
            None,
            b"AAAAAATAGCTCTGTAGGCACCATCAATCATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIIII"
        ))
    }

    #[test]
    fn test_barcode_trimming_left_single() {
        let res = process_record(
            get_record_with_barcode_at_start_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "LEFT",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 28);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 72);
        assert_eq!(res.get_sequence(0), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    #[test]
    fn test_barcode_trimming_any_right_single() {
        let res = process_record(
            get_record_with_barcode_at_end_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "ANY",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 0);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 44);
    }

    #[test]
    fn test_barcode_trimming_any_left_single() {
        let res = process_record(
            get_record_with_barcode_at_start_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "ANY",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 28);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 72);
    }

    #[test]
    fn test_barcode_trimming_ltail_single() {
        let res = process_record(
            get_record_with_barcode_at_start_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "LTAIL",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 28);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 72);
    }

    #[test]
    fn test_barcode_trimming_rtail_single() {
        let res = process_record(
            get_record_with_barcode_at_end_single(),
            Args::parse_from(&[
                "flexbar_rs", 
                "--reads", "test.fq", 
                "--barcode-trim-end", "RTAIL",
                "--barcode-min-overlap", "10"
            ])
        );
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.start, 0);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range.end, 44);
    }

    #[test]
    fn test_with_last_barcode_base_missing_in_sequence_single() {
        let r = SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:28608:1038",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII"
        );
        let res = process_record(
            RecordDto::single(r),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(0), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    #[test]
    fn test_record_with_with_13_bases_barcode_overlap_single() {
        let r = SingeRecordDto::with_attrs(
            "generic",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCT",
            b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIII"
        );
        let res = process_record(
            RecordDto::single(r),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags",
            ])
        );
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(0), b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGG");
    }

    #[test]
    fn test_record_error_no_tolerance_single() {
        let r = SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:37924:1000",
            None,
            b"GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCTCCCTGGTAGCTCTGAGGCACCATCAATAGATCGGAA",
            b"99I999-99-9III--I9IIIIIIIIIIIIIII9IIIIIIIIIIII9IIIIIIIII9III-IIIII9II9I"
        );
        let res = process_record(
            RecordDto::single(r),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_record_error_with_tolerance_single() {
        let r = SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:37924:1000",
            None,
            b"GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCTCCCTGGTAGCTCTGAGGCACCATCAATAGATCGGAA",
            b"99I999-99-9III--I9IIIIIIIIIIIIIII9IIIIIIIIIIII9IIIIIIIII9III-IIIII9II9I"
        );
        let res = process_record(
            RecordDto::single(r),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--barcode-error-rate", "0.1",
                "--umi-tags"
            ])
        );
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.get_sequence(0), b"GACGATCACGGGGCGACGCCTGTCTGAGCGTCGCT");
    }

    #[test]
    fn test_partial_overlap_record_no_tolerance_single() {
        let r = SingeRecordDto::with_attrs(
            "VH00748:119:AAGJ22MM5:1:1101:39591:1000",
            Some("1:N:0:TACGGTCT+GCAATGGA"),
            b"GATAGCAACTCTAGCATCACCAGTATTAGAGGCTTTCTGTAGCTCTGTAGGCACATCAATAGATCGGAAGA",
            b"I9I9III9IIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        );
        let res = process_record(
            RecordDto::single(r),
            Args::parse_from(&[
                "flexbar_rs",
                "--reads", "test.fq",
                "--barcode-trim-end", "RIGHT",
                "--barcode-min-overlap", "11",
                "--umi-tags"
            ])
        );
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_two_pass_alignment_fallback_single() {
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
            run_type: RunType::Single,
            ..Context::empty()
        };
        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();
        
        let rec = RecordDto::single(SingeRecordDto::with_attrs("read1", None, b"AAAAATTTTTAAATTT", b"IIIIIIIIIIIIIIII"));
        let res = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
    }

    #[test]
    fn test_barcode_splitting_paired() {
        let mut args = Args::parse_from(["flexbar", "--barcode-trim-end", "Right"]);
        args.barcode_trim_end = TrimEnd::Right;
        args.barcode_min_overlap = Some(1);

        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);

        let rec_in = RecordInTransformation::from_dto(rec);
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..12);
    }

    #[test]
    fn test_barcode_splitting_paired_only_read1() {
        let mut args = Args::parse_from(["flexbar", "--barcode-trim-end", "Right"]);
        args.barcode_trim_end = TrimEnd::Right;
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // Barcode in R1, should be detected
        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..12);

        // Barcode in R2, but not in R1. Should NOT be detected because context.barcodes2 is empty.
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        // It should NOT be bc1 because that's only for R1.
        assert!(result.alignment_results.barcode_tag.is_none());
        assert!(result.alignment_results.per_read[0].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_barcode_splitting_paired_only_read2() {
        let mut args = Args::parse_from(["flexbar", "--barcode-trim-end", "Right"]);
        args.barcode_trim_end = TrimEnd::Right;
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.barcodes2.push(SequenceWithId {
            id: "bc2".to_string(),
            sequence: "CCCC".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // Barcode in R2, should be detected
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTCCCC", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc2".to_string()));
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..12);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..8);

        // Barcode in R1, but context.barcodes is empty. Should NOT be detected.
        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTCCCC", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert!(result.alignment_results.barcode_tag.is_none());
        assert!(result.alignment_results.per_read[1].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_barcode_splitting_paired_both_reads_different_barcodes() {
        let mut args = Args::parse_from(["flexbar", "--barcode-trim-end", "Right"]);
        args.barcode_trim_end = TrimEnd::Right;
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });
        context.barcodes2.push(SequenceWithId {
            id: "bc2".to_string(),
            sequence: "CCCC".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // bc1 in R1, bc2 in R2. Should detect both, but BarcodeSplitter picks the BEST one (they likely score same).
        // Let's make them clearly different.
        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII"); // Score for AAAA in R1
        let r2 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTCCCC", b"IIIIIIIIIIII"); // Score for CCCC in R2
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        
        // It will pick one of them.
        assert!(result.alignment_results.barcode_tag == Some("barcode_bc1".to_string()) || 
                result.alignment_results.barcode_tag == Some("barcode_bc2".to_string()));

        // Case: bc1 in R2 (WRONG READ). Should NOT be detected because bc1 is in context.barcodes (for R1).
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        // It will NOT find bc1 in R2 because bc1 is not in context.barcodes2.
        // It will NOT find anything in R1.
        assert!(result.alignment_results.barcode_tag.is_none() || result.alignment_results.barcode_tag != Some("barcode_bc1".to_string()));
        
        // Case: bc2 in R1 (WRONG READ). Should NOT be detected.
        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTCCCC", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert!(result.alignment_results.barcode_tag.is_none() || result.alignment_results.barcode_tag != Some("barcode_bc2".to_string()));
    }

    #[test]
    fn test_barcode_splitting_separate_barcode_read_single() {
        let mut args = Args::parse_from(["flexbar", "--barcode-reads", "bar.fq"]);
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::SingleBarcoded;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // Barcode in separate barcode read
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let bar = SingeRecordDto::with_attrs("id", None, b"AAAA", b"IIII");
        let rec = RecordDto::single_barcoded(r1, bar);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);

        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // No removal for separate barcode read
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..12);
    }

    #[test]
    fn test_barcode_splitting_separate_barcode_read_paired() {
        let mut args = Args::parse_from(["flexbar", "--barcode-reads", "bar.fq"]);
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::PairedBarcoded;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // Barcode in separate barcode read
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGGGGG", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"HHHHHHHHHHHH", b"IIIIIIIIIIII");
        let bar = SingeRecordDto::with_attrs("id", None, b"AAAA", b"IIII");
        let rec = RecordDto::paired_barcoded(r1, r2, bar);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);

        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // No removal
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..12);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..12);
    }

    #[test]
    fn test_barcode_splitting_separate_barcode_read_ignores_read_barcodes() {
        let mut args = Args::parse_from(["flexbar", "--barcode-reads", "bar.fq", "--barcode-trim-end", "Right"]);
        args.barcode_trim_end = TrimEnd::Right;
        args.barcode_min_overlap = Some(4);

        let mut context = Context::empty();
        context.run_type = RunType::SingleBarcoded;
        context.barcodes.push(SequenceWithId {
            id: "bc1".to_string(),
            sequence: "AAAA".to_string(),
        });
        context.barcodes.push(SequenceWithId {
            id: "bc2".to_string(),
            sequence: "CCCC".to_string(),
        });

        let step = BarcodeSplitter::new(&args, &context);
        let mut state = step.create_state();

        // R1 contains bc1 at the end, but barcode read contains bc2.
        // It should detect bc2 and NOT trim bc1 from R1.
        let r1 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGAAAA", b"IIIIIIIIIIII");
        let bar = SingeRecordDto::with_attrs("id", None, b"CCCC", b"IIII");
        let rec = RecordDto::single_barcoded(r1, bar);
        let result = step.apply(RecordInTransformation::from_dto(rec), &mut state);

        assert_eq!(result.alignment_results.barcode_tag, Some("barcode_bc2".to_string()));
        // No removal because separate barcode read was used
        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..12);
        assert_eq!(result.get_sequence(0), b"GGGGGGGGAAAA");
    }
}
