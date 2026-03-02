use log::trace;
use crate::context::Context;
use crate::args::Args;
use crate::types::{SingeRecordDto, QualTrimType};
use crate::steps::single_adapter_remover::{SingleAdapterRemover, SingleAdapterRemoverState};
use crate::steps::single_uncalled_sequence_filter::{SingleUncalledSequenceFilter, SingleUncalledSequenceFilterState};
use crate::steps::single_barcode_splitter::{SingleBarcodeSplitter, SingleBarcodeSplitterState};
use crate::steps::single_quality_trim::{SingleQualityTrim, SingleQualityTrimState};
use crate::steps::single_filter_length::{SingleFilterLength, SingleFilterLengthState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable};

pub struct SingleWorker {
    pub adapter_remover: Option<SingleAdapterRemover>,
    pub barcode_splitter: Option<SingleBarcodeSplitter>,
    pub quality_trim: Option<SingleQualityTrim>,
    pub filter_length: SingleFilterLength,
    pub uncalled_sequence_filter: Option<SingleUncalledSequenceFilter>,
    pub qtrim_post_removal: bool,
}

pub struct SingleWorkerState {
    pub adapter_remover: Option<SingleAdapterRemoverState>,
    pub barcode_splitter: Option<SingleBarcodeSplitterState>,
    pub quality_trim: Option<SingleQualityTrimState>,
    pub filter_length: SingleFilterLengthState,
    pub uncalled_sequence_filter: Option<SingleUncalledSequenceFilterState>,
}

impl SingleWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let uncalled_sequence_filter = match args.max_uncalled {
            0 => None,
            _ => Some(SingleUncalledSequenceFilter::new(args, context))
        };
        let quality_trim = match args.qtrim {
            QualTrimType::Off => None,
            _ => Some(SingleQualityTrim::new(args, context))
        };
        let barcode_splitter = match context.barcodes.is_empty() {
            true => None,
            _ => Some(SingleBarcodeSplitter::new(args, context))
        };
        let adapter_remover = match context.adapters.is_empty() {
            true => None,
            _ => Some(SingleAdapterRemover::new(args, context))
        };
        let filter_length = SingleFilterLength::new(args, context);
        
        let qtrim_post_removal = args.qtrim_post_removal;
        if args.trace {
            trace!("Single worker pipeline:");
            if let Some(s) = &uncalled_sequence_filter {
                trace!("- Step SingleUncalledSequenceFilter (max_uncalled: {})",
                    s.max_uncalled
                );
            }
            if !qtrim_post_removal && let Some(s) = &quality_trim {
                trace!("- Step SingleQualityTrim (threshold: {}, window size: {})",
                    s.threshold, s.win_size
                );
            }
            if let Some(s) = &barcode_splitter {
                trace!("- Step SingleBarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})",
                    s.trim_end, s.barcode_keep, s.barcode_unassigned
                );
            }
            if let Some(s) = &adapter_remover {
                trace!("- Step SingleAdapterRemover (trim end {:?})", s.trim_end);
            }
            if qtrim_post_removal && let Some(s) = &quality_trim {
                trace!("- Step SingleQualityTrim (threshold: {}, window size: {})",
                    s.threshold, s.win_size
                );
            }
            trace!("- Step SingleFilterLength (min length: {})",
                filter_length.min_read_length
            );
        }

        Self {
            quality_trim,
            barcode_splitter,
            adapter_remover,
            filter_length,
            uncalled_sequence_filter,
            qtrim_post_removal,
        }
    }

    pub fn create_state(&self) -> SingleWorkerState {
        SingleWorkerState {
            adapter_remover: self.adapter_remover.as_ref().map(|s| s.create_state()),
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            quality_trim: self.quality_trim.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
            uncalled_sequence_filter: self.uncalled_sequence_filter.as_ref().map(|s| s.create_state()),
        }
    }
    
    pub fn process_record(
        &self,
        record: &SingeRecordDto,
        state: &mut SingleWorkerState
    ) -> FinalRecord<SingeRecordDto> {
        let mut rec_in = record.get_record_in_transformation();

        if let Some(step) = &self.uncalled_sequence_filter {
            rec_in = step.apply(rec_in, state.uncalled_sequence_filter.as_mut().unwrap());
        }
        if !self.qtrim_post_removal && let Some(step) = &self.quality_trim {
            rec_in = step.apply(rec_in, state.quality_trim.as_mut().unwrap());
        }
        if let Some(step) = &self.barcode_splitter {
            rec_in = step.apply(rec_in, state.barcode_splitter.as_mut().unwrap());
        }
        if let Some(step) = &self.adapter_remover {
            rec_in = step.apply(rec_in, state.adapter_remover.as_mut().unwrap());
        }
        if self.qtrim_post_removal && let Some(step) = &self.quality_trim {
            rec_in = step.apply(rec_in, state.quality_trim.as_mut().unwrap());
        }
        rec_in = self.filter_length.apply(rec_in, &mut state.filter_length);
        
        finalize(rec_in)
    }
}

fn finalize(rec: RecordInTransformation<SingeRecordDto>) -> FinalRecord<SingeRecordDto> {
    if rec.alignment_results.record_dropped {
        return FinalRecord {
            record: None,
            output_tag: None,
            stats: rec.stats,
        }
    }

    let mut desc_parts = Vec::new();
    match {rec.record.desc()} {
        Some(desc) => desc_parts.push(desc),
        None => (),
    }

    for tag in &rec.alignment_results.removal_tags {
        desc_parts.push(tag);
    }

    for tag in &rec.alignment_results.umi_tags {
        desc_parts.push(tag);
    }

    let new_desc: Option<String> = if desc_parts.is_empty() {
        None
    } else {
        Some(desc_parts.join("_"))
    };

    let final_record = SingeRecordDto::with_attrs(
        rec.record.id(),
        new_desc.as_deref(),
        rec.get_sequence(),
        rec.get_qual(),
    );

    FinalRecord {
        record: Some(final_record),
        output_tag: rec.alignment_results.barcode_tag,
        stats: rec.stats,
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{Context, SequenceWithId};
    use crate::args::Args;
    use crate::types::SingeRecordDto;
    use clap::Parser;
    
    fn process_record(args: &Args, record: &SingeRecordDto) -> FinalRecord<SingeRecordDto> {
        let context = Context::new(args);
        let worker = SingleWorker::new(args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_left() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        // Barcode at start
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(result.output_tag, Some("barcode_bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
        // Stats are now in writer, so worker.process doesn't see them yet.
        // We'll trust GenericWriter to handle it.
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);

        // 2 'N's - should keep
        let record = SingeRecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let result = process_record(&args, &record);
        assert_eq!(result.stats.total_discarded_uncalled, 0);

        // 3 'N's - should discard
        let record = SingeRecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let result = process_record(&args, &record);
        assert_eq!(result.stats.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal"));
    }

    #[test]
    fn test_removal_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "adapter".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_0"));
    }

    #[test]
    fn test_removal_tags_both() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", 
            "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4",
            "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5",
            "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "adapter".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGTCTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        // barcode removal appends _Flexbar_removal
        // adapter removal appends _Flexbar_removal_0
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_Flexbar_removal_0"));
    }

    #[test]
    fn test_quality_trim_tail() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "TAIL", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        
        // Threshold 20 + 33 (sanger) = 53 ('5')
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_quality_trim_win() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "WIN", "--qtrim-threshold", "20", "--qtrim-win-size", "5", "--min-read-length", "0"]);
        
        // Threshold 53. Window size 5.
        // "IIIII" sum = 5 * 73 = 365. Threshold sum = 5 * 53 = 265.
        // "44444" sum = 5 * 52 = 260.
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGTACGT", b"IIIII44444III");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_quality_trim_bwa() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "BWA", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "NNGGGGNN".to_string()}
        ];

        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        // Barcode NNGGGGNN matches AAGGGGTT
        let record = SingeRecordDto::with_attrs("read1", None, b"AAGGGGTTACGT", b"IIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(result.output_tag, Some("barcode_bc1".to_string()));
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("AA_TT"));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_umi_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATNNNN", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "adapter".to_string(), sequence: "AGATNNNN".to_string() }];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCCCC", b"IIIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("CCCC"));
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_min_read_length() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);

        // Length 10 - should keep
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789", b"IIIIIIIIII");
        let result = process_record(&args, &record);
        assert!(result.record.is_some(), "Should keep read of length 10 when min is 10");
        assert_eq!(result.stats.total_discarded_length, 0);

        // Length 9 - should discard
        let record = SingeRecordDto::with_attrs("read2", None, b"012345678", b"IIIIIIIII");
        let result = process_record(&args, &record);
        assert!(result.record.is_none(), "Should discard read of length 9 when min is 10");
        assert_eq!(result.stats.total_discarded_length, 1);
    }
}