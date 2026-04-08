use log::debug;
use crate::context::Context;
use crate::args::Args;
use crate::types::{SingeRecordDto, QualTrimType};
use crate::steps::single::pre_trim::{SinglePreTrim, SinglePreTrimState};
use crate::steps::single::uncalled_filter::{SingleUncalledFilter, SingleUncalledFilterState};
use crate::steps::single::adapter_remover::{SingleAdapterRemover, SingleAdapterRemoverState};
use crate::steps::single::barcode_splitter::{SingleBarcodeSplitter, SingleBarcodeSplitterState};
use crate::steps::single::quality_trim::{SingleQualityTrim, SingleQualityTrimState};
use crate::steps::single::filter_length::{SingleFilterLength, SingleFilterLengthState};
use crate::steps::single::post_trim::{SinglePostTrim, SinglePostTrimState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::workers::util::finalize_read_stats;

pub struct SingleWorker {
    pub pre_trim: Option<SinglePreTrim>,
    pub adapter_remover: Option<SingleAdapterRemover>,
    pub barcode_splitter: Option<SingleBarcodeSplitter>,
    pub quality_trim: Option<SingleQualityTrim>,
    pub filter_length: SingleFilterLength,
    pub post_trim: SinglePostTrim,
    pub uncalled_filter: Option<SingleUncalledFilter>,
    pub qtrim_post_removal: bool,
}

pub struct SingleWorkerState {
    pub pre_trim: Option<SinglePreTrimState>,
    pub adapter_remover: Option<SingleAdapterRemoverState>,
    pub barcode_splitter: Option<SingleBarcodeSplitterState>,
    pub quality_trim: Option<SingleQualityTrimState>,
    pub filter_length: SingleFilterLengthState,
    pub post_trim: SinglePostTrimState,
    pub uncalled_filter: Option<SingleUncalledFilterState>,
}

impl SingleWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let pre_trim = if args.pre_trim_left > 0 || args.pre_trim_right > 0 {
            Some(SinglePreTrim::new(args, context))
        } else {
            None
        };
        let uncalled_filter = if args.max_uncalled > 0 {
            Some(SingleUncalledFilter::new(args, context))
        } else {
            None
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
        let post_trim = SinglePostTrim::new(args, context);
        
        let qtrim_post_removal = args.qtrim_post_removal;
        if args.trace {
            debug!("Single worker pipeline:");
            if let Some(s) = &pre_trim {
                debug!("- Step SinglePreTrim (left: {}, right: {})", s.pre_trim_left, s.pre_trim_right);
            }
            if let Some(s) = &uncalled_filter {
                debug!("- Step SingleUncalledFilter (max_uncalled: {})",
                    s.max_uncalled
                );
            }
            if !qtrim_post_removal && let Some(s) = &quality_trim {
                debug!("- Step SingleQualityTrim (threshold: {}, window size: {})",
                    s.inner.threshold, s.inner.win_size
                );
            }
            if let Some(s) = &barcode_splitter {
                debug!("- Step SingleBarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})",
                    s.trim_end, s.barcode_keep, s.barcode_unassigned
                );
            }
            if let Some(s) = &adapter_remover {
                debug!("- Step SingleAdapterRemover (trim end {:?})", s.inner.trim_end);
            }
            if qtrim_post_removal && let Some(s) = &quality_trim {
                debug!("- Step SingleQualityTrim (threshold: {}, window size: {})",
                    s.inner.threshold, s.inner.win_size
                );
            }
            debug!("- Step SingleFilterLength (min length: {})",
                filter_length.min_read_length
            );
            if post_trim.post_trim_length > 0 {
                debug!("- Step PostTrim (length: {})", post_trim.post_trim_length);
            }
        }

        Self {
            pre_trim,
            quality_trim,
            barcode_splitter,
            adapter_remover,
            filter_length,
            post_trim,
            uncalled_filter,
            qtrim_post_removal,
        }
    }

    pub fn create_state(&self) -> SingleWorkerState {
        SingleWorkerState {
            pre_trim: self.pre_trim.as_ref().map(|s| s.create_state()),
            adapter_remover: self.adapter_remover.as_ref().map(|s| s.create_state()),
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            quality_trim: self.quality_trim.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
            post_trim: self.post_trim.create_state(),
            uncalled_filter: self.uncalled_filter.as_ref().map(|s| s.create_state()),
        }
    }
    
    pub fn process_record(
        &self,
        record: &SingeRecordDto,
        state: &mut SingleWorkerState
    ) -> FinalRecord<SingeRecordDto> {
        let mut rec_in: RecordInTransformation<SingeRecordDto, SingleAlignmentResults> = record.get_record_in_transformation();
        rec_in.stats.total_bases_read = record.seq().len();

        if let Some(step) = &self.pre_trim {
            rec_in = step.apply(rec_in, state.pre_trim.as_mut().unwrap());
        }
        if let Some(step) = &self.uncalled_filter {
            rec_in = step.apply(rec_in, state.uncalled_filter.as_mut().unwrap());
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
        
        if self.post_trim.post_trim_length > 0 {
            rec_in = self.post_trim.apply(rec_in, &mut state.post_trim);
        }
        finalize(rec_in)
    }
}

fn finalize(mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>) -> FinalRecord<SingeRecordDto> {
    let original_len = rec.record.seq().len();
    finalize_read_stats(&mut rec.stats, &rec.alignment_results.read, original_len);

    if rec.alignment_results.read.dropped {
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

    for tag in &rec.alignment_results.read.removal_tags {
        desc_parts.push(tag);
    }

    for tag in &rec.alignment_results.read.umi_tags {
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
        context.adapters = vec![SequenceWithId { id: "ad1".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        let processed = result.record.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_ad1"));
    }

    #[test]
    fn test_removal_tags_both() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", 
            "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4",
            "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5",
            "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "ad1".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
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
        assert_eq!(processed.desc(), Some("Flexbar_removal_Flexbar_removal_ad1"));
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

    #[test]
    fn test_pre_trim() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--pre-trim-left", "5", "--pre-trim-right", "3", "--min-read-length", "0"]);
        let context = Context::empty();
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();

        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let res = worker.process_record(&record, &mut state);
        
        // FinalRecord structure has record: Option<SingeRecordDto>
        assert!(res.record.is_some(), "Record should not be discarded");
        let _rec = res.record.unwrap();
        
        assert_eq!(res.stats.bases_pre_trim, 8);
        assert_eq!(res.stats.total_bases_read, 15);
    }

    #[test]
    fn test_worker_statistics() {
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--adapter-seq", "AAAA", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4",
            "--qtrim", "TAIL", "--qtrim-threshold", "20",
            "--min-read-length", "5"
        ]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        context.adapters = vec![SequenceWithId{id: "ad1".to_string(), sequence: "AAAA".to_string()}];
        
        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        
        // Input length 13: GGGG (4) + 12345 (5) + AAAA (4) = 13
        // After barcode removal (LTAIL default): 12345 (5) + AAAA (4) = 9. Removed 4.
        // After adapter removal (RIGHT): 12345 (5). Removed 4.
        
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGG12345AAAA", b"IIIIIIIIIIIII");
        let result = worker.process_record(&record, &mut state);
        
        assert_eq!(result.stats.total_bases_read, 13);
        assert_eq!(result.stats.bases_barcode_removal, 4);
        assert_eq!(result.stats.bases_adapter_removal, 4);
        assert_eq!(result.stats.bases_quality_trim, 0);
        assert_eq!(result.stats.bases_discarded, 0);
        
        // Test with quality trim
        let record2 = SingeRecordDto::with_attrs("read2", None, b"GGGG12345AAAA!!", b"IIIIIIIIIIIII!! "); // !! are low qual
        // 4 (GGGG) + 5 (12345) + 4 (AAAA) + 2 (!!) = 15
        let result2 = worker.process_record(&record2, &mut state);
        assert_eq!(result2.stats.total_bases_read, 15);
        assert_eq!(result2.stats.bases_barcode_removal, 4);
        assert_eq!(result2.stats.bases_adapter_removal, 4);
        assert_eq!(result2.stats.bases_quality_trim, 2);
        assert_eq!(result2.stats.bases_discarded, 0);

        // Test discarded read
        // min-read-length is 5.
        // GGGG (4) + 12 (2) + AAAA (4) = 10.
        // After barcode (4) and adapter (4) removal, length is 2. 
        // 2 < 5, so dropped.
        let record3 = SingeRecordDto::with_attrs("read3", None, b"GGGG12AAAA", b"IIIIIIIIII");
        let result3 = worker.process_record(&record3, &mut state);
        assert!(result3.record.is_none());
        assert_eq!(result3.stats.total_bases_read, 10);
        assert_eq!(result3.stats.bases_discarded, 10);
        assert_eq!(result3.stats.bases_barcode_removal, 0);
        assert_eq!(result3.stats.bases_adapter_removal, 0);
        assert_eq!(result3.stats.bases_quality_trim, 0);
        assert_eq!(result3.stats.bases_post_trim, 0);
    }

    #[test]
    fn test_post_trim() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10", "--min-read-length", "0"]);
        
        // Read of length 15 should be trimmed to 10
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let result = process_record(&args, &record);
        let processed = result.record.expect("Should keep read");
        
        assert_eq!(processed.seq().len(), 10);
        assert_eq!(processed.seq(), b"0123456789");
        assert_eq!(result.stats.bases_post_trim, 5);
        assert_eq!(result.stats.total_bases_read, 15);
        
        // Read of length 8 should not be trimmed
        let record2 = SingeRecordDto::with_attrs("read2", None, b"01234567", b"IIIIIIII");
        let result2 = process_record(&args, &record2);
        let processed2 = result2.record.expect("Should keep read");
        
        assert_eq!(processed2.seq().len(), 8);
        assert_eq!(result2.stats.bases_post_trim, 0);
    }

    fn process_with_reference_setup(record: &SingeRecordDto) -> FinalRecord<SingeRecordDto> {
        let args = Args::parse_from(&[
            "flexbar_rs",
            "--removal-tags",
            "--umi-tags",
            "--min-read-length", "18",
            "--max-uncalled", "1",
            "--qtrim-threshold", "10",
            "--barcode-unassigned",
            "--post-trim-length", "248",
            "--adapter-min-overlap", "6",
            "--qtrim-format", "sanger",
        ]);

        let barcodes = vec![
            SequenceWithId { id: "SmartSeqUMI".to_string(), sequence: "ATTGCGCAATGNNNNNNNN".to_string() },
        ];

        let adapters = vec![
            SequenceWithId { id: "Trans1".to_string(), sequence: "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans1_rc".to_string(), sequence: "CTGTCTCTTATACACATCTGACGCTGCCGACGA".to_string() },
            SequenceWithId { id: "Trans2".to_string(), sequence: "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans2_rc".to_string(), sequence: "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC".to_string() },
        ];

        let mut context = Context::empty();
        context.barcodes = barcodes;
        context.adapters = adapters;

        let worker = SingleWorker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_adapter_with_errors() {
        // Read:
        // ATTGCGCAATGCTTGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGCTGTCTCTTATACACAAC CCGAGCCCACGAGAATGCCGGTCAGATCTCGTATG
        // ATTGCGCAATGNNNNNNNN (Barcode: LTAIL - perfect match)            \ ,              \
        // (Adapter: RIGHT - poor match but OK)            CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        let r = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1221:21866:34444",
            Some("N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"ATTGCGCAATGCTTGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGCTGTCTCTTATACACAACCCGAGCCCACGAGAATGCCGGTCAGATCTCGTATG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,,,,,,,,,,FFF:F,F,F,::F,F,FF:FF,F,F,F,FF:F:,:,FF,,,:,,F:",
        );
        let res = process_with_reference_setup(&r);
        assert!(res.record.is_some(), "Should keep read");
        assert_eq!(res.record.unwrap().seq(), b"AAAAAAAAAAAAAAAAAAAAAAAAAATCG");
    }

    // #[test]
    // fn test_adapter_with_errors_in_overlap() {
    //     // Read:
    //     // GATAAACAGAGAGCAAGGCCTCTCCCAGGAGGTCACCAAAACTGCCCACCTACCTCCTCCCTAACCCACTGGCTTCTGCTTTCCTCTCTCCCTGT TCTTAT
    //     // (Adapter: RIGHT - poor match but OK)                                                       CTGTCTCTTATACACATCTGACGCTGCCGACGA
    //     let r = SingeRecordDto::with_attrs(
    //         "A00620:147:H225CDMXY:2:1101:17183:24831",
    //         Some("N:0:TGCCGGTCAG+TTGTATCAGG"),
    //         b"GATAAACAGAGAGCAAGGCCTCTCCCAGGAGGTCACCAAAACTGCCCACCTACCTCCTCCCTAACCCACTGGCTTCTGCTTTCCTCTCTCCCTGTTCTTAT",
    //         b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F",
    //     );
    //     // It seems that we found here BioRust bug:
    //     // If just extract the relevant parts from the sequences, BioRust finds this alignment
    //     // with score 3, and we discard it as invalid:
    //     //
    //     // CTCCCTGTTCTTAT
    //     //    |\\\\||||||
    //     //    CTGTCTCTTAT...
    //     //
    //     // The legacy seqan library used by the old flexbar found this alignment:
    //     // CTCCCTGT TCTTAT
    //     //     |||| ||||||
    //     //     CTGTCTCTTAT...
    //     //
    //     // This alignment has 10 matches, minus 6 for the opening gab, results in score 4.
    //     // See also: https://github.com/rust-bio/rust-bio/issues/656
    //     let res = process_with_reference_setup(&r);
    //     assert!(res.record.is_some(), "Should keep read");
    //     assert_eq!(res.record.unwrap().seq(), b"GATAAACAGAGAGCAAGGCCTCTCCCAGGAGGTCACCAAAACTGCCCACCTACCTCCTCCCTAACCCACTGGCTTCTGCTTTCCTCTCTCC");
    // }
}