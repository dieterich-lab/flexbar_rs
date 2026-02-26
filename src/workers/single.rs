use crate::context::Context;
use crate::args::Args;
use crate::types::{SingeRecordDto, QualTrimType};
use crate::steps::single_adapter_remover::{SingleAdapterRemover, SingleAdapterRemoverSetup};
use crate::steps::single_uncalled_sequence_filter::{SingleUncalledSequenceFilter, SingleUncalledSequenceFilterSetup};
use crate::steps::single_barcode_splitter::{SingleBarcodeSplitter, SingleBarcodeSplitterSetup};
use crate::steps::single_quality_trim::{SingleQualityTrim, SingleQualityTrimSetup};
use crate::steps::single_filter_length::{SingleFilterLength, SingleFilterLengthSetup};
use crate::steps::single_finalize::SingleFinalizeStep;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::workers::generic::{GenericWorker, WorkerResult};
use std::sync::Arc;
use log::{debug, trace};

use std::any::Any;
use crate::statistics::WorkerStatistics;

struct SingleWorkerSetup {
    adapter_remover: Option<SingleAdapterRemoverSetup>,
    barcode_splitter: Option<SingleBarcodeSplitterSetup>,
    quality_trim: Option<SingleQualityTrimSetup>,
    filter_length: SingleFilterLengthSetup,
    uncalled_sequence_filter: Option<SingleUncalledSequenceFilterSetup>,
    qtrim_post_removal: bool,
    trace: bool
}

struct SingleWorker {
    adapter_remover: Option<SingleAdapterRemover>,
    barcode_splitter: Option<SingleBarcodeSplitter>,
    quality_trim: Option<SingleQualityTrim>,
    filter_length: SingleFilterLength,
    uncalled_sequence_filter: Option<SingleUncalledSequenceFilter>,
}

impl SingleWorkerSetup {
    pub fn new(args: &Args, context: &Context) -> Self {
        // Setups
        let uncalled_sequence_filter = match args.max_uncalled {
            0 => None,
            _ =>Some(SingleUncalledSequenceFilter::create_setup(args, context))
        };
        let quality_trim = match args.qtrim {
            QualTrimType::Off => None,
            _ => Some(SingleQualityTrim::create_setup(args,context))
        };
        let barcode_splitter = match context.barcodes.is_empty() {
            true => None,
            _ => Some(SingleBarcodeSplitter::create_setup(args, context))
        };
        let adapter_remover = match context.adapters.is_empty() {
            true => None,
            _ => Some(SingleAdapterRemover::create_setup(args, context))
        };
        let filter_length = SingleFilterLength::create_setup(args, context);
        
        let qtrim_post_removal = args.qtrim_post_removal;
        let trace = args.trace;
        if trace {
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
            trace
        }
    }
}

impl SingleWorker {
    fn new(setup: &SingleWorkerSetup) -> Self {
        Self {
            adapter_remover: match &setup.adapter_remover {
                Some(s) => Some(SingleAdapterRemover::new(s)),
                None => None,
            },
            barcode_splitter: match &setup.barcode_splitter {
                Some(s) => Some(SingleBarcodeSplitter::new(s)),
                None => None,
            },

            quality_trim: match &setup.quality_trim {
                Some(s) => Some(SingleQualityTrim::new(s)),
                None => None,
            },
            uncalled_sequence_filter: match &setup.uncalled_sequence_filter { 
                Some(s) => Some(SingleUncalledSequenceFilter::new(s)),
                None => None,
            },
            filter_length: SingleFilterLength::new(&setup.filter_length),
        }
    }
    fn process_record(
        &mut self,
        record: &SingeRecordDto,
        setup: &SingleWorkerSetup
    ) -> WorkerResult<SingeRecordDto> {
        let mut rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: Default::default(),
            stats: WorkerStatistics::default(),
        };

        if let Some(s) = &setup.uncalled_sequence_filter {
            rec_in = self.uncalled_sequence_filter.as_mut().unwrap().apply(rec_in, &s);
        }
        if !setup.qtrim_post_removal && let Some(s) = &setup.quality_trim {
            rec_in = self.quality_trim.as_mut().unwrap().apply(rec_in, &s);
        }
        if let Some(s) = &setup.barcode_splitter {
            rec_in = self.barcode_splitter.as_mut().unwrap().apply(rec_in, &s);
        }
        if let Some(s) = &setup.adapter_remover {
            rec_in = self.adapter_remover.as_mut().unwrap().apply(rec_in, &s);
        }
        if setup.qtrim_post_removal && let Some(s) = &setup.quality_trim {
            rec_in = self.quality_trim.as_mut().unwrap().apply(rec_in, &s);
        }
        rec_in = self.filter_length.apply(rec_in, &setup.filter_length);

        let output_record = match rec_in.alignment_results.record_dropped {
            true => None,
            false => Some(
                SingeRecordDto::with_attrs(
                    rec_in.record.id(),
                    rec_in.record.desc(),
                    rec_in.get_sequence(),
                    rec_in.get_qual(),
                )
            )
        };
        WorkerResult {
            record: output_record,
            output_tag: rec_in.alignment_results.barcode_tag,
            stats: rec_in.stats,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{Context, SequenceWithId};
    use crate::args::Args;
    use crate::types::SingeRecordDto;
    use clap::Parser;
    
    fn process_record(args: &Args, record: &SingeRecordDto) -> WorkerResult<SingeRecordDto> {
        let context = Context::new(args);
        let setup = SingleWorkerSetup::new(args, &context);
        let mut worker = SingleWorker::new(&setup);
        worker.process_record(record, &setup)
    }

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let context = Context::new(&args);
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let (processed_opt, _, _stats) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_left() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let context = Context::new(&args);
         
         
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let (processed_opt, _, _stats) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        // Barcode at start
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (processed_opt, tag_opt, _stats) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(tag_opt, Some("barcode_bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
        // Stats are now in writer, so worker.process doesn't see them yet.
        // We'll trust GenericWriter to handle it.
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let context = Context::new(&args);
         
        let mut worker = SingleWorker::new(&args, &context);

        // 2 'N's - should keep
        let record = SingeRecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let (_, _, stats) = worker.process(&record);
        assert_eq!(stats.total_discarded_uncalled, 0);

        // 3 'N's - should discard
        let record = SingeRecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let (_, _, stats) = worker.process(&record);
        assert_eq!(stats.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal"));
    }

    #[test]
    fn test_removal_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec!["AGATCGGAAGAG".to_string()];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
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
        context.adapters = vec!["AGATCGGAAGAG".to_string()];
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGTCTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIIIIIIII");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        // barcode removal appends _Flexbar_removal
        // adapter removal appends _Flexbar_removal_0
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_Flexbar_removal_0"));
    }

    #[test]
    fn test_quality_trim_tail() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "TAIL", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let context = Context::new(&args);
         
        let mut worker = SingleWorker::new(&args, &context);

        // Threshold 20 + 33 (sanger) = 53 ('5')
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_quality_trim_win() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "WIN", "--qtrim-threshold", "20", "--qtrim-win-size", "5", "--min-read-length", "0"]);
        let context = Context::new(&args);
         
        let mut worker = SingleWorker::new(&args, &context);

        // Threshold 53. Window size 5.
        // "IIIII" sum = 5 * 73 = 365. Threshold sum = 5 * 53 = 265.
        // "44444" sum = 5 * 52 = 260.
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGTACGT", b"IIIII44444III");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_quality_trim_bwa() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "BWA", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let context = Context::new(&args);
         
        let mut worker = SingleWorker::new(&args, &context);

        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "NNGGGGNN".to_string()}
        ];

        let mut worker = SingleWorker::new(&args, &context);
        
        // Barcode NNGGGGNN matches AAGGGGTT
        let record = SingeRecordDto::with_attrs("read1", None, b"AAGGGGTTACGT", b"IIIIIIIIIIII");
        let (processed_opt, tag_opt, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(tag_opt, Some("barcode_bc1".to_string()));
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("AA_TT"));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_umi_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATNNNN", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec!["AGATNNNN".to_string()];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCCCC", b"IIIIIIIIIIIII");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("CCCC"));
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_min_read_length() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let context = Context::new(&args);
         
        let mut worker = SingleWorker::new(&args, &context);

        // Length 10 - should keep
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789", b"IIIIIIIIII");
        let (processed_opt, _, stats) = worker.process(&record);
        assert!(processed_opt.is_some(), "Should keep read of length 10 when min is 10");
        assert_eq!(stats.total_discarded_length, 0);

        // Length 9 - should discard
        let record = SingeRecordDto::with_attrs("read2", None, b"012345678", b"IIIIIIIII");
        let (processed_opt, _, stats) = worker.process(&record);
        assert!(processed_opt.is_none(), "Should discard read of length 9 when min is 10");
        assert_eq!(stats.total_discarded_length, 1);
    }
}