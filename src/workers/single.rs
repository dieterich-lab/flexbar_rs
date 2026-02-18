use crate::context::Context;
use crate::args::Args;
use crate::types::{SingeRecordDto, QualTrimType};
use crate::steps::single_adapter_remover::SingleAdapterRemover;
use crate::steps::single_uncalled_sequence_filter::SingleUncalledSequenceFilter;
use crate::steps::single_barcode_splitter::SingleBarcodeSplitter;
use crate::steps::single_quality_trim::SingleQualityTrim;
use crate::steps::single_filter_length::SingleFilterLength;
use crate::steps::single_finalize::SingleFinalizeStep;
use crate::steps::generic::GenericStep;
use crate::workers::generic::GenericWorker;
use std::sync::Arc;
use log::debug;

pub struct SingleWorker (GenericWorker<SingeRecordDto> );

impl SingleWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let mut steps: Vec<Arc<dyn GenericStep<SingeRecordDto>>> = Vec::new();
        debug!("Creating new single worker pipeline:");

        if args.max_uncalled > 0 {
            debug!(
                "- Step SingleUncalledSequenceFilter (max_uncalled: {})", args.max_uncalled
            );
            steps.push(
                Arc::new(SingleUncalledSequenceFilter::new(args))
            );
        }
        if args.qtrim != QualTrimType::Off && !args.qtrim_post_removal {
            debug!(
                "- Step SingleQualityTrim (mode {:?}, format {:?}, window size {})",
                args.qtrim,
                context.format,
                args.qtrim_win_size
            );
            steps.push(
                Arc::new(SingleQualityTrim::new(args, context.format))
            );
        }
        if !context.barcodes.is_empty() {
            debug!(
                "- Step SingleBarcodeSplitter (min overlap {:?}, err rate {}, trim end {:?}, keep: {})",
                args.barcode_min_overlap,
                args.barcode_error_rate,
                args.barcode_trim_end,
                args.barcode_keep
            );
            steps.push(
                Arc::new(SingleBarcodeSplitter::new(args, context.barcodes.clone()))
            );
        }
        if !context.adapters.is_empty() {
            debug!(
                "- Step SingleAdapterRemover (min overlap {}, err rate {}, trim_end {:?}, , revcomp {:?})",
                args.adapter_min_overlap,
                args.adapter_error_rate,
                args.adapter_trim_end,
                args.adapter_revcomp,
            );
            steps.push(
                Arc::new(SingleAdapterRemover::new(args, &context.adapters))
            );
        }

        if args.qtrim != QualTrimType::Off && args.qtrim_post_removal {
            debug!(
                "- Step SingleQualityTrim (mode {:?}, format {:?}, window size {})",
                args.qtrim,
                context.format,
                args.qtrim_win_size
            );
            steps.push(
                Arc::new(SingleQualityTrim::new(args, context.format))
            )
        }

        steps.push(Arc::new(SingleFilterLength::new(args)));
        debug!("Single worker pipeline created.");
        Self (GenericWorker::new(steps, Arc::new(SingleFinalizeStep)))
    }
    
    pub fn process(&self, record: &SingeRecordDto) -> (Option<SingeRecordDto>, Option<String>, crate::statistics::WorkerStatistics) {
        self.0.process(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{Context, SequenceWithId};
    use crate::args::Args;
    use crate::types::SingeRecordDto;
    use clap::Parser;

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let context = Context::new(&args);
        
        let worker = SingleWorker::new(&args, &context);
        
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
         
         
        
        let worker = SingleWorker::new(&args, &context);
        
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
        context.barcodes = Arc::new(vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}]);
        
        let worker = SingleWorker::new(&args, &context);
        
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
         
        let worker = SingleWorker::new(&args, &context);

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
        context.barcodes = Arc::new(vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ]);
        
        let worker = SingleWorker::new(&args, &context);
        
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
        context.adapters = Arc::new(vec!["AGATCGGAAGAG".to_string()]);
        
        let worker = SingleWorker::new(&args, &context);
        
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
        context.adapters = Arc::new(vec!["AGATCGGAAGAG".to_string()]);
        context.barcodes = Arc::new(vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ]);
        
        let worker = SingleWorker::new(&args, &context);
        
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
         
        let worker = SingleWorker::new(&args, &context);

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
         
        let worker = SingleWorker::new(&args, &context);

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
         
        let worker = SingleWorker::new(&args, &context);

        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = Arc::new(vec![
            SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ]);

        let worker = SingleWorker::new(&args, &context);
        
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
        context.adapters = Arc::new(vec!["AGATNNNN".to_string()]);
        
        let worker = SingleWorker::new(&args, &context);
        
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
         
        let worker = SingleWorker::new(&args, &context);

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