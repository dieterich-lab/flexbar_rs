use crate::context::Context;
use crate::args::Args;
use crate::types::{SingeRecordDto, QualTrimType};
use crate::steps::single_adapter_remover::{SingleAdapterRemover, SingleAdapterRemoverSetup};
use crate::steps::single_uncalled_sequence_filter::{SingleUncalledSequenceFilter, SingleUncalledSequenceFilterSetup};
use crate::steps::single_barcode_splitter::{SingleBarcodeSplitter, SingleBarcodeSplitterSetup};
use crate::steps::single_quality_trim::{SingleQualityTrim, SingleQualityTrimSetup};
use crate::steps::single_filter_length::{SingleFilterLength, SingleFilterLengthSetup};
use crate::steps::single_finalize::SingleFinalizeStep;
use crate::steps::generic::{GenericStep, FinalizeStep, Transformable};
use crate::workers::generic::GenericWorker;
use crate::statistics::WorkerStatistics;
use log::debug;

pub struct SingleWorker {
    uncalled_filter: Option<(SingleUncalledSequenceFilter, SingleUncalledSequenceFilterSetup)>,
    pre_quality_trim: Option<(SingleQualityTrim, SingleQualityTrimSetup)>,
    barcode_splitter: Option<(SingleBarcodeSplitter, SingleBarcodeSplitterSetup)>,
    adapter_remover: Option<(SingleAdapterRemover, SingleAdapterRemoverSetup)>,
    post_quality_trim: Option<(SingleQualityTrim, SingleQualityTrimSetup)>,
    length_filter: (SingleFilterLength, SingleFilterLengthSetup),
    finalize: SingleFinalizeStep,
}

impl SingleWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        debug!("Creating new single worker pipeline:");

        let uncalled_filter = if args.max_uncalled > 0 {
            debug!("- Step SingleUncalledSequenceFilter (max_uncalled: {})", args.max_uncalled);
            let setup = SingleUncalledSequenceFilter::create_setup(args, context);
            Some((SingleUncalledSequenceFilter::new(setup.clone()), setup))
        } else {
            None
        };

        let (pre_quality_trim, post_quality_trim) = if args.qtrim != QualTrimType::Off {
            let setup = SingleQualityTrim::create_setup(args, context);
            if args.qtrim_post_removal {
                debug!("- Step SingleQualityTrim (post-removal)");
                (None, Some((SingleQualityTrim::new(setup.clone()), setup)))
            } else {
                debug!("- Step SingleQualityTrim (pre-removal)");
                (Some((SingleQualityTrim::new(setup.clone()), setup)), None)
            }
        } else {
            (None, None)
        };

        let barcode_splitter = if !context.barcodes.is_empty() {
            debug!("- Step SingleBarcodeSplitter");
            let setup = SingleBarcodeSplitter::create_setup(args, context);
            Some((SingleBarcodeSplitter::new(setup.clone()), setup))
        } else {
            None
        };

        let adapter_remover = if !context.adapters.is_empty() {
            debug!("- Step SingleAdapterRemover");
            let setup = SingleAdapterRemover::create_setup(args, context);
            Some((SingleAdapterRemover::new(setup.clone()), setup))
        } else {
            None
        };

        let length_filter_setup = SingleFilterLength::create_setup(args, context);
        let length_filter = (SingleFilterLength::new(length_filter_setup.clone()), length_filter_setup);

        debug!("Single worker pipeline created.");

        Self {
            uncalled_filter,
            pre_quality_trim,
            barcode_splitter,
            adapter_remover,
            post_quality_trim,
            length_filter,
            finalize: SingleFinalizeStep,
        }
    }

    pub fn process(&mut self, record: &SingeRecordDto) -> (Option<SingeRecordDto>, Option<String>, WorkerStatistics) {
        let mut rec_in = record.get_record_in_transformation();

        if let Some((step, setup)) = &mut self.uncalled_filter {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped { return self.finalize.apply(rec_in); }
        }

        if let Some((step, setup)) = &mut self.pre_quality_trim {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped { return self.finalize.apply(rec_in); }
        }

        if let Some((step, setup)) = &mut self.barcode_splitter {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped { return self.finalize.apply(rec_in); }
        }

        if let Some((step, setup)) = &mut self.adapter_remover {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped { return self.finalize.apply(rec_in); }
        }

        if let Some((step, setup)) = &mut self.post_quality_trim {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped { return self.finalize.apply(rec_in); }
        }

        rec_in = self.length_filter.0.apply(rec_in, &self.length_filter.1);
        
        self.finalize.apply(rec_in)
    }
}

impl GenericWorker<SingeRecordDto> for SingleWorker {
    fn process(&mut self, record: &SingeRecordDto) -> (Option<SingeRecordDto>, Option<String>, WorkerStatistics) {
        self.process(record)
    }
}

impl Clone for SingleWorker {
    fn clone(&self) -> Self {
        Self {
            uncalled_filter: self.uncalled_filter.as_ref().map(|(_step, setup)| (SingleUncalledSequenceFilter::new(setup.clone()), setup.clone())),
            pre_quality_trim: self.pre_quality_trim.as_ref().map(|(_step, setup)| (SingleQualityTrim::new(setup.clone()), setup.clone())),
            barcode_splitter: self.barcode_splitter.as_ref().map(|(_step, setup)| (SingleBarcodeSplitter::new(setup.clone()), setup.clone())),
            adapter_remover: self.adapter_remover.as_ref().map(|(_step, setup)| (SingleAdapterRemover::new(setup.clone()), setup.clone())),
            post_quality_trim: self.post_quality_trim.as_ref().map(|(_step, setup)| (SingleQualityTrim::new(setup.clone()), setup.clone())),
            length_filter: (SingleFilterLength::new(self.length_filter.1.clone()), self.length_filter.1.clone()),
            finalize: SingleFinalizeStep,
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

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let mut context = Context::default();
        context.adapters = vec!["AGATCGGAAGAG".to_string()];
        
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
        let mut context = Context::default();
        context.adapters = vec!["AGATCGGAAGAG".to_string()];
        
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
        let mut context = Context::default();
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        
        let mut worker = SingleWorker::new(&args, &context);
        
        // Barcode at start
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (processed_opt, tag_opt, _stats) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");
        
        assert_eq!(tag_opt, Some("barcode_bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let mut context = Context::default();
         
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
        let mut context = Context::default();
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
        let mut context = Context::default();
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
        let mut context = Context::default();
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
        let mut context = Context::default();
        context.format = crate::types::FileFormat::Fastq;
         
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
        let mut context = Context::default();
        context.format = crate::types::FileFormat::Fastq;
         
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
        let mut context = Context::default();
        context.format = crate::types::FileFormat::Fastq;
         
        let mut worker = SingleWorker::new(&args, &context);

        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (processed_opt, _, _) = worker.process(&record);
        let processed = processed_opt.expect("Should not be dropped by length filter");

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::default();
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
        let mut context = Context::default();
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
        let mut context = Context::default();
         
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