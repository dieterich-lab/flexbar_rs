use crate::context::Context;
use crate::config::Config;
use crate::types::{RecordDto, QualTrimType};
use crate::steps::single_adapter_remover::SingleAdapterRemover;
use crate::steps::single_uncalled_sequence_filter::SingleUncalledSequenceFilter;
use crate::steps::single_barcode_splitter::SingleBarcodeSplitter;
use crate::steps::single_quality_trim::SingleQualityTrim;
use crate::steps::add_umi_tags::AddUmiTags;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::workers::generic::GenericWorker;
use std::sync::Arc;

pub struct SingleWorker ( GenericWorker<RecordDto> );

impl SingleWorker {
    pub fn new(config: &Config, context: &Context) -> Self {
        let mut steps: Vec<Arc<dyn GenericStep<RecordDto>>> = Vec::new();

        if config.args.max_uncalled > 0 {
            steps.push(
                Arc::new(SingleUncalledSequenceFilter::new(config))
            );
        }
        if config.args.qtrim != QualTrimType::Off && !config.args.qtrim_post_removal {
            steps.push(
                Arc::new(SingleQualityTrim::new(config))
            );
        }
        if !context.barcodes.is_empty() {
            steps.push(
                Arc::new(SingleBarcodeSplitter::new(config, &context.barcodes))
            );
        }
        if !context.adapters.is_empty() {
            steps.push(
                Arc::new(SingleAdapterRemover::new(config, &context.adapters))
            );
        }

        if config.args.qtrim != QualTrimType::Off && config.args.qtrim_post_removal {
            steps.push(
                Arc::new(SingleQualityTrim::new(config))
            )
        }

        if config.args.umi_tags {
            steps.push(
                Arc::new(AddUmiTags::new())
            );
        }

        Self (GenericWorker::new(steps))
    }
    
    pub fn process(&self, record: &RecordDto) -> RecordInTransformation<RecordDto> {
        self.0.process(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::Context;
    use crate::args::Args;
    use crate::config::Config;
    use crate::types::RecordDto;
    use clap::Parser;

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_left() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.barcodes.push(("bc1".to_string(), "GGGG".to_string()));
        
        let worker = SingleWorker::new(&config, &context);
        
        // Barcode at start
        let record = RecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
        assert_eq!(*res.stats.barcode_counts.get("bc1").unwrap(), 1);
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        // 2 'N's - should keep
        let record = RecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let res = worker.process(&record);
        assert_eq!(res.stats.total_discarded_uncalled, 0);

        // 3 'N's - should discard
        let record = RecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let res = worker.process(&record);
        assert_eq!(res.stats.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.barcodes.push(("bc1".to_string(), "GGGG".to_string()));
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(processed.id(), "read1_Flexbar_removal");
    }

    #[test]
    fn test_removal_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--removal-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(processed.id(), "read1_Flexbar_removal_0");
    }

    #[test]
    fn test_removal_tags_both() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", 
            "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4",
            "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5",
            "--removal-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.barcodes.push(("bc1".to_string(), "GGGG".to_string()));
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"GGGGACGTCTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        // barcode removal appends _Flexbar_removal
        // adapter removal appends _Flexbar_removal_0
        assert_eq!(processed.id(), "read1_Flexbar_removal_Flexbar_removal_0");
    }

    #[test]
    fn test_quality_trim_tail() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "TAIL", "--qtrim-threshold", "20"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        // Threshold 20 + 33 (sanger) = 53 ('5')
        let record = RecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let res = worker.process(&record);
        let processed = res.record.unwrap();

        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_quality_trim_win() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "WIN", "--qtrim-threshold", "20", "--qtrim-win-size", "5"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        // Threshold 53. Window size 5.
        // "IIIII" sum = 5 * 73 = 365. Threshold sum = 5 * 53 = 265.
        // "44444" sum = 5 * 52 = 260.
        let record = RecordDto::with_attrs("read1", None, b"ACTGTACGTACGT", b"IIIII44444III");
        let res = worker.process(&record);
        let processed = res.record.unwrap();

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_quality_trim_bwa() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "BWA", "--qtrim-threshold", "20"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        let record = RecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let res = worker.process(&record);
        let processed = res.record.unwrap();

        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.barcodes.push(("bc1".to_string(), "NNGGGGNN".to_string()));
        
        let worker = SingleWorker::new(&config, &context);
        
        // Barcode NNGGGGNN matches AAGGGGTT
        let record = RecordDto::with_attrs("read1", None, b"AAGGGGTTACGT", b"IIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert_eq!(res.alignment_results.umi_tags, vec!["AA".to_string(), "TT".to_string()]);
        assert_eq!(processed.id(), "read1_AA_TT");
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_umi_tags_adapter() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATNNNN", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4", "--umi-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATNNNN".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"ACTGTAGATCCCC", b"IIIIIIIIIIIII");
        let res = worker.process(&record);
        let processed = res.record.unwrap();
        
        assert_eq!(res.alignment_results.umi_tags, vec!["CCCC".to_string()]);
        assert_eq!(processed.id(), "read1_CCCC");
        assert_eq!(processed.seq(), b"ACTGT");
    }
}