use crate::context::Context;
use crate::config::Config;
use crate::types::RecordDto;
use crate::steps::single_adapter_remover::SingleAdapterRemover;
use crate::steps::single_uncalled_sequence_filter::SingleUncalledSequenceFilter;
use crate::steps::single_barcode_splitter::SingleBarcodeSplitter;

pub struct SingleWorker<'a> {
    adapter_remover: SingleAdapterRemover<'a>,
    uncalled_filter: SingleUncalledSequenceFilter,
    barcode_splitter: SingleBarcodeSplitter<'a>,
}

impl<'a> SingleWorker<'a> {
    pub fn new(config: &'a Config, context: &'a Context) -> Self {
        Self {
            adapter_remover: SingleAdapterRemover::new(config, &context.adapters),
            uncalled_filter: SingleUncalledSequenceFilter::new(config),
            barcode_splitter: SingleBarcodeSplitter::new(config, &context.barcodes),
        }
    }

    pub fn process(&self, record: &RecordDto) -> (Option<RecordDto>, crate::statistics::WorkerStatistics, Option<String>) {
        let mut worker_stats = crate::statistics::WorkerStatistics::default();
        if self.uncalled_filter.filter(record, &mut worker_stats) {
            return (None, worker_stats, None);
        }
        let (processed_record, barcode_tag) = self.barcode_splitter.split_barcodes(record, &mut worker_stats);
        let processed_record = self.adapter_remover.remove_adapters(&processed_record, &mut worker_stats);
        (Some(processed_record), worker_stats, barcode_tag)
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
        let (processed_opt, _worker_stats, _tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
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
        let (processed_opt, _worker_stats, _tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
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
        let (processed_opt, worker_stats, tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
        assert_eq!(tag, Some("bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
        assert_eq!(*worker_stats.barcode_counts.get("bc1").unwrap(), 1);
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        // 2 'N's - should keep
        let record = RecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let (_processed, worker_stats, _tag) = worker.process(&record);
        assert_eq!(worker_stats.total_discarded_uncalled, 0);

        // 3 'N's - should discard
        let record = RecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let (_processed, worker_stats, _tag) = worker.process(&record);
        assert_eq!(worker_stats.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.barcodes.push(("bc1".to_string(), "GGGG".to_string()));
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = RecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (processed_opt, _worker_stats, _tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
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
        let (processed_opt, _worker_stats, _tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
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
        let (processed_opt, _worker_stats, _tag) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
        // barcode removal appends _Flexbar_removal
        // adapter removal appends _Flexbar_removal_0
        assert_eq!(processed.id(), "read1_Flexbar_removal_Flexbar_removal_0");
    }
}