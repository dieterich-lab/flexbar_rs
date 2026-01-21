use crate::context::Context;
use crate::config::Config;
use crate::types::RecordDto;
use crate::steps::single_adapter_remover::SingleAdapterRemover;
use crate::steps::single_uncalled_sequence_filter::SingleUncalledSequenceFilter;

pub struct SingleWorker<'a> {
    adapter_remover: SingleAdapterRemover<'a>,
    uncalled_filter: SingleUncalledSequenceFilter,
}

impl<'a> SingleWorker<'a> {
    pub fn new(config: &'a Config, context: &'a Context) -> Self {
        Self {
            adapter_remover: SingleAdapterRemover::new(config, &context.adapters),
            uncalled_filter: SingleUncalledSequenceFilter::new(config),
        }
    }

    pub fn process(&self, record: &RecordDto) -> (Option<RecordDto>, crate::statistics::WorkerStatistics) {
        let mut worker_stats = crate::statistics::WorkerStatistics::default();
        if self.uncalled_filter.filter(record, &mut worker_stats) {
            return (None, worker_stats);
        }
        let processed_record = self.adapter_remover.remove_adapters(record, &mut worker_stats);
        (Some(processed_record), worker_stats)
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
        let (processed_opt, _worker_stats) = worker.process(&record);
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
        let (processed_opt, _worker_stats) = worker.process(&record);
        let processed = processed_opt.unwrap();
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_max_uncalled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let config = Config::new(args);
        let context = Context::new();
        let worker = SingleWorker::new(&config, &context);

        // 2 'N's - should keep
        let record = RecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let (_processed, worker_stats) = worker.process(&record);
        assert_eq!(worker_stats.total_discarded_uncalled, 0);

        // 3 'N's - should discard
        let record = RecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let (_processed, worker_stats) = worker.process(&record);
        assert_eq!(worker_stats.total_discarded_uncalled, 1);
    }
}