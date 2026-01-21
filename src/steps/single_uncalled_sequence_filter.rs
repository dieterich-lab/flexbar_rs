use crate::config::Config;
use crate::types::RecordDto;
use crate::statistics::WorkerStatistics;

pub struct SingleUncalledSequenceFilter {
    max_uncalled: usize,
}

impl SingleUncalledSequenceFilter {
    pub fn new(config: &Config) -> Self {
        Self {
            max_uncalled: config.args.max_uncalled,
        }
    }

    pub fn filter(&self, record: &RecordDto, stats: &mut WorkerStatistics) -> bool {
        let uncalled_count = record.seq().iter().filter(|&&b| b == b'N' || b == b'n').count();
        if uncalled_count > self.max_uncalled {
            stats.total_discarded_uncalled = 1;
            return true;
        }
        false
    }
}
