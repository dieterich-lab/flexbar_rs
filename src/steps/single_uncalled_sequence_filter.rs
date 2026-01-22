use crate::config::Config;
use crate::types::RecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};

pub struct SingleUncalledSequenceFilter {
    max_uncalled: usize,
}

impl SingleUncalledSequenceFilter {
    pub fn new(config: &Config) -> Self {
        Self {
            max_uncalled: config.args.max_uncalled,
        }
    }
}

impl GenericStep<RecordDto> for SingleUncalledSequenceFilter {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if let Some(ref record) = rec.record {
            let uncalled_count = record.seq().iter().filter(|&&b| b == b'N' || b == b'n').count();
            if uncalled_count > self.max_uncalled {
                rec.record = None;
                rec.stats.total_discarded_uncalled += 1;
            }
        }
        rec
    }
}
