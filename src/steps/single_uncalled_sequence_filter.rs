use crate::config::Config;
use crate::types::RecordDto;
use log::trace;
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
        if !rec.alignment_results.record_dropped {
            trace!("Uncalled filter processing record: {}", rec.record.id());
            let seq = &rec.record.seq()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
            let uncalled_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
            if uncalled_count > self.max_uncalled {
                trace!("Uncalled filter dropped record: {}", rec.record.id());
                rec.alignment_results.record_dropped = true;
                rec.stats.total_discarded_uncalled += 1;
            }
        }
        rec
    }
}
