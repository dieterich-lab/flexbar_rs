use crate::args::Args;
use crate::types::SingeRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};

use std::any::Any;

pub struct SingleUncalledSequenceFilter {
    max_uncalled: usize,
}

impl SingleUncalledSequenceFilter {
    pub fn new(args: &Args) -> Self {
        Self {
            max_uncalled: args.max_uncalled,
        }
    }
}

pub struct SingleUncalledSequenceFilterSetup;

impl GenericStep<SingeRecordDto> for SingleUncalledSequenceFilter {
    type Setup = Box<dyn Any + Send + Sync>;
    fn create_setup(&self) -> Box<dyn Any + Send + Sync> {
        Box::new(SingleUncalledSequenceFilterSetup)
    }

    fn clone_setup(&self, _setup: &Box<dyn Any + Send + Sync>) -> Box<dyn Any + Send + Sync> {
        Box::new(SingleUncalledSequenceFilterSetup)
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>, _setup: &mut Box<dyn Any + Send + Sync>) -> RecordInTransformation<SingeRecordDto> {
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
