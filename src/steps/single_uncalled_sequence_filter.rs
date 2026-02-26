use crate::args::Args;
use crate::types::SingeRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;

use std::any::Any;

pub struct SingleUncalledSequenceFilter {}

#[derive(Clone)]
pub struct SingleUncalledSequenceFilterSetup {
    pub max_uncalled: usize,
}

impl GenericStep<SingeRecordDto> for SingleUncalledSequenceFilter {
    type Setup = SingleUncalledSequenceFilterSetup;

    fn create_setup(args: &Args, _context: &Context) -> SingleUncalledSequenceFilterSetup {
        SingleUncalledSequenceFilterSetup {
            max_uncalled: args.max_uncalled,
        }
    }

    fn new(_setup: &SingleUncalledSequenceFilterSetup) -> Self {
        Self {}
    }

    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &SingleUncalledSequenceFilterSetup) -> RecordInTransformation<SingeRecordDto> {
        if !rec.alignment_results.record_dropped {
            trace!("Uncalled filter processing record: {}", rec.record.id());
            let seq = &rec.record.seq()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
            let uncalled_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
            if uncalled_count > setup.max_uncalled {
                trace!("Uncalled filter dropped record: {}", rec.record.id());
                rec.alignment_results.record_dropped = true;
                rec.stats.total_discarded_uncalled += 1;
            }
        }
        rec
    }
}
