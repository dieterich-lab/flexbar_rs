use crate::args::Args;
use crate::types::SingeRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, SingleAlignmentResults};
use crate::context::Context;

pub struct SingleUncalledSequenceFilter {
    pub max_uncalled: usize,
}

pub struct SingleUncalledSequenceFilterState {}

impl GenericStep<SingeRecordDto> for SingleUncalledSequenceFilter {
    type State = SingleUncalledSequenceFilterState;
    type AlignmentResult = SingleAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            max_uncalled: args.max_uncalled,
        }
    }

    fn create_state(&self) -> Self::State {
        SingleUncalledSequenceFilterState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        if !rec.alignment_results.record_dropped {
            trace!("Uncalled filter processing record: {}", rec.record.id());
            let seq = &rec.record.seq()[rec.alignment_results.read.read_range.clone()];
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
