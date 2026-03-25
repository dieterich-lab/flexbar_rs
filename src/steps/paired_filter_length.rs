use crate::args::Args;
use crate::types::PairedRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, PairedAlignmentResults};
use crate::context::Context;

pub struct PairedFilterLength {
    pub min_read_length: usize,
}

pub struct PairedFilterLengthState {}

impl GenericStep<PairedRecordDto> for PairedFilterLength {
    type State = PairedFilterLengthState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            min_read_length: args.min_read_length,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedFilterLengthState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        if !rec.alignment_results.record_dropped && self.min_read_length > 0 {
            let len1 = rec.alignment_results.per_read[0].read_range.end - rec.alignment_results.per_read[0].read_range.start;
            let len2 = rec.alignment_results.per_read[1].read_range.end - rec.alignment_results.per_read[1].read_range.start;

            if len1 < self.min_read_length || len2 < self.min_read_length {
                trace!("Length filter dropped paired record: {} (lengths: {}, {})", rec.record.id(), len1, len2);
                rec.alignment_results.record_dropped = true;
                rec.stats.total_discarded_length += 1;
            }
        }
        rec
    }
}
