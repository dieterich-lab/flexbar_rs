use crate::args::Args;
use crate::types::PairedRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;

pub struct PairedFilterLength {
    pub min_read_length: usize,
    pub single_reads: bool,
    pub single_reads_paired: bool,
}

pub struct PairedFilterLengthState {}

impl GenericStep<PairedRecordDto> for PairedFilterLength {
    type State = PairedFilterLengthState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            min_read_length: args.min_read_length,
            single_reads: args.single_reads,
            single_reads_paired: args.single_reads_paired,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedFilterLengthState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        if self.min_read_length == 0 {
            return rec;
        }

        let len1 = if rec.alignment_results.per_read[0].dropped { 0 } else {
            rec.alignment_results.per_read[0].read_range.end - rec.alignment_results.per_read[0].read_range.start
        };
        let len2 = if rec.alignment_results.per_read[1].dropped { 0 } else {
            rec.alignment_results.per_read[1].read_range.end - rec.alignment_results.per_read[1].read_range.start
        };

        let short1 = !rec.alignment_results.per_read[0].dropped && len1 < self.min_read_length;
        let short2 = !rec.alignment_results.per_read[1].dropped && len2 < self.min_read_length;

        if short1 || short2 {
            if self.single_reads_paired {
                if short1 {
                    rec.alignment_results.per_read[0].read_range = 0..1;
                    rec.stats.total_single_reads_replaced += 1;
                    trace!("Length filter replaced R1 with N: {} (length: {})", rec.record.id(), len1);
                }
                if short2 {
                    rec.alignment_results.per_read[1].read_range = 0..1;
                    rec.stats.total_single_reads_replaced += 1;
                    trace!("Length filter replaced R2 with N: {} (length: {})", rec.record.id(), len2);
                }
            } else if self.single_reads {
                if short1 {
                    rec.alignment_results.per_read[0].dropped = true;
                    rec.stats.total_single_reads_written += 1;
                    trace!("Length filter dropped R1, keeping R2: {} (length: {})", rec.record.id(), len1);
                }
                if short2 {
                    rec.alignment_results.per_read[1].dropped = true;
                    rec.stats.total_single_reads_written += 1;
                    trace!("Length filter dropped R2, keeping R1: {} (length: {})", rec.record.id(), len2);
                }
            } else {
                rec.alignment_results.per_read[0].dropped = true;
                rec.alignment_results.per_read[1].dropped = true;
                rec.stats.total_discarded_length += 1;
                rec.stats.total_single_reads_dropped += 2;
                trace!("Length filter dropped paired record: {} (lengths: {}, {})", rec.record.id(), len1, len2);
            }
        }
        rec
    }
}
