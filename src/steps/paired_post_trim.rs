use crate::args::Args;
use crate::types::PairedRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;
use crate::steps::utils::post_trim::apply_post_trim;

pub struct PairedPostTrim {
    pub post_trim_length: usize,
}

pub struct PairedPostTrimState {}

impl GenericStep<PairedRecordDto> for PairedPostTrim {
    type State = PairedPostTrimState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            post_trim_length: args.post_trim_length,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedPostTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        for i in 0..2 {
            apply_post_trim(&mut rec.alignment_results.per_read[i], self.post_trim_length);
        }
        rec
    }
}
