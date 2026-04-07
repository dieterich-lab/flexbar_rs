use crate::args::Args;
use crate::types::PairedRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;
use crate::steps::utils::quality_trim::QualityTrim;

pub struct PairedQualityTrim {
    pub inner: QualityTrim,
}

pub struct PairedQualityTrimState {}

impl GenericStep<PairedRecordDto> for PairedQualityTrim {
    type State = PairedQualityTrimState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        Self {
            inner: QualityTrim::new(args, context),
        }
    }

    fn create_state(&self) -> Self::State {
        PairedQualityTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        for i in 0..2 {
            self.inner.trim_read(
                &rec.record.reads[i], Some(i),
                &mut rec.alignment_results.per_read[i],
                &mut rec.stats
            );
        }
        rec
    }
}
