use crate::args::Args;
use crate::types::{SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::context::{Context};
use crate::multi_aligner::{MultiAlignerState};
use crate::steps::utils::adapter_remover::AdapterRemover;

pub struct SingleAdapterRemover {
    pub inner: AdapterRemover,
}

pub struct SingleAdapterRemoverState {
    pub multi_aligner_state: MultiAlignerState,
}

impl GenericStep<SingeRecordDto> for SingleAdapterRemover {
    type State = SingleAdapterRemoverState;
    type AlignmentResult = SingleAlignmentResults;
    
    fn new(args: &Args, context: &Context) -> Self {
        Self {
            inner: AdapterRemover::new(args, &context.adapters),
        }
    }

    fn create_state(&self) -> Self::State {
        Self::State {
            multi_aligner_state: self.inner.create_state(),
        }
    }

    fn apply(
        &self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, state: &mut Self::State
    ) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        if rec.alignment_results.read.dropped { return rec }
        
        if self.inner.trace {
            trace!("Adapter remover processing record: {}", rec.record.id());
        }
        
        self.inner.remove_adapters(
            &mut state.multi_aligner_state,
            &rec.record,
            None,
            &mut rec.alignment_results.read,
            &mut rec.stats
        );
        
        rec
    }
}
