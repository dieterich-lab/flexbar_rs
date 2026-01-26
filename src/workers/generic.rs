use crate::steps::generic::{GenericStep, FinalizeStep, Transformable};
use crate::statistics::WorkerStatistics;
use std::sync::Arc;

pub struct GenericWorker<RecordType: Transformable> {
    pub steps: Vec<Arc<dyn GenericStep<RecordType>>>,
    pub finalize_step: Arc<dyn FinalizeStep<RecordType>>,
}

impl<RecordType: Clone + Transformable> GenericWorker<RecordType> {
    pub fn new(steps: Vec<Arc<dyn GenericStep<RecordType>>>, finalize_step: Arc<dyn FinalizeStep<RecordType>>) -> Self {
        Self {steps, finalize_step}
    }

    pub fn process(&self, record: &RecordType) -> (Option<RecordType>, Option<String>, WorkerStatistics) {
        let mut rec_in = record.get_record_in_transformation();

        for step in &self.steps {
            rec_in = step.apply(rec_in);
            if rec_in.alignment_results.record_dropped {
                break;
            }
        }

        self.finalize_step.apply(rec_in)
    }
}
