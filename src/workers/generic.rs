use crate::steps::generic::{GenericStep, RecordInTransformation, AlignmentResults};
use crate::statistics::WorkerStatistics;
use std::sync::Arc;

pub struct GenericWorker<RecordType> {
    pub steps: Vec<Arc<dyn GenericStep<RecordType>>>,
}

impl<RecordType: Clone> GenericWorker<RecordType> {
    pub fn new(steps: Vec<Arc<dyn GenericStep<RecordType>>>) -> Self {
        Self {steps}
    }

    pub fn process(&self, record: &RecordType) -> RecordInTransformation<RecordType> {
        let mut rec_in = RecordInTransformation::<RecordType> {
            record: Some(record.clone()),
            alignment_results: AlignmentResults::default(),
            stats: WorkerStatistics::default(),
        };

        for step in &self.steps {
            rec_in = step.apply(rec_in);
        }

        rec_in
    }
}
