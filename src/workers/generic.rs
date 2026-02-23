use crate::steps::generic::{GenericStep, FinalizeStep, Transformable};
use crate::statistics::WorkerStatistics;
use std::sync::Arc;
use std::any::Any;

pub struct GenericWorker<RecordType: Transformable> {
    pub steps: Vec<Arc<dyn GenericStep<RecordType, Setup = Box<dyn Any + Send + Sync>>>>,
    pub finalize_step: Arc<dyn FinalizeStep<RecordType>>,
    pub setups: Vec<Box<dyn Any + Send + Sync>>,
}

impl<RecordType: Transformable> Clone for GenericWorker<RecordType> {
    fn clone(&self) -> Self {
        let setups = self.steps.iter().zip(self.setups.iter())
            .map(|(step, setup)| step.clone_setup(setup))
            .collect();
        Self {
            steps: self.steps.clone(),
            finalize_step: self.finalize_step.clone(),
            setups,
        }
    }
}

impl<RecordType: Clone + Transformable + Send + Sync + 'static> GenericWorker<RecordType> {
    pub fn new(steps: Vec<Arc<dyn GenericStep<RecordType, Setup = Box<dyn Any + Send + Sync>>>>, finalize_step: Arc<dyn FinalizeStep<RecordType>>) -> Self {
        let setups = steps.iter().map(|s| s.create_setup()).collect();
        Self {steps, finalize_step, setups}
    }

    pub fn process(&mut self, record: &RecordType) -> (Option<RecordType>, Option<String>, WorkerStatistics) {
        let mut rec_in = record.get_record_in_transformation();

        for (step, setup) in self.steps.iter().zip(self.setups.iter_mut()) {
            rec_in = step.apply(rec_in, setup);
            if rec_in.alignment_results.record_dropped {
                break;
            }
        }

        self.finalize_step.apply(rec_in)
    }
}
