use crate::steps::generic::{Transformable, RecordInTransformation};
use crate::statistics::WorkerStatistics;

pub struct WorkerResult<Record> {
    pub(crate) record: Option<Record>,
    pub(crate) output_tag: Option<String>,
    pub(crate) stats: WorkerStatistics,
}
pub trait GenericWorker<RecordType: Transformable, StepId, Setup: Sized> {
    fn apply_step(
        &mut self,
        step_id: &StepId,
        mut record: RecordInTransformation<RecordType>,
        setup: &Setup
    ) -> RecordInTransformation<RecordType>;
    fn finalize(
        &self,
        record: &mut RecordInTransformation<RecordType>,
        setup: &Setup
    ) -> (Option<RecordType>, Option<String>, WorkerStatistics);
}

pub trait GenericWorkerSetup<RecordType: Transformable, StepId, Setup> {
    type Worker: GenericWorker<RecordType, StepId, Setup>;
    fn get_step_sequence(&self) -> &[StepId];
    fn create_worker(&self) -> Box<Self::Worker>;
}


impl<RecordType: Transformable, StepId, Setup: GenericWorkerSetup<RecordType, StepId, Setup>>
    dyn GenericWorker<RecordType, StepId, Setup>
{
    fn process(&mut self, record: RecordType, setup: &Setup)
        -> (Option<RecordType>, Option<String>, WorkerStatistics) {
        let mut rec_in = record.get_record_in_transformation();
        for step_id in setup.get_step_sequence() {
            rec_in = self.apply_step(step_id, &mut rec_in, setup);
            if rec_in.alignment_results.record_dropped {
                break;
            }
        }
        self.finalize(&mut rec_in, setup)
    }
}

struct WorkerFactory<
    's,
    RecordType: Transformable,
    StepId,
    Setup: GenericWorkerSetup<RecordType, StepId, Setup>+Sized,
> {
    cache: Option<Box<dyn GenericWorker<RecordType, StepId, Setup>>>,
    setup: &'s dyn GenericWorkerSetup<RecordType, StepId, Setup, Worker=Setup::Worker>,
}

impl <
    's,
    RecordType: Transformable, 
    StepId,
    Setup: GenericWorkerSetup<RecordType, StepId, Setup>+Sized
>
    WorkerFactory<'s, RecordType, StepId, Setup>
{
    pub fn new(setup: &'s Setup) -> Self {
        WorkerFactory {
            cache: None,
            setup,
        }
    }
    
    pub fn get_worker(&mut self) -> &mut dyn GenericWorker<RecordType, StepId, Setup> where <Setup as GenericWorkerSetup<RecordType, StepId, Setup>>::Worker: 'static {
        if self.cache.is_none() {
            self.cache = Some(self.setup.create_worker());
        }
        self.cache.as_mut().unwrap().as_mut()
    }
}