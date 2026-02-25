use crate::steps::generic::{Transformable, RecordInTransformation};
use crate::statistics::WorkerStatistics;



pub trait GenericWorker<RecordType: Transformable> {
    fn process(&mut self, record: &RecordType) -> (Option<RecordType>, Option<String>, WorkerStatistics);
}