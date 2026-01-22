use crate::statistics::WorkerStatistics;

pub struct RecordInTransformation<RecordType> {
    pub record: Option<RecordType>,
    pub barcode_tag: Option<String>,
    pub stats: WorkerStatistics,
}

pub trait GenericStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> RecordInTransformation<RecordType>;
}

