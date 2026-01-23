use crate::statistics::WorkerStatistics;

#[derive(Default, Clone)]
pub struct AlignmentResults {
    pub barcode_tag: Option<String>,
    pub umi_tags: Vec<String>,
    pub removal_tags: Vec<String>,
    pub seq_start: usize,
    pub seq_end: usize,
    pub record_dropped: bool,
}

pub struct RecordInTransformation<RecordType> {
    pub record: RecordType,
    pub alignment_results: AlignmentResults,
    pub stats: WorkerStatistics,
}

pub trait GenericStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> RecordInTransformation<RecordType>;
}

pub trait FinalizeStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> (Option<RecordType>, Option<String>, WorkerStatistics);
}

