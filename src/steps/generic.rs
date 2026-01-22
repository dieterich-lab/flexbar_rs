use crate::statistics::WorkerStatistics;

#[derive(Default, Clone)]
pub struct AlignmentResults {
    pub barcode_tag: Option<String>,
    pub umi_tags: Vec<String>,
}

pub struct RecordInTransformation<RecordType> {
    pub record: Option<RecordType>,
    pub alignment_results: AlignmentResults,
    pub stats: WorkerStatistics,
}

pub trait GenericStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> RecordInTransformation<RecordType>;
}

