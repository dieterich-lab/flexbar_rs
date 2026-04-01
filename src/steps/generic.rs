use crate::statistics::WorkerStatistics;
use crate::args::Args;
use crate::context::Context;
use std::ops::Range;


#[derive(Default, Clone)]
pub struct PerReadAlignmentResults {
    pub umi_tags: Vec<String>,
    pub removal_tags: Vec<String>,
    pub read_range: Range<usize>,
    pub dropped: bool,
    pub bases_quality_trim: usize,
    pub bases_adapter_removal: usize,
    pub bases_barcode_removal: usize,
}

pub struct RecordInTransformation<RecordType, AR> {
    pub record: RecordType,
    pub alignment_results: AR,
    pub stats: WorkerStatistics,
}

pub trait RecordId {
    fn record_id(&self) -> &str;
}

impl<RecordType, AR> RecordInTransformation<RecordType, AR> where RecordType: RecordId {
    pub fn record_id(&self) -> &str {
        self.record.record_id()
    }
}

pub trait Transformable {
    type AlignmentResults;
    fn get_record_in_transformation(&self) -> RecordInTransformation<Self, Self::AlignmentResults> where Self: Sized;
}

pub trait GenericStep<RecordType>: Send + Sync {
    type State: Send;
    type AlignmentResult;
    fn new(args: &Args, context: &Context) -> Self;
    fn create_state(&self) -> Self::State;
    fn apply(
        &self, rec: RecordInTransformation<RecordType, Self::AlignmentResult>, state: &mut Self::State
    ) -> RecordInTransformation<RecordType, Self::AlignmentResult>;
}

pub struct FinalRecord<RecordType> {
    pub record: Option<RecordType>,
    pub output_tag: Option<String>,
    pub stats: WorkerStatistics,
}
