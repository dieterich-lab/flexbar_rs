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
    pub bases_post_trim: usize,
    pub bases_pre_trim: usize,
}

#[derive(Default, Clone)]
pub struct AlignmentResults {
    pub barcode_tag: Option<String>,
    pub per_read: [Option<PerReadAlignmentResults>; 2],
}

pub struct RecordInTransformation<RecordType> {
    pub record: RecordType,
    pub alignment_results: AlignmentResults,
    pub stats: WorkerStatistics,
}

pub trait RecordId {
    fn record_id(&self) -> &str;
}

impl<RecordType> RecordInTransformation<RecordType> where RecordType: RecordId {
    pub fn record_id(&self) -> &str {
        self.record.record_id()
    }
}

pub trait Transformable {
    fn get_record_in_transformation(&self) -> RecordInTransformation<Self> where Self: Sized;
}

pub trait GenericStep<RecordType>: Send + Sync {
    type State: Send;
    fn new(args: &Args, context: &Context) -> Self;
    fn create_state(&self) -> Self::State;
    fn apply(
        &self, rec: RecordInTransformation<RecordType>, state: &mut Self::State
    ) -> RecordInTransformation<RecordType>;
}

pub struct FinalRecord<RecordType> {
    pub record: Option<RecordType>,
    pub output_tag: Option<String>,
    pub stats: WorkerStatistics,
}
