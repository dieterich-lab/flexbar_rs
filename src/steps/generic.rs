use crate::statistics::WorkerStatistics;
use crate::types::SingeRecordDto;

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

impl RecordInTransformation<SingeRecordDto> {
    pub fn get_sequence_as_str(&self) -> String {
        String::from_utf8_lossy(&self.record.seq()[
            self.alignment_results.seq_start..self.alignment_results.seq_end
        ]).to_string()
    }
}

pub trait Transformable {
    fn get_record_in_transformation(&self) -> RecordInTransformation<Self> where Self: Sized;
}

impl Transformable for SingeRecordDto {
    fn get_record_in_transformation(&self) -> RecordInTransformation<SingeRecordDto> {
        RecordInTransformation {
            record: self.clone(),
            alignment_results: AlignmentResults {
                barcode_tag: None,
                umi_tags: Vec::new(),
                removal_tags: Vec::new(),
                seq_start: 0,
                seq_end: self.seq().len(),
                record_dropped: false,
            },
            stats: WorkerStatistics::default(),
        }
    }   
}


pub trait GenericStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> RecordInTransformation<RecordType>;
}

pub trait FinalizeStep<RecordType>: Send + Sync {
    fn apply(&self, rec: RecordInTransformation<RecordType>) -> (Option<RecordType>, Option<String>, WorkerStatistics);
}

