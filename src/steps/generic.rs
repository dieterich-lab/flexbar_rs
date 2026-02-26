use crate::statistics::WorkerStatistics;
use crate::types::{SingeRecordDto, TrimEnd};
use crate::args::Args;
use crate::context::Context;

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
    pub fn get_sequence(&self) -> &[u8] {
        &self.record.seq()[self.alignment_results.seq_start..self.alignment_results.seq_end]
    }
    
    pub fn get_qual(&self) -> &[u8] {
        &self.record.qual()[self.alignment_results.seq_start..self.alignment_results.seq_end]
    }
    
    pub fn get_sequence_as_str(&self) -> String {
        String::from_utf8_lossy(self.get_sequence()).to_string()
    }

    pub fn trim(&mut self, trim_end: TrimEnd, start: usize, end: usize) {
        match trim_end {
            TrimEnd::Right | TrimEnd::RTail => self.alignment_results.seq_end =
                self.alignment_results.seq_start + start,
            TrimEnd::Left | TrimEnd::LTail => self.alignment_results.seq_start += end,
            TrimEnd::Any => {
                let right_length = self.get_sequence().len() - end;
                if start > right_length {
                    self.alignment_results.seq_end = self.alignment_results.seq_start + start;
                } else {
                    self.alignment_results.seq_start += end;
                }
            }
        }
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
    type Setup;
    fn create_setup(args: &Args, context: &Context) -> Self::Setup;
    fn new(setup: &Self::Setup) -> Self;
    fn apply(&mut self, rec: RecordInTransformation<RecordType>, setup: &Self::Setup) -> RecordInTransformation<RecordType>;
}

pub struct FinalRecord<RecordType> {
    pub record: Option<RecordType>,
    pub output_tag: Option<String>,
    pub stats: WorkerStatistics,
}

