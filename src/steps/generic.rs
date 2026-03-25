use crate::statistics::WorkerStatistics;
use crate::types::{SingeRecordDto, PairedRecordDto, TrimEnd};
use crate::args::Args;
use crate::context::Context;
use std::ops::Range;

#[derive(Default, Clone)]
pub struct PerReadAlignmentResults {
    pub umi_tags: Vec<String>,
    pub removal_tags: Vec<String>,
    pub read_range: Range<usize>,
    pub dropped: bool,
}

#[derive(Default, Clone)]
pub struct SingleAlignmentResults {
    pub barcode_tag: Option<String>,
    pub read: PerReadAlignmentResults,
}

#[derive(Default, Clone)]
pub struct PairedAlignmentResults {
    pub barcode_tag: Option<String>,
    pub per_read: [PerReadAlignmentResults; 2],
}

pub struct RecordInTransformation<RecordType, AR> {
    pub record: RecordType,
    pub alignment_results: AR,
    pub stats: WorkerStatistics,
}

impl RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
    pub fn get_sequence(&self) -> &[u8] {
        &self.record.seq()[self.alignment_results.read.read_range.clone()]
    }
    
    pub fn get_qual(&self) -> &[u8] {
        &self.record.qual()[self.alignment_results.read.read_range.clone()]
    }
    
    pub fn get_sequence_as_str(&self) -> String {
        String::from_utf8_lossy(self.get_sequence()).to_string()
    }

    pub fn trim(&mut self, trim_end: TrimEnd, start: usize, end: usize) {
        match trim_end {
            TrimEnd::Right | TrimEnd::RTail => {
                self.alignment_results.read.read_range.end = self.alignment_results.read.read_range.start + start;
            }
            TrimEnd::Left | TrimEnd::LTail => {
                self.alignment_results.read.read_range.start += end;
            }
            TrimEnd::Any => {
                let current_len = self.alignment_results.read.read_range.end - self.alignment_results.read.read_range.start;
                let right_length = current_len - end;
                if start > right_length {
                    self.alignment_results.read.read_range.end = self.alignment_results.read.read_range.start + start;
                } else {
                    self.alignment_results.read.read_range.start += end;
                }
            }
        }
    }
}

pub trait Transformable {
    type AR;
    fn get_record_in_transformation(&self) -> RecordInTransformation<Self, Self::AR> where Self: Sized;
}

impl Transformable for SingeRecordDto {
    type AR = SingleAlignmentResults;
    fn get_record_in_transformation(&self) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        RecordInTransformation {
            record: self.clone(),
            alignment_results: SingleAlignmentResults {
                barcode_tag: None,
                read: PerReadAlignmentResults {
                    umi_tags: Vec::new(),
                    removal_tags: Vec::new(),
                    read_range: 0..self.seq().len(),
                    dropped: false,
                },
            },
            stats: WorkerStatistics::default(),
        }
    }   
}

impl Transformable for PairedRecordDto {
    type AR = PairedAlignmentResults;
    fn get_record_in_transformation(&self) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        RecordInTransformation {
            record: self.clone(),
            alignment_results: PairedAlignmentResults {
                barcode_tag: None,
                per_read: [
                    PerReadAlignmentResults {
                        umi_tags: Vec::new(),
                        removal_tags: Vec::new(),
                        read_range: 0..self.records[0].seq().len(),
                        dropped: false,
                    },
                    PerReadAlignmentResults {
                        umi_tags: Vec::new(),
                        removal_tags: Vec::new(),
                        read_range: 0..self.records[1].seq().len(),
                        dropped: false,
                    },
                ],
            },
            stats: WorkerStatistics::default(),
        }
    }
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

