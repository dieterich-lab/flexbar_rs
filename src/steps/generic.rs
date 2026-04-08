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

impl AlignmentResults {
    pub fn set_barcode_tag(&mut self, tag: Option<String>) {
        self.barcode_tag = tag;
    }
}

use crate::types::{RecordDto, TrimEnd};

pub struct RecordInTransformation {
    pub record: RecordDto,
    pub alignment_results: AlignmentResults,
    pub stats: WorkerStatistics,
}

impl RecordInTransformation {
    pub fn from_dto(record: RecordDto) -> Self {
        let per_read = match &record {
            RecordDto::Single(s) => {
                let len = s.seq().len();
                [
                    Some(PerReadAlignmentResults {
                        read_range: 0..len,
                        ..Default::default()
                    }),
                    None,
                ]
            }
            RecordDto::Paired(p) => {
                let len0 = p.reads[0].seq().len();
                let len1 = p.reads[1].seq().len();
                [
                    Some(PerReadAlignmentResults {
                        read_range: 0..len0,
                        ..Default::default()
                    }),
                    Some(PerReadAlignmentResults {
                        read_range: 0..len1,
                        ..Default::default()
                    }),
                ]
            }
        };

        RecordInTransformation {
            record,
            alignment_results: AlignmentResults {
                barcode_tag: None,
                per_read,
            },
            stats: WorkerStatistics::default(),
        }
    }

    pub fn record_id(&self) -> &str {
        // Prefer unified helper on RecordDto
        self.record.id(None)
    }

    pub fn get_sequence(&self, read_index: usize) -> &[u8] {
        match &self.record {
            RecordDto::Single(s) => {
                if read_index != 0 { panic!("Single read only has index 0"); }
                &s.seq()[self.alignment_results.per_read[0].as_ref().unwrap().read_range.clone()]
            }
            RecordDto::Paired(p) => {
                &p.reads[read_index].seq()[self.alignment_results.per_read[read_index].as_ref().unwrap().read_range.clone()]
            }
        }
    }

    pub fn get_qual(&self, read_index: usize) -> &[u8] {
        match &self.record {
            RecordDto::Single(s) => {
                if read_index != 0 { panic!("Single read only has index 0"); }
                &s.qual()[self.alignment_results.per_read[0].as_ref().unwrap().read_range.clone()]
            }
            RecordDto::Paired(p) => {
                &p.reads[read_index].qual()[self.alignment_results.per_read[read_index].as_ref().unwrap().read_range.clone()]
            }
        }
    }

    pub fn trim(&mut self, read_index: usize, trim_end: TrimEnd, start: usize, end: usize) {
        let range = &mut self.alignment_results.per_read[read_index].as_mut().unwrap().read_range;
        match trim_end {
            TrimEnd::Right | TrimEnd::RTail => {
                range.end = range.start + start;
            }
            TrimEnd::Left | TrimEnd::LTail => {
                range.start += end;
            }
            TrimEnd::Any => {
                let current_len = range.end - range.start;
                let right_length = current_len - end;
                if start > right_length {
                    range.end = range.start + start;
                } else {
                    range.start += end;
                }
            }
        }
    }
}

pub trait GenericStep: Send + Sync {
    type State: Send;
    fn new(args: &Args, context: &Context) -> Self;
    fn create_state(&self) -> Self::State;
    fn apply(
        &self, rec: RecordInTransformation, state: &mut Self::State
    ) -> RecordInTransformation;
}

pub struct FinalRecord {
    pub record: Option<RecordDto>,
    pub output_tag: Option<String>,
    pub stats: WorkerStatistics,
}
