use crate::records::unified::{GenericRawRead, UnifiedRecord};
use std::ops::Range;
use crate::statistics::WorkerStatistics;
use crate::args::Args;
use crate::context::Context;

#[derive(Default, Clone)]
pub struct PerReadAlignmentResults {
    pub umi_tags: Vec<String>,
    pub removal_tags: Vec<String>,
    pub read_range: Range<usize>,
    pub dropped: bool,
    // per-read stats - deleted bases
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

use crate::types::TrimEnd;

pub struct RecordInTransformation {
    pub record: UnifiedRecord,
    pub alignment_results: AlignmentResults,
    pub stats: WorkerStatistics,
}

impl RecordInTransformation {
    pub fn from_dto(record: UnifiedRecord) -> Self {
        let mut per_read: [Option<PerReadAlignmentResults>; 2] = [None, None];
        if let Some(read) = &record.reads[0] {
            per_read[0] = Some(PerReadAlignmentResults {
                read_range: 0..read.seq().len(),
                ..Default::default()
            });
        }
        if let Some(read) = &record.reads[1] {
            per_read[1] = Some(PerReadAlignmentResults {
                read_range: 0..read.seq().len(),
                ..Default::default()
            });
        }

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
        self.record.id(0)
    }

    pub fn read_id(&self, read_index: usize) -> String {
        if self.record.reads[1].is_some() {
            let id_str = self.record.reads[read_index].as_ref().unwrap().id();
            format!("{} (R{})", id_str, read_index + 1).to_string()
        } else {
            self.record.reads[read_index].as_ref().unwrap().id().to_string()
        }
    }

    pub fn get_sequence(&self, read_index: usize) -> &[u8] {
        if let Some(read) = &self.record.reads[read_index] {
            &read.seq()[self.alignment_results.per_read[read_index].as_ref().unwrap().read_range.clone()]
        } else {
            panic!("Invalid read index {}", read_index);
        }
    }

    pub fn get_qual(&self, read_index: usize) -> Option<&[u8]> {
        if let Some(read) = &self.record.reads[read_index] {
            match &read.qual() {
                Some(q) => Some(&q[self.alignment_results.per_read[read_index].as_ref().unwrap().read_range.clone()]),
                None => None
            }
        } else {
            panic!("Invalid read index {}", read_index);
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

    #[inline]
    pub fn get_read_length(&self, read_index: usize) -> usize {
        match self.alignment_results.per_read.get(read_index) {
            Some(Some(x)) => x.read_range.end - x.read_range.start,
            _ => panic!("Invalid read index {}", read_index),
        }
    }

    #[inline]
    pub fn drop_read(&mut self, read_index: usize) {
        self.alignment_results.per_read[read_index].as_mut().unwrap().dropped = true;
    }

    #[inline]
    pub fn is_dropped_read(&self, read_index: usize) -> bool {
        self.alignment_results.per_read[read_index].as_ref().unwrap().dropped
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
