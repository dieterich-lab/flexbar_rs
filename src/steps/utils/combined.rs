use crate::types::{RecordDto, TrimEnd};
use crate::steps::generic::{RecordInTransformation, PerReadAlignmentResults, Transformable, RecordId, AlignmentResults};
use crate::statistics::WorkerStatistics;
use crate::steps::utils::BarcodeAlignmentResults;

impl BarcodeAlignmentResults for AlignmentResults {
    fn set_barcode_tag(&mut self, tag: Option<String>) {
        self.barcode_tag = tag;
    }
}

impl AsMut<AlignmentResults> for AlignmentResults {
    fn as_mut(&mut self) -> &mut AlignmentResults {
        self
    }
}

impl RecordId for RecordDto {
    fn record_id(&self) -> &str {
        match self {
            RecordDto::Single(s) => s.id(),
            RecordDto::Paired(p) => p.id(),
        }
    }
}

impl RecordInTransformation<RecordDto> {
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

impl Transformable for RecordDto {
    fn get_record_in_transformation(&self) -> RecordInTransformation<RecordDto> {
        match self {
            RecordDto::Single(s) => {
                RecordInTransformation {
                    record: self.clone(),
                    alignment_results: AlignmentResults {
                        barcode_tag: None,
                        per_read: [
                            Some(PerReadAlignmentResults {
                                read_range: 0..s.seq().len(),
                                ..Default::default()
                            }),
                            None,
                        ],
                    },
                    stats: WorkerStatistics::default(),
                }
            }
            RecordDto::Paired(p) => {
                RecordInTransformation {
                    record: self.clone(),
                    alignment_results: AlignmentResults {
                        barcode_tag: None,
                        per_read: [
                            Some(PerReadAlignmentResults {
                                read_range: 0..p.reads[0].seq().len(),
                                ..Default::default()
                            }),
                            Some(PerReadAlignmentResults {
                                read_range: 0..p.reads[1].seq().len(),
                                ..Default::default()
                            }),
                        ],
                    },
                    stats: WorkerStatistics::default(),
                }
            }
        }
    }
}
