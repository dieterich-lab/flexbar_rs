use crate::types::{SingeRecordDto, TrimEnd};
use crate::steps::generic::{RecordInTransformation, PerReadAlignmentResults, Transformable, RecordId};
use crate::statistics::WorkerStatistics;
use crate::steps::utils::BarcodeAlignmentResults;

#[derive(Default, Clone)]
pub struct SingleAlignmentResults {
    pub barcode_tag: Option<String>,
    pub read: PerReadAlignmentResults,
}

impl BarcodeAlignmentResults for SingleAlignmentResults {
    fn set_barcode_tag(&mut self, tag: Option<String>) {
        self.barcode_tag = tag;
    }
}

impl AsMut<SingleAlignmentResults> for SingleAlignmentResults {
    fn as_mut(&mut self) -> &mut SingleAlignmentResults {
        self
    }
}

impl RecordId for SingeRecordDto {
    fn record_id(&self) -> &str {
        self.id()
    }
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

impl Transformable for SingeRecordDto {
    type AlignmentResults = SingleAlignmentResults;
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
