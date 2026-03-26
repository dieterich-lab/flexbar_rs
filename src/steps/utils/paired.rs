use crate::types::{PairedRecordDto, TrimEnd};
use crate::steps::generic::{RecordInTransformation, PerReadAlignmentResults, Transformable, RecordId};
use crate::statistics::WorkerStatistics;
use crate::steps::utils::BarcodeAlignmentResults;

#[derive(Default, Clone)]
pub struct PairedAlignmentResults {
    pub barcode_tag: Option<String>,
    pub per_read: [PerReadAlignmentResults; 2],
}

impl BarcodeAlignmentResults for PairedAlignmentResults {
    fn set_barcode_tag(&mut self, tag: Option<String>) {
        self.barcode_tag = tag;
    }
}

impl AsMut<PairedAlignmentResults> for PairedAlignmentResults {
    fn as_mut(&mut self) -> &mut PairedAlignmentResults {
        self
    }
}

impl RecordId for PairedRecordDto {
    fn record_id(&self) -> &str {
        self.id()
    }
}

impl RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
    pub fn get_sequence(&self, read_index: usize) -> &[u8] {
        &self.record.records[read_index].seq()[self.alignment_results.per_read[read_index].read_range.clone()]
    }

    pub fn trim(&mut self, read_index: usize, trim_end: TrimEnd, start: usize, end: usize) {
        let range = &mut self.alignment_results.per_read[read_index].read_range;
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

impl Transformable for PairedRecordDto {
    type AlignmentResults = PairedAlignmentResults;
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
