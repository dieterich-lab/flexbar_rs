use crate::types::RecordDto;
use crate::steps::generic::{RecordInTransformation, PerReadAlignmentResults, AlignmentResults};
use crate::statistics::WorkerStatistics;



impl RecordDto {
    pub fn get_record_in_transformation(&self) -> RecordInTransformation {
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
