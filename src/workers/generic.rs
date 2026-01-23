use crate::steps::generic::{GenericStep, FinalizeStep, RecordInTransformation, AlignmentResults};
use crate::statistics::WorkerStatistics;
use std::sync::Arc;

pub struct GenericWorker<RecordType> {
    pub steps: Vec<Arc<dyn GenericStep<RecordType>>>,
    pub finalize_step: Arc<dyn FinalizeStep<RecordType>>,
}

pub trait HasSequence {
    fn seq(&self) -> &[u8];
}

impl HasSequence for crate::types::RecordDto {
    fn seq(&self) -> &[u8] {
        self.seq()
    }
}

impl<RecordType: Clone + HasSequence> GenericWorker<RecordType> {
    pub fn new(steps: Vec<Arc<dyn GenericStep<RecordType>>>, finalize_step: Arc<dyn FinalizeStep<RecordType>>) -> Self {
        Self {steps, finalize_step}
    }

    pub fn process(&self, record: &RecordType) -> (Option<RecordType>, Option<String>, WorkerStatistics) {
        let mut rec_in = RecordInTransformation::<RecordType> {
            record: record.clone(),
            alignment_results: AlignmentResults {
                barcode_tag: None,
                umi_tags: Vec::new(),
                removal_tags: Vec::new(),
                seq_start: 0,
                seq_end: record.seq().len(),
                record_dropped: false,
            },
            stats: WorkerStatistics::default(),
        };

        for step in &self.steps {
            rec_in = step.apply(rec_in);
            if rec_in.alignment_results.record_dropped {
                break;
            }
        }

        self.finalize_step.apply(rec_in)
    }
}
