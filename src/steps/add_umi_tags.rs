use crate::types::RecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};

pub struct AddUmiTags;

impl AddUmiTags {
    pub fn new() -> Self {
        Self
    }
}

impl GenericStep<RecordDto> for AddUmiTags {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if let Some(record) = rec.record {
            trace!("AddUmiTags processing record: {}", record.id());
            if !rec.alignment_results.umi_tags.is_empty() {
                trace!("AddUmiTags adding UMIs: {:?} to record: {}", rec.alignment_results.umi_tags, record.id());
                let mut new_id = record.id().to_string();
                for tag in &rec.alignment_results.umi_tags {
                    new_id.push('_');
                    new_id.push_str(tag);
                }
                rec.record = Some(RecordDto::with_attrs(
                    &new_id,
                    record.desc(),
                    record.seq(),
                    record.qual(),
                ));
            } else {
                rec.record = Some(record);
            }
        }
        rec
    }
}
