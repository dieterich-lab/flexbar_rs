use crate::statistics::WorkerStatistics;
use crate::steps::generic::{FinalizeStep, RecordInTransformation};
use crate::types::RecordDto;

pub struct SingleFinalizeStep;

impl FinalizeStep<RecordDto> for SingleFinalizeStep {
    fn apply(&self, rec: RecordInTransformation<RecordDto>) -> (Option<RecordDto>, Option<String>, WorkerStatistics) {
        if rec.alignment_results.record_dropped {
            return (None, None, rec.stats);
        }

        let mut desc_parts = Vec::new();
        match {rec.record.desc()} {
            Some(desc) => desc_parts.push(desc),
            None => (),
        }

        for tag in &rec.alignment_results.removal_tags {
            desc_parts.push(tag);
        }

        for tag in &rec.alignment_results.umi_tags {
            desc_parts.push(tag);
        }

        let new_desc: Option<&str> = match desc_parts.is_empty() {
            true => None,
            false => Some(&desc_parts.join("_")),
        };

        let new_seq = rec.record.seq()[rec.alignment_results.seq_start..rec.alignment_results.seq_end].to_vec();
        let new_qual = rec.record.qual()[rec.alignment_results.seq_start..rec.alignment_results.seq_end].to_vec();

        let final_record = RecordDto::with_attrs(
            rec.record.id(),
            new_desc,
            &new_seq,
            &new_qual,
        );

        (Some(final_record), rec.alignment_results.barcode_tag, rec.stats)
    }
}
