use crate::statistics::WorkerStatistics;
use crate::steps::generic::PerReadAlignmentResults;

pub fn finalize_read_stats(stats: &mut WorkerStatistics, per_read: &PerReadAlignmentResults, original_len: usize) {
    if per_read.dropped {
        stats.bases_discarded += original_len;
    } else {
        stats.bases_quality_trim += per_read.bases_quality_trim;
        stats.bases_adapter_removal += per_read.bases_adapter_removal;
        stats.bases_barcode_removal += per_read.bases_barcode_removal;
    }
}
