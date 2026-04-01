use crate::steps::generic::PerReadAlignmentResults;

pub fn apply_post_trim(per_read: &mut PerReadAlignmentResults, post_trim_length: usize) {
    if per_read.dropped || post_trim_length == 0 {
        return;
    }
    let current_length = per_read.read_range.end - per_read.read_range.start;
    if current_length > post_trim_length {
        per_read.bases_post_trim += current_length - post_trim_length;
        per_read.read_range.end = per_read.read_range.start + post_trim_length;
    }
}
