use crate::steps::generic::PerReadAlignmentResults;

pub fn apply_pre_trim(per_read: &mut PerReadAlignmentResults, pre_trim_left: usize, pre_trim_right: usize) {
    if per_read.dropped || (pre_trim_left == 0 && pre_trim_right == 0) {
        return;
    }
    
    let mut current_length = per_read.read_range.end - per_read.read_range.start;
    
    if pre_trim_left > 0 {
        let trim_left = pre_trim_left.min(current_length);
        per_read.read_range.start += trim_left;
        per_read.bases_pre_trim += trim_left;
        current_length -= trim_left;
    }
    
    if pre_trim_right > 0 {
        let trim_right = pre_trim_right.min(current_length);
        per_read.read_range.end -= trim_right;
        per_read.bases_pre_trim += trim_right;
    }
}
