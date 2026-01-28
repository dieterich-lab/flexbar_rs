use bio::alignment::{Alignment, AlignmentOperation};
use crate::types::{TrimEnd};

pub fn extract_umi_tags(
    alignment: &Alignment,
    reference: &[u8],
    read: &[u8],
) -> Vec<String> {
    let mut umi_tags = Vec::new();
    let mut umi = String::new();
    let mut x = alignment.xstart;
    let mut y = alignment.ystart;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                if reference[x] == b'N' {
                    umi.push(read[y] as char);
                } else if !umi.is_empty() {
                    umi_tags.push(umi);
                    umi = String::new();
                }
                x += 1;
                y += 1;
            }
            AlignmentOperation::Del => {
                x += 1;
            }
            AlignmentOperation::Ins => {
                y += 1;
            }
            _ => {}
        }
    }
    if !umi.is_empty() {
        umi_tags.push(umi);
    }
    umi_tags
}

pub fn get_trim_range(
    len: usize,
    alignment: &Alignment,
    trim_end: TrimEnd,
) -> (usize, usize) {
    match trim_end {
        TrimEnd::Right => {
            (0, alignment.ystart)
        }
        TrimEnd::RTail => {
            (0, alignment.yend)
        }
        TrimEnd::Left => {
            (alignment.yend, len)
        }
        TrimEnd::LTail => {
            (alignment.ystart, len)
        }
        TrimEnd::Any => {
            let left_dist = alignment.ystart;
            let right_dist = len - alignment.yend;
            if left_dist <= right_dist {
                (alignment.yend, len)
            } else {
                (0, alignment.ystart)
            }
        }
    }
}
