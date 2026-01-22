use bio::alignment::{Alignment, AlignmentOperation};
use crate::types::{RecordDto, TrimEnd};

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

pub fn trim_sequence(
    record: &RecordDto,
    alignment: &Alignment,
    trim_end: TrimEnd,
) -> (Vec<u8>, Vec<u8>) {
    let mut new_seq = record.seq().to_vec();
    let mut new_qual = record.qual().to_vec();

    match trim_end {
        TrimEnd::Right | TrimEnd::RTail => {
            new_seq.truncate(alignment.ystart);
            new_qual.truncate(alignment.ystart);
        }
        TrimEnd::Left | TrimEnd::LTail => {
            new_seq = new_seq[alignment.yend..].to_vec();
            new_qual = new_qual[alignment.yend..].to_vec();
        }
        TrimEnd::Any => {
            if alignment.ystart < record.seq().len() - alignment.yend {
                new_seq = new_seq[alignment.yend..].to_vec();
                new_qual = new_qual[alignment.yend..].to_vec();
            } else {
                new_seq.truncate(alignment.ystart);
                new_qual.truncate(alignment.ystart);
            }
        }
    }
    (new_seq, new_qual)
}
