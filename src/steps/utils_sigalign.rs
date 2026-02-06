use log::trace;
use sigalign::{Aligner, algorithms::Local, algorithms::SemiGlobal, Reference};
use sigalign::results::{Alignment, AlignmentOperation as SigOp, QueryAlignment};
use crate::logging::log_fatal;

pub fn extract_umi_tags_sig(
    alignment: &Alignment,
    reference: &[u8],
    read: &[u8],
) -> Vec<String> {
    let mut umi_tags = Vec::new();
    let mut umi = String::new();
    let mut x = alignment.position.query.0 as usize;
    let mut y = alignment.position.target.0 as usize;

    for op_count in &alignment.operations {
        match op_count.operation {
            SigOp::Match | SigOp::Subst => {
                for _ in 0..op_count.count {
                    if reference[y] == b'N' {
                        umi.push(read[x] as char);
                    } else if !umi.is_empty() {
                        umi_tags.push(umi);
                        umi = String::new();
                    }
                    x += 1;
                    y += 1;
                }
            }
            SigOp::Insertion => {
                x += op_count.count as usize;
            }
            SigOp::Deletion => {
                y += op_count.count as usize;
            }
        }
    }
    if !umi.is_empty() {
        umi_tags.push(umi);
    }
    umi_tags
}

pub fn format_alignment_sig(alignment: &Alignment, query: &[u8], target: &[u8]) -> (String, String, String) {
    let mut line1 = String::new();
    let mut line2 = String::new();
    let mut line3 = String::new();

    let mut q_pos = alignment.position.query.0 as usize;
    let mut t_pos = alignment.position.target.0 as usize;

    for op_count in &alignment.operations {
        match op_count.operation {
            SigOp::Match => {
                for _ in 0..op_count.count {
                    let q_char = query[q_pos];
                    let t_char = target[t_pos];
                    line1.push(q_char as char);
                    line3.push(t_char as char);
                    if q_char == b'N' || t_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('|');
                    }
                    q_pos += 1;
                    t_pos += 1;
                }
            }
            SigOp::Subst => {
                for _ in 0..op_count.count {
                    let q_char = query[q_pos];
                    let t_char = target[t_pos];
                    line1.push(q_char as char);
                    line3.push(t_char as char);
                    if q_char == b'N' || t_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('\\');
                    }
                    q_pos += 1;
                    t_pos += 1;
                }
            }
            SigOp::Insertion => {
                for _ in 0..op_count.count {
                    line1.push(query[q_pos] as char);
                    line2.push(' ');
                    line3.push(' ');
                    q_pos += 1;
                }
            }
            SigOp::Deletion => {
                for _ in 0..op_count.count {
                    line1.push(' ');
                    line2.push(' ');
                    line3.push(target[t_pos] as char);
                    t_pos += 1;
                }
            }
        }
    }
    (line1, line2, line3)
}

pub fn trace_alignment_sig(alignment: &Alignment, query: &[u8], target: &[u8]) {
    let (line1, line2, line3) = format_alignment_sig(alignment, query, target);

    trace!(
        "Alignment with penalty {}:\n  {}\n  {}\n  {}",
        alignment.penalty, line1, line2, line3
    );
}

pub struct AlignerAlgorithmParameters {
    pub mismatch: u32,
    pub gap_open: u32,
    pub gap_extend: u32,
    pub min_length: u32,
    pub max_penalty_per_length: f32,
}

pub fn get_local_aligner_from_parameters(
    params: &AlignerAlgorithmParameters) -> Aligner<Local>
{
    let algorithm = Local::new(
        params.mismatch,
        params.gap_open,
        params.gap_extend,
        params.min_length,
        params.max_penalty_per_length,
    );
    Aligner::new(algorithm.unwrap())
}

pub fn get_semiglobal_aligner_from_parameters(params: &AlignerAlgorithmParameters)
    -> Aligner<SemiGlobal>
{
    let algorithm = SemiGlobal::new(
        params.mismatch,
        params.gap_open,
        params.gap_extend,
        params.min_length,
        params.max_penalty_per_length,
    );
    Aligner::new(algorithm.unwrap())
}

pub fn usize_to_u32(x: usize, desc: &str) -> u32 {
    if x > u32::MAX as usize {
        log_fatal(&format!("Parameter {} too large for 32 bit: {}", desc, x));
    }
    x as u32
}