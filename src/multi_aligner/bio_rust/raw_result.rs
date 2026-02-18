use bio::alignment::{Alignment, AlignmentOperation};
use crate::multi_aligner::generic::{AlignerMode, GenericRawResult};

impl GenericRawResult for Alignment {
    fn validate(&self, mode: AlignerMode, y: &[u8], error_rate: f32, min_overlap: usize) -> bool {
        let required_overlap = match mode {
            AlignerMode::Local => y.len(),
            _ => min_overlap,
        };
        let overlap = self.xend - self.xstart;
        let allowed_errors = (overlap as f32 * error_rate) as usize;
        check_match(self, y, required_overlap, allowed_errors)
    }

    fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String> {
        let mut umi_tags = Vec::new();
        let mut umi = String::new();
        let mut x_pos = self.xstart;
        let mut y_pos = self.ystart;

        for op in &self.operations {
            match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    if y[y_pos] == b'N' {
                        umi.push(x[x_pos] as char);
                    } else if !umi.is_empty() {
                        umi_tags.push(umi);
                        umi = String::new();
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Del => {
                    x_pos += 1;
                }
                AlignmentOperation::Ins => {
                    y_pos += 1;
                }
                _ => {}
            }
        }
        if !umi.is_empty() {
            umi_tags.push(umi);
        }
        umi_tags
    }

    fn format(&self, x: &[u8], y: &[u8]) -> String {
        let mut line1 = String::new();
        let mut line2 = String::new();
        let mut line3 = String::new();

        let mut x_pos = 0;
        let mut y_pos = 0;

        for op in &self.operations {
            match op {
                AlignmentOperation::Match => {
                    let x_char = x[x_pos];
                    let y_char = y[y_pos];
                    line1.push(x_char as char);
                    line3.push(y_char as char);
                    if x_char == b'N' || y_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('|');
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Subst => {
                    let x_char = x[x_pos];
                    let y_char = y[y_pos];
                    line1.push(x_char as char);
                    line3.push(y_char as char);
                    if x_char == b'N' || y_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('\\');
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Ins => {
                    line1.push(x[x_pos] as char);
                    line2.push(' ');
                    line3.push(' ');
                    x_pos += 1;
                }
                AlignmentOperation::Del => {
                    line1.push(' ');
                    line2.push(' ');
                    line3.push(y[y_pos] as char);
                    y_pos += 1;
                }
                AlignmentOperation::Xclip(len) => {
                    for _ in 0..*len {
                        line1.push(x[x_pos] as char);
                        line2.push(' ');
                        line3.push(' ');
                        x_pos += 1;
                    }
                }
                AlignmentOperation::Yclip(len) => {
                    for _ in 0..*len {
                        line1.push(' ');
                        line2.push(' ');
                        line3.push(y[y_pos] as char);
                        y_pos += 1;
                    }
                }
            }
        }
        format!("\
            Alignment with score {}:\n  {}\n  {}\n  {}",
                self.score, line1, line2, line3
        )
    }
}
fn check_match(
    alignment: &Alignment,
    y: &[u8],
    required_overlap: usize,
    allowed_errors: usize
) -> bool {
    let mut y_pos = alignment.ystart;
    let mut bases_to_check = required_overlap;
    let mut free_errors: i32 = allowed_errors as i32;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {
                y_pos += 1;
                bases_to_check -= 1;
            }
            AlignmentOperation::Subst => {
                if y[y_pos] != b'N' {
                    free_errors -= 1;
                }
                y_pos += 1;
                bases_to_check -= 1;
            }
            AlignmentOperation::Ins => {
                bases_to_check -= 1;
                free_errors -= 1;
            }
            AlignmentOperation::Del => {
                y_pos += 1;
                free_errors -= 1;
            }
            _ => {}
        }
        if free_errors < 0 {
            return false
        }
        if bases_to_check <= 0 {
            return true
        }
    }
    bases_to_check <= free_errors as usize
}
