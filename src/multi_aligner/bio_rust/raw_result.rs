use std::ops::Range;
use bio::alignment::{Alignment, AlignmentOperation};
use crate::multi_aligner::generic::GenericRawResult;

impl GenericRawResult for Alignment {
    fn validate(
        &self, x: &[u8], y: &[u8], error_rate: f32, min_overlap: Option<usize>, x_offset: usize, y_range: Range<usize>
    ) -> bool {
        let y_part = &y[y_range];
        let check_length = if self.xstart + x_offset > self.ystart {
            std::cmp::min(y_part.len(), x.len() - self.xstart - x_offset)
        } else {
            std::cmp::min(y_part.len(), y_part.len() - (self.ystart - self.xstart - x_offset))
        };
        match min_overlap {
            Some(x) => if x > check_length {return false},
            None => ()
        }
        let allowed_errors = (check_length as f32 * error_rate) as usize;
        check_match(self, y_part, check_length, allowed_errors)
    }

    fn get_umi_tags(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> Vec<String> {
        let y_part = &y[y_range];
        let mut umi_tags = Vec::new();
        let mut umi = String::new();
        let mut x_pos = self.xstart + x_offset;
        let mut y_pos = self.ystart;

        for op in &self.operations {
            match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    if y_pos < y_part.len() && y_part[y_pos] == b'N' {
                        if x_pos < x.len() {
                            umi.push(x[x_pos] as char);
                        }
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

    fn format(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> String {
        let y_part = &y[y_range];
        let mut line1 = String::new();
        let mut line2 = String::new();
        let mut line3 = String::new();

        let mut x_pos = 0;
        let mut y_pos = 0;

        for _ in 0..x_offset {
            if x_pos < x.len() {
                line1.push(x[x_pos] as char);
                line2.push('´');
                line3.push(' ');
                x_pos += 1;
            }
        }
        for op in &self.operations {
            match op {
                AlignmentOperation::Match => {
                    let x_char = if x_pos < x.len() { x[x_pos] } else { b' ' };
                    let y_char = if y_pos < y_part.len() { y_part[y_pos] } else { b' ' };
                    line1.push(if x_pos < x.len() { x_char as char } else { '-' });
                    line3.push(if y_pos < y_part.len() { y_char as char } else { '-' });
                    if x_char == b'N' || y_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('|');
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Subst => {
                    let x_char = if x_pos < x.len() { x[x_pos] } else { b' ' };
                    let y_char = if y_pos < y_part.len() { y_part[y_pos] } else { b' ' };
                    line1.push(if x_pos < x.len() { x_char as char } else { '-' });
                    line3.push(if y_pos < y_part.len() { y_char as char } else { '-' });
                    if x_char == b'N' || y_char == b'N' {
                        line2.push('.');
                    } else {
                        line2.push('\\');
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Ins => {
                    line1.push(if x_pos < x.len() { x[x_pos] as char } else { '-' });
                    line2.push(' ');
                    line3.push(' ');
                    x_pos += 1;
                }
                AlignmentOperation::Del => {
                    line1.push(' ');
                    line2.push(' ');
                    line3.push(if y_pos < y_part.len() { y_part[y_pos] as char } else { '-' });
                    y_pos += 1;
                }
                AlignmentOperation::Xclip(len) => {
                    for _ in 0..*len {
                        line1.push(if x_pos < x.len() { x[x_pos] as char } else { '-' });
                        line2.push('´');
                        line3.push(' ');
                        x_pos += 1;
                    }
                }
                AlignmentOperation::Yclip(len) => {
                    for _ in 0..*len {
                        line1.push(' ');
                        line2.push(',');
                        line3.push(if y_pos < y_part.len() { y_part[y_pos] as char } else { '-' });
                        y_pos += 1;
                    }
                }
            }
        }
        // Append remaining x if any
        while x_pos < x.len() {
            line1.push(x[x_pos] as char);
            line2.push(' ');
            line3.push(' ');
            x_pos += 1;
        }

        format!("\
            Alignment with score {}:\n  {}\n  {}\n  {}\n",
                self.score, line1, line2, line3
        )
    }
}
fn check_match(
    alignment: &Alignment,
    y: &[u8],
    check_length: usize,
    allowed_errors: usize
) -> bool {
    let mut y_pos = alignment.ystart;
    let mut bases_to_check = check_length;
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
                bases_to_check -= 1;
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
