use rust_htslib::bam::record::{Record, Aux, CigarString, Cigar};
use crate::records::generic_read::GenericRawRead;
use std::ops::Range;
use std::cell::OnceCell;

#[derive(Clone, Debug)]
pub struct BamRecordWrapper {
    pub record: Record,
    pub seq_cache: OnceCell<Vec<u8>>,
}

impl BamRecordWrapper {
    pub fn new(record: Record) -> Self {
        Self {
            record,
            seq_cache: OnceCell::new(),
        }
    }
}

impl GenericRawRead for BamRecordWrapper {
    #[inline(always)]
    fn seq(&self) -> &[u8] {
        self.seq_cache.get_or_init(|| self.record.seq().as_bytes())
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        Some(self.record.qual())
    }

    #[inline(always)]
    fn id(&self) -> &str {
        std::str::from_utf8(self.record.qname()).unwrap_or("unknown")
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        let mut new_record = self.record.clone();
        let cigar = self.record.cigar().take();
        new_record.set(self.record.qname(), Some(&cigar), seq, self.record.qual());
        Self::new(new_record)
    }

    fn finalize(
        &self,
        range: &Range<usize>,
        removal_tags: &Vec<String>,
        umi_tags: &Vec<String>,
        paired_read_dropped: bool,
        bam_desc_tag: &str,
    ) -> Self {
        let mut new_record = Record::new();
        let seq = self.record.seq().as_bytes();
        let qual = self.record.qual();
        
        let new_seq = &seq[range.clone()];
        let new_qual = &qual[range.clone()];
        
        let (new_pos, new_cigar) = adapt_cigar(&self.record, range);

        new_record.set(self.record.qname(), Some(&new_cigar), new_seq, new_qual);
        new_record.set_pos(new_pos);

        // Copy other fields if necessary, or let adapt_cigar handle it.
        // For now we set basic fields.

        let mut desc_parts: Vec<String> = Vec::new();
        for tag in removal_tags {
            desc_parts.push(tag.clone());
        }
        for tag in umi_tags {
            desc_parts.push(tag.clone());
        }
        if paired_read_dropped {
            desc_parts.push("single".to_string());
        }

        if !desc_parts.is_empty() {
            let desc = desc_parts.join("_");
            new_record.push_aux(bam_desc_tag.as_bytes().try_into().unwrap(), Aux::String(&desc)).unwrap();
        }

        Self::new(new_record)
    }
}

fn adapt_cigar(old_record: &Record, range: &Range<usize>) -> (i64, CigarString) {
    let old_pos = old_record.pos() as usize;
    let start = range.start;
    let end = range.end;

    let new_pos = if old_pos < start {
        0
    } else if old_pos > end {
        end - start
    } else {
        old_pos - start
    };
    let mut new_cigar: Vec<Cigar> = old_record.cigar().iter().cloned().collect();
    if start > old_pos {
        new_cigar = cut_cigar_at_start(new_cigar, (start - old_pos) as u32);
    }
    let covered_length = get_cigar_covered_length(&new_cigar);
    if new_pos + covered_length as usize > end - start {
        let new_covered_length = (end - start) as u32 - new_pos as u32;
        new_cigar = shorten_cigar_to_length(new_cigar, new_covered_length);
    }
    (new_pos as i64, CigarString(new_cigar))
}

fn cut_cigar_at_start(cigar: Vec<Cigar>, delta: u32) -> Vec<Cigar> {
    let mut new_cigar = Vec::new();
    let mut need_to_delete_length = 0;
    let mut remaining_delta = delta;
    for op in cigar {
        let (op_length, del_length) = get_op_covered_and_delete_length(op);
        if remaining_delta > 0 {
            need_to_delete_length += del_length;
            if op_length <= remaining_delta {
                remaining_delta -= op_length;
                need_to_delete_length += op_length;
                continue;
            }
        }
        // We reached an operation, which at least partially overlaps with the new cigar
        if remaining_delta <= 0 {  // If we are done, just copy the remaining operation
            new_cigar.push(op);
            continue;
        }

        // Partial operation - may only happen with op_length > 0
        let length_left = op_length - remaining_delta;
        need_to_delete_length += remaining_delta;
        remaining_delta = 0;

        // Insert initial delete if we still need to
        if need_to_delete_length > 0 {
            new_cigar.push(Cigar::Del(need_to_delete_length));
            need_to_delete_length = 0;
        }

        new_cigar.push(match op {
            Cigar::Match(_) => Cigar::Match(length_left),
            Cigar::Ins(_) => Cigar::Ins(length_left),
            Cigar::SoftClip(_) => Cigar::SoftClip(length_left),
            Cigar::Equal(_) => Cigar::Equal(length_left),
            Cigar::Diff(_) => Cigar::Diff(length_left),
            _ => unreachable!(),
        });
    }
    if need_to_delete_length > 0 {
        new_cigar.push(Cigar::Del(need_to_delete_length));
    }
    new_cigar
}

fn get_op_covered_and_delete_length(op: Cigar) -> (u32, u32) {
    match op {
        Cigar::Match(x) |
        Cigar::Ins(x) |
        Cigar::SoftClip(x) |
        Cigar::Equal(x) |
        Cigar::Diff(x) => (x, 0),
        Cigar::Del(x) => (0, x),
        _ => (0, 0),
    }
}

fn get_cigar_covered_length(cigar: &[Cigar]) -> u32 {
    cigar.iter()
        .filter_map(|op| match op {
            Cigar::Match(x) |
            Cigar::Ins(x) |
            Cigar::SoftClip(x) |
            Cigar::Equal(x) |
            Cigar::Diff(x) => Some(x),
            _ => None,
        })
        .sum()
}

fn shorten_cigar_to_length(cigar: Vec<Cigar>, length: u32) -> Vec<Cigar> {
    let mut new_cigar = Vec::new();
    let mut covered_length: u32 = 0;
    let mut length_to_delete: u32 = 0;
    for op in cigar {
        let (op_length, del_length) =get_op_covered_and_delete_length(op);
        if covered_length + op_length < length
            || del_length == 0 && covered_length + op_length == length
        {
            covered_length += op_length;
            new_cigar.push(op);
            continue;
        }
        if covered_length >= length {
            length_to_delete += op_length + del_length;
            continue;
        }
        // We are in an operation overlaping with the new end
        let remaining_length = length - covered_length;
        match op {
            Cigar::Del(_) => length_to_delete += remaining_length,
            Cigar::Match(_) => new_cigar.push(Cigar::Match(remaining_length)),
            Cigar::Ins(_) => new_cigar.push(Cigar::Ins(remaining_length)),
            Cigar::SoftClip(_) => new_cigar.push(Cigar::SoftClip(remaining_length)),
            Cigar::Equal(_) => new_cigar.push(Cigar::Equal(remaining_length)),
            Cigar::Diff(_) => new_cigar.push(Cigar::Diff(remaining_length)),
            _ => unreachable!(),
        }
        length_to_delete += op_length - remaining_length;
    }
    new_cigar.push(Cigar::Del(length_to_delete));
    new_cigar
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Record, CigarString};

    fn create_test_record(pos: i64, cigar: &str) -> Record {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mut record = Record::new();
        let cigar_obj: CigarString = cigar.as_bytes().try_into().unwrap();
        record.set(b"test", Some(&cigar_obj), seq, qual);
        record.set_pos(pos);
        record
    }

    #[test]
    fn test_adapt_cigar_remaining_sequence_before_alignment() {
        let record = create_test_record(12, "2D3M2I4M");
        let range = 1..7;
        let (new_pos, new_cigar) = adapt_cigar(&record, &range);
        assert_eq!(new_pos, 6);
        assert_eq!(format!("{}", new_cigar), "11D");
    }

    #[test]
    fn test_adapt_cigar_remaining_sequence_after_alignment() {
        let record = create_test_record(12, "2D3M2I4M");
        let range = 23..40;
        let (new_pos, new_cigar) = adapt_cigar(&record, &range);
        assert_eq!(new_pos, 0);
        assert_eq!(format!("{}", new_cigar), "11D");
    }

    #[test]
    fn test_adapt_cigar_remaining_sequence_includes_alignment() {
        let record = create_test_record(12, "2D3M2I4M");
        let range = 10..30;
        let (new_pos, new_cigar) = adapt_cigar(&record, &range);
        assert_eq!(new_pos, 2);
        assert_eq!(format!("{}", new_cigar), "2D3M2I4M");
    }

    #[test]
    fn test_adapt_cigar_remaining_sequence_skips_start_of_alignment() {
        let record = create_test_record(12, "2D3M2I4M");
        let range = 13..30;
        let (new_pos, new_cigar) = adapt_cigar(&record, &range);
        assert_eq!(new_pos, 0);
        assert_eq!(format!("{}", new_cigar), "3D2M2I4M");
    }

    #[test]
    fn test_adapt_cigar_remaing_sequence_truncates_alignment() {
        let record = create_test_record(12, "2D3M2I4M");
        let range = 0..19;
        let (new_pos, new_cigar) = adapt_cigar(&record, &range);
        assert_eq!(new_pos, 12);
        assert_eq!(format!("{}", new_cigar), "2D3M2I2M2D");
    }
}
