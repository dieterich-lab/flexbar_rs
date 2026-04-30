use bio::utils::TextSlice;
use serde::{Deserialize, Serialize};

pub type SingeRecordDto = bio::io::fastq::Record;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnifiedRecord {
    pub reads: [Option<SingeRecordDto>; 2],
    pub barcode_read: Option<SingeRecordDto>
}

macro_rules! implement_by_record_index {
    ($name: ident, $ReturnType: ty) => {
        #[inline(always)]
        pub fn $name(&self, read_index: usize) -> $ReturnType {
            match read_index {
                0 => {
                    if let Some(read) = &self.reads[0] {
                        read.$name()
                    } else {
                        panic!("Invalid read index {:?}", read_index)
                    }
                }
                1 => {
                    if let Some(read) = &self.reads[1] {
                        read.$name()
                    } else {
                        panic!("Invalid read index {:?}", read_index)
                    }
                }
                _ => panic!("Invalid read index {:?}", read_index),
            }
        }
    };
}


impl UnifiedRecord {
    pub fn single(read: SingeRecordDto) -> Self {
        Self { reads: [Some(read), None], barcode_read: None }
    }

    pub fn single_barcoded(read: SingeRecordDto, barcode: SingeRecordDto) -> Self {
        Self { reads: [Some(read), None], barcode_read: Some(barcode) }
    }

    pub fn paired(r1: SingeRecordDto, r2: SingeRecordDto) -> Self {
        Self { reads: [Some(r1), Some(r2)], barcode_read: None }
    }

    pub fn paired_barcoded(r1: SingeRecordDto, r2: SingeRecordDto, barcode: SingeRecordDto) -> Self {
        Self { reads: [Some(r1), Some(r2)], barcode_read: Some(barcode) }
    }

    pub fn is_paired(&self) -> bool {
        self.reads[1].is_some()
    }

    implement_by_record_index!(id, &str);
    implement_by_record_index!(seq, TextSlice);
    implement_by_record_index!(desc, Option<&str>);
    implement_by_record_index!(qual, &[u8]);

    #[inline(always)]
    pub fn record_id(&self) -> &str {
        if let Some(read) = &self.reads[0] {
            read.id()
        } else {
            panic!("RecordDto has no reads")
        }
    }

    #[inline(always)]
    pub fn patch_read_sequence(&mut self, read_index: usize, new_sequence: &[u8]) {
        let old_read = self.reads[read_index].take().unwrap();
        let new_read = SingeRecordDto::with_attrs(
            old_read.id(),
            old_read.desc(),
            new_sequence,
            old_read.qual(),
        );
        self.reads[read_index] = Some(new_read);
    }

    #[inline(always)]
    pub fn patch_barcode_read(&mut self, new_sequence: &[u8]) {
        let old_barcode = self.barcode_read.take().unwrap();
        let new_barcode = SingeRecordDto::with_attrs(
            old_barcode.id(),
            old_barcode.desc(),
            new_sequence,
            old_barcode.qual(),
        );
        self.barcode_read = Some(new_barcode);
    }
}