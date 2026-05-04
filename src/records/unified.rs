use serde::{Deserialize, Serialize};

pub use crate::records::unified_read::UnifiedRead;
pub use crate::records::generic_read::GenericRawRead;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnifiedRecord {
    pub reads: [Option<UnifiedRead>; 2],
    pub barcode_read: Option<UnifiedRead>
}

macro_rules! implement_by_read_index {
    ($name: ident, $ReturnType: ty) => {
        #[inline(always)]
        pub fn $name(&self, read_index: usize) -> $ReturnType {
            match read_index {
                0 => {
                    if let Some(read) = &self.reads[0] {
                        read.$name()
                    } else {
                        panic!("Read index {:?} not present", read_index)
                    }
                }
                1 => {
                    if let Some(read) = &self.reads[1] {
                        read.$name()
                    } else {
                        panic!("Read index {:?} not present", read_index)
                    }
                }
                _ => panic!("Invalid read index {:?}", read_index),
            }
        }
    };
}


impl UnifiedRecord {
    pub fn single(read: UnifiedRead) -> Self {
        Self { reads: [Some(read), None], barcode_read: None }
    }

    pub fn single_barcoded(read: UnifiedRead, barcode: UnifiedRead) -> Self {
        Self { reads: [Some(read), None], barcode_read: Some(barcode) }
    }

    pub fn paired(r1: UnifiedRead, r2: UnifiedRead) -> Self {
        Self { reads: [Some(r1), Some(r2)], barcode_read: None }
    }

    pub fn paired_barcoded(r1: UnifiedRead, r2: UnifiedRead, barcode: UnifiedRead) -> Self {
        Self { reads: [Some(r1), Some(r2)], barcode_read: Some(barcode) }
    }

    pub fn is_paired(&self) -> bool {
        self.reads[1].is_some()
    }

    implement_by_read_index!(id, &str);
    implement_by_read_index!(seq, &[u8]);
    implement_by_read_index!(desc, Option<&str>);
    implement_by_read_index!(qual, Option<&[u8]>);

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
        let new_read = UnifiedRead::new_fastq(
            old_read.id(),
            old_read.desc_str(),
            new_sequence,
            old_read.qual().unwrap_or(&[]),
        );
        self.reads[read_index] = Some(new_read);
    }

    #[inline(always)]
    pub fn patch_barcode_read(&mut self, new_sequence: &[u8]) {
        let old_barcode = self.barcode_read.take().unwrap();
        let new_barcode = UnifiedRead::new_fastq(
            old_barcode.id(),
            old_barcode.desc_str(),
            new_sequence,
            old_barcode.qual().unwrap_or(&[]),
        );
        self.barcode_read = Some(new_barcode);
    }
}