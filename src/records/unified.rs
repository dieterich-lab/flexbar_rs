
pub use crate::records::unified_read::UnifiedRead;
pub use crate::records::generic_read::GenericRawRead;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::AlignmentResults;

#[derive(Clone, Debug)]
pub struct UnifiedRecord {
    pub reads: [Option<UnifiedRead>; 2],
    pub barcode_read: Option<UnifiedRead>
}

pub struct AnnotatedRecord {
    pub record: UnifiedRecord,
    pub barcode_id: Option<String>,
    pub removal_tags: [Vec<String>; 2],
    pub umi_tags: Vec<String>,
    pub paired_read_dropped: bool,
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
    // implement_by_read_index!(seq, &[u8]);
    // implement_by_read_index!(desc, Option<&str>);
    // implement_by_read_index!(qual, Option<&[u8]>);

    #[inline(always)]
    pub fn record_id(&self) -> &str {
        if let Some(read) = &self.reads[0] {
            read.id()
        } else {
            panic!("Record has no reads")
        }
    }

    #[inline(always)]
    pub fn patch_read_sequence(&mut self, read_index: usize, new_sequence: &[u8]) {
        let old_read = self.reads[read_index].as_ref().unwrap();
        self.reads[read_index] = Some(old_read.set_seq(new_sequence));
    }

    #[inline(always)]
    pub fn patch_barcode_read(&mut self, new_sequence: &[u8]) {
        let old_barcode = self.barcode_read.as_ref().unwrap();
        self.barcode_read = Some(old_barcode.set_seq(new_sequence));
    }

    pub fn finalize(
        self,
        alignment_results: AlignmentResults,
        umi_tags: &Vec<String>,
        mut stats: WorkerStatistics,
    ) -> (AnnotatedRecord, WorkerStatistics) {
        let mut reads: [Option<UnifiedRead>; 2] = [None, None];
        let mut removal_tags: [Vec<String>; 2] = [vec![], vec![]];

        for i in 0..2 {
            if let Some(old_read) = &self.reads[i] {
                let original_len = old_read.seq().len();
                let results = alignment_results.per_read[i].as_ref().unwrap();

                if results.dropped {
                    stats.bases_discarded += original_len;
                    continue;
                } else {
                    stats.bases_quality_trim += results.bases_quality_trim;
                    stats.bases_adapter_removal += results.bases_adapter_removal;
                    stats.bases_barcode_removal += results.bases_barcode_removal;
                    stats.bases_post_trim += results.bases_post_trim;
                    stats.bases_pre_trim += results.bases_pre_trim;
                }

                removal_tags[i] = results.removal_tags.clone();
                reads[i] = Some(old_read.trim(&results.read_range));
            }
        }

        let paired_read_dropped = self.is_paired() && (
            alignment_results.per_read[0].as_ref().map_or(false, |r| r.dropped) ||
            alignment_results.per_read[1].as_ref().map_or(false, |r| r.dropped)
        );
        
        let result = AnnotatedRecord {
            record: UnifiedRecord { reads, barcode_read: None },
            barcode_id: alignment_results.barcode_id,
            removal_tags,
            umi_tags: umi_tags.clone(),
            paired_read_dropped,
        };

        (result, stats)
    }
}