use std::ops::{Add, AddAssign};

#[derive(Debug, Default, Clone)]
pub struct WorkerStatistics {
    pub total_adapters_removed: usize,
    pub total_discarded_uncalled: usize,
    pub total_discarded_barcode: usize,
    pub total_discarded_length: usize,
    pub total_single_reads_written: usize,
    pub total_single_reads_replaced: usize,
    pub total_single_reads_dropped: usize,
    pub total_bases_read: usize,
    pub bases_quality_trim: usize,
    pub bases_discarded: usize,
    pub bases_adapter_removal: usize,
    pub bases_barcode_removal: usize,
    pub bases_post_trim: usize,
    pub bases_pre_trim: usize,
}

impl WorkerStatistics {
    pub fn is_discarded(&self) -> bool {
        self.total_discarded_uncalled > 0 || self.total_discarded_barcode > 0 || self.total_discarded_length > 0
    }
}

impl Add for WorkerStatistics {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self.total_adapters_removed += other.total_adapters_removed;
        self.total_discarded_uncalled += other.total_discarded_uncalled;
        self.total_discarded_barcode += other.total_discarded_barcode;
        self.total_discarded_length += other.total_discarded_length;
        self.total_single_reads_written += other.total_single_reads_written;
        self.total_single_reads_replaced += other.total_single_reads_replaced;
        self.total_single_reads_dropped += other.total_single_reads_dropped;
        self.total_bases_read += other.total_bases_read;
        self.bases_quality_trim += other.bases_quality_trim;
        self.bases_discarded += other.bases_discarded;
        self.bases_adapter_removal += other.bases_adapter_removal;
        self.bases_barcode_removal += other.bases_barcode_removal;
        self.bases_post_trim += other.bases_post_trim;
        self.bases_pre_trim += other.bases_pre_trim;
        self
    }
}

impl AddAssign for WorkerStatistics {
    fn add_assign(&mut self, other: Self) {
        self.total_adapters_removed += other.total_adapters_removed;
        self.total_discarded_uncalled += other.total_discarded_uncalled;
        self.total_discarded_barcode += other.total_discarded_barcode;
        self.total_discarded_length += other.total_discarded_length;
        self.total_single_reads_written += other.total_single_reads_written;
        self.total_single_reads_replaced += other.total_single_reads_replaced;
        self.total_single_reads_dropped += other.total_single_reads_dropped;
        self.total_bases_read += other.total_bases_read;
        self.bases_quality_trim += other.bases_quality_trim;
        self.bases_discarded += other.bases_discarded;
        self.bases_adapter_removal += other.bases_adapter_removal;
        self.bases_barcode_removal += other.bases_barcode_removal;
        self.bases_post_trim += other.bases_post_trim;
        self.bases_pre_trim += other.bases_pre_trim;
    }
}
