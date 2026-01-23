use std::ops::Add;
use std::collections::HashMap;

#[derive(Debug, Default, Clone)]
pub struct WriterStatistics {
    pub total_records_written: usize,
    pub total_bases_written: usize,
    pub barcode_counts: HashMap<String, usize>,
}

impl Add for WriterStatistics {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self.total_records_written += other.total_records_written;
        self.total_bases_written += other.total_bases_written;
        for (id, count) in other.barcode_counts {
            *self.barcode_counts.entry(id).or_insert(0) += count;
        }
        self
    }
}
