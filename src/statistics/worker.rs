use std::collections::HashMap;
use std::ops::{Add, AddAssign};

#[derive(Debug, Default, Clone)]
pub struct WorkerStatistics {
    pub total_adapters_removed: usize,
    pub total_discarded_uncalled: usize,
    pub total_discarded_barcode: usize,
    pub total_discarded_length: usize,
    pub barcode_counts: HashMap<String, usize>,
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
        for (id, count) in other.barcode_counts {
            *self.barcode_counts.entry(id).or_insert(0) += count;
        }
        self
    }
}

impl AddAssign for WorkerStatistics {
    fn add_assign(&mut self, other: Self) {
        self.total_adapters_removed += other.total_adapters_removed;
        self.total_discarded_uncalled += other.total_discarded_uncalled;
        self.total_discarded_barcode += other.total_discarded_barcode;
        self.total_discarded_length += other.total_discarded_length;
        for (id, count) in other.barcode_counts {
            *self.barcode_counts.entry(id).or_insert(0) += count;
        }
    }
}
