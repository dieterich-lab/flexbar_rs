use std::collections::HashMap;
use std::ops::{Add, AddAssign};

#[derive(Debug, Default, Clone)]
pub struct WorkerStatistics {
    pub total_adapters_removed: u64,
    pub total_discarded_uncalled: usize,
    pub barcode_counts: HashMap<String, u64>,
}

impl WorkerStatistics {
    pub fn is_discarded(&self) -> bool {
        self.total_discarded_uncalled > 0
    }
}

impl Add for WorkerStatistics {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self.total_adapters_removed += other.total_adapters_removed;
        self.total_discarded_uncalled += other.total_discarded_uncalled;
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
        for (id, count) in other.barcode_counts {
            *self.barcode_counts.entry(id).or_insert(0) += count;
        }
    }
}
