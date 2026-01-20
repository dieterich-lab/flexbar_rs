use std::ops::{Add, AddAssign};

#[derive(Debug, Default, Clone, Copy)]
pub struct WorkerStatistics {
    pub total_adapters_removed: u64,
    pub total_discarded_uncalled: usize,
}

impl WorkerStatistics {
    pub fn is_discarded(&self) -> bool {
        self.total_discarded_uncalled > 0
    }
}

impl Add for WorkerStatistics {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            total_adapters_removed: self.total_adapters_removed + other.total_adapters_removed,
            total_discarded_uncalled: self.total_discarded_uncalled + other.total_discarded_uncalled,
        }
    }
}

impl AddAssign for WorkerStatistics {
    fn add_assign(&mut self, other: Self) {
        self.total_adapters_removed += other.total_adapters_removed;
        self.total_discarded_uncalled += other.total_discarded_uncalled;
    }
}
