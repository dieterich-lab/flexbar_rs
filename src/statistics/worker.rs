use std::ops::Add;

#[derive(Debug, Default, Clone, Copy)]
pub struct WorkerStatistics {
    pub total_adapters_removed: u64,
}

impl Add for WorkerStatistics {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            total_adapters_removed: self.total_adapters_removed + other.total_adapters_removed,
        }
    }
}
