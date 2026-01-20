use std::ops::Add;

#[derive(Debug, Default, Clone, Copy)]
pub struct WriterStatistics {
    pub total_records_written: u64,
}

impl Add for WriterStatistics {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            total_records_written: self.total_records_written + other.total_records_written,
        }
    }
}
