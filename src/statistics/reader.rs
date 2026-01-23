use std::ops::Add;

#[derive(Debug, Default, Clone, Copy)]
pub struct ReaderStatistics {
    pub total_records_read: usize,
    pub total_bases_read: usize,
}

impl Add for ReaderStatistics {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            total_records_read: self.total_records_read + other.total_records_read,
            total_bases_read: self.total_bases_read + other.total_bases_read,
        }
    }
}
