use std::ops::Add;

mod reader;
mod worker;
mod writer;

pub use reader::ReaderStatistics;
pub use worker::WorkerStatistics;
pub use writer::WriterStatistics;

#[derive(Debug, Default, Clone, Copy)]
pub struct Statistics {
    pub reader: ReaderStatistics,
    pub worker: WorkerStatistics,
    pub writer: WriterStatistics,
}

impl Add for Statistics {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            reader: self.reader + other.reader,
            worker: self.worker + other.worker,
            writer: self.writer + other.writer,
        }
    }
}

impl Statistics {
    pub fn print(&self) {
        eprintln!("Flexbar statistics:");
        eprintln!("  Total records read:     {}", self.reader.total_records_read);
        eprintln!("  Total adapters removed: {}", self.worker.total_adapters_removed);
        eprintln!("  Total records written:  {}", self.writer.total_records_written);
    }
}
