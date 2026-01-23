use std::ops::Add;

mod reader;
mod worker;
mod writer;

pub use reader::ReaderStatistics;
pub use worker::WorkerStatistics;
pub use writer::WriterStatistics;

#[derive(Debug, Default, Clone)]
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

use log::info;

impl Statistics {
    pub fn print(&self) {
        info!("Flexbar statistics:");
        info!("  Total records read:           {}", self.reader.total_records_read);
        info!("  Total discarded (uncalled):   {}", self.worker.total_discarded_uncalled);
        info!("  Total discarded (barcode):    {}", self.worker.total_discarded_barcode);
        info!("  Total adapters removed:       {}", self.worker.total_adapters_removed);
        if !self.worker.barcode_counts.is_empty() {
            let mut sorted_barcodes: Vec<_> = self.worker.barcode_counts.iter().collect();
            sorted_barcodes.sort_by_key(|(id, _)| *id);
            for (id, count) in sorted_barcodes {
                info!("  Barcode {}: {}", id, count);
            }
        }
        info!("  Total records written:        {}", self.writer.total_records_written);
    }
}
