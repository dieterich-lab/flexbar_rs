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
        info!("  Total bases read:             {}", self.reader.total_bases_read);
        info!("");

        let total_discards = self.worker.total_discarded_uncalled
            + self.worker.total_discarded_barcode
            + self.worker.total_discarded_length;

        info!("Discarded reads:");
        info!("* Uncalled:                     {}", self.worker.total_discarded_uncalled);
        info!("* Barcode:                      {}", self.worker.total_discarded_barcode);
        info!("* Length:                       {}", self.worker.total_discarded_length);
        info!("Total discards:                 {}", total_discards);

        let remaining = self.writer.total_records_written;
        let remaining_pct = if self.reader.total_records_read > 0 {
            (remaining as f64 / self.reader.total_records_read as f64) * 100.0
        } else {
            0.0
        };

        info!("Remaining reads:                {} ({:.1}%)", remaining, remaining_pct);
        info!("");

        info!("Total adapters removed:         {}", self.worker.total_adapters_removed);

        if !self.writer.barcode_counts.is_empty() {
            info!("");
            info!("Processed barcodes:");
            let mut sorted_barcodes: Vec<_> = self.writer.barcode_counts.iter().collect();
            sorted_barcodes.sort_by_key(|(id, _)| *id);
            for (id, count) in sorted_barcodes {
                info!("* Barcode {}: {}", id, count);
            }
        }

        let bases_written_pct = if self.reader.total_bases_read > 0 {
            (self.writer.total_bases_written as f64 / self.reader.total_bases_read as f64) * 100.0
        } else {
            0.0
        };

        info!("");
        info!("Total bases read:               {}", self.reader.total_bases_read);
        info!("Total bases written:            {} ({:.1}%)", self.writer.total_bases_written, bases_written_pct);
    }
}
