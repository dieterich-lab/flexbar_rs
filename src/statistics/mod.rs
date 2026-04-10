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

        info!("Discarded reads:                {}", total_discards);
        info!("* Uncalled:                     {}", self.worker.total_discarded_uncalled);
        info!("* Barcode:                      {}", self.worker.total_discarded_barcode);
        info!("* Length:                       {}", self.worker.total_discarded_length);
        if self.worker.total_single_reads_written > 0 {
            info!("  ** Single reads written:         {}", self.worker.total_single_reads_written);
        }
        if self.worker.total_single_reads_replaced > 0 {
            info!("  ** Reads replaced by N:          {}", self.worker.total_single_reads_replaced);
        }
        if self.worker.total_single_reads_dropped > 0 {
            info!("  ** Single reads dropped:         {}", self.worker.total_single_reads_dropped);
        }
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

        info!("");
        info!("Total bases read:               {}", self.worker.total_bases_read);
        info!("- bases in discarded reads:     {}", self.worker.bases_discarded);
        info!("- removed by pre trim:          {}", self.worker.bases_pre_trim);
        info!("- removed by quality trim:      {}", self.worker.bases_quality_trim);
        info!("- removed by barcode removal:   {}", self.worker.bases_barcode_removal);
        info!("- removed by adapter removal:   {}", self.worker.bases_adapter_removal);
        info!("- removed by post trim:         {}", self.worker.bases_post_trim);

        let remaining_bases = self.writer.total_bases_written;
        let bases_written_pct = if self.worker.total_bases_read > 0 {
            (remaining_bases as f64 / self.worker.total_bases_read as f64) * 100.0
        } else {
            0.0
        };

        info!("- remaining bases written:      {} ({:.1}%)", remaining_bases, bases_written_pct);
    }
}
