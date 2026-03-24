use log::trace;
use crate::types::{RecordDto};

pub mod fasta;
pub mod fastq;
mod generic;

pub trait SequenceReader: Iterator<Item=RecordDto> {
    fn compression(&self) -> crate::types::CompressionType;
    fn statistics(&self) -> crate::statistics::ReaderStatistics;

    fn next_batch(&mut self, batch_size: usize) -> Option<Vec<RecordDto>> {
        trace!("Starting to read new batch ...");
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            if let Some(record) = self.next() {
                batch.push(record);
            } else {
                break;
            }
        }
        if batch.is_empty() {
            trace!("Fetched empty batch. No records available.");
            None
        } else {
            trace!("Fetched bach of {} record(s).", batch.len());
            Some(batch)
        }
    }
}

pub struct SequenceReaderBatchIterator<'a> {
    reader: &'a mut dyn SequenceReader,
    batch_size: usize,
}

impl<'a> Iterator for SequenceReaderBatchIterator<'a> {
    type Item = Vec<RecordDto>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next_batch(self.batch_size)
    }
}

pub fn get_batch_iter(reader: &mut dyn SequenceReader, batch_size: usize) -> SequenceReaderBatchIterator<'_> {
    SequenceReaderBatchIterator {
        reader,
        batch_size,
    }
}