use std::io::{self, BufRead};
use log::{debug, trace};
use crate::types::{CompressionType, FileFormat};
use crate::logging::log_fatal;
use crate::args::Args;
use crate::context::Context;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::readers::backends::ReaderBackend;
use crate::readers::backends::fasta::FastaBackend;
use crate::readers::backends::fastq::FastqBackend;
use crate::readers::backends::bam::BamBackend;
use crate::records::unified::{UnifiedRecord, UnifiedRead, GenericRawRead};
use crate::statistics::ReaderStatistics;
use crate::utils::io::get_raw_reader;

pub struct SequenceReader<Reader: BufRead + Send> {
    backends: [Option<ReaderBackend<Reader>>; 2],
    barcode_backend: Option<ReaderBackend<Reader>>,
    compression: CompressionType,
    last_id: String,
    statistics: ReaderStatistics,
    skip: usize,
    count: Option<usize>,
    trace: bool,
    _phantom: std::marker::PhantomData<(Reader, UnifiedRecord)>,
}

impl<Reader: BufRead + Send> SequenceReader<Reader> {
    pub fn new(
        readers: Vec<Reader>,
        sources: Vec<String>,
        file_format: FileFormat,
        compression: CompressionType,
        skip: usize,
        count: Option<usize>,
        trace: bool,
        barcode_reader: Option<Reader>,
        barcode_source: Option<String>,
    ) -> Self {
        debug!("Creating new sequence reader with source:");
        for (i, source) in sources.iter().enumerate() {
            debug!("- Source {}: {}", i+1, source);
        }
        if let Some(source) = &barcode_source {
            debug!("- Barcode Source: {}", source);
        }
        debug!("- Format: {:?}", file_format);
        debug!("- Compression: {:?}", compression);
        debug!("- Skip {}/Count {:?}", skip, count);
        let mut backends: [Option<ReaderBackend<Reader>>; 2] = [None, None];
        if readers.len() > 2 {
            log_fatal("GenericReader only supports up to two input files.");
        }
        for (i, reader) in readers.into_iter().enumerate() {
            backends[i] = match file_format {
                FileFormat::Fasta => Some(
                    ReaderBackend::Fasta(FastaBackend::new(reader, &sources[i]))
                ),
                FileFormat::Fastq => Some(
                    ReaderBackend::Fastq(FastqBackend::new(reader, &sources[i]))
                ),
                FileFormat::Bam => Some(
                    ReaderBackend::Bam(BamBackend::new(reader, &sources[i]))
                ),
            };
        }

        let barcode_backend = if let (Some(reader), Some(source)) = (barcode_reader, barcode_source) {
            match file_format {
                FileFormat::Fasta => Some(
                    ReaderBackend::Fasta(FastaBackend::new(reader, &source))
                ),
                FileFormat::Fastq => Some(
                    ReaderBackend::Fastq(FastqBackend::new(reader, &source))
                ),
                FileFormat::Bam => Some(
                    ReaderBackend::Bam(BamBackend::new(reader, &source))
                ),
            }
        } else {
            None
        };

        Self {
            backends,
            barcode_backend,
            compression,
            last_id: String::new(),
            statistics: crate::statistics::ReaderStatistics::default(),
            skip,
            count,
            trace,
            _phantom: std::marker::PhantomData,
        }
    }

    fn read_record(&mut self) -> Option<UnifiedRecord> {
        let barcode = if let Some(backend) = &mut self.barcode_backend {
            let b_rec = backend.read_record();
            if b_rec.is_none() {
                return None;
            }
            b_rec
        } else {
            None
        };

        if self.backends[1].is_none() {
            let rec = self.backends[0].as_mut().unwrap().read_record();
            if rec.is_none() {
                return None;
            }
            if let Some(barcode) = barcode {
                Some(UnifiedRecord::single_barcoded(rec.unwrap(), barcode))
            } else {
                Some(UnifiedRecord::single(rec.unwrap()))
            }
        } else {
            let mut records: [Option<UnifiedRead>; 2] = [None, None];
            for i in 0..2 {
                records[i] = self.backends[i].as_mut().unwrap().read_record();
            }
            match (records, barcode) {
                ([None, None], None) => None,
                ([Some(r1), Some(r2)], None) => {
                    Some(UnifiedRecord::paired(r1, r2))
                },
                ([Some(r1), Some(r2)], Some(barcode)) => {
                    Some(UnifiedRecord::paired_barcoded(r1, r2, barcode))
                },
                _ => {
                    log_fatal("Input files have different number of records.");
                }
            }
        }
    }

    pub fn compression(&self) -> CompressionType {
        self.compression
    }

    pub fn statistics(&self) -> crate::statistics::ReaderStatistics {
        self.statistics
    }

    fn try_get_record(&mut self) -> Option<UnifiedRecord> {
        match self.read_record() {
            Some(record) => {
                self.last_id = Self::truncate_id(record.id(0));
                Some(record)
            },
            None => None,
        }
    }

    fn truncate_id(id: &str) -> String {
        if id.len() > 64 {
            id[..64].to_string() + "..."
        } else {
            id.to_string()
        }
    }

    pub fn count_read(&mut self, record: &UnifiedRead, read_index: usize) {
        self.statistics.total_records_read += 1;
        self.statistics.total_bases_read += record.seq().len();
        if self.trace {
            let cooked_id = if self.backends[1].is_none() {
                record.id()
            } else {
                &*format!("{} Read {}", record.id(), read_index + 1)
            };
            trace!("Reader processed record:\nID: {}\nSeq: {}\nQual: {}",
                cooked_id,
                String::from_utf8_lossy(record.seq()),
                String::from_utf8_lossy(record.qual().unwrap_or(&[]))
            );
        }
    }
}

impl SequenceReader<io::BufReader<Box<dyn io::Read + Send>>>
{
    pub fn create_form_args(args: &Args, context: &Context) -> Self {
        let mut paths = Vec::new();
        paths.push(args.reads.as_ref().unwrap_or_else(|| {
            log_fatal("No input reads provided");
        }));
        if args.reads2.is_some() {
            paths.push(args.reads2.as_ref().unwrap());
        }

        let mut barcode_path = None;
        if args.barcode_reads.is_some() {
            barcode_path = args.barcode_reads.as_ref();
        }

        let mut readers = Vec::with_capacity(paths.len());
        let mut sources = Vec::with_capacity(paths.len());
        let mut compression = CompressionType::Uncompressed;

        for path in paths {
            let (reader, source, comp) = get_raw_reader(path, args.read_buffer);
            readers.push(reader);
            sources.push(source);
            compression = comp;
        }

        let mut barcode_reader = None;
        let mut barcode_source = None;

        if let Some(path) = barcode_path {
            let (reader, source, _) = get_raw_reader(path, args.read_buffer);
            barcode_reader = Some(reader);
            barcode_source = Some(source);
        }

        Self::new(readers, sources, context.format, compression, args.skip, args.count, args.trace, barcode_reader, barcode_source)
    }
}

impl<Reader: BufRead + Send> Iterator for SequenceReader<Reader> {
    type Item = UnifiedRecord;

    fn next(&mut self) -> Option<UnifiedRecord> {
        loop {
            return match self.try_get_record() {
                Some(record) => {
                    if self.skip > 0 {
                        self.skip -= 1;
                        if self.skip == 0 {
                            trace!("Done skipped records according to -- skip parameter.");
                        }
                        continue;
                    }
                    if let Some(count) = self.count {
                        if count > 0 {
                            self.count = Some(count - 1);
                        } else {
                            trace!("Stopped reading records due to --count parameter.");
                            return None;
                        }
                    }
                    for i in 0..2 {
                        if let Some(read) = &record.reads[i] {
                            self.count_read(read, i);
                        }
                    }
                    Some(record)
                },
                None => None
            }
        }
    }
}

pub struct SequenceReaderBatchIterator<'a, Reader: BufRead + Send> {
    reader: &'a mut SequenceReader<Reader>,
    batch_size: usize,
}

impl<'a, Reader: BufRead + Send> SequenceReaderBatchIterator<'a, Reader> {
    pub fn new(reader: &'a mut SequenceReader<Reader>, batch_size: usize) -> Self {
        Self { reader, batch_size }
    }
}

impl<'a, Reader: BufRead + Send> Iterator
    for SequenceReaderBatchIterator<'a, Reader>
{
    type Item = Vec<UnifiedRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        trace!("Starting to read new batch ...");
        let mut batch = Vec::with_capacity(self.batch_size);
        for _ in 0..self.batch_size {
            if let Some(record) = self.reader.next() {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_generic_reader_statistics() {
        let data = b"@id1\nACGT\n+\n!!!!\n";
        let reader = BufReader::new(&data[..]);
        let mut gen_reader: SequenceReader<_> =
            SequenceReader::new(
                vec![reader],
                vec!["test".to_string()],
                FileFormat::Fastq,
                CompressionType::Uncompressed,
                0,
                None,
                true,
                None,
                None
            );
        
        let _ = gen_reader.next();
        let stats = gen_reader.statistics();
        assert_eq!(stats.total_records_read, 1);
        assert_eq!(stats.total_bases_read, 4);
    }

    #[test]
    fn test_generic_reader_statistics_paired() {
        let data1 = b"@id1/1\nACGT\n+\n!!!!\n";
        let data2 = b"@id1/2\nTGCA\n+\n????\n";
        let reader1 = BufReader::new(&data1[..]);
        let reader2 = BufReader::new(&data2[..]);

        let mut gen_reader: SequenceReader<_> =
            SequenceReader::new(
                vec![reader1, reader2],
                vec!["test1".to_string(), "test2".to_string()],
                FileFormat::Fastq,
                CompressionType::Uncompressed,
                0,
                None,
                true,
                None,
                None
            );

        let _ = gen_reader.next();
        let stats = gen_reader.statistics();
        assert_eq!(stats.total_records_read, 2, "Should count both reads in a pair");
        assert_eq!(stats.total_bases_read, 8, "Should count bases from both reads in a pair");
    }
}
