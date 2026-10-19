use std::cell::RefCell;
use std::io::BufRead;
use rayon::prelude::*;
use crate::args::Args;
use crate::writers::sequence::SequenceWriter;
use crate::types::{CompressionType, FileFormat};
use crate::statistics::{Statistics, WorkerStatistics};
use log::debug;
use crate::context;
use crate::records::header::SequenceHeader;
use crate::records::unified::AnnotatedRecord;
use crate::worker::{Worker, WorkerState};
use crossbeam_channel::{bounded, Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;
use crate::readers::sequence::{SequenceReader, SequenceReaderBatchIterator};
use crate::records::unified::UnifiedRecord;

pub fn do_it(args: Args) {
    let context = context::Context::new(&args);

    let output_prefix = args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| args.target.clone());
    let reader = SequenceReader::create_form_args(&args, &context);
    let header = reader.get_header();
    let compression = get_compression(&reader, &args);
    let worker = Worker::new(&args, &context);
    let writer = get_writer(&args, &output_prefix, compression, header, &reader);

    let (batch_in_tx, batch_in_rx) = bounded::<Vec<UnifiedRecord>>(16);
    let (batch_out_tx, batch_out_rx) = bounded::<ProcessedBatch>(16);

    let reader_handle = spawn_reader_thread(reader, batch_in_tx, args.batch_size);
    let writer_handle = spawn_writer_thread(writer, batch_out_rx);

    process_batch_orchestrator(batch_in_rx, batch_out_tx, &worker);

    let reader = reader_handle.join().expect("Reader thread panicked");
    let (writer, worker_stats) = writer_handle.join().expect("Writer thread panicked");

    writer.write_length_distribution().expect("Failed to write length distribution");

    let stats = Statistics {
        reader: reader.statistics(),
        worker: worker_stats,
        writer: writer.statistics(),
    };
    stats.print();
}

fn spawn_reader_thread<Reader: BufRead + Send + 'static>(
    mut reader: SequenceReader<Reader>,
    batch_in_tx: Sender<Vec<UnifiedRecord>>,
    batch_size: usize,
) -> JoinHandle<SequenceReader<Reader>> {
    thread::spawn(move || {
        for batch in SequenceReaderBatchIterator::new(&mut reader, batch_size) {
            if batch_in_tx.send(batch).is_err() {
                break;
            }
        }
        reader
    })
}

fn spawn_writer_thread(
    mut writer: SequenceWriter<'static>,
    batch_out_rx: Receiver<ProcessedBatch>,
) -> JoinHandle<(SequenceWriter<'static>, WorkerStatistics)> {
    thread::spawn(move || {
        let mut worker_stats = WorkerStatistics::default();
        while let Ok(batch_out) = batch_out_rx.recv() {
            worker_stats += batch_out.stats;
            for rec in batch_out.records {
                writer.write(&rec)
                    .expect("Failed to write record");
            }
        }
        (writer, worker_stats)
    })
}

fn process_batch_orchestrator(
    batch_in_rx: Receiver<Vec<UnifiedRecord>>,
    batch_out_tx: Sender<ProcessedBatch>,
    worker: &Worker,
) {
    let mut batch_num = 0;
    while let Ok(batch_in) = batch_in_rx.recv() {
        batch_num += 1;
        debug!("Processing batch {}", batch_num);

        let batch_out: Vec<(AnnotatedRecord, WorkerStatistics)> = batch_in
            .into_par_iter()
            .map(|rec| process_record(&rec, worker))
            .collect();

        let mut stats = WorkerStatistics::default();
        let mut records = Vec::with_capacity(batch_out.len());

        for (r, s) in batch_out {
            stats += s;
            records.push(r);
        }

        if batch_out_tx.send(ProcessedBatch { records, stats }).is_err() {
            break;
        }
    }
}

struct ProcessedBatch {
    records: Vec<AnnotatedRecord>,
    stats: WorkerStatistics,
}

fn get_compression<Reader: BufRead + Send>(reader: &SequenceReader<Reader>, args: &Args)
    -> CompressionType
{
    match args.output_compression {
        Some(CompressionType::Keep) => reader.compression(),
        Some(c) => c,
        None => CompressionType::Uncompressed,
    }
}

use crate::utils::io::get_writer_factory;

fn get_writer<Reader: BufRead + Send>(args: &Args, output_prefix: &str, compression: CompressionType, header: SequenceHeader, reader: &SequenceReader<Reader>) -> SequenceWriter<'static> {
    let p = match args.single_reads {
        true => "-".to_string(),
        false => output_prefix.to_string(),
    };
    let factory = get_writer_factory(p, compression, args.write_buffer);
    let format = match args.output_format {
        Some(x) => x,
        None => if reader.format() == FileFormat::Fasta {FileFormat::Fastq} else {reader.format()},
    };
    SequenceWriter::new(format, factory, args.length_dist, header, args.barcode_no_split, args.barcode_mark)
}

thread_local! {
    static WORKER_STATE: RefCell<Option<WorkerState>> = RefCell::new(None);
}

fn process_record(
    record: &UnifiedRecord, worker: &Worker
) -> (AnnotatedRecord, WorkerStatistics) {
    WORKER_STATE.with(|state| {
        let mut state = state.borrow_mut();
        if state.is_none() {
            *state = Some(worker.create_state());
        }
        worker.process_record(record, state.as_mut().unwrap())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::records::generic_read::GenericRawRead;
    use crate::utils::io::WriterFactory;

    #[test]
    fn test_bam_quality_roundtrip() {
        let fastq_data: Vec<u8> = b"@read1\nACGT\n+\nIIII\n".to_vec();
        let cursor = Cursor::new(fastq_data);
        let mut reader = SequenceReader::new(
            vec![cursor],
            vec!["test.fastq".to_string()],
            FileFormat::Fastq,
            CompressionType::Uncompressed,
            0,
            None,
            false,
            None,
            None,
        );

        let record = reader.next().expect("Failed to read FASTQ record");
        let read = record.reads[0].as_ref().expect("Missing read 0");
        assert_eq!(read.qual(), Some(&b"IIII"[..]));

        let temp_file = tempfile::NamedTempFile::new().expect("Failed to create temp file");
        let temp_path = temp_file.path().to_str().expect("Invalid path").to_string();
        let temp_path_clone = temp_path.clone();

        let factory: WriterFactory = Box::new(move |_tag, _read_index, _format, _output_type| {
            let cursor = Box::new(Cursor::new(Vec::new()));
            Ok((cursor, temp_path_clone.clone()))
        });

        {
            let mut writer = SequenceWriter::new(
                FileFormat::Bam,
                factory,
                false,
                SequenceHeader::Default,
                false,
                false,
            );

            let annotated = AnnotatedRecord {
                record: record.clone(),
                barcode_id: None,
                removal_tags: [vec![], vec![]],
                umi_tags: vec![],
                paired_read_dropped: false,
            };

            writer.write(&annotated).expect("Failed to write BAM record");
        }

        let bam_bytes = std::fs::read(&temp_path).expect("Failed to read temp BAM file");
        let bam_cursor = Cursor::new(bam_bytes);
        let mut bam_reader = SequenceReader::new(
            vec![bam_cursor],
            vec![temp_path.clone()],
            FileFormat::Bam,
            CompressionType::Uncompressed,
            0,
            None,
            false,
            None,
            None,
        );

        let roundtrip_record = bam_reader.next().expect("Failed to read BAM record");
        let roundtrip_read = roundtrip_record.reads[0].as_ref().expect("Missing read 0");
        assert_eq!(roundtrip_read.qual(), read.qual());
    }

    #[test]
    fn test_samtools_bam_sequence_and_quality() {
        use std::io::Write;

        let bam_data: &[u8] = &[
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
            0xa8, 0x00, 0x45, 0x8d, 0xcb, 0x0a, 0x82, 0x40, 0x14, 0x40, 0x6f, 0xcb, 0xeb, 0x57, 0xdc, 0x1f,
            0x70, 0xc0, 0x82, 0x16, 0xb3, 0xf2, 0x05, 0x16, 0xa4, 0x23, 0x05, 0xed, 0x07, 0xbc, 0x92, 0xa0,
            0x8e, 0xce, 0x8c, 0x3d, 0xbe, 0xa6, 0x5f, 0x0d, 0x84, 0x6a, 0x77, 0x38, 0x1c, 0x38, 0x69, 0x52,
            0x6e, 0xde, 0x00, 0x10, 0x1f, 0x72, 0xbc, 0x56, 0x32, 0x12, 0x7b, 0xbc, 0x28, 0xb9, 0x8c, 0xce,
            0x58, 0xcf, 0x0d, 0x16, 0x4a, 0xce, 0x0b, 0xdb, 0x57, 0x10, 0x67, 0x0a, 0xcf, 0x7c, 0x67, 0xeb,
            0x98, 0x1e, 0x9d, 0xbf, 0x49, 0x72, 0x7a, 0xf0, 0xc6, 0xf4, 0x8e, 0x5a, 0xed, 0xfc, 0x4c, 0xe1,
            0x48, 0xa1, 0xa1, 0x49, 0x77, 0x96, 0x1b, 0xb1, 0xaa, 0x20, 0xae, 0x0b, 0x3c, 0xe6, 0xf2, 0x1b,
            0x62, 0x5d, 0xfd, 0x79, 0x7d, 0x6d, 0x77, 0x22, 0xc2, 0xec, 0xf4, 0xb3, 0xd4, 0x0d, 0x93, 0xb1,
            0x9e, 0x42, 0x45, 0x69, 0x52, 0x12, 0x3f, 0xf5, 0x30, 0xf5, 0x2c, 0xda, 0x39, 0x00, 0x00, 0xf8,
            0x00, 0xbc, 0x3a, 0xec, 0xc2, 0xab, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00,
            0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x3c, 0x00, 0xd3, 0x61, 0x60, 0x60, 0xf8,
            0x0f, 0x05, 0x6c, 0x0c, 0x1e, 0x42, 0x0c, 0x0c, 0x2c, 0x0c, 0x2c, 0x48, 0x62, 0x0c, 0x0c, 0x0c,
            0x0c, 0x45, 0xa9, 0x89, 0x29, 0x86, 0x0c, 0x42, 0x1e, 0x1a, 0x1a, 0x1a, 0x1a, 0x00, 0x94, 0x43,
            0x07, 0x4e, 0x30, 0x00, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
            0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00,
        ];

        let mut temp_file = tempfile::NamedTempFile::new().expect("Failed to create temp file");
        temp_file.write_all(bam_data).expect("Failed to write BAM data");
        let temp_path = temp_file.path().to_str().expect("Invalid path").to_string();

        let cursor = Cursor::new(Vec::new());
        let mut reader = SequenceReader::new(
            vec![cursor],
            vec![temp_path],
            FileFormat::Bam,
            CompressionType::Uncompressed,
            0,
            None,
            false,
            None,
            None,
        );

        let record = reader.next().expect("Failed to read BAM record");
        let read = record.reads[0].as_ref().expect("Missing read 0");
        assert_eq!(read.seq(), b"ACGT");
        assert_eq!(read.qual(), Some(&b"IIII"[..]));
    }
}