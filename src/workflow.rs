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