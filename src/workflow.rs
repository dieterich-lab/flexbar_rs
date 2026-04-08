use std::cell::{RefCell};
use rayon::prelude::*;
use crate::args::Args;
use crate::readers::fasta::{SingleReaderFasta, PairedReaderFasta};
use crate::readers::fastq::{SingleReaderFastq, PairedReaderFastq};
use crate::readers::{SequenceReader, get_batch_iter};
use crate::writers::fasta::{SingleWriterFasta};
use crate::writers::fastq::{SingleWriterFastq};
use crate::types::{FileFormat, CompressionType, RecordDto};
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};
use crate::logging::log_fatal;
use log::{debug};
use crate::context;
use crate::steps::generic::FinalRecord;
use crate::worker::{Worker, WorkerState};
use crossbeam_channel::{bounded, Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;

pub fn do_it(args: Args) {
    let context = context::Context::new(&args);
    let reads_path = args.reads.as_ref().unwrap_or_else(|| {
        log_fatal("No input reads provided");
    });
    let output_prefix = args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| args.target.clone());
    let reader = get_reader(reads_path, args.reads2.as_ref(), context.format, &args);
    let compression = get_compression(&reader, &args);
    let worker = Worker::new(&args, &context);
    let writer = get_writer(&args, &output_prefix, compression);

    let (batch_in_tx, batch_in_rx) = bounded::<Vec<RecordDto>>(16);
    let (batch_out_tx, batch_out_rx) = bounded::<ProcessedBatch>(16);

    let reader_handle = spawn_reader_thread(reader, batch_in_tx, args.batch_size);
    let writer_handle = spawn_writer_thread(writer, batch_out_rx);

    process_batch_orchestrator(batch_in_rx, batch_out_tx, &worker);

    let reader = reader_handle.join().expect("Reader thread panicked");
    let (writer, worker_stats) = writer_handle.join().expect("Writer thread panicked");

    let stats = Statistics {
        reader: reader.statistics(),
        worker: worker_stats,
        writer: writer.statistics(),
    };
    stats.print();
}

fn spawn_reader_thread(
    mut reader: Box<dyn SequenceReader + Send>,
    batch_in_tx: Sender<Vec<RecordDto>>,
    batch_size: usize,
) -> JoinHandle<Box<dyn SequenceReader + Send>> {
    thread::spawn(move || {
        for batch in get_batch_iter(reader.as_mut(), batch_size) {
            if batch_in_tx.send(batch).is_err() {
                break;
            }
        }
        reader
    })
}

fn spawn_writer_thread(
    mut writer: Box<dyn Writer>,
    batch_out_rx: Receiver<ProcessedBatch>,
) -> JoinHandle<(Box<dyn Writer>, WorkerStatistics)> {
    thread::spawn(move || {
        let mut worker_stats = WorkerStatistics::default();
        while let Ok(batch_out) = batch_out_rx.recv() {
            worker_stats += batch_out.stats;
            for rec in batch_out.records {
                writer.write(&rec.0, rec.1.as_deref())
                    .expect("Failed to write record");
            }
        }
        (writer, worker_stats)
    })
}

fn process_batch_orchestrator(
    batch_in_rx: Receiver<Vec<RecordDto>>,
    batch_out_tx: Sender<ProcessedBatch>,
    worker: &Worker,
) {
    let mut batch_num = 0;
    while let Ok(batch_in) = batch_in_rx.recv() {
        batch_num += 1;
        debug!("Processing batch {}", batch_num);

        let batch_out: Vec<FinalRecord<RecordDto>> = batch_in
            .into_par_iter()
            .map(|rec| process_record(&rec, worker))
            .collect();

        let mut stats = WorkerStatistics::default();
        let mut records = Vec::with_capacity(batch_out.len());

        for res in batch_out {
            stats += res.stats;
            if let Some(r) = res.record {
                records.push((r, res.output_tag));
            }
        }

        if batch_out_tx.send(ProcessedBatch { records, stats }).is_err() {
            break;
        }
    }
}

struct ProcessedBatch {
    records: Vec<(RecordDto, Option<String>)>,
    stats: WorkerStatistics,
}

fn get_reader(
    reads_path: &Path, reads_path2: Option<&std::path::PathBuf>, format: FileFormat, args: &Args
) -> Box<dyn SequenceReader + Send> {
    let mut paths = vec![reads_path];
    if let Some(p) = reads_path2 {
        paths.push(p);
    }
    match (format, reads_path2.is_some()) {
        (FileFormat::Fasta, false) => Box::new(SingleReaderFasta::create_form_args(paths, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTA reader: {}", e)))),
        (FileFormat::Fasta, true) => Box::new(PairedReaderFasta::create_form_args(paths, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create paired FASTA reader: {}", e)))),
        (FileFormat::Fastq, false) => Box::new(SingleReaderFastq::create_form_args(paths, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTQ reader: {}", e)))),
        (FileFormat::Fastq, true) => Box::new(PairedReaderFastq::create_form_args(paths, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create paired FASTQ reader: {}", e)))),
    }
}

fn get_compression(reader: &Box<dyn SequenceReader + Send>, args: &Args) -> CompressionType {
    match args.output_compression {
        Some(CompressionType::Keep) => reader.compression(),
        Some(c) => c,
        None => CompressionType::Uncompressed,
    }
}

use crate::utils::io::{get_writer_factory};

fn get_writer(args: &Args, output_prefix: &str, compression: CompressionType) -> Box<dyn Writer> {
    let p = match args.single_reads {
        true => "-".to_string(),
        false => output_prefix.to_string(),
    };
    let factory = get_writer_factory(p, compression, args.write_buffer);
    if args.fasta_output {
        Box::new(SingleWriterFasta::new(factory))
    } else {
        Box::new(SingleWriterFastq::new(factory))
    }
}

thread_local! {
    static WORKER_STATE: RefCell<Option<WorkerState>> = RefCell::new(None);
}

fn process_record(
    record: &RecordDto, worker: &Worker
) -> FinalRecord<RecordDto> {
    WORKER_STATE.with(|state| {
        let mut state = state.borrow_mut();
        if state.is_none() {
            *state = Some(worker.create_state());
        }
        worker.process_record(record, state.as_mut().unwrap())
    })
}