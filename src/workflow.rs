use std::cell::{RefCell};
use rayon::prelude::*;
use crate::args::Args;
use crate::readers::fasta::SingleReaderFasta;
use crate::readers::fastq::SingleReaderFastq;
use crate::readers::{SequenceReader, get_batch_iter};
use crate::writers::fasta::SingleWriterFasta;
use crate::writers::fastq::SingleWriterFastq;
use crate::types::{FileFormat, CompressionType, RecordDto};
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};
use crate::logging::log_fatal;
use log::{debug};
use crate::context;
use crate::steps::generic::FinalRecord;
use crate::workers::{Worker, WorkerState};

pub fn do_it(args: Args) {
    let context = context::Context::new(&args);
    let reads_path = args.reads.as_ref().unwrap_or_else(|| {
        log_fatal("No input reads provided");
    });
    let output_prefix = args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| args.target.clone());
    let mut reader = get_single_reader(reads_path, context.format, &args);
    let compression = get_compression(&reader, &args);
    let worker = Worker::new(&args, &context);
    let mut writer = get_single_writer(&args, &output_prefix, compression);

    let worker_stats = process_batch(
        &mut reader, &worker, &mut writer, args.batch_size
    );

    let stats = Statistics {
        reader: reader.statistics(),
        worker: worker_stats,
        writer: writer.statistics(),
    };
    stats.print();
}

fn process_batch(
    reader: &mut Box<dyn SequenceReader + Send>,
    worker: &Worker,
    writer: &mut Box<dyn Writer>,
    batch_size: usize,
) -> WorkerStatistics {
    let mut worker_stats = WorkerStatistics::default();
    let mut batch_num = 0;

    for batch_in in get_batch_iter(reader.as_mut(), batch_size) {
        batch_num += 1;
        debug!("Processing batch {}", batch_num);
        let batch_out: Vec<FinalRecord<RecordDto>> = batch_in
            .into_par_iter()
            .map(|rec| process_record(&rec, worker) )
            .collect();

        for rec in batch_out {
            worker_stats += rec.stats;
            if let Some(r) = rec.record {
                writer.write(&r, rec.output_tag.as_deref())
                    .expect("Failed to write record");
            }
        }
    }

    worker_stats
}

fn get_single_reader(reads_path: &Path, format: FileFormat, args: &Args) -> Box<dyn SequenceReader + Send> {
    match format {
        FileFormat::Fasta => Box::new(SingleReaderFasta::create_form_args(reads_path, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTA reader: {}", e)))),
        FileFormat::Fastq => Box::new(SingleReaderFastq::create_form_args(reads_path, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTQ reader: {}", e)))),
    }
}

fn get_compression(reader: &Box<dyn SequenceReader + Send>, args: &Args) -> CompressionType {
    match args.output_compression {
        Some(CompressionType::Keep) => reader.compression(),
        Some(c) => c,
        None => CompressionType::Uncompressed,
    }
}

fn get_single_writer(args: &Args, output_prefix: &str, compression: CompressionType) -> Box<dyn Writer> {
    let p = match args.single_reads {
        true => "-",
        false => output_prefix,
    };
    if args.fasta_output {
        Box::new(SingleWriterFasta::new(p, compression)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTA writer: {}", e))))
    } else {
        Box::new(SingleWriterFastq::new(p, compression)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTQ writer: {}", e))))
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