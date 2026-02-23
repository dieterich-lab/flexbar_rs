use rayon::prelude::*;
use crate::args::Args;
use crate::context::Context;
use crate::readers::fasta::SingleReaderFasta;
use crate::readers::fastq::SingleReaderFastq;
use crate::readers::{SequenceReader, get_batch_iter};
use crate::writers::fasta::SingleWriterFasta;
use crate::writers::fastq::SingleWriterFastq;
use crate::types::{FileFormat, CompressionType, SingeRecordDto};
use crate::workers::single::SingleWorker;
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};
use crate::logging::log_fatal;
use log::{debug};
use crate::context;

pub fn do_it(args: Args) {
    let context = context::Context::new(&args);
    let reads_path = args.reads.as_ref().unwrap_or_else(|| {
        log_fatal("No input reads provided");
    });
    let output_prefix = args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| args.target.clone());

    let stats = process_single(&args, &context, reads_path, &output_prefix);
    stats.print();
}

fn process_single(args: &Args, context: &Context, reads_path: &Path, output_prefix: &str) -> Statistics {
    let mut reader = get_single_reader(reads_path, context.format, args);
    let compression = get_compression(&reader, args);
    let worker = SingleWorker::new(args, context);
    let mut writer = get_single_writer(args, output_prefix, compression);
    let worker_stats = process_batch(&mut reader, &worker, &mut writer, args.batch_size, args.count);

    Statistics {
        reader: reader.statistics(),
        worker: worker_stats,
        writer: writer.statistics(),
    }
}


fn process_batch(
    reader: &mut Box<dyn SequenceReader + Send>,
    worker_shared: &SingleWorker,
    writer: &mut Box<dyn Writer>,
    batch_size: usize,
    count: Option<usize>,
) -> WorkerStatistics {
    let mut worker_stats = WorkerStatistics::default();
    let mut batch_num = 0;
    let mut processed_total = 0;

    for mut batch_in in get_batch_iter(reader.as_mut(), batch_size) {
        batch_num += 1;
        debug!("Processing batch {}", batch_num);

        if let Some(c) = count {
            if processed_total + batch_in.len() > c {
                let remaining = c - processed_total;
                batch_in.truncate(remaining);
            }
        }

        let batch_len = batch_in.len();
        
        let batch_out: Vec<(Option<SingeRecordDto>, Option<String>, WorkerStatistics)> = batch_in
            .into_par_iter()
            .map(|rec| {
                let mut worker = worker_shared.clone();
                worker.process(&rec)
            })
            .collect();

        for (rec_opt, tag_opt, stats) in batch_out {
            worker_stats += stats;
            if let Some(rec) = rec_opt {
                writer.write(&rec, tag_opt.as_deref()).expect("Failed to write record");
            }
        }

        processed_total += batch_len;
        if let Some(c) = count {
            if processed_total >= c {
                break;
            }
        }
    }

    worker_stats
}

fn get_single_reader(reads_path: &Path, format: FileFormat, args: &crate::args::Args) -> Box<dyn SequenceReader + Send> {
    match format {
        FileFormat::Fasta => Box::new(SingleReaderFasta::create_form_args(reads_path, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTA reader: {}", e)))),
        FileFormat::Fastq => Box::new(SingleReaderFastq::create_form_args(reads_path, args)
            .unwrap_or_else(|e| log_fatal(&format!("Failed to create FASTQ reader: {}", e)))),
    }
}

fn get_compression(reader: &Box<dyn SequenceReader + Send>, args: &crate::args::Args) -> CompressionType {
    match args.output_compression {
        Some(CompressionType::Keep) => reader.compression(),
        Some(c) => c,
        None => CompressionType::Uncompressed,
    }
}

fn get_single_writer(args: &crate::args::Args, output_prefix: &str, compression: CompressionType) -> Box<dyn Writer> {
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
