use rayon::prelude::*;
use crate::config::Config;
use crate::context::Context;
use crate::readers::fasta::SingleReaderFasta;
use crate::readers::fastq::SingleReaderFastq;
use crate::readers::{SequenceReader, get_batch_iter};
use crate::writers::fasta::SingleWriterFasta;
use crate::writers::fastq::SingleWriterFastq;
use crate::types::{FileFormat, CompressionType, RecordDto};
use crate::workers::single::SingleWorker;
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};
use log::info;

pub fn do_it(config: Config) {
    let mut context = Context::new();
    context.load(&config);

    if !context.adapters.is_empty() {
        for (i, adapter) in context.adapters.iter().enumerate() {
            info!("Adapter 1-{}: {}", i + 1, adapter);
        }
    }
    if !context.adapters2.is_empty() {
        for (i, adapter) in context.adapters2.iter().enumerate() {
            info!("Adapter 2-{}: {}", i + 1, adapter);
        }
    }
    if !context.barcodes.is_empty() {
        for (i, (id, barcode)) in context.barcodes.iter().enumerate() {
            info!("Barcode 1-{}: {} {}", i + 1, id, barcode);
        }
    }
    if !context.barcodes2.is_empty() {
        for (i, (id, barcode)) in context.barcodes2.iter().enumerate() {
            info!("Barcode 2-{}: {} {}", i + 1, id, barcode);
        }
    }

    let reads_path = config.args.reads.as_ref().expect("No input reads provided");
    let output_prefix = config.args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| config.args.target.clone());

    let stats = process_single(&config, &context, reads_path, &output_prefix);
    stats.print();
}

fn process_single(config: &Config, context: &Context, reads_path: &Path, output_prefix: &str) -> Statistics {
    let mut reader = get_single_reader(reads_path, config.format);
    let compression = get_compression(&reader, config);
    let worker = SingleWorker::new(config, context);
    let mut writer = get_single_writer(config, output_prefix, compression);

    let worker_stats = process_batch(&mut reader, &worker, &mut writer, config.batch_size);

    Statistics {
        reader: reader.statistics(),
        worker: worker_stats,
        writer: writer.statistics(),
    }
}

fn process_batch(
    reader: &mut Box<dyn SequenceReader + Send>,
    worker: &SingleWorker,
    writer: &mut Box<dyn Writer>,
    batch_size: usize,
) -> WorkerStatistics {
    let mut worker_stats = WorkerStatistics::default();

    for batch_in in get_batch_iter(reader.as_mut(), batch_size) {
        let batch_out: Vec<(Option<RecordDto>, WorkerStatistics, Option<String>)> = batch_in
            .into_par_iter()
            .map(|rec| worker.process(&rec))
            .collect();

        for (rec_opt, stats, tag) in batch_out {
            worker_stats += stats;
            if let Some(rec) = rec_opt {
                writer.write(&rec, tag.as_deref()).expect("Failed to write record");
            }
        }
    }

    worker_stats
}

fn get_single_reader(reads_path: &Path, format: FileFormat) -> Box<dyn SequenceReader + Send> {
    match format {
        FileFormat::Fasta => Box::new(SingleReaderFasta::from_path(reads_path)
            .expect("Failed to create FASTA reader")),
        FileFormat::Fastq => Box::new(SingleReaderFastq::from_path(reads_path)
            .expect("Failed to create FASTQ reader")),
    }
}

fn get_compression(reader: &Box<dyn SequenceReader + Send>, config: &Config) -> CompressionType {
    match config.args.output_compression {
        Some(CompressionType::Keep) => reader.compression(),
        Some(c) => c,
        None => CompressionType::Uncompressed,
    }
}

fn get_single_writer(config: &Config, output_prefix: &str, compression: CompressionType) -> Box<dyn Writer> {
    let p = match config.args.single_reads {
        true => "-",
        false => output_prefix,
    };
    if config.args.fasta_output {
        Box::new(SingleWriterFasta::new(p, compression).expect("Failed to create FASTA writer"))
    } else {
        Box::new(SingleWriterFastq::new(p, compression).expect("Failed to create FASTQ writer"))
    }
}
