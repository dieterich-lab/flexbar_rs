use rayon::prelude::*;
use crate::config::Config;
use crate::context::Context;
use crate::readers::fasta::SingleReaderFasta;
use crate::readers::fastq::SingleReaderFastq;
use crate::readers::SequenceReader;
use crate::writers::fasta::SingleWriterFasta;
use crate::writers::fastq::SingleWriterFastq;
use crate::types::{FileFormat, CompressionType, RecordDto};
use crate::workers::single::SingleWorker;
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};

pub fn do_it(config: Config) {
    let mut context = Context::new();
    context.load(&config);

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

    let mut last_batch: bool = false;

    while !last_batch {
        let mut batch_in: Vec<RecordDto> = Vec::with_capacity(batch_size);
        let mut actual_batch_size = batch_size;
        for i in 0..batch_size {
            match reader.next() {
                Some(rec) => { batch_in.push(rec) },
                None => { last_batch = true; actual_batch_size = i; },
            }
        }

        if actual_batch_size == 0 { continue }

        let batch_out: Vec<Option<(RecordDto, WorkerStatistics)>> = batch_in
            .into_par_iter()
            .map(|rec| worker.process(&rec))
            .collect();
        for result in batch_out {
            if let Some((processed, w_stats)) = result {
                worker_stats += w_stats;
                if !w_stats.is_discarded() {
                    writer.write(&processed).expect("Failed to write record");
                }
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
