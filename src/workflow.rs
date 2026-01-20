use rayon::prelude::*;
use crate::config::Config;
use crate::context::Context;
use crate::readers::fasta::SingleReaderFasta;
use crate::readers::fastq::SingleReaderFastq;
use crate::readers::SequenceReader;
use crate::writers::fasta::SingleWriterFasta;
use crate::writers::fastq::SingleWriterFastq;
use crate::types::{FileFormat, CompressionType};
use crate::workers::single::SingleWorker;
use std::path::Path;
use crate::writers::Writer;
use crate::statistics::{Statistics, WorkerStatistics};

use std::sync::mpsc::sync_channel;
use std::thread;

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
    let reader = get_single_reader(reads_path, config.format);
    let compression = get_compression(&reader, config);
    let worker = SingleWorker::new(config, context);
    let mut writer = get_single_writer(config, output_prefix, compression);

    let (tx, rx) = sync_channel(100);
    let (stats_tx, stats_rx) = std::sync::mpsc::channel();

    let writer_stats = thread::scope(|s| {
        s.spawn(move || {
            let reader_stats_val = reader.statistics();
            let (worker_stats, reader_stats) = reader
                .par_bridge()
                .map(|record| {
                    let result = worker.process(record.clone());
                    (record, result)
                })
                .fold(|| (WorkerStatistics::default(), crate::statistics::ReaderStatistics::default()), |(mut w_acc, r_acc), (_record, result)| {
                    if let Some((processed, w_stats)) = result {
                        w_acc += w_stats;
                        
                        if !w_stats.is_discarded() {
                            tx.send(processed).expect("Failed to send record to writer");
                        }
                    }
                    (w_acc, r_acc)
                })
                .reduce(|| (WorkerStatistics::default(), crate::statistics::ReaderStatistics::default()), |(w1, r1), (w2, r2)| (w1 + w2, r1 + r2));

            stats_tx.send((reader_stats_val + reader_stats, worker_stats)).ok();
        });

        for record in rx {
            writer.write(&record).expect("Failed to write record");
        }
        writer.statistics()
    });

    let (reader_stats, worker_stats) = stats_rx.recv().expect("Failed to receive statistics from worker thread");
    
    Statistics {
        reader: reader_stats,
        worker: worker_stats,
        writer: writer_stats,
    }
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

fn get_single_writer(config: &Config, output_prefix: &str, compression: CompressionType) -> Box<dyn Writer<bio::io::fastq::Record>> {
    if config.args.stdout_reads {
        if config.args.fasta_output {
            Box::new(SingleWriterFasta::new("-", compression).expect("Failed to create FASTA writer"))
        } else {
            Box::new(SingleWriterFastq::new("-", compression).expect("Failed to create FASTQ writer"))
        }
    } else if config.args.fasta_output {
        Box::new(SingleWriterFasta::new(output_prefix, compression).expect("Failed to create FASTA writer"))
    } else {
        Box::new(SingleWriterFastq::new(output_prefix, compression).expect("Failed to create FASTQ writer"))
    }
}
