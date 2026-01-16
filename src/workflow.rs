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

use std::sync::mpsc::sync_channel;
use std::thread;

pub fn do_it(config: Config) {
    let mut context = Context::new();
    context.load(&config);

    let reads_path = config.args.reads.as_ref().expect("No input reads provided");
    let output_prefix = config.args.output_reads.as_ref()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| config.args.target.clone());

    process_single(&config, &context, reads_path, &output_prefix);
}

fn process_single(config: &Config, context: &Context, reads_path: &Path, output_prefix: &str) {
    let reader: Box<dyn SequenceReader + Send> = match config.format {
        FileFormat::Fasta => Box::new(SingleReaderFasta::from_path(reads_path).expect("Failed to create FASTA reader")),
        FileFormat::Fastq => Box::new(SingleReaderFastq::from_path(reads_path).expect("Failed to create FASTQ reader")),
    };

    let compression = match config.args.output_compression {
        Some(CompressionType::Keep) | None => reader.compression(),
        Some(c) => c,
    };

    let mut writer: Box<dyn Writer<bio::io::fastq::Record>> = if config.args.stdout_reads {
        if config.args.fasta_output {
            Box::new(SingleWriterFasta::new("-", compression).expect("Failed to create FASTA writer"))
        } else {
            Box::new(SingleWriterFastq::new("-", compression).expect("Failed to create FASTQ writer"))
        }
    } else if config.args.fasta_output {
        Box::new(SingleWriterFasta::new(output_prefix, compression).expect("Failed to create FASTA writer"))
    } else {
        Box::new(SingleWriterFastq::new(output_prefix, compression).expect("Failed to create FASTQ writer"))
    };
    
    let worker = SingleWorker::new(context);

    let (tx, rx) = sync_channel(100);

    thread::scope(|s| {
        s.spawn(move || {
            reader
                .par_bridge()
                .filter_map(|record| worker.process(record))
                .for_each_with(tx, |tx, record| {
                    tx.send(record).expect("Failed to send record to writer");
                });
        });

        for record in rx {
            writer.write(&record).expect("Failed to write record");
        }
    });
}
