pub mod args;
pub mod config;
pub mod context;
pub mod types;
pub mod utils;
pub mod statistics;
mod readers;
mod writers;
mod workers;
mod workflow;
mod steps;

use std::env;
use std::io::Write;
use log::info;
use chrono::Local;

fn main() {
    let args = args::Args::get_args();
    let config = config::Config::new(args.clone());

    setup_logging(&config.args);

    info!("Command line: {}", env::args().collect::<Vec<_>>().join(" "));

    let threads = config.args.threads;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    workflow::do_it(config);
}

fn setup_logging(args: &args::Args) {
    let log_level = if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };

    if args.stderr_log {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level.to_string()))
            .format(|buf, record| {
                writeln!(buf, "[{}] {}", Local::now().format("%Y-%m-%d %H:%M:%S"), record.args())
            })
            .init();
    } else {
        let log_path = get_log_path(args);

        if let Some(parent) = log_path.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                std::fs::create_dir_all(parent).unwrap_or_else(|e| {
                    panic!("Failed to create log directory {}: {}", parent.display(), e);
                });
            }
        }

        let log_file = std::fs::File::options()
            .append(true)
            .create(true)
            .open(&log_path)
            .unwrap_or_else(|e| {
                panic!("Failed to open log file {}: {}", log_path.display(), e);
            });

        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level.to_string()))
            .target(env_logger::Target::Pipe(Box::new(log_file)))
            .format(|buf, record| {
                writeln!(buf, "[{}] {}", Local::now().format("%Y-%m-%d %H:%M:%S"), record.args())
            })
            .init();
    }
}

fn get_log_path(args: &args::Args) -> std::path::PathBuf {
    if let Some(path) = &args.output_log {
        path.clone()
    } else {
        let mut path = if let Some(out) = &args.output_reads {
            out.clone()
        } else {
            std::path::PathBuf::from(&args.target)
        };
        let mut filename = path.file_name().unwrap_or_default().to_os_string();
        filename.push(".log");
        path.set_file_name(filename);
        path
    }
}
