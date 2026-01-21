use std::io::Write;
use log::{error, LevelFilter};
use chrono::Local;
use std::process::exit;
use crate::args::Args;

use std::sync::atomic::{AtomicBool, Ordering};

static LOG_TO_STDERR: AtomicBool = AtomicBool::new(false);

pub fn setup_logging(args: &Args) {
    let log_level = if args.debug {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    if args.stderr_log {
        LOG_TO_STDERR.store(true, Ordering::SeqCst);
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
                    eprintln!("Failed to create log directory {}: {}", parent.display(), e);
                    exit(1);
                });
            }
        }

        let log_file = std::fs::File::options()
            .append(true)
            .create(true)
            .open(&log_path)
            .unwrap_or_else(|e| {
                eprintln!("Failed to open log file {}: {}", log_path.display(), e);
                exit(1);
            });

        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level.to_string()))
            .target(env_logger::Target::Pipe(Box::new(log_file)))
            .format(|buf, record| {
                writeln!(buf, "[{}] {}", Local::now().format("%Y-%m-%d %H:%M:%S"), record.args())
            })
            .init();
    }
}

pub fn get_log_path(args: &Args) -> std::path::PathBuf {
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

pub fn log_fatal(msg: &str) -> ! {
    error!("{}", msg);
    if !LOG_TO_STDERR.load(Ordering::SeqCst) {
        eprintln!("[{}] {}", Local::now().format("%Y-%m-%d %H:%M:%S"), msg);
    }
    exit(1);
}
