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

fn main() {
    let args = args::Args::get_args();
    setup_logging(&args);
    let threads = args.threads;
    let config = config::Config::new(args);
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    workflow::do_it(config);
}

fn setup_logging(args: &args::Args) {
    let mut builder = env_logger::Builder::new();
    builder.target(env_logger::Target::Stderr);
    
    if args.debug {
        builder.filter_level(log::LevelFilter::Debug);
    } else {
        builder.filter_level(log::LevelFilter::Info);
    }
    
    if let Ok(rust_log) = std::env::var("RUST_LOG") {
        builder.parse_filters(&rust_log);
    }
    
    builder.init();
}
