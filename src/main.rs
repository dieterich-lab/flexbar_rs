pub mod args;
pub mod config;
pub mod context;
pub mod types;
pub mod utils;
pub mod statistics;
pub mod logging;
mod readers;
mod writers;
mod workers;
mod workflow;
mod steps;

use std::env;
use log::info;

fn main() {
    let args = args::Args::get_args();
    let config = config::Config::new(args.clone());

    logging::setup_logging(&config.args);

    info!("Command line: {}", env::args().collect::<Vec<_>>().join(" "));

    let threads = config.args.threads;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    workflow::do_it(config);
}
