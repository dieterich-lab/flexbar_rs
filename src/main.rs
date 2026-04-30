pub mod args;
pub mod context;
pub mod types;
pub mod utils;
pub mod statistics;
pub mod logging;
mod readers;
mod writers;
mod workflow;
mod steps;
mod multi_aligner;
mod worker;
mod records;

extern crate jemallocator;

#[global_allocator]
static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

use std::env;
use log::info;

fn main() {
    let args = args::Args::get_args();
    logging::setup_logging(&args);
    info!("Command line: {}", env::args().collect::<Vec<_>>().join(" "));
    let threads = args.threads;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    workflow::do_it(args);
}
