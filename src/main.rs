pub mod args;
pub mod config;
pub mod context;
pub mod types;
pub mod utils;
mod readers;
mod writers;
mod workers;
mod workflow;

fn main() {
    let args = args::Args::get_args();
    setup_logging(&args);
    let threads = args.threads;
    let config = config::Config::new(args);
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    workflow::do_it(config);
}

fn setup_logging(args: &args::Args) {
    if args.debug {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::init();
    }

}
