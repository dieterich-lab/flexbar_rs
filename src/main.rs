pub mod args;
pub mod config;
pub mod types;

fn main() {
    let args = args::Args::get_args();
    setup_logging(&args);
    let _ = config::Config::new(args);
}

fn setup_logging(args: &args::Args) {
    if args.debug {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::init();
    }

}
