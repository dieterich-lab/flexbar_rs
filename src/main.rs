pub mod args;

fn main() {
    let _args = args::Args::get_args();
    println!("Hello, flexbar_rs!");
}
