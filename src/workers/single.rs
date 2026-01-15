pub struct SingleWorker;

impl SingleWorker {
    pub fn new() -> Self {
        Self
    }

    pub fn process<T>(&self, input: T) -> Option<T> {
        Some(input)
    }
}