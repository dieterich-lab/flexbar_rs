use crate::context::Context;

pub struct SingleWorker<'a> {
    context: &'a Context,
}

impl<'a> SingleWorker<'a> {
    pub fn new(context: &'a Context) -> Self {
        Self { context }
    }

    pub fn process<T>(&self, input: T) -> Option<T> {
        // Here we will eventually use self.context for adapter removal
        Some(input)
    }
}