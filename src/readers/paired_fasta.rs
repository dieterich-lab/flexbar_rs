use std::io::{self, BufRead};
use crate::types::{RecordDto, PairedRecordDto};
use crate::readers::generic::GenericReaderBackend;
use crate::readers::fasta::FastaBackend;
use crate::logging::log_fatal;

pub struct PairedFastaBackend<R: BufRead + Send> {
    backend1: FastaBackend<R>,
    backend2: FastaBackend<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R> for PairedFastaBackend<R> {
    type Error = io::Error;

    fn new(mut readers: Vec<R>) -> Self {
        if readers.len() != 2 {
            log_fatal("PairedFastaBackend requires exactly two input files.");
        }
        let r2 = readers.pop().unwrap();
        let r1 = readers.pop().unwrap();
        Self {
            backend1: FastaBackend::new(vec![r1]),
            backend2: FastaBackend::new(vec![r2]),
        }
    }

    fn read_record(&mut self) -> Result<Option<RecordDto>, Self::Error> {
        let rec1 = self.backend1.read_record()?;
        let rec2 = self.backend2.read_record()?;

        match (rec1, rec2) {
            (Some(RecordDto::Single(r1)), Some(RecordDto::Single(r2))) => {
                let paired = PairedRecordDto::new(r1, r2);
                if !paired.validate() {
                    log_fatal(&format!("Validation failed for paired records: {} and {}", paired.records[0].id(), paired.records[1].id()));
                }
                Ok(Some(RecordDto::Paired(paired)))
            },
            (None, None) => Ok(None),
            _ => {
                log_fatal("Paired files have different number of records.");
            }
        }
    }
}
