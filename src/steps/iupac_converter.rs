use crate::args::Args;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use crate::records::generic_read::GenericRawRead;

pub struct IupacConverter {}

pub struct IupacConverterState {}

impl GenericStep for IupacConverter {
    type State = IupacConverterState;

    fn new(_args: &Args, _context: &Context) -> Self {
        Self {}
    }

    fn create_state(&self) -> Self::State {
        IupacConverterState {}
    }

    fn apply(&self, mut rec: RecordInTransformation, _state: &mut Self::State) -> RecordInTransformation {
        for i in 0..2 {
            if rec.record.reads[i].is_none() || rec.is_dropped_read(i) {
                continue;
            }
            if rec.record.reads[i].is_some() {
                let (converted, count) = get_iupac_converted_read(
                    rec.get_sequence(i)
                );
                if count != 0 {
                    rec.record.patch_read_sequence(i, &converted);
                }
            }
        }
        if rec.record.barcode_read.is_some() {
            let (converted, count) = get_iupac_converted_read(rec.record.barcode_read.as_ref().unwrap().seq());
            if count != 0 {
                rec.record.patch_barcode_read(&converted);
            }
        }
        rec
    }
}

fn get_iupac_converted_read(raw: &[u8]) -> (Vec<u8>, usize) {
    let mut result: Vec<u8> = Vec::with_capacity(raw.len());
    let mut count: usize = 0;
    for base in raw.iter() {
        match *base {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' =>  {
                result.push(*base);
            },
            _ => {
                result.push(b'N');
                count += 1;
            }
        }
    }
    (result, count)
}

#[cfg(test)]
mod tests {
    use crate::records::unified::{UnifiedRecord, UnifiedRead, GenericRawRead};
    use super::*;
    #[test]
    fn test_iupac_conversion() {
        let args = Args::default();
        let context = Context::empty();
        let step = IupacConverter::new(&args, &context);
        let inp_seq = b"ATCGNRYSWKMBDHVatcgnryswkmbdhv";
        let out_seq = b"ATCGNNNNNNNNNNNatcgNNNNNNNNNNN";

        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("read1", None, inp_seq, b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
        let record = UnifiedRecord::single(r1);
        
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
        let output_seq = res.record.reads[0].as_ref().unwrap().seq();

        assert_eq!(output_seq, out_seq);
    }

    #[test]
    fn test_iupac_conversion_paired() {
        let args = Args::default();
        let context = Context::empty();
        let step = IupacConverter::new(&args, &context);
        let inp_seq1 = b"ATCGN";
        let out_seq1 = b"ATCGN";
        let inp_seq2 = b"RYSwk";

        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("read1", None, inp_seq1, b"IIIII");
        let r2 = UnifiedRead::new_fastq("read2", None, inp_seq2, b"IIIII");
        let record = UnifiedRecord::paired(r1, r2);
        
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
        assert_eq!(res.record.reads[0].as_ref().unwrap().seq(), out_seq1);
        assert_eq!(res.record.reads[1].as_ref().unwrap().seq(), b"NNNNN");
    }

    #[test]
    fn test_iupac_conversion_with_barcode() {
        let args = Args::default();
        let context = Context::empty();
        let step = IupacConverter::new(&args, &context);
        
        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("read1", None, b"ATCG", b"IIII");
        let bc = UnifiedRead::new_fastq("bc1", None, b"NRY", b"III");
        let record = UnifiedRecord::single_barcoded(r1, bc);
        
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
        assert_eq!(res.record.barcode_read.as_ref().unwrap().seq(), b"NNN");
    }
}
