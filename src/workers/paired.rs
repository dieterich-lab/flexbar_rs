use log::debug;
use crate::context::Context;
use crate::args::Args;
use crate::types::{PairedRecordDto, SingeRecordDto};
use crate::steps::paired_barcode_splitter::{PairedBarcodeSplitter, PairedBarcodeSplitterState};
use crate::steps::paired_filter_length::{PairedFilterLength, PairedFilterLengthState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable};
use crate::steps::utils::paired::PairedAlignmentResults;

pub struct PairedWorker {
    pub barcode_splitter: Option<PairedBarcodeSplitter>,
    pub filter_length: PairedFilterLength,
}

pub struct PairedWorkerState {
    pub barcode_splitter: Option<PairedBarcodeSplitterState>,
    pub filter_length: PairedFilterLengthState,
}

impl PairedWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let barcode_splitter = match context.barcodes.is_empty() {
            true => None,
            _ => Some(PairedBarcodeSplitter::new(args, context))
        };
        let filter_length = PairedFilterLength::new(args, context);
        
        if args.trace {
            debug!("Paired worker pipeline:");
            if let Some(s) = &barcode_splitter {
                debug!("- Step PairedBarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})",
                    s.trim_end, s.barcode_keep, s.barcode_unassigned
                );
            }
            debug!("- Step PairedFilterLength (min length: {})",
                filter_length.min_read_length
            );
        }

        Self {
            barcode_splitter,
            filter_length,
        }
    }

    pub fn create_state(&self) -> PairedWorkerState {
        PairedWorkerState {
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
        }
    }
    
    pub fn process_record(
        &self,
        record: &PairedRecordDto,
        state: &mut PairedWorkerState
    ) -> FinalRecord<PairedRecordDto> {
        let mut rec_in: RecordInTransformation<PairedRecordDto, PairedAlignmentResults> = record.get_record_in_transformation();
        rec_in.stats.total_bases_read = record.records[0].seq().len() + record.records[1].seq().len();

        if let Some(step) = &self.barcode_splitter {
            rec_in = step.apply(rec_in, state.barcode_splitter.as_mut().unwrap());
        }

        rec_in = self.filter_length.apply(rec_in, &mut state.filter_length);

        finalize(rec_in)
    }
}

fn finalize(rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>) -> FinalRecord<PairedRecordDto> {
    if rec.alignment_results.per_read[0].dropped && rec.alignment_results.per_read[1].dropped {
        return FinalRecord {
            record: None,
            output_tag: None,
            stats: rec.stats,
        }
    }

    let mut records = rec.record.records;
    let mut good_reads: u32 = 0;
    for i in 0..2 {
        let results = &rec.alignment_results.per_read[i];
        let old_rec = &records[i];
        
        if results.dropped {
            records[i] = SingeRecordDto::with_attrs(
                old_rec.id(),
                None,
                &[],
                &[],
            );
            continue;
        }
        good_reads += 1;

        let range = &results.read_range;
        
        let mut desc_parts = Vec::new();
        if let Some(desc) = old_rec.desc() {
            desc_parts.push(desc);
        }
        for tag in &results.removal_tags {
            desc_parts.push(tag);
        }
        for tag in &results.umi_tags {
            desc_parts.push(tag);
        }
        
        if rec.alignment_results.per_read[1-i].dropped {
            desc_parts.push("single");
        }

        let new_desc = if desc_parts.is_empty() {
            None
        } else {
            Some(desc_parts.join("_"))
        };

        records[i] = SingeRecordDto::with_attrs(
            old_rec.id(),
            new_desc.as_deref(),
            &old_rec.seq()[range.clone()],
            &old_rec.qual()[range.clone()],
        );
    }
    
    let record = if good_reads > 0 {
        Some(PairedRecordDto { records })
    } else {
        None
    };

    FinalRecord {
        record,
        output_tag: rec.alignment_results.barcode_tag,
        stats: rec.stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use crate::context::SequenceWithId;
    use crate::types::{SingeRecordDto, PairedRecordDto, TrimEnd};

    fn process_record(args: Args, record: &PairedRecordDto, barcodes: Vec<SequenceWithId>, barcodes2: Vec<SequenceWithId>) -> FinalRecord<PairedRecordDto> {
        let mut context = Context::empty();
        context.barcodes = barcodes;
        context.barcodes2 = barcodes2;
        let worker = PairedWorker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_barcode_splitting_better_match_first() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        let barcodes2 = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        // R1 starts with AAAAA (perfect match), R2 starts with TTTTT (no match)
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAAGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, barcodes2);
        
        assert_eq!(res.output_tag, Some("barcode_B1".to_string()));
        let records = res.record.unwrap().records;
        // B1 (AAAAA) should be removed from R1
        assert_eq!(records[0].seq(), b"GGGGG");
        // R2 should be unchanged
        assert_eq!(records[1].seq(), b"TTTTTGGGGG");
    }

    #[test]
    fn test_barcode_splitting_better_match_second() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        
        let barcodes = vec![
            SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() },
        ];
        let barcodes2 = vec![
            SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() },
        ];
        
        // R1 starts with TTTTT, R2 starts with CCCCC
        let r1 = SingeRecordDto::with_attrs("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"CCCCCGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, barcodes2);
        
        assert_eq!(res.output_tag, Some("barcode_B2".to_string()));
        let records = res.record.unwrap().records;
        assert_eq!(records[0].seq(), b"TTTTTGGGGG");
        assert_eq!(records[1].seq(), b"GGGGG");
    }

    #[test]
    fn test_barcode_splitting_no_match_unassigned() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-unassigned", "--min-read-length", "0"]);
        
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        let barcodes2 = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        let r1 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, barcodes2);
        
        assert_eq!(res.output_tag, Some("barcode_unassigned_".to_string()));
        assert!(res.record.is_some());
    }

    #[test]
    fn test_barcode_splitting_no_match_no_unassigned() {
        let args = Args::parse_from(&["flexbar_rs", "--min-read-length", "0"]);
        
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        let barcodes2 = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        let r1 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, barcodes2);
        
        assert_eq!(res.output_tag, None);
        assert!(res.record.is_none());
    }

    #[test]
    fn test_length_filter_drop_default() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        
        // R1 is 14 bases, B1 (AAAAA) at start, R2 is 10 bases.
        // After B1 removal: R1 becomes 9 bases (short).
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII"); // 14 bases
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");    // 10 bases
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, vec![]);
        
        // Default: both dropped if one is short
        assert!(res.record.is_none());
        assert_eq!(res.stats.total_single_reads_dropped, 2); // R2 is dropped as counterpart of short R1, and R1 also counted here?
        assert_eq!(res.stats.total_discarded_length, 2);    // Both R1 and R2 are discarded
    }

    #[test]
    fn test_length_filter_single_reads() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads = true;
        
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII"); // 14 -> 9 (short)
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");    // 10 (good)
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, vec![]);
        
        assert!(res.record.is_some());
        let records = res.record.unwrap().records;
        assert_eq!(records[0].seq(), b""); // R1 dropped
        assert_eq!(records[1].seq(), b"0123456789"); // R2 kept
        assert_eq!(res.stats.total_single_reads_written, 1);
    }

    #[test]
    fn test_length_filter_single_reads_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads_paired = true;
        
        let barcodes = vec![
            SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() },
        ];
        
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII"); // 14 -> 9 (short)
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");    // 10 (good)
        let rec = PairedRecordDto::new(r1, r2);
        
        let res = process_record(args, &rec, barcodes, vec![]);
        
        assert!(res.record.is_some());
        let records = res.record.unwrap().records;
        assert_eq!(records[0].seq(), b"N");
        assert_eq!(records[1].seq(), b"0123456789");
        assert_eq!(res.stats.total_single_reads_replaced, 1);
    }
}
