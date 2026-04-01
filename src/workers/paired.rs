use log::debug;
use crate::context::Context;
use crate::args::Args;
use crate::types::{PairedRecordDto, SingeRecordDto, QualTrimType};
use crate::steps::paired_adapter_remover::{PairedAdapterRemover, PairedAdapterRemoverState};
use crate::steps::paired_barcode_splitter::{PairedBarcodeSplitter, PairedBarcodeSplitterState};
use crate::steps::paired_quality_trim::{PairedQualityTrim, PairedQualityTrimState};
use crate::steps::paired_filter_length::{PairedFilterLength, PairedFilterLengthState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable};
use crate::steps::utils::paired::PairedAlignmentResults;

pub struct PairedWorker {
    pub barcode_splitter: Option<PairedBarcodeSplitter>,
    pub adapter_remover: Option<PairedAdapterRemover>,
    pub quality_trim: Option<PairedQualityTrim>,
    pub filter_length: PairedFilterLength,
    pub qtrim_post_removal: bool,
}

pub struct PairedWorkerState {
    pub barcode_splitter: Option<PairedBarcodeSplitterState>,
    pub adapter_remover: Option<PairedAdapterRemoverState>,
    pub quality_trim: Option<PairedQualityTrimState>,
    pub filter_length: PairedFilterLengthState,
}

impl PairedWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let quality_trim = match args.qtrim {
            QualTrimType::Off => None,
            _ => Some(PairedQualityTrim::new(args, context))
        };
        let barcode_splitter = match context.barcodes.is_empty() {
            true => None,
            _ => Some(PairedBarcodeSplitter::new(args, context))
        };
        let adapter_remover = match context.adapters.is_empty() && context.adapters2.is_empty() {
            true => None,
            _ => Some(PairedAdapterRemover::new(args, context))
        };
        let filter_length = PairedFilterLength::new(args, context);
        
        let qtrim_post_removal = args.qtrim_post_removal;
        if args.trace {
            debug!("Paired worker pipeline:");
            if !qtrim_post_removal && let Some(s) = &quality_trim {
                debug!("- Step PairedQualityTrim (threshold: {}, window size: {})",
                    s.inner.threshold, s.inner.win_size
                );
            }
            if let Some(s) = &barcode_splitter {
                debug!("- Step PairedBarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})",
                    s.trim_end, s.barcode_keep, s.barcode_unassigned
                );
            }
            if let Some(_) = &adapter_remover {
                debug!("- Step PairedAdapterRemover (trim end {:?})", args.adapter_trim_end);
            }
            if qtrim_post_removal && let Some(s) = &quality_trim {
                debug!("- Step PairedQualityTrim (threshold: {}, window size: {})",
                    s.inner.threshold, s.inner.win_size
                );
            }
            debug!("- Step PairedFilterLength (min length: {})",
                filter_length.min_read_length
            );
        }

        Self {
            barcode_splitter,
            adapter_remover,
            quality_trim,
            filter_length,
            qtrim_post_removal,
        }
    }

    pub fn create_state(&self) -> PairedWorkerState {
        PairedWorkerState {
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            adapter_remover: self.adapter_remover.as_ref().map(|s| s.create_state()),
            quality_trim: self.quality_trim.as_ref().map(|s| s.create_state()),
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

        if !self.qtrim_post_removal && let Some(step) = &self.quality_trim {
            rec_in = step.apply(rec_in, state.quality_trim.as_mut().unwrap());
        }

        if let Some(step) = &self.barcode_splitter {
            rec_in = step.apply(rec_in, state.barcode_splitter.as_mut().unwrap());
        }

        if let Some(step) = &self.adapter_remover {
            rec_in = step.apply(rec_in, state.adapter_remover.as_mut().unwrap());
        }

        if self.qtrim_post_removal && let Some(step) = &self.quality_trim {
            rec_in = step.apply(rec_in, state.quality_trim.as_mut().unwrap());
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
    use crate::types::{SingeRecordDto, PairedRecordDto, TrimEnd, FileFormat};

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
        
        assert_eq!(res.output_tag, Some("barcode_unassigned".to_string()));
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

    #[test]
    fn test_quality_trimming() {
        let mut args = Args::default();
        args.qtrim = QualTrimType::Tail;
        args.qtrim_threshold = 20; // 20 + 33 = 53 ('5')
        args.min_read_length = 0;
        
        // R1: 5555544444 (qual: 8888877777 - all above 53) -> No trim
        // R2: 5555544444 (qual: 8888833333 - last 5 below 53) -> Trimmed to 5 bases
        // Note: '8' is 56, '7' is 55, '3' is 51. 51 < 53.
        let r1 = SingeRecordDto::with_attrs("R1", None, b"5555544444", b"8888877777");
        let r2 = SingeRecordDto::with_attrs("R2", None, b"5555544444", b"8888833333");
        let rec = PairedRecordDto::new(r1, r2);
        
        let mut context = Context::empty();
        context.format = FileFormat::Fastq;
        let worker = PairedWorker::new(&args, &context);
        let mut state = worker.create_state();
        let res = worker.process_record(&rec, &mut state);
        
        assert!(res.record.is_some());
        let records = res.record.unwrap().records;
        assert_eq!(records[0].seq(), b"5555544444");
        assert_eq!(records[1].seq(), b"55555");
        assert_eq!(res.stats.bases_quality_trim, 5);
    }

    fn process_with_reference_setup(record: &PairedRecordDto) -> FinalRecord<PairedRecordDto> {
        let args = Args::parse_from(&[
            "flexbar_rs",
            "--min-read-length", "18",
            "--max-uncalled", "1",
            "--qtrim", "TAIL",
            "--qtrim-threshold", "10",
            "--barcode-unassigned",
            "--post-trim-length", "248",
            "--adapter-min-overlap", "6",
            "--qtrim-format", "sanger",
        ]);

        let barcodes = vec![
            SequenceWithId { id: "SmartSeqUMI".to_string(), sequence: "ATTGCGCAATGNNNNNNNN".to_string() },
        ];

        let adapters = vec![
            SequenceWithId { id: "Trans1".to_string(), sequence: "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans1_rc".to_string(), sequence: "CTGTCTCTTATACACATCTGACGCTGCCGACGA".to_string() },
            SequenceWithId { id: "Trans2".to_string(), sequence: "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans2_rc".to_string(), sequence: "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC".to_string() },
        ];

        let mut context = Context::empty();
        context.barcodes = barcodes;
        context.adapters = adapters.clone();
        context.adapters2 = adapters;
        context.format = FileFormat::Fastq;

        let worker = PairedWorker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_near_miss_with_adapter_trim() {
        let r1 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:17562:1031 1:N:0:TGCCGGTCAG+TTGTATCAGG",
            None,
            b"ATTGCGCAATGGTAGACTGGGGCTAGACTAGGCTGTCTTTGAATAACACTTAAAAATACATTAAAAGGAAAACCCAACTAATAACCTTGTTATTCTCTTCC",
            b"FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
        );
        let r2 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:17562:1031 2:N:0:TGCCGGTCAG+TTGTATCAGG",
            None,
            b"GAATACAGAGATGTGGTCCCCACCATCAAGGGACCAGAAACTTATTGCACAGTAAGTTGAATGTTGTGGAACTTTGATAGGTATTTCCTGAAAAGCTGACA",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFF,FFFFFFFFFF",
        );

        let rec = PairedRecordDto::new(r1, r2);
        let res = process_with_reference_setup(&rec);

        assert_eq!(res.output_tag, Some("barcode_SmartSeqUMI".to_string()));
        
        assert_eq!(res.record.as_ref().unwrap().records[0].seq(),
                   b"GGGCTAGACTAGGCTGTCTTTGAATAACACTTAAAAATACATTAAAAGGAAAACCCAACTAATAACCTTGTTATTCTCTTCC"
        );
        
        let stats = res.stats;
        assert_eq!(stats.total_adapters_removed, 0, "Expected no adapters to be found");
    }

    #[test]
    fn test_short_read_discard() {
        let r1 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:5936:2253 1:N:0:TGCCGGTCAG+GTGTATCAGG",
            None,
            b"ATTGCGCAATGGCTGTCTCTTATACACATCTCCGAGCCCACGAGACTGCCGGTCAGATCTCGTATGCCGGCTTCTGCTTGAAAAGGGGGGGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF,FF::FF,FF,:,F,,F,::F,F,,,:,F:FFFFFFFFFFF!",
        );
        let r2 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:5936:2253 2:N:0:TGCCGGTCAG+GTGTATCAGG",
            None,
            b"CCATTGCGCAATCTGTCTCTTATACACATCTGACGCTGCCGACGATTGTATCAGGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAATGGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFF::,,,FFFFFFFFFF!",
        );

        let rec = PairedRecordDto::new(r1, r2);
        let res = process_with_reference_setup(&rec);

        assert!(res.record.is_none(), "Record should be discarded");
    }

    #[test]
    fn test_zero_length_after_adapter_removal() {
        // This one is an actual change in behaviour to the old flexbar.
        // First read: ATTGCGCAATGACAGGACCTGTCTCTTATACACATCTCCGAGCCCAC ...
        // Barcode:    ATTGCGCAATGNNNNNNNN    - found and removed by both programs
        // Adapter:                      CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        //
        // Please note:
        // - The adapter overlaps with the barcode in the first base.
        // - flexbar_rs will detect the adapter due to the default --adapter-error-rate of 0.1.
        // - flexbar fails to detect/remove the adapter (why?).
        // - If detected, the default for --adapter-trim-end "RIGHT" will be used, truncating the
        //   read to 0 (!) bases, leading to a discard.
        let r1 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:16468:4993",
            Some("1:N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"ATTGCGCAATGACAGGACCTGTCTCTTATACACATCTCCGAGCCCACGAGACTGCCGGTCAGATCTCGTATGCCGTCTTCTGCTTGAAAAGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,::,,F:,,,F,,:,,F,,F::,,,::F:FF:F",
        );
        let r2 = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1101:16468:4993",
            Some("2:N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"GTCCTGTCATTGCGCAATCTGTCTCTTATACACATCTGACGCTGCCGACGATTGTATCAGGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFF,:,,,FF:",
        );

        let rec = PairedRecordDto::new(r1, r2);
        let res = process_with_reference_setup(&rec);

        assert!(res.record.is_none(), "Record should be discarded");
    }
}
