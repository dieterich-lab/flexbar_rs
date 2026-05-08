use log::debug;
use crate::args::Args;
use crate::context::Context;
use crate::records::unified::{UnifiedRecord, AnnotatedRecord};
use crate::records::generic_read::GenericRawRead;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::{adapter_remover::AdapterRemoverStep, barcode_splitter::BarcodeSplitter, filter_length::FilterLength, iupac_converter::IupacConverter, post_trim::PostTrim, pre_trim::PreTrim, quality_trim::QualityTrim, uncalled_filter::UncalledFilter};
use crate::steps::pre_trim::PreTrimState;
use crate::steps::filter_length::FilterLengthState;
use crate::steps::post_trim::PostTrimState;
use crate::steps::barcode_splitter::BarcodeSplitterState;
use crate::steps::adapter_remover::AdapterRemoverStepState;
use crate::steps::iupac_converter::IupacConverterState;
use crate::types::QualTrimType;

pub struct Worker {
    pub iupac_converter: Option<IupacConverter>,
    pub pre_trim: Option<PreTrim>,
    pub barcode_splitter: Option<BarcodeSplitter>,
    pub adapter_remover: Option<AdapterRemoverStep>,
    pub quality_trim: Option<QualityTrim>,
    pub filter_length: FilterLength,
    pub post_trim: PostTrim,
    pub uncalled_filter: Option<UncalledFilter>,
    pub qtrim_post_removal: bool,
}

pub struct WorkerState {
    pub iupac_converter: Option<IupacConverterState>,
    pub pre_trim: Option<PreTrimState>,
    pub barcode_splitter: Option<BarcodeSplitterState>,
    pub adapter_remover: Option<AdapterRemoverStepState>,
    pub quality_trim: Option<()>,
    pub filter_length: FilterLengthState,
    pub post_trim: PostTrimState,
    pub uncalled_filter: Option<()>,
}

impl Worker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let iupac_converter = if args.iupac {
            Some(IupacConverter::new(args, context))
        } else { None };

        let pre_trim = if args.pre_trim_left > 0 || args.pre_trim_right > 0 {
            Some(PreTrim::new(args, context))
        } else { None };

        let uncalled_filter = if args.max_uncalled > 0 {
            Some(UncalledFilter::new(args, context))
        } else { None };

        let quality_trim = match args.qtrim {
            QualTrimType::Off => None,
            _ => Some(QualityTrim::new(args, context)),
        };

        let barcode_splitter = match context.barcodes.is_empty() && context.barcodes2.is_empty() {
            true => None,
            _ => Some(BarcodeSplitter::new(args, context)),
        };

        let adapter_remover = match context.adapters.is_empty() && context.adapters2.is_empty() {
            true => None,
            _ => Some(AdapterRemoverStep::new(args, context)),
        };

        let filter_length = FilterLength::new(args, context);
        let post_trim = PostTrim::new(args, context);

        let qtrim_post_removal = args.qtrim_post_removal;

        debug!("Unified worker pipeline:");
        if iupac_converter.is_some() {
            debug!("- Step IupacConverter");
        }
        if let Some(s) = &pre_trim {
            debug!("- Step PreTrim (left: {}, right: {})", s.pre_trim_left, s.pre_trim_right);
        }
        if let Some(s) = &uncalled_filter {
            debug!("- Step UncalledFilter (max_uncalled: {})", s.max_uncalled);
        }
        if !qtrim_post_removal { if let Some(s) = &quality_trim {
            debug!("- Step QualityTrim (threshold: {}, window size: {})", s.threshold, s.win_size);
        }}
        if let Some(s) = &barcode_splitter {
            debug!("- Step BarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})", s.trim_end, s.barcode_keep, s.barcode_unassigned);
        }
        if adapter_remover.is_some() {
            debug!("- Step AdapterRemoverStep");
        }
        if qtrim_post_removal { if let Some(s) = &quality_trim {
            debug!("- Step QualityTrim (threshold: {}, window size: {})", s.threshold, s.win_size);
        }}
        debug!("- Step FilterLength (min length: {})", filter_length.min_read_length);
        if post_trim.post_trim_length > 0 {
            debug!("- Step PostTrim (length: {})", post_trim.post_trim_length);
        }

        Self {
            iupac_converter, pre_trim, barcode_splitter, adapter_remover, quality_trim,
            filter_length, post_trim, uncalled_filter, qtrim_post_removal,
        }
    }

    pub fn create_state(&self) -> WorkerState {
        WorkerState {
            iupac_converter: self.iupac_converter.as_ref().map(|s| s.create_state()),
            pre_trim: self.pre_trim.as_ref().map(|s| s.create_state()),
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            adapter_remover: self.adapter_remover.as_ref().map(|s| s.create_state()),
            quality_trim: self.quality_trim.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
            post_trim: self.post_trim.create_state(),
            uncalled_filter: self.uncalled_filter.as_ref().map(|s| s.create_state()),
        }
    }

    pub fn process_record(&self, record: &UnifiedRecord, state: &mut WorkerState)
        -> (AnnotatedRecord, WorkerStatistics)
    {
        let mut rec_in = RecordInTransformation::from_dto(record.clone());
        let mut total_bases = 0;
        for read in &record.reads {
            if let Some(r) = read {
                total_bases += r.seq().len();
            }
        }
        rec_in.stats.total_bases_read = total_bases;

        if let Some(step) = &self.iupac_converter {
            rec_in = step.apply(rec_in, state.iupac_converter.as_mut().unwrap());
        }
        if let Some(step) = &self.pre_trim {
            rec_in = step.apply(rec_in, state.pre_trim.as_mut().unwrap());
        }
        if let Some(step) = &self.uncalled_filter {
            rec_in = step.apply(rec_in, state.uncalled_filter.as_mut().unwrap());
        }
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
        if self.post_trim.post_trim_length > 0 { rec_in = self.post_trim.apply(rec_in, &mut state.post_trim); }
        self.finalize(rec_in)
    }

    fn finalize(&self, rec: RecordInTransformation) -> (AnnotatedRecord, WorkerStatistics) {
        let umi_tags = if rec.alignment_results.per_read[0].is_some() && !rec.alignment_results.per_read[0].as_ref().unwrap().umi_tags.is_empty() {
            rec.alignment_results.per_read[0].as_ref().unwrap().umi_tags.clone()
        } else if rec.alignment_results.per_read[1].is_some() {
            rec.alignment_results.per_read[1].as_ref().unwrap().umi_tags.clone()
        } else {
            vec![]
        };

        rec.record.finalize(rec.alignment_results, &umi_tags, rec.stats)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use crate::context::{Context, SequenceWithId};
    use crate::records::unified::{UnifiedRecord, UnifiedRead};
    use crate::records::generic_read::GenericRawRead;
    use crate::types::{FileFormat, TrimEnd};

    fn process_record_unified(args: &Args, record: &UnifiedRecord) -> (AnnotatedRecord, WorkerStatistics) {
        let context = Context::new(args);
        let worker = Worker::new(args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    fn process_record_paired_helper(args: Args, record: &UnifiedRecord, barcodes: Vec<SequenceWithId>, barcodes2: Vec<SequenceWithId>) -> (AnnotatedRecord, WorkerStatistics) {
        let mut context = Context::empty();
        context.barcodes = barcodes;
        context.barcodes2 = barcodes2;
        context.run_type = crate::types::RunType::Paired;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(&record, &mut state)
    }

    // Copied and adapted from src/workers/paired.rs
    #[test]
    fn test_barcode_splitting_better_match_first_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"AAAAAGGGGG", b"IIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.barcode_id, Some("B1".to_string()));
        let records = match res.record.reads { [Some(r1), Some(r2)] => [r1, r2], _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"GGGGG");
        assert_eq!(records[1].seq(), b"TTTTTGGGGG");
    }

    #[test]
    fn test_barcode_splitting_better_match_second_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"CCCCCGGGGG", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.barcode_id, Some("B2".to_string()));
        let records = match res.record.reads { [Some(r1), Some(r2)] => [r1, r2], _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"TTTTTGGGGG");
        assert_eq!(records[1].seq(), b"GGGGG");
    }

    #[test]
    fn test_barcode_splitting_no_match_unassigned_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-unassigned", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.barcode_id, Some("unassigned".to_string()));
        assert!(res.record.reads[0].is_some());
        assert!(res.record.reads[1].is_some());
    }

    #[test]
    fn test_barcode_splitting_no_match_no_unassigned_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.barcode_id, None);
        assert!(res.record.reads[0].is_none());
        assert!(res.record.reads[1].is_none());
    }

    #[test]
    fn test_length_filter_drop_default_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, vec![]);
        assert!(res.record.reads[0].is_none());
        assert!(res.record.reads[1].is_none());
    }

    #[test]
    fn test_length_filter_single_reads_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads = true;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, vec![]);
        assert!(res.record.reads[0].is_none());
        if let Some(r2) = res.record.reads[1].as_ref() {
            assert_eq!(r2.seq(), b"0123456789");
        } else {
            panic!("expected second read");
        }
    }

    #[test]
    fn test_length_filter_single_reads_paired_tagged_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads_paired = true;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = UnifiedRead::new_fastq("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_record_paired_helper(args, &rec, barcodes, vec![]);
        let records = match res.record.reads { [Some(r1), Some(r2)] => [r1, r2], _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"N");
        assert_eq!(records[1].seq(), b"0123456789");
    }

    #[test]
    fn test_quality_trimming_paired() {
        let mut args = Args::default();
        args.qtrim = QualTrimType::Tail;
        args.qtrim_threshold = 20;
        args.min_read_length = 0;
        let r1 = UnifiedRead::new_fastq("R1", None, b"5555544444", b"8888877777");
        let r2 = UnifiedRead::new_fastq("R2", None, b"5555544444", b"8888833333");
        let rec = UnifiedRecord::paired(r1, r2);
        let mut context = Context::empty();
        context.format = FileFormat::Fastq;
        context.run_type = crate::types::RunType::Paired;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let (res, _) = worker.process_record(&rec, &mut state);
        let records = match res.record.reads { [Some(r1), Some(r2)] => [r1, r2], _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"5555544444");
        assert_eq!(records[1].seq(), b"55555");
    }

    // Single worker tests adapted
    #[test]
    fn test_process_right_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let (result, _) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), Some(&b"IIIII"[..]));
    }

    #[test]
    fn test_process_left_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let (result, _) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), Some(&b"IIIII"[..]));
    }

    #[test]
    fn test_process_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(result.barcode_id, Some("bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_max_uncalled_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let record = UnifiedRead::new_fastq("read1", None, b"ACNNGT", b"IIIIII");
        let (_, stats) = process_record_unified(&args, &UnifiedRecord::single(record));
        assert_eq!(stats.total_discarded_uncalled, 0);
        let record2 = UnifiedRead::new_fastq("read2", None, b"ACNNGN", b"IIIIII");
        let (_, stats2) = process_record_unified(&args, &UnifiedRecord::single(record2));
        assert_eq!(stats2.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"GGGGACGT", b"IIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(result.removal_tags[0].join("_"), "Flexbar_removal");
    }

    #[test]
    fn test_removal_tags_adapter_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "ad1".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(result.removal_tags[0].join("_"), "Flexbar_removal_ad1");
    }

    #[test]
    fn test_removal_tags_both_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq",
            "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4",
            "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5",
            "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "ad1".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"GGGGACGTCTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(result.removal_tags[0].join("_"), "Flexbar_removal_Flexbar_removal_ad1");
    }

    #[test]
    fn test_quality_trim_tail_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "TAIL", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (result, _) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), Some(&b"IIIII"[..]));
    }

    #[test]
    fn test_quality_trim_win_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "WIN", "--qtrim-threshold", "20", "--qtrim-win-size", "5", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTACGTACGT", b"IIIII44444III");
        let (result, _) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_quality_trim_bwa_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "BWA", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTACGT", b"IIIII4444");
        let (result, _) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "NNGGGGNN".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"AAGGGGTTACGT", b"IIIIIIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(result.barcode_id, Some("bc1".to_string()));
        assert_eq!(processed.id(), "read1");
        assert_eq!(result.umi_tags.join("_"), "AA_TT");
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_umi_tags_adapter_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATNNNN", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "adapter".to_string(), sequence: "AGATNNNN".to_string() }];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"ACTGTAGATCCCC", b"IIIIIIIIIIIII");
        let (result, _) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(result.umi_tags.join("_"), "CCCC");
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_min_read_length_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let record = UnifiedRead::new_fastq("read1", None, b"0123456789", b"IIIIIIIIII");
        let (_, stats) = process_record_unified(&args, &UnifiedRecord::single(record));
        assert_eq!(stats.total_bases_read, 10);
        let record2 = UnifiedRead::new_fastq("read2", None, b"012345678", b"IIIIIIIII");
        let (result2, stats2) = process_record_unified(&args, &UnifiedRecord::single(record2));
        assert!(result2.record.reads[0].is_none());
        assert!(result2.record.reads[1].is_none());
        assert_eq!(stats2.total_bases_read, 9);
    }

    #[test]
    fn test_pre_trim_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--pre-trim-left", "5", "--pre-trim-right", "3", "--min-read-length", "0"]);
        let context = Context::empty();
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = UnifiedRead::new_fastq("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let (_, stats) = worker.process_record(&UnifiedRecord::single(record), &mut state);
        assert_eq!(stats.bases_pre_trim, 8);
        assert_eq!(stats.total_bases_read, 15);
    }

    #[test]
    fn test_post_trim_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10", "--min-read-length", "0"]);
        let record = UnifiedRead::new_fastq("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let (result, stats) = process_record_unified(&args, &UnifiedRecord::single(record));
        let processed = match result.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed.seq().len(), 10);
        assert_eq!(processed.seq(), b"0123456789");
        assert_eq!(stats.bases_post_trim, 5);
        assert_eq!(stats.total_bases_read, 15);
        let record2 = UnifiedRead::new_fastq("read2", None, b"01234567", b"IIIIIIII");
        let (result2, stats2) = process_record_unified(&args, &UnifiedRecord::single(record2));
        let processed2 = match result2.record.reads { [Some(r), None] => r, _ => unreachable!() };
        assert_eq!(processed2.seq().len(), 8);
        assert_eq!(stats2.bases_post_trim, 0);
    }

    fn process_with_reference_setup(record: &UnifiedRecord) -> (AnnotatedRecord, WorkerStatistics) {
        let args = Args::parse_from(&[
            "flexbar_rs",
            "--removal-tags",
            "--umi-tags",
            "--min-read-length", "18",
            "--max-uncalled", "1",
            "--qtrim-threshold", "10",
            "--barcode-unassigned",
            "--post-trim-length", "248",
            "--adapter-min-overlap", "6",
            "--qtrim-format", "sanger",
        ]);
        let barcodes = vec![SequenceWithId { id: "SmartSeqUMI".to_string(), sequence: "ATTGCGCAATGNNNNNNNN".to_string() }];
        let adapters = vec![
            SequenceWithId { id: "Trans1".to_string(), sequence: "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans1_rc".to_string(), sequence: "CTGTCTCTTATACACATCTGACGCTGCCGACGA".to_string() },
            SequenceWithId { id: "Trans2".to_string(), sequence: "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".to_string() },
            SequenceWithId { id: "Trans2_rc".to_string(), sequence: "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC".to_string() },
        ];
        let mut context = Context::empty();
        if record.reads[1].is_none() {
            context.barcodes = barcodes.clone();
            context.adapters = adapters.clone();
            context.run_type = crate::types::RunType::Single;
        } else {
            context.barcodes = barcodes.clone();
            context.adapters = adapters.clone();
            context.adapters2 = adapters.clone();
            context.run_type = crate::types::RunType::Paired;
        }
        context.format = FileFormat::Fastq;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_near_miss_with_adapter_trim_paired() {
        let r1 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:17562:1031 1:N:0:TGCCGGTCAG+TTGTATCAGG",
            None,
            b"ATTGCGCAATGGTAGACTGGGGCTAGACTAGGCTGTCTTTGAATAACACTTAAAAATACATTAAAAGGAAAACCCAACTAATAACCTTGTTATTCTCTTCC",
            b"FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
        );
        let r2 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:17562:1031 2:N:0:TGCCGGTCAG+TTGTATCAGG",
            None,
            b"GAATACAGAGATGTGGTCCCCACCATCAAGGGACCAGAAACTTATTGCACAGTAAGTTGAATGTTGTGGAACTTTGATAGGTATTTCCTGAAAAGCTGACA",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFF,FFFFFFFFFF",
        );
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, stats) = process_with_reference_setup(&rec);
        assert_eq!(res.barcode_id, Some("SmartSeqUMI".to_string()));
        let reads = match res.record.reads {
            [Some(x1), Some(x2)] => [x1, x2],
            _ => panic!("expected paired")
        };
        assert_eq!(reads[0].seq(), b"GGGCTAGACTAGGCTGTCTTTGAATAACACTTAAAAATACATTAAAAGGAAAACCCAACTAATAACCTTGTTATTCTCTTCC");
        assert_eq!(stats.total_adapters_removed, 0, "Expected no adapters to be found");
    }

    #[test]
    fn test_short_read_discard_paired() {
        let r1 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:5936:2253 1:N:0:TGCCGGTCAG+GTGTATCAGG",
            None,
            b"ATTGCGCAATGGCTGTCTCTTATACACATCTCCGAGCCCACGAGACTGCCGGTCAGATCTCGTATGCCGGCTTCTGCTTGAAAAGGGGGGGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF,FF::FF,FF,:,F,,F,::F,F,,,:,F:FFFFFFFFFFF!",
        );
        let r2 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:5936:2253 2:N:0:TGCCGGTCAG+GTGTATCAGG",
            None,
            b"CCATTGCGCAATCTGTCTCTTATACACATCTGACGCTGCCGACGATTGTATCAGGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAATGGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFF::,,,FFFFFFFFFF!",
        );
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_with_reference_setup(&rec);
        assert!(res.record.reads[0].is_none(), "Record should be discarded");
        assert!(res.record.reads[1].is_none(), "Record should be discarded");
    }

    #[test]
    fn test_zero_length_after_adapter_removal_paired() {
        let r1 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:16468:4993",
            Some("1:N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"ATTGCGCAATGACAGGACCTGTCTCTTATACACATCTCCGAGCCCACGAGACTGCCGGTCAGATCTCGTATGCCGTCTTCTGCTTGAAAAGGGGGGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,::,,F:,,,F,,:,,F,,F::,,,::F:FF:F",
        );
        let r2 = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1101:16468:4993",
            Some("2:N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"GTCCTGTCATTGCGCAATCTGTCTCTTATACACATCTGACGCTGCCGACGATTGTATCAGGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAGGGGGG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFF,:,,,FF:",
        );
        let rec = UnifiedRecord::paired(r1, r2);
        let (res, _) = process_with_reference_setup(&rec);
        assert!(res.record.reads[0].is_none(), "Record should be discarded");
        assert!(res.record.reads[1].is_none(), "Record should be discarded");
    }

    #[test]
    fn test_adapter_with_errors_single() {
        let r = UnifiedRead::new_fastq(
            "A00620:147:H225CDMXY:2:1221:21866:34444",
            Some("N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"ATTGCGCAATGCTTGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGCTGTCTCTTATACACAACCCGAGCCCACGAGAATGCCGGTCAGATCTCGTATG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,,,,,,,,,,FFF:F,F,F,::F,F,FF:FF,F,F,F,FF:F:,:,FF,,,:,,F:",
        );
        let (res, _) = process_with_reference_setup(&UnifiedRecord::single(r));
        let reads = match res.record.reads {
            [Some(r1), None] => [r1],
            _ => panic!()
        };
        assert_eq!(reads[0].seq(), b"AAAAAAAAAAAAAAAAAAAAAAAAAATCG");
    }

    // Included commented test from single worker
    // #[test]
    // fn test_adapter_with_errors_in_overlap_single() {
    //     let r = SingeRecordDto::with_attrs(
    //         "A00620:147:H225CDMXY:2:1101:17183:24831",
    //         Some("N:0:TGCCGGTCAG+TTGTATCAGG"),
    //         b"GATAAACAGAGAGCAAGGCCTCTCCCAGGAGGTCACCAAAACTGCCCACCTACCTCCTCCCTAACCCACTGGCTTCTGCTTTCCTCTCTCCCTGTTCTTAT",
    //         b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F",
    //     );
    //     let res = process_with_reference_setup(&RecordDto::Single(r));
    //     assert!(res.record.is_some(), "Should keep read");
    //     if let Some(RecordDto::Single(s)) = res.record { assert_eq!(s.seq(), b"GATAAACAGAGAGCAAGGCCTCTCCCAGGAGGTCACCAAAACTGCCCACCTACCTCCTCCCTAACCCACTGGCTTCTGCTTTCCTCTCTCC"); } else { panic!(); }
    // }
}
