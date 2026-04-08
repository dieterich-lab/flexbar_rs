use log::debug;
use crate::args::Args;
use crate::context::Context;
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation};
use crate::steps::{pre_trim::PreTrim, uncalled_filter::UncalledFilter, adapter_remover::AdapterRemoverStep, barcode_splitter::BarcodeSplitter, quality_trim::QualityTrimStep, filter_length::FilterLength, post_trim::PostTrim};
use crate::steps::pre_trim::PreTrimState;
use crate::steps::filter_length::FilterLengthState;
use crate::steps::post_trim::PostTrimState;
use crate::steps::barcode_splitter::BarcodeSplitterState;
use crate::steps::adapter_remover::AdapterRemoverStepState;
use crate::types::{RecordDto, SingeRecordDto, PairedRecordDto, QualTrimType};

pub struct Worker {
    pub pre_trim: Option<PreTrim>,
    pub barcode_splitter: Option<BarcodeSplitter>,
    pub adapter_remover: Option<AdapterRemoverStep>,
    pub quality_trim: Option<QualityTrimStep>,
    pub filter_length: FilterLength,
    pub post_trim: PostTrim,
    pub uncalled_filter: Option<UncalledFilter>,
    pub qtrim_post_removal: bool,
}

pub struct WorkerState {
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
        let pre_trim = if args.pre_trim_left > 0 || args.pre_trim_right > 0 {
            Some(PreTrim::new(args, context))
        } else { None };

        let uncalled_filter = if args.max_uncalled > 0 {
            Some(UncalledFilter::new(args, context))
        } else { None };

        let quality_trim = match args.qtrim {
            QualTrimType::Off => None,
            _ => Some(QualityTrimStep::new(args, context)),
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
        if args.trace {
            debug!("Unified worker pipeline:");
            if let Some(s) = &pre_trim {
                debug!("- Step PreTrim (left: {}, right: {})", s.pre_trim_left, s.pre_trim_right);
            }
            if let Some(s) = &uncalled_filter {
                debug!("- Step UncalledFilter (max_uncalled: {})", s.max_uncalled);
            }
            if !qtrim_post_removal { if let Some(s) = &quality_trim {
                debug!("- Step QualityTrim (threshold: {}, window size: {})", s.inner.threshold, s.inner.win_size);
            }}
            if let Some(s) = &barcode_splitter {
                debug!("- Step BarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})", s.trim_end, s.barcode_keep, s.barcode_unassigned);
            }
            if adapter_remover.is_some() {
                debug!("- Step AdapterRemoverStep");
            }
            if qtrim_post_removal { if let Some(s) = &quality_trim {
                debug!("- Step QualityTrim (threshold: {}, window size: {})", s.inner.threshold, s.inner.win_size);
            }}
            debug!("- Step FilterLength (min length: {})", filter_length.min_read_length);
            if post_trim.post_trim_length > 0 {
                debug!("- Step PostTrim (length: {})", post_trim.post_trim_length);
            }
        }

        Self { pre_trim, barcode_splitter, adapter_remover, quality_trim, filter_length, post_trim, uncalled_filter, qtrim_post_removal }
    }

    pub fn create_state(&self) -> WorkerState {
        WorkerState {
            pre_trim: self.pre_trim.as_ref().map(|s| s.create_state()),
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            adapter_remover: self.adapter_remover.as_ref().map(|s| s.create_state()),
            quality_trim: self.quality_trim.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
            post_trim: self.post_trim.create_state(),
            uncalled_filter: self.uncalled_filter.as_ref().map(|s| s.create_state()),
        }
    }

    pub fn process_record(&self, record: &RecordDto, state: &mut WorkerState) -> FinalRecord {
        let mut rec_in = record.get_record_in_transformation();
        match record {
            RecordDto::Single(r) => { rec_in.stats.total_bases_read = r.seq().len(); },
            RecordDto::Paired(p) => { rec_in.stats.total_bases_read = p.reads[0].seq().len() + p.reads[1].seq().len(); },
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
        finalize(rec_in)
    }
}

fn finalize_read_stats(stats: &mut crate::statistics::WorkerStatistics, per_read: &crate::steps::generic::PerReadAlignmentResults, original_len: usize) {
    if per_read.dropped {
        stats.bases_discarded += original_len;
    } else {
        stats.bases_quality_trim += per_read.bases_quality_trim;
        stats.bases_adapter_removal += per_read.bases_adapter_removal;
        stats.bases_barcode_removal += per_read.bases_barcode_removal;
        stats.bases_post_trim += per_read.bases_post_trim;
        stats.bases_pre_trim += per_read.bases_pre_trim;
    }
}

fn finalize(mut rec: RecordInTransformation) -> FinalRecord {
    match &rec.record {
        RecordDto::Single(r) => {
            let original_len = r.seq().len();
            let per = rec.alignment_results.per_read[0].as_ref().unwrap();
            finalize_read_stats(&mut rec.stats, per, original_len);
            if per.dropped {
                return FinalRecord { record: None, output_tag: None, stats: rec.stats };
            }
            let mut desc_parts = Vec::new();
            if let Some(desc) = r.desc() { desc_parts.push(desc); }
            for tag in &per.removal_tags { desc_parts.push(tag); }
            for tag in &per.umi_tags { desc_parts.push(tag); }
            let new_desc: Option<String> = if desc_parts.is_empty() { None } else { Some(desc_parts.join("_")) };
            let final_record = SingeRecordDto::with_attrs(r.id(), new_desc.as_deref(), &r.seq()[per.read_range.clone()], &r.qual()[per.read_range.clone()]);
            FinalRecord { record: Some(RecordDto::Single(final_record)), output_tag: rec.alignment_results.barcode_tag, stats: rec.stats }
        }
        RecordDto::Paired(p) => {
            for i in 0..2 {
                let original_len = p.reads[i].seq().len();
                finalize_read_stats(&mut rec.stats, rec.alignment_results.per_read[i].as_ref().unwrap(), original_len);
            }
            if rec.alignment_results.per_read[0].as_ref().unwrap().dropped && rec.alignment_results.per_read[1].as_ref().unwrap().dropped {
                return FinalRecord { record: None, output_tag: None, stats: rec.stats };
            }
            let mut records = p.reads.clone();
            let umi_tags = if rec.alignment_results.per_read[0].as_ref().unwrap().umi_tags.is_empty() {
                rec.alignment_results.per_read[1].as_ref().unwrap().umi_tags.clone()
            } else { rec.alignment_results.per_read[0].as_ref().unwrap().umi_tags.clone() };
            let mut good_reads: u32 = 0;
            for i in 0..2 {
                let results = rec.alignment_results.per_read[i].as_ref().unwrap();
                let old_rec = &p.reads[i];
                if results.dropped {
                    records[i] = SingeRecordDto::with_attrs(old_rec.id(), None, &[], &[]);
                    continue;
                }
                good_reads += 1;
                let range = &results.read_range;
                let mut desc_parts = Vec::new();
                if let Some(desc) = old_rec.desc() { desc_parts.push(desc); }
                for tag in &results.removal_tags { desc_parts.push(tag); }
                for tag in umi_tags.iter() { desc_parts.push(tag); }
                if rec.alignment_results.per_read[1 - i].as_ref().unwrap().dropped { desc_parts.push("single"); }
                let new_desc = if desc_parts.is_empty() { None } else { Some(desc_parts.join("_")) };
                records[i] = SingeRecordDto::with_attrs(old_rec.id(), new_desc.as_deref(), &old_rec.seq()[range.clone()], &old_rec.qual()[range.clone()]);
            }
            let record = if good_reads > 0 { Some(RecordDto::Paired(PairedRecordDto { reads: records })) } else { None };
            FinalRecord { record, output_tag: rec.alignment_results.barcode_tag, stats: rec.stats }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use crate::context::{Context, SequenceWithId};
    use crate::types::{SingeRecordDto, PairedRecordDto, TrimEnd, FileFormat};

    fn process_record_unified(args: &Args, record: &RecordDto) -> FinalRecord {
        let context = Context::new(args);
        let worker = Worker::new(args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    fn process_record_paired_helper(args: Args, record: &PairedRecordDto, barcodes: Vec<SequenceWithId>, barcodes2: Vec<SequenceWithId>) -> FinalRecord {
        let mut context = Context::empty();
        context.barcodes = barcodes;
        context.barcodes2 = barcodes2;
        context.run_type = crate::types::RunType::Paired;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let rec = RecordDto::Paired(record.clone());
        worker.process_record(&rec, &mut state)
    }

    // Copied and adapted from src/workers/paired.rs
    #[test]
    fn test_barcode_splitting_better_match_first_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAAGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.output_tag, Some("barcode_B1".to_string()));
        let records = match res.record.unwrap() { RecordDto::Paired(p) => p.reads, _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"GGGGG");
        assert_eq!(records[1].seq(), b"TTTTTGGGGG");
    }

    #[test]
    fn test_barcode_splitting_better_match_second_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "5", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B2".to_string(), sequence: "CCCCC".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"TTTTTGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"CCCCCGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.output_tag, Some("barcode_B2".to_string()));
        let records = match res.record.unwrap() { RecordDto::Paired(p) => p.reads, _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"TTTTTGGGGG");
        assert_eq!(records[1].seq(), b"GGGGG");
    }

    #[test]
    fn test_barcode_splitting_no_match_unassigned_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--barcode-unassigned", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.output_tag, Some("barcode_unassigned".to_string()));
        assert!(res.record.is_some());
    }

    #[test]
    fn test_barcode_splitting_no_match_no_unassigned_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--min-read-length", "0"]);
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let barcodes2 = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"GGGGGGGGGG", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, barcodes2);
        assert_eq!(res.output_tag, None);
        assert!(res.record.is_none());
    }

    #[test]
    fn test_length_filter_drop_default_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, vec![]);
        assert!(res.record.is_none());
    }

    #[test]
    fn test_length_filter_single_reads_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads = true;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, vec![]);
        assert!(res.record.is_some());
        let records = match res.record.unwrap() { RecordDto::Paired(p) => p.reads, _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"");
        assert_eq!(records[1].seq(), b"0123456789");
    }

    #[test]
    fn test_length_filter_single_reads_paired_tagged_paired() {
        let mut args = Args::default();
        args.barcode_trim_end = TrimEnd::Left;
        args.barcode_min_overlap = Some(5);
        args.min_read_length = 10;
        args.single_reads_paired = true;
        let barcodes = vec![SequenceWithId { id: "B1".to_string(), sequence: "AAAAA".to_string() }];
        let r1 = SingeRecordDto::with_attrs("R1", None, b"AAAAA987654321", b"IIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("R1", None, b"0123456789", b"IIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);
        let res = process_record_paired_helper(args, &rec, barcodes, vec![]);
        assert!(res.record.is_some());
        let records = match res.record.unwrap() { RecordDto::Paired(p) => p.reads, _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"N");
        assert_eq!(records[1].seq(), b"0123456789");
    }

    #[test]
    fn test_quality_trimming_paired() {
        let mut args = Args::default();
        args.qtrim = QualTrimType::Tail;
        args.qtrim_threshold = 20;
        args.min_read_length = 0;
        let r1 = SingeRecordDto::with_attrs("R1", None, b"5555544444", b"8888877777");
        let r2 = SingeRecordDto::with_attrs("R2", None, b"5555544444", b"8888833333");
        let rec = PairedRecordDto::new(r1, r2);
        let mut context = Context::empty();
        context.format = FileFormat::Fastq;
        context.run_type = crate::types::RunType::Paired;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let res = worker.process_record(&RecordDto::Paired(rec), &mut state);
        assert!(res.record.is_some());
        let records = match res.record.unwrap() { RecordDto::Paired(p) => p.reads, _ => panic!("expected paired") };
        assert_eq!(records[0].seq(), b"5555544444");
        assert_eq!(records[1].seq(), b"55555");
    }

    // Single worker tests adapted
    #[test]
    fn test_process_right_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_left_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(result.output_tag, Some("barcode_bc1".to_string()));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_max_uncalled_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--max-uncalled", "2"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACNNGT", b"IIIIII");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        assert_eq!(result.stats.total_discarded_uncalled, 0);
        let record2 = SingeRecordDto::with_attrs("read2", None, b"ACNNGN", b"IIIIII");
        let result2 = process_record_unified(&args, &RecordDto::Single(record2));
        assert_eq!(result2.stats.total_discarded_uncalled, 1);
    }

    #[test]
    fn test_removal_tags_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGT", b"IIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal"));
    }

    #[test]
    fn test_removal_tags_adapter_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5", "--removal-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "ad1".to_string(), sequence: "AGATCGGAAGAG".to_string() }];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_ad1"));
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
        let record = SingeRecordDto::with_attrs("read1", None, b"GGGGACGTCTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("Flexbar_removal_Flexbar_removal_ad1"));
    }

    #[test]
    fn test_quality_trim_tail_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "TAIL", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_quality_trim_win_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "WIN", "--qtrim-threshold", "20", "--qtrim-win-size", "5", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGTACGT", b"IIIII44444III");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_quality_trim_bwa_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--qtrim", "BWA", "--qtrim-threshold", "20", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTACGT", b"IIIII4444");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_umi_tags_barcode_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.barcodes = vec![SequenceWithId{id: "bc1".to_string(), sequence: "NNGGGGNN".to_string()}];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"AAGGGGTTACGT", b"IIIIIIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(result.output_tag, Some("barcode_bc1".to_string()));
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("AA_TT"));
        assert_eq!(processed.seq(), b"ACGT");
    }

    #[test]
    fn test_umi_tags_adapter_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATNNNN", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "4", "--umi-tags", "--min-read-length", "0"]);
        let mut context = Context::new(&args);
        context.adapters = vec![SequenceWithId { id: "adapter".to_string(), sequence: "AGATNNNN".to_string() }];
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"ACTGTAGATCCCC", b"IIIIIIIIIIIII");
        let result = worker.process_record(&RecordDto::Single(record), &mut state);
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.id(), "read1");
        assert_eq!(processed.desc(), Some("CCCC"));
        assert_eq!(processed.seq(), b"ACTGT");
    }

    #[test]
    fn test_min_read_length_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789", b"IIIIIIIIII");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        assert!(result.record.is_some());
        assert_eq!(result.stats.total_bases_read, 10);
        let record2 = SingeRecordDto::with_attrs("read2", None, b"012345678", b"IIIIIIIII");
        let result2 = process_record_unified(&args, &RecordDto::Single(record2));
        assert!(result2.record.is_none());
        assert_eq!(result2.stats.total_bases_read, 9);
    }

    #[test]
    fn test_pre_trim_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--pre-trim-left", "5", "--pre-trim-right", "3", "--min-read-length", "0"]);
        let context = Context::empty();
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let res = worker.process_record(&RecordDto::Single(record), &mut state);
        assert!(res.record.is_some());
        assert_eq!(res.stats.bases_pre_trim, 8);
        assert_eq!(res.stats.total_bases_read, 15);
    }

    #[test]
    fn test_post_trim_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10", "--min-read-length", "0"]);
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let result = process_record_unified(&args, &RecordDto::Single(record));
        let processed = match result.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed.seq().len(), 10);
        assert_eq!(processed.seq(), b"0123456789");
        assert_eq!(result.stats.bases_post_trim, 5);
        assert_eq!(result.stats.total_bases_read, 15);
        let record2 = SingeRecordDto::with_attrs("read2", None, b"01234567", b"IIIIIIII");
        let result2 = process_record_unified(&args, &RecordDto::Single(record2));
        let processed2 = match result2.record.expect("kept") { RecordDto::Single(r) => r, _ => unreachable!() };
        assert_eq!(processed2.seq().len(), 8);
        assert_eq!(result2.stats.bases_post_trim, 0);
    }

    fn process_with_reference_setup(record: &RecordDto) -> FinalRecord {
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
        match record {
            RecordDto::Single(_) => {
                context.barcodes = barcodes.clone();
                context.adapters = adapters.clone();
                context.run_type = crate::types::RunType::Single;
            }
            RecordDto::Paired(_) => {
                context.barcodes = barcodes.clone();
                context.adapters = adapters.clone();
                context.adapters2 = adapters.clone();
                context.run_type = crate::types::RunType::Paired;
            }
        }
        context.format = FileFormat::Fastq;
        let worker = Worker::new(&args, &context);
        let mut state = worker.create_state();
        worker.process_record(record, &mut state)
    }

    #[test]
    fn test_near_miss_with_adapter_trim_paired() {
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
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));
        let res = process_with_reference_setup(&rec);
        assert_eq!(res.output_tag, Some("barcode_SmartSeqUMI".to_string()));
        if let Some(RecordDto::Paired(p)) = res.record {
            assert_eq!(p.reads[0].seq(), b"GGGCTAGACTAGGCTGTCTTTGAATAACACTTAAAAATACATTAAAAGGAAAACCCAACTAATAACCTTGTTATTCTCTTCC");
        } else { panic!("expected paired"); }
        let stats = res.stats;
        assert_eq!(stats.total_adapters_removed, 0, "Expected no adapters to be found");
    }

    #[test]
    fn test_short_read_discard_paired() {
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
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));
        let res = process_with_reference_setup(&rec);
        assert!(res.record.is_none(), "Record should be discarded");
    }

    #[test]
    fn test_zero_length_after_adapter_removal_paired() {
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
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));
        let res = process_with_reference_setup(&rec);
        assert!(res.record.is_none(), "Record should be discarded");
    }

    #[test]
    fn test_adapter_with_errors_single() {
        let r = SingeRecordDto::with_attrs(
            "A00620:147:H225CDMXY:2:1221:21866:34444",
            Some("N:0:TGCCGGTCAG+TTGTATCAGG"),
            b"ATTGCGCAATGCTTGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGCTGTCTCTTATACACAACCCGAGCCCACGAGAATGCCGGTCAGATCTCGTATG",
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,,,,,,,,,,FFF:F,F,F,::F,F,FF:FF,F,F,F,FF:F:,:,FF,,,:,,F:",
        );
        let res = process_with_reference_setup(&RecordDto::Single(r));
        assert!(res.record.is_some(), "Should keep read");
        if let Some(RecordDto::Single(s)) = res.record { assert_eq!(s.seq(), b"AAAAAAAAAAAAAAAAAAAAAAAAAATCG"); } else { panic!(); }
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
