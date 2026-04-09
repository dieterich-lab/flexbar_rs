use crate::args::Args;
use crate::types::{AdapterReadset, RunType};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use crate::multi_aligner::MultiAlignerState;
use crate::types::{SingeRecordDto, TrimEnd};
use crate::context::SequenceWithId;
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, AlignerConfig};
use crate::multi_aligner::result::MultiAlignerResult;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::PerReadAlignmentResults;

pub struct AdapterRemover {
    pub adapters: Vec<SequenceWithId>,
    pub multi_aligner: MultiAligner,
    pub multi_aligner_setup: MultiAlignerSetup,
    pub trim_end: TrimEnd,
    pub adapter_cycles: usize,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

impl AdapterRemover {
    pub fn new(args: &Args, adapters: &[SequenceWithId]) -> Self {
        let config = AlignerConfig {
            match_score: args.adapter_match,
            mismatch_score: args.adapter_mismatch,
            gap_open: args.adapter_gap_open,
            gap_extend: args.adapter_gap_extend,
            error_rate: args.adapter_error_rate,
            min_overlap: Some(args.adapter_min_overlap),
            cycles: args.adapter_cycles,
        };
        let multi_aligner_setup = MultiAlignerSetup::new(
            "Adapter Remover".to_string(),
            args.adapter_trim_end,
            adapters,
            config,
            args.adapter_kmer,
        );
        let multi_aligner = multi_aligner_setup.create_aligner();
        Self {
            adapters: adapters.to_vec(),
            multi_aligner,
            multi_aligner_setup,
            trim_end: args.adapter_trim_end,
            adapter_cycles: args.adapter_cycles,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    pub fn create_state(&self) -> MultiAlignerState {
        self.multi_aligner_setup.create_state(&self.multi_aligner)
    }

    pub fn remove_adapters(
        &self,
        state: &mut MultiAlignerState,
        read: &SingeRecordDto,
        read_index: Option<usize>,
        per_read: &mut PerReadAlignmentResults,
        stats: &mut WorkerStatistics,
    ) {
        let mut cycle_count = 0;
        loop {
            let current_seq = &read.seq()[per_read.read_range.clone()];
            let best_alignment = self.multi_aligner.align(state, &self.multi_aligner_setup, current_seq);
            if let Some(alignment) = best_alignment {
                let read_seq = &read.seq()[per_read.read_range.clone()];
                self.trace_alignment(read.id(), read_index, read_seq, &alignment);
                let adapter = &self.adapters[alignment.reference_index].sequence;
                self.handle_umi_tags(read, read_index, per_read, &alignment, adapter);
                self.trim_record(read, read_index, per_read, stats, &alignment);
                stats.total_adapters_removed += 1;
                cycle_count += 1;
                if self.adapter_cycles > 0 && cycle_count >= self.adapter_cycles {
                    break;
                }
            } else {
                break;
            }
        }
    }

    fn trace_alignment(
        &self,
        read_id: &str,
        read_index: Option<usize>,
        read_seq: &[u8],
        alignment: &MultiAlignerResult
    ) {
        if self.trace {
            let adapter_seq = self.adapters[alignment.reference_index].sequence.as_bytes();
            let read_str = match read_index {
                Some(i) => format!(" with read {}", i + 1),
                None => "".to_string(),
            };
            trace!(
                    "Adapter remover on {} selected the following alignment{}:",
                    read_id, read_str
                );
            trace!("{}", alignment.format(read_seq, adapter_seq));
        }
    }

    fn handle_umi_tags(
        &self,
        read: &SingeRecordDto,
        read_index: Option<usize>,
        per_read: &mut PerReadAlignmentResults,
        alignment: &MultiAlignerResult,
        adapter: &str,
    ) {
        if !self.umi_tags { return; }
        let current_seq = &read.seq()[per_read.read_range.clone()];
        let tags = alignment.get_umi_tags(current_seq, adapter.as_bytes());
        if self.trace && !tags.is_empty() {
            let read_index_str = match read_index {
                Some(i) => format!(" (R{})", i + 1),
                None => "".to_string(),
            };
            trace!("Adapter remover on {}{}: Extracted UMI tags {:?}", read.id(), read_index_str, tags);
        }
        per_read.umi_tags.extend(tags);
    }

    fn trim_record(
        &self,
        _read: &SingeRecordDto,
        _read_index: Option<usize>,
        per_read: &mut PerReadAlignmentResults,
        _stats: &mut WorkerStatistics,
        alignment: &MultiAlignerResult,
    ) {
        let old_len = per_read.read_range.end - per_read.read_range.start;
        
        let start = alignment.get_start_index();
        let end = alignment.get_end_index();

        match self.trim_end {
            TrimEnd::Right | TrimEnd::RTail => {
                per_read.read_range.end = per_read.read_range.start + start;
            }
            TrimEnd::Left | TrimEnd::LTail => {
                per_read.read_range.start += end;
            }
            TrimEnd::Any => {
                let current_len = per_read.read_range.end - per_read.read_range.start;
                let right_length = current_len - end;
                if start > right_length {
                    per_read.read_range.end = per_read.read_range.start + start;
                } else {
                    per_read.read_range.start += end;
                }
            }
        }

        let new_len = per_read.read_range.end - per_read.read_range.start;
        per_read.bases_adapter_removal += old_len - new_len;
        if self.removal_tags {
            let id = &self.adapters[alignment.reference_index].id;
            per_read.removal_tags.push(
                format!("Flexbar_removal_{}", id))
            ;
        }
    }
}

pub struct AdapterRemoverStep {
    pub inner: [Option<AdapterRemover>; 2],
}

pub struct AdapterRemoverStepState {
    pub states: [Option<MultiAlignerState>; 2],
}

impl GenericStep for AdapterRemoverStep {
    type State = AdapterRemoverStepState;

    fn new(args: &Args, context: &Context) -> Self {
        let inner  = match context.run_type {
            RunType::Single => {
                [Some(AdapterRemover::new(args, &context.adapters)), None]
            },
            RunType::Paired => {
                match args.adapter_read_set {
                    AdapterReadset::One => [
                        Some(AdapterRemover::new(args, &context.adapters)),
                        None,
                    ],
                    AdapterReadset::Two => [
                        None,
                        Some(AdapterRemover::new(args, &context.adapters2)),
                    ],
                    AdapterReadset::None => [
                        Some(AdapterRemover::new(args, &context.adapters)),
                        Some(AdapterRemover::new(args, &context.adapters2)),
                    ],
                }
            }
            _ => {
                panic!("Barcoded run types not implemented yet");
            }
        };

        Self { inner }
    }

    fn create_state(&self) -> Self::State {
        Self::State {
            states: [
                self.inner[0].as_ref().map(|s| s.create_state()),
                self.inner[1].as_ref().map(|s| s.create_state()),
            ],
        }
    }

    fn apply(&self, mut rec: RecordInTransformation, state: &mut Self::State) -> RecordInTransformation {
        for i in 0..2 {
            if let Some(remover) = &self.inner[i] {
                if rec.alignment_results.per_read[i].as_ref().unwrap().dropped { continue; }
                if remover.trace {
                    let read_id = rec.read_id(i);
                    trace!("Adapter remover processing read {}", read_id);
                }
                remover.remove_adapters(
                    state.states[i].as_mut().unwrap(),
                    &rec.record.reads[i].as_ref().unwrap(),
                    Some(i),
                    rec.alignment_results.per_read[i].as_mut().unwrap(),
                    &mut rec.stats
                );
            }
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use crate::context::{Context, SequenceWithId};
    use crate::types::{RecordDto, SingeRecordDto};
    use crate::steps::generic::GenericStep;
    use clap::Parser;

    #[test]
    fn test_adapter_removal_paired() {
        let args = Args::parse_from(["flexbar", "--adapter-trim-end", "Right"]);
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.adapters.push(SequenceWithId {
            id: "ad1".to_string(),
            sequence: "AAAA".to_string(),
        });
        context.adapters2.push(SequenceWithId {
            id: "ad2".to_string(),
            sequence: "CCCC".to_string(),
        });

        let step = AdapterRemoverStep::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGCCCC", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);

        let rec_in = RecordInTransformation::from_dto(rec);
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.stats.total_adapters_removed, 2);
    }

    #[test]
    fn test_adapter_removal_paired_read_set_one() {
        let args = Args::parse_from(["flexbar", "--adapter-trim-end", "Right", "--adapter-read-set", "1"]);

        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        context.adapters.push(SequenceWithId {
            id: "ad1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = AdapterRemoverStep::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGAAAA", b"IIIIIIIIIIII");
        let rec = RecordDto::paired(r1, r2);

        let rec_in = RecordInTransformation::from_dto(rec);
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..12);
        assert_eq!(result.stats.total_adapters_removed, 1);
    }

    #[test]
    fn test_adapter_removal_single() {
        let args = Args::parse_from(["flexbar", "--adapter-trim-end", "Right"]);

        let mut context = Context::empty();
        context.run_type = RunType::Single;
        context.adapters.push(SequenceWithId {
            id: "ad1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = AdapterRemoverStep::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let rec = RecordDto::single(r1);

        let rec_in = RecordInTransformation::from_dto(rec);
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.stats.total_adapters_removed, 1);
    }
}
