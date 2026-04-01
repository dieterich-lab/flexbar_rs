use log::trace;
use crate::args::Args;
use crate::types::{SingeRecordDto, TrimEnd};
use crate::context::SequenceWithId;
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, AlignerConfig, MultiAlignerState};
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
