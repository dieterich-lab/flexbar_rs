use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, SingleAlignmentResults};
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, AlignerConfig, MultiAlignerState};
use crate::multi_aligner::result::MultiAlignerResult;

pub struct SingleAdapterRemover {
    pub adapters: Vec<SequenceWithId>,
    pub multi_aligner: MultiAligner,
    pub multi_aligner_setup: MultiAlignerSetup,
    pub trim_end: TrimEnd,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

pub struct SingleAdapterRemoverState {
    pub multi_aligner_state: MultiAlignerState,
}

impl SingleAdapterRemover {
    pub fn remove_adapters(&self, state: &mut SingleAdapterRemoverState, rec: &mut RecordInTransformation<SingeRecordDto, SingleAlignmentResults>) {
        let current_seq = &rec.get_sequence();
        let best_alignment = self.multi_aligner.align(&mut state.multi_aligner_state, &self.multi_aligner_setup, current_seq);
        if let Some(alignment) = best_alignment {
            let adapter = &self.adapters[alignment.reference_index].sequence;
            self.handle_umi_tags(rec, &alignment, adapter);
            self.trim_record(rec, &alignment);
            rec.stats.total_adapters_removed += 1;
        }
    }

    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto, SingleAlignmentResults>,
        alignment: &MultiAlignerResult,
        adapter: &str,
    ) {
        if !self.umi_tags { return; }
        let tags = alignment.get_umi_tags(rec.get_sequence(), adapter.as_bytes());
        if self.trace && !tags.is_empty() {
            trace!("Adapter remover on {}: Extracted UMI tags {:?}", rec.record.id(), tags);
        }
        rec.alignment_results.read.umi_tags.extend(tags);
    }

    fn trim_record(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto, SingleAlignmentResults>,
        alignment: &MultiAlignerResult,
    ) {
        let old_len = rec.alignment_results.read.read_range.end - rec.alignment_results.read.read_range.start;
        rec.trim(self.trim_end, alignment.get_start_index(), alignment.get_end_index());
        let new_len = rec.alignment_results.read.read_range.end - rec.alignment_results.read.read_range.start;
        rec.stats.bases_adapter_removal += old_len - new_len;
        if self.removal_tags {
            rec.alignment_results.read.removal_tags.push(
                format!("Flexbar_removal_{}", alignment.reference_index))
            ;
        }
    }
}

impl GenericStep<SingeRecordDto> for SingleAdapterRemover {
    type State = SingleAdapterRemoverState;
    type AlignmentResult = SingleAlignmentResults;
    
    fn new(args: &Args, context: &Context) -> Self {
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
            &context.adapters,
            config,
            args.adapter_kmer,
        );
        let multi_aligner = multi_aligner_setup.create_aligner();
        Self {
            adapters: context.adapters.clone(),
            multi_aligner,
            multi_aligner_setup,
            trim_end: args.adapter_trim_end,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn create_state(&self) -> Self::State {
        Self::State {
            multi_aligner_state: self.multi_aligner_setup.create_state(&self.multi_aligner),
        }
    }

    fn apply(
        &self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, state: &mut Self::State
    ) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        if rec.alignment_results.record_dropped { return rec }
        
        if self.trace {
            trace!("Adapter remover processing record: {}", rec.record.id());
        }
        self.remove_adapters(state, &mut rec);
        rec
    }
}
