use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, AlignerConfig};
use crate::multi_aligner::result::MultiAlignerResult;

#[derive(Clone)]
pub struct SingleAdapterRemoverSetup {
    pub adapters: Vec<String>,
    pub multi_aligner_setup: MultiAlignerSetup,
    pub trim_end: TrimEnd,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

pub struct SingleAdapterRemover {
    multi_aligner: MultiAligner,
}

impl SingleAdapterRemover {
    pub fn remove_adapters(&mut self, setup: &SingleAdapterRemoverSetup, rec: &mut RecordInTransformation<SingeRecordDto>) {
        let current_seq = &rec.get_sequence();
        let best_alignment = self.multi_aligner.align(&setup.multi_aligner_setup, current_seq);
        if let Some(alignment) = best_alignment {
            let adapter = &setup.adapters[alignment.reference_index];
            self.handle_umi_tags(rec, &alignment, adapter, setup);
            self.trim_record(rec, &alignment, setup);
            rec.stats.total_adapters_removed += 1;
        }
    }

    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        alignment: &MultiAlignerResult,
        adapter: &str,
        setup: &SingleAdapterRemoverSetup,
    ) {
        if !setup.umi_tags { return; }
        let tags = alignment.get_umi_tags(rec.get_sequence(), adapter.as_bytes());
        if setup.trace && !tags.is_empty() {
            trace!("Adapter remover on {}: Extracted UMI tags {:?}", rec.record.id(), tags);
        }
        rec.alignment_results.umi_tags.extend(tags);
    }

    fn trim_record(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        alignment: &MultiAlignerResult,
        setup: &SingleAdapterRemoverSetup,
    ) {
        rec.trim(setup.trim_end, alignment.get_start_index(), alignment.get_end_index());
        if setup.removal_tags {
            rec.alignment_results.removal_tags.push(
                format!("Flexbar_removal_{}", alignment.reference_index))
            ;
        }
    }
}

impl GenericStep<SingeRecordDto> for SingleAdapterRemover {
    type Setup = SingleAdapterRemoverSetup;
    
    fn create_setup(args: &Args, context: &Context) -> SingleAdapterRemoverSetup {
        let adapters = context.adapters.clone();
        let references: Vec<SequenceWithId> = adapters.iter().enumerate()
            .map(|(index, seq)|
                SequenceWithId{id: format!("adapter_{}", index), sequence: seq.to_string()}
            )
            .collect();
        let config = AlignerConfig {
            match_score: args.adapter_match,
            mismatch_score: args.adapter_mismatch,
            gap_open: args.adapter_gap_open,
            gap_extend: args.adapter_gap_extend,
            error_rate: args.adapter_error_rate,
            min_overlap: Some(args.adapter_min_overlap),
        };
        let multi_aligner_setup = MultiAlignerSetup::new(
            "Adapter Remover".to_string(),
            args.adapter_trim_end,
            references,
            config,
            args.adapter_kmer
        );
        SingleAdapterRemoverSetup {
            adapters,
            multi_aligner_setup,
            trim_end: args.adapter_trim_end,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn new(setup: &SingleAdapterRemoverSetup) -> Self {
        Self {
            multi_aligner: setup.multi_aligner_setup.create_aligner(),
        }
    }

    fn apply(
        &mut self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &SingleAdapterRemoverSetup
    ) -> RecordInTransformation<SingeRecordDto> {
        if rec.alignment_results.record_dropped { return rec }
        
        if setup.trace {
            trace!("Adapter remover processing record: {}", rec.record.id());
        }
        self.remove_adapters(setup, &mut rec);
        rec
    }
}
