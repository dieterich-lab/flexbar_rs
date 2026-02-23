use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::SequenceWithId;
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup};
use crate::multi_aligner::result::MultiAlignerResult;

use std::any::Any;

pub struct SingleAdapterRemover {
    adapters: Vec<String>,
    multi_aligner: MultiAligner,
    adapter_trim_end: TrimEnd,
    removal_tags: bool,
    umi_tags: bool,
    trace: bool,
}

pub struct SingleAdapterRemoverSetup {
    pub multi_aligner_setup: MultiAlignerSetup,
}

impl SingleAdapterRemover {
    pub fn new(args: &Args, adapters: &[String]) -> Self {
        let references: Vec<SequenceWithId> = adapters.iter().enumerate()
                .map(|(index, seq)|
                    SequenceWithId{id: format!("adapter_{}", index), sequence: seq.to_string()}
                )
                .collect();
        let multi_aligner = MultiAligner::new(
            "Adapter Remover".to_string(),
            references,
            args.adapter_trim_end,
            args.adapter_error_rate,
            Some(args.adapter_min_overlap),
            args.adapter_match,
            args.adapter_mismatch,
            args.adapter_gap_open,
            args.adapter_gap_extend,
            args.adapter_kmer,
        );

        Self {
            adapters: adapters.to_vec(),
            multi_aligner,
            adapter_trim_end: args.adapter_trim_end,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    pub fn remove_adapters(&self, setup: &mut SingleAdapterRemoverSetup, rec: &mut RecordInTransformation<SingeRecordDto>) {
        let current_seq = &rec.get_sequence();
        let best_alignment = self.multi_aligner.align(&mut setup.multi_aligner_setup, current_seq);
        if let Some(alignment) = best_alignment {
            let adapter = &self.adapters[alignment.reference_index];
            self.handle_umi_tags(rec, &alignment, adapter);
            self.trim_record(rec, &alignment);
            rec.stats.total_adapters_removed += 1;
        }
    }

    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        alignment: &MultiAlignerResult,
        adapter: &str
    ) {
        if !self.umi_tags { return; }
        let tags = alignment.get_umi_tags(rec.get_sequence(), adapter.as_bytes());
        if self.trace && !tags.is_empty() {
            trace!("Adapter remover on {}: Extracted UMI tags {:?}", rec.record.id(), tags);
        }
        rec.alignment_results.umi_tags.extend(tags);
    }

    fn trim_record(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        alignment: &MultiAlignerResult
    ) {
        rec.trim(self.adapter_trim_end, alignment.get_start_index(), alignment.get_end_index());
        if self.removal_tags {
            rec.alignment_results.removal_tags.push(
                format!("Flexbar_removal_{}", alignment.reference_index))
            ;
        }
    }
}

impl GenericStep<SingeRecordDto> for SingleAdapterRemover {
    type Setup = Box<dyn Any + Send + Sync>;
    fn create_setup(&self) -> Box<dyn Any + Send + Sync> {
        Box::new(SingleAdapterRemoverSetup {
            multi_aligner_setup: self.multi_aligner.create_setup(),
        })
    }

    fn clone_setup(&self, _setup: &Box<dyn Any + Send + Sync>) -> Box<dyn Any + Send + Sync> {
        self.create_setup()
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &mut Box<dyn Any + Send + Sync>) -> RecordInTransformation<SingeRecordDto> {
        if rec.alignment_results.record_dropped { return rec }
        let setup = setup.downcast_mut::<SingleAdapterRemoverSetup>().unwrap();
        if self.trace {
            trace!("Adapter remover processing record: {}", rec.record.id());
        }
        self.remove_adapters(setup, &mut rec);
        rec
    }
}
