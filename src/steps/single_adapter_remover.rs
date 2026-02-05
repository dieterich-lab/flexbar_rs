use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::Alignment;
use crate::utils::aligner::AlignerFactory;

pub struct SingleAdapterRemover {
    adapters: Vec<String>,
    adapter_trim_end: TrimEnd,
    adapter_min_overlap: usize,
    adapter_error_rate: f32,
    removal_tags: bool,
    aligner_factory: AlignerFactory,
    umi_tags: bool,
    max_adapter_len: usize,
    adapter_match: i32,
    adapter_mismatch: i32,
}

impl SingleAdapterRemover {
    pub fn new(args: &Args, adapters: &[String]) -> Self {
        let max_adapter_len = adapters.iter().map(|a| a.len()).max().unwrap_or(0);

        let aligner_factory = AlignerFactory::new(
            args.adapter_trim_end,
            args.adapter_match,
            args.adapter_mismatch,
            args.adapter_gap_open,
            args.adapter_gap_extend,
        );

        Self {
            adapters: adapters.to_vec(),
            adapter_trim_end: args.adapter_trim_end,
            adapter_min_overlap: args.adapter_min_overlap,
            adapter_error_rate: args.adapter_error_rate,
            removal_tags: args.removal_tags,
            aligner_factory,
            umi_tags: args.umi_tags,
            max_adapter_len,
            adapter_match: args.adapter_match,
            adapter_mismatch: args.adapter_mismatch,
        }
    }

    pub fn remove_adapters(&self, rec: &mut RecordInTransformation<SingeRecordDto>) {
        let current_seq = &rec.record.seq()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
        let mut aligner = self.aligner_factory.get(current_seq.len(), self.max_adapter_len);

        let mut best_alignment: Option<(Alignment, usize)> = None;

        for (i, adapter) in self.adapters.iter().enumerate() {
            let alignment = aligner.custom(current_seq, adapter.as_bytes());
            let overlap = alignment.xend - alignment.xstart;
            let adapter_overlap = alignment.yend - alignment.ystart;
            
            if overlap >= self.adapter_min_overlap && adapter_overlap >= self.adapter_min_overlap {
                let matches = (alignment.score - overlap as i32 * self.adapter_mismatch) / (self.adapter_match - self.adapter_mismatch);
                let errors = overlap as i32 - matches;
                let error_rate = errors as f32 / overlap as f32;

                if error_rate <= self.adapter_error_rate {
                    if best_alignment.as_ref().map_or(true, |(best, _)| alignment.score > best.score) {
                        best_alignment = Some((alignment, i));
                    }
                }
            }
        }

        if let Some((alignment, i)) = best_alignment {
            if self.umi_tags {
                rec.alignment_results.umi_tags.extend(extract_umi_tags(&alignment, self.adapters[i].as_bytes(), current_seq));
            }

            let (rel_start, rel_end) = get_trim_range(current_seq.len(), &alignment, self.adapter_trim_end);
            let old_start = rec.alignment_results.seq_start;
            rec.alignment_results.seq_start = old_start + rel_start;
            rec.alignment_results.seq_end = old_start + rel_end;
            
            rec.stats.total_adapters_removed += 1;
            
            if self.removal_tags {
                rec.alignment_results.removal_tags.push(format!("Flexbar_removal_{}", i));
            }
        }
    }
}

impl GenericStep<SingeRecordDto> for SingleAdapterRemover {
    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        if !rec.alignment_results.record_dropped {
            trace!("Adapter remover processing record: {}", rec.record.id());
            self.remove_adapters(&mut rec);
        }
        rec
    }
}
