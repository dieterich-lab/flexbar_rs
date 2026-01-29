use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;

pub struct SingleAdapterRemover {
    adapters: Vec<String>,
    adapter_trim_end: TrimEnd,
    adapter_min_overlap: usize,
    adapter_error_rate: f32,
    removal_tags: bool,
    adapter_match: i32,
    adapter_mismatch: i32,
    adapter_gap_open: i32,
    adapter_gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    umi_tags: bool,
    max_adapter_len: usize,
}

impl SingleAdapterRemover {
    pub fn new(args: &Args, adapters: &[String]) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match args.adapter_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (MIN_SCORE, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, MIN_SCORE, 0, 0),
        };

        let max_adapter_len = adapters.iter().map(|a| a.len()).max().unwrap_or(0);

        Self {
            adapters: adapters.to_vec(),
            adapter_trim_end: args.adapter_trim_end,
            adapter_min_overlap: args.adapter_min_overlap,
            adapter_error_rate: args.adapter_error_rate,
            removal_tags: args.removal_tags,
            adapter_match: args.adapter_match,
            adapter_mismatch: args.adapter_mismatch,
            adapter_gap_open: args.adapter_gap_open,
            adapter_gap_extend: args.adapter_gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            umi_tags: args.umi_tags,
            max_adapter_len,
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.adapter_match;
        let mm = self.adapter_mismatch;
        Scoring {
            match_fn: move |a, b| if a == b || a == b'N' { m } else { mm },
            match_scores: None,
            gap_open: self.adapter_gap_open,
            gap_extend: self.adapter_gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn remove_adapters(&self, rec: &mut RecordInTransformation<SingeRecordDto>) {
        let current_seq = &rec.record.seq()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
        let mut aligner = Aligner::with_capacity_and_scoring(
            current_seq.len(),
            self.max_adapter_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, usize)> = None;

        for (i, adapter) in self.adapters.iter().enumerate() {
            let alignment = aligner.semiglobal(adapter.as_bytes(), current_seq);
            let overlap = alignment.yend - alignment.ystart;
            let adapter_overlap = alignment.xend - alignment.xstart;
            
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
    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        if !rec.alignment_results.record_dropped {
            trace!("Adapter remover processing record: {}", rec.record.id());
            self.remove_adapters(&mut rec);
        }
        rec
    }
}
