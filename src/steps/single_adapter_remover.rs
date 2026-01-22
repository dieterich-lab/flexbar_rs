use crate::config::Config;
use crate::types::{TrimEnd, RecordDto};
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, trim_sequence};
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
    pub fn new(config: &Config, adapters: &[String]) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match config.args.adapter_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (0, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, 0, 0, 0),
        };

        let max_adapter_len = adapters.iter().map(|a| a.len()).max().unwrap_or(0);

        Self {
            adapters: adapters.to_vec(),
            adapter_trim_end: config.args.adapter_trim_end,
            adapter_min_overlap: config.args.adapter_min_overlap,
            adapter_error_rate: config.args.adapter_error_rate,
            removal_tags: config.args.removal_tags,
            adapter_match: config.args.adapter_match,
            adapter_mismatch: config.args.adapter_mismatch,
            adapter_gap_open: config.args.adapter_gap_open,
            adapter_gap_extend: config.args.adapter_gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            umi_tags: config.args.umi_tags,
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

    pub fn remove_adapters(&self, record: &RecordDto, stats: &mut WorkerStatistics) -> (RecordDto, Vec<String>) {
        let mut aligner = Aligner::with_capacity_and_scoring(
            record.seq().len(),
            self.max_adapter_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, usize)> = None;

        for (i, adapter) in self.adapters.iter().enumerate() {
            let alignment = aligner.semiglobal(adapter.as_bytes(), record.seq());
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
            let umi_tags = if self.umi_tags {
                extract_umi_tags(&alignment, self.adapters[i].as_bytes(), record.seq())
            } else {
                Vec::new()
            };

            let (new_seq, new_qual) = trim_sequence(record, &alignment, self.adapter_trim_end);
            
            stats.total_adapters_removed += 1;
            
            let mut id_str = record.id().to_string();
            if self.removal_tags {
                id_str.push_str(&format!("_Flexbar_removal_{}", i));
            }
            
            (RecordDto::with_attrs(&id_str, record.desc(), &new_seq, &new_qual), umi_tags)
        } else {
            (record.clone(), Vec::new())
        }
    }
}

impl GenericStep<RecordDto> for SingleAdapterRemover {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if let Some(record) = rec.record {
            let (new_record, umi_tags) = self.remove_adapters(&record, &mut rec.stats);
            rec.record = Some(new_record);
            rec.alignment_results.umi_tags.extend(umi_tags);
        }
        rec
    }
}
