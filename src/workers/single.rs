use crate::context::Context;
use crate::config::Config;
use crate::types::TrimEnd;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::io::fastq::Record;

pub struct SingleWorker<'a> {
    config: &'a Config,
    context: &'a Context,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    max_adapter_len: usize,
}

impl<'a> SingleWorker<'a> {
    pub fn new(config: &'a Config, context: &'a Context) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match config.args.adapter_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (0, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, 0, 0, 0),
        };

        let max_adapter_len = context.adapters.iter().map(|a| a.len()).max().unwrap_or(0);

        Self {
            config,
            context,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            max_adapter_len,
        }
    }

    fn get_scoring(&self) -> Scoring<fn(u8, u8) -> i32> {
        Scoring {
            match_fn: |a, b| if a == b { 1 } else { -1 },
            match_scores: None,
            gap_open: self.config.args.adapter_gap,
            gap_extend: self.config.args.adapter_gap,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn process(&self, mut record: Record) -> Option<Record> {
        if self.context.adapters.is_empty() {
            return Some(record);
        }

        let mut aligner = Aligner::with_capacity_and_scoring(
            record.seq().len(),
            self.max_adapter_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<Alignment> = None;

        for adapter in &self.context.adapters {
            let alignment = aligner.semiglobal(adapter.as_bytes(), record.seq());
            let overlap = alignment.yend - alignment.ystart;
            let adapter_overlap = alignment.xend - alignment.xstart;
            
            if overlap >= self.config.args.adapter_min_overlap && adapter_overlap >= self.config.args.adapter_min_overlap {
                let matches = (alignment.score + overlap as i32) / 2;
                let errors = overlap as i32 - matches;
                let error_rate = errors as f32 / overlap as f32;

                if error_rate <= self.config.args.adapter_error_rate {
                    if best_alignment.as_ref().map_or(true, |best| alignment.score > best.score) {
                        best_alignment = Some(alignment);
                    }
                }
            }
        }

        if let Some(alignment) = best_alignment {
            let mut new_seq = record.seq().to_vec();
            let mut new_qual = record.qual().to_vec();

            match self.config.args.adapter_trim_end {
                TrimEnd::Right | TrimEnd::RTail => {
                    new_seq.truncate(alignment.ystart);
                    new_qual.truncate(alignment.ystart);
                }
                TrimEnd::Left | TrimEnd::LTail => {
                    new_seq = new_seq[alignment.yend..].to_vec();
                    new_qual = new_qual[alignment.yend..].to_vec();
                }
                TrimEnd::Any => {
                    if alignment.ystart < record.seq().len() - alignment.yend {
                         new_seq = new_seq[alignment.yend..].to_vec();
                         new_qual = new_qual[alignment.yend..].to_vec();
                    } else {
                         new_seq.truncate(alignment.ystart);
                         new_qual.truncate(alignment.ystart);
                    }
                }
            }
            
            record = Record::with_attrs(record.id(), record.desc(), &new_seq, &new_qual);
        }

        Some(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::Context;
    use crate::args::Args;
    use crate::config::Config;
    use bio::io::fastq::Record;
    use clap::Parser;

    #[test]
    fn test_process_right() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "RIGHT", "--adapter-min-overlap", "5"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = Record::with_attrs("read1", None, b"ACTGTAGATCGGAAGAG", b"IIIIIIIIIIIIIIIII");
        let processed = worker.process(record).unwrap();
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }

    #[test]
    fn test_process_left() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--adapter-seq", "AGATCGGAAGAG", "--adapter-trim-end", "LEFT", "--adapter-min-overlap", "5"]);
        let config = Config::new(args);
        let mut context = Context::new();
        context.adapters.push("AGATCGGAAGAG".to_string());
        
        let worker = SingleWorker::new(&config, &context);
        
        let record = Record::with_attrs("read1", None, b"AGATCGGAAGAGACTGT", b"IIIIIIIIIIIIIIIII");
        let processed = worker.process(record).unwrap();
        
        assert_eq!(processed.seq(), b"ACTGT");
        assert_eq!(processed.qual(), b"IIIII");
    }
}