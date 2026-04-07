use crate::args::Args;
use crate::types::{PairedRecordDto, AdapterReadset};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;
use crate::multi_aligner::MultiAlignerState;
use crate::steps::utils::adapter_remover::AdapterRemover;

pub struct PairedAdapterRemover {
    pub inner: [Option<AdapterRemover>; 2],
}

pub struct PairedAdapterRemoverState {
    pub states: [Option<MultiAlignerState>; 2],
}

impl GenericStep<PairedRecordDto> for PairedAdapterRemover {
    type State = PairedAdapterRemoverState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        let (apply1, apply2) = match args.adapter_read_set {
            AdapterReadset::One => (true, false),
            AdapterReadset::Two => (false, true),
            AdapterReadset::None => (true, true),
        };

        let inner = [
            if apply1 && !context.adapters.is_empty() {
                Some(AdapterRemover::new(args, &context.adapters))
            } else { None },
            if apply2 && !context.adapters2.is_empty() {
                Some(AdapterRemover::new(args, &context.adapters2))
            } else { None },
        ];

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

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        for i in 0..2 {
            if rec.alignment_results.per_read[i].dropped { continue; }
            
            if let Some(remover) = &self.inner[i] {
                if remover.trace {
                    trace!("Adapter remover processing record: {} (R{})", rec.record.id(), i + 1);
                }
                remover.remove_adapters(
                    state.states[i].as_mut().unwrap(),
                    &rec.record.reads[i],
                    Some(i),
                    &mut rec.alignment_results.per_read[i],
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
    use crate::types::{PairedRecordDto, SingeRecordDto, AdapterReadset, TrimEnd};
    use crate::steps::generic::{GenericStep, Transformable};
    use clap::Parser;

    #[test]
    fn test_paired_adapter_removal() {
        let mut args = Args::parse_from(["flexbar", "--adapter-trim-end", "Right"]);
        args.adapter_trim_end = TrimEnd::Right;
        args.adapter_match = 1;
        args.adapter_mismatch = -1;
        args.adapter_gap_open = -6;
        args.adapter_gap_extend = -1;
        args.adapter_min_overlap = 1;

        let mut context = Context::empty();
        context.adapters.push(SequenceWithId {
            id: "ad1".to_string(),
            sequence: "AAAA".to_string(),
        });
        context.adapters2.push(SequenceWithId {
            id: "ad2".to_string(),
            sequence: "CCCC".to_string(),
        });

        let step = PairedAdapterRemover::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGCCCC", b"IIIIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);

        let rec_in = rec.get_record_in_transformation();
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].read_range, 0..8);
        assert_eq!(result.stats.total_adapters_removed, 2);
    }

    #[test]
    fn test_paired_adapter_removal_read_set_one() {
        let mut args = Args::parse_from(["flexbar", "--adapter-trim-end", "Right", "--adapter-read-set", "1"]);
        args.adapter_trim_end = TrimEnd::Right;
        args.adapter_read_set = AdapterReadset::One;
        args.adapter_min_overlap = 1;

        let mut context = Context::empty();
        context.adapters.push(SequenceWithId {
            id: "ad1".to_string(),
            sequence: "AAAA".to_string(),
        });

        let step = PairedAdapterRemover::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"TTTTTTTTAAAA", b"IIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"GGGGGGGGAAAA", b"IIIIIIIIIIII");
        let rec = PairedRecordDto::new(r1, r2);

        let rec_in = rec.get_record_in_transformation();
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].read_range, 0..8);
        assert_eq!(result.alignment_results.per_read[1].read_range, 0..12);
        assert_eq!(result.stats.total_adapters_removed, 1);
    }
}
