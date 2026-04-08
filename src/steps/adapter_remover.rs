use crate::args::Args;
use crate::types::{RecordDto, AdapterReadset, RunType};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, CombinedAligmentResults};
use crate::context::Context;
use crate::multi_aligner::MultiAlignerState;
use crate::steps::utils::adapter_remover::AdapterRemover;

pub struct AdapterRemoverStep {
    pub inner: [Option<AdapterRemover>; 2],
}

pub struct AdapterRemoverStepState {
    pub states: [Option<MultiAlignerState>; 2],
}

impl GenericStep<RecordDto> for AdapterRemoverStep {
    type State = AdapterRemoverStepState;
    type AlignmentResult = CombinedAligmentResults;

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

    fn apply(&self, mut rec: RecordInTransformation<RecordDto, CombinedAligmentResults>, state: &mut Self::State) -> RecordInTransformation<RecordDto, CombinedAligmentResults> {
        match &rec.record {
            RecordDto::Single(s) => {
                if rec.alignment_results.per_read[0].as_ref().unwrap().dropped { return rec; }
                if let Some(remover) = &self.inner[0] {
                    if remover.trace {
                        trace!("Adapter remover processing record: {}", s.id());
                    }
                    remover.remove_adapters(
                        state.states[0].as_mut().unwrap(),
                        s,
                        None,
                        rec.alignment_results.per_read[0].as_mut().unwrap(),
                        &mut rec.stats
                    );
                }
            }
            RecordDto::Paired(p) => {
                for i in 0..2 {
                    if rec.alignment_results.per_read[i].as_ref().unwrap().dropped { continue; }
                    
                    if let Some(remover) = &self.inner[i] {
                        if remover.trace {
                            trace!("Adapter remover processing record: {} (R{})", p.id(), i + 1);
                        }
                        remover.remove_adapters(
                            state.states[i].as_mut().unwrap(),
                            &p.reads[i],
                            Some(i),
                            rec.alignment_results.per_read[i].as_mut().unwrap(),
                            &mut rec.stats
                        );
                    }
                }
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
    use crate::types::{RecordDto, SingeRecordDto, PairedRecordDto};
    use crate::steps::generic::{GenericStep, Transformable};
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
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));

        let rec_in = rec.get_record_in_transformation();
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
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));

        let rec_in = rec.get_record_in_transformation();
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
        let rec = RecordDto::Single(r1);

        let rec_in = rec.get_record_in_transformation();
        let result = step.apply(rec_in, &mut state);

        assert_eq!(result.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..8);
        assert_eq!(result.stats.total_adapters_removed, 1);
    }
}
