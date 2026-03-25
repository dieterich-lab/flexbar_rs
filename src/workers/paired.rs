use log::debug;
use crate::context::Context;
use crate::args::Args;
use crate::types::{PairedRecordDto, SingeRecordDto};
use crate::steps::paired_barcode_splitter::{PairedBarcodeSplitter, PairedBarcodeSplitterState};
use crate::steps::paired_filter_length::{PairedFilterLength, PairedFilterLengthState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable, PairedAlignmentResults};

pub struct PairedWorker {
    pub barcode_splitter: Option<PairedBarcodeSplitter>,
    pub filter_length: PairedFilterLength,
}

pub struct PairedWorkerState {
    pub barcode_splitter: Option<PairedBarcodeSplitterState>,
    pub filter_length: PairedFilterLengthState,
}

impl PairedWorker {
    pub fn new(args: &Args, context: &Context) -> Self {
        let barcode_splitter = match context.barcodes.is_empty() {
            true => None,
            _ => Some(PairedBarcodeSplitter::new(args, context))
        };
        let filter_length = PairedFilterLength::new(args, context);
        
        if args.trace {
            debug!("Paired worker pipeline:");
            if let Some(s) = &barcode_splitter {
                debug!("- Step PairedBarcodeSplitter (trim end {:?}, keep: {}, unassigned: {})",
                    s.trim_end, s.barcode_keep, s.barcode_unassigned
                );
            }
            debug!("- Step PairedFilterLength (min length: {})",
                filter_length.min_read_length
            );
        }

        Self {
            barcode_splitter,
            filter_length,
        }
    }

    pub fn create_state(&self) -> PairedWorkerState {
        PairedWorkerState {
            barcode_splitter: self.barcode_splitter.as_ref().map(|s| s.create_state()),
            filter_length: self.filter_length.create_state(),
        }
    }
    
    pub fn process_record(
        &self,
        record: &PairedRecordDto,
        state: &mut PairedWorkerState
    ) -> FinalRecord<PairedRecordDto> {
        let mut rec_in: RecordInTransformation<PairedRecordDto, PairedAlignmentResults> = record.get_record_in_transformation();
        rec_in.stats.total_bases_read = record.records[0].seq().len() + record.records[1].seq().len();

        if let Some(step) = &self.barcode_splitter {
            rec_in = step.apply(rec_in, state.barcode_splitter.as_mut().unwrap());
        }

        rec_in = self.filter_length.apply(rec_in, &mut state.filter_length);

        finalize(rec_in)
    }
}

fn finalize(rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>) -> FinalRecord<PairedRecordDto> {
    if rec.alignment_results.record_dropped {
        return FinalRecord {
            record: None,
            output_tag: None,
            stats: rec.stats,
        }
    }

    let mut records = rec.record.records;
    
    // Apply trimming results from alignment_results to both records
    for i in 0..2 {
        let results = &rec.alignment_results.per_read[i];
        let range = &results.read_range;
        let old_rec = &records[i];
        
        let mut desc_parts = Vec::new();
        if let Some(desc) = old_rec.desc() {
            desc_parts.push(desc);
        }
        for tag in &results.removal_tags {
            desc_parts.push(tag);
        }
        for tag in &results.umi_tags {
            desc_parts.push(tag);
        }

        let new_desc = if desc_parts.is_empty() {
            None
        } else {
            Some(desc_parts.join("_"))
        };

        records[i] = SingeRecordDto::with_attrs(
            old_rec.id(),
            new_desc.as_deref(),
            &old_rec.seq()[range.clone()],
            &old_rec.qual()[range.clone()],
        );
    }

    FinalRecord {
        record: Some(PairedRecordDto { records }),
        output_tag: rec.alignment_results.barcode_tag,
        stats: rec.stats,
    }
}
