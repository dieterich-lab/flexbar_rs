use log::debug;
use crate::context::Context;
use crate::args::Args;
use crate::types::{PairedRecordDto, SingeRecordDto};
use crate::steps::paired_barcode_splitter::{PairedBarcodeSplitter, PairedBarcodeSplitterState};
use crate::steps::paired_filter_length::{PairedFilterLength, PairedFilterLengthState};
use crate::steps::generic::{FinalRecord, GenericStep, RecordInTransformation, Transformable};
use crate::steps::utils::paired::PairedAlignmentResults;

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
    if rec.alignment_results.per_read[0].dropped && rec.alignment_results.per_read[1].dropped {
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
        let old_rec = &records[i];
        
        if results.dropped {
            records[i] = SingeRecordDto::with_attrs(
                old_rec.id(),
                None,
                &[],
                &[],
            );
            continue;
        }

        let range = &results.read_range;
        
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
        
        if rec.alignment_results.per_read[1-i].dropped {
            desc_parts.push("single");
        }

        let new_desc = if desc_parts.is_empty() {
            None
        } else {
            Some(desc_parts.join("_"))
        };

        let seq = if range.start == 0 && range.end == 1 && old_rec.seq().len() > 0 && !results.dropped && 
                    (range.end - range.start < old_rec.seq().len()) && old_rec.seq()[0] != b'N' && range.clone().len() == 1 {
             // This is a bit hacky to detect the 'N' replacement from paired_filter_length
             // Actually, let's just use the range. If range is 0..1 and it's not dropped, and original was longer...
             // Wait, paired_filter_length sets range to 0..1.
             // If we want to replace with 'N', we should do it here or there.
             // The requirement says: replace short reads with the sequence 'N', using the quality of the first base in the original sequence.
             b"N".as_slice()
        } else {
             &old_rec.seq()[range.clone()]
        };

        let qual = if range.start == 0 && range.end == 1 && old_rec.seq().len() > 0 && seq == b"N" {
            &old_rec.qual()[0..1]
        } else {
            &old_rec.qual()[range.clone()]
        };

        records[i] = SingeRecordDto::with_attrs(
            old_rec.id(),
            new_desc.as_deref(),
            seq,
            qual,
        );
    }

    FinalRecord {
        record: Some(PairedRecordDto { records }),
        output_tag: rec.alignment_results.barcode_tag,
        stats: rec.stats,
    }
}
