use log::trace;
use crate::steps::generic::{RecordInTransformation, RecordId};

pub mod combined;

pub fn assign_barcode_tag<RecordType>(
    rec: &mut RecordInTransformation<RecordType>,
    tag: String,
    trace: bool
) where RecordType: RecordId {
    if trace {
        trace!("Barcode splitter on {}: Assigned tag {}", rec.record_id(), tag);
    }
    rec.alignment_results.set_barcode_tag(Some(tag));
}

pub trait BarcodeAlignmentResults {
    fn set_barcode_tag(&mut self, tag: Option<String>);
}
