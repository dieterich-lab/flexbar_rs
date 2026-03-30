use log::trace;
use crate::steps::generic::{RecordInTransformation, RecordId};

pub mod single;
pub mod paired;
pub mod barcode;
pub mod quality_trim;
pub mod adapter_remover;

pub fn assign_barcode_tag<RecordType, AlignmentResults, T: BarcodeAlignmentResults>(
    rec: &mut RecordInTransformation<RecordType, AlignmentResults>,
    tag: String,
    trace: bool
) where AlignmentResults: AsMut<T>, RecordType: RecordId {
    if trace {
        trace!("Barcode splitter on {}: Assigned tag {}", rec.record_id(), tag);
    }
    rec.alignment_results.as_mut().set_barcode_tag(Some(tag));
}

pub trait BarcodeAlignmentResults {
    fn set_barcode_tag(&mut self, tag: Option<String>);
}
