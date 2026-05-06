
use std::ops::Range;

pub trait GenericRawRead {
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
    fn id(&self) -> &str;
    fn set_seq(&self, seq: &[u8]) -> Self;

    fn trim(&self, range: &Range<usize>) -> Self;
}
