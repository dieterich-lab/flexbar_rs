use std::ops::Range;

pub trait GenericRawRead {
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
    fn id(&self) -> &str;
    fn desc(&self) -> &[u8];
    fn transform(self, range: Range<usize>, extra_desc_parts: Vec<String>, replace_with_dummy: bool) -> Self;
}
