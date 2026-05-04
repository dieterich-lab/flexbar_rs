
pub trait GenericRawRead {
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
    fn id(&self) -> &str;
    fn desc(&self) -> Option<&str>;
    fn desc_str(&self) -> Option<&str>;
}
