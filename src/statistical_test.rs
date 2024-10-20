mod one_sample;

pub use one_sample::gauss_test;

pub enum TestKind {
    UpperOneSided,
    LowerOneSided,
    TwoSided,
}

pub struct TestOutput<S, P> {
    pub statistics: S,
    pub pvalue: Option<P>,
}
