mod one_sample;
mod two_sample;
pub use one_sample::*;
pub use two_sample::*;

/// Which side to test.
#[derive(Debug)]
pub enum TestTSide {
    /// - For one sample tests:
    /// H<sub>0</sub> : μ ≤ μ<sub>0</sub> vs H<sub>1</sub> : μ > μ<sub>0</sub>
    UpperOneSided,
    /// - For one sample tests:
    /// H<sub>0</sub> : μ ≥ μ<sub>0</sub> vs H<sub>1</sub> : μ < μ<sub>0</sub>
    LowerOneSided,
    /// - For one sample tests:
    /// H<sub>0</sub> : μ = μ<sub>0</sub> vs H<sub>1</sub> : μ != μ<sub>0</sub>
    TwoSided,
}

/// Output values of a test.
#[derive(Debug)]
pub struct TestOutput<S, P> {
    pub statistics: S,
    pub pvalue: Option<P>,
}
