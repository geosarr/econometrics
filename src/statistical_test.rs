mod one_sample;

pub use one_sample::*;

/// Which side to test.
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
pub struct TestOutput<S, P> {
    pub statistics: S,
    pub pvalue: Option<P>,
}
