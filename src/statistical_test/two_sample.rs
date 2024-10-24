use num_traits::{Float, FloatConst};

use crate::cdf_t;

use super::{TestOutput, TestTSide};

/// Two samples T-test assuming equal population variances.
///
/// Assuming having independent populations X<sub>1,1</sub>, . . . ,
/// X<sub>1,n<sub>1</sub></sub> and X<sub>2,1</sub>, . . . ,
/// X<sub>2,n<sub>2</sub></sub> such that:
/// - for all j=1..n<sub>1</sub>, X<sub>1,j</sub> ∼ N(μ<sub>1</sub>,
///   σ<sup>2</sup><sub>1</sub>) are i.i.d.,
/// - for all j=1..n<sub>2</sub>, X<sub>2,j</sub> ∼ N(μ<sub>2</sub>,
///   σ<sup>2</sup><sub>2</sub>) are i.i.d.,
///
/// where N(μ, σ<sup>2</sup>) designates the gaussian law with mean μ and
/// variance σ<sup>2</sup>.
///
/// Let ∆ = μ<sub>1</sub> - μ<sub>2</sub> be the spread of the samples means.
/// This function can test the following:
///
/// (a) H<sub>0</sub> : ∆ ≤ ∆<sub>0</sub> vs H<sub>1</sub> : ∆ > ∆<sub>0</sub>
/// (upper one-sided test).
///
/// (b) H<sub>0</sub> : ∆ ≥ ∆<sub>0</sub> vs H<sub>1</sub> : ∆ < ∆<sub>0</sub>
/// (lower one-sided test).
///
/// (c) H<sub>0</sub> : ∆ = ∆<sub>0</sub> vs H<sub>1</sub> : ∆ != ∆<sub>0</sub>
/// (two-sided test).
pub fn two_sample_homoscedastic_ttest<T: Float + FloatConst>(
    delta0: T,
    sample_means: (T, T),
    sample_sizes: (T, T),
    sample_vars: (T, T),
    test_type: TestTSide,
) -> Option<TestOutput<T, T>> {
    let (m1, m2) = sample_means;
    let (n1, n2) = sample_sizes;
    let (v1, v2) = sample_vars;
    let one = T::one();
    let two = one + one;
    let n = n1 + n2;
    if (n1 < one) | (n2 < one) | (n < two) {
        return None;
    }
    let var = (one / (n - two)) * ((n1 - one) * v1 + (n2 - one) * v2);
    let stat = if var > T::zero() {
        (m1 - m2 - delta0) / (var.sqrt() * (one / n1 + one / n2).sqrt())
    } else {
        return None;
    };
    let pvalue = match test_type {
        TestTSide::UpperOneSided => cdf_t(stat, n - two).map(|cdf_stat| one - cdf_stat),
        TestTSide::LowerOneSided => cdf_t(stat, n - two),
        TestTSide::TwoSided => {
            cdf_t(stat.abs(), n - two).map(|cdf_abs_stat| two * (one - cdf_abs_stat))
        }
    };
    Some(TestOutput {
        statistics: stat,
        pvalue,
    })
}
