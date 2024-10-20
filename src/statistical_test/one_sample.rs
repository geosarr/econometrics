use num_traits::Float;

use crate::{cdf_n01, mean};

use super::{TestKind, TestOutput};

/// Gauss-Test.
///
/// Assuming the samples X<sub>1</sub>, . . . , X<sub>n</sub> ∼ N(μ,
/// σ<sup>2</sup>) are i.i.d., where N(μ, σ<sup>2</sup>) is the gaussian law
/// with mean μ and **known variance** σ<sup>2</sup>. This function can test the
/// following:
///
/// (a) H<sub>0</sub> : μ ≤ μ<sub>0</sub> vs H<sub>1</sub> : μ > μ<sub>0</sub>
/// (upper one-sided test).
///
/// (b) H<sub>0</sub> : μ ≥ μ<sub>0</sub> vs H<sub>1</sub> : μ < μ<sub>0</sub>
/// (lower one-sided test).
///
/// (c) H<sub>0</sub> : μ = μ<sub>0</sub> vs H<sub>1</sub> : μ != μ<sub>0</sub>
/// (two-sided test).
pub fn gauss_test<S, T: Float>(
    sample: &S,
    mu0: T,
    sigma: T,
    test_type: TestKind,
) -> Option<TestOutput<T, T>>
where
    for<'a> &'a S: IntoIterator<Item = &'a T>,
{
    let stat = if let Some((mean, n_sample)) = mean(sample) {
        (mean - mu0) * n_sample.sqrt() / sigma
    } else {
        return None;
    };
    let one = T::one();
    let two = one + one;
    let pvalue = match test_type {
        TestKind::UpperOneSided => f32::from(stat)
            .map(cdf_n01)
            .map(T::from)
            .map(|cdf_stat| one - cdf_stat),
        TestKind::LowerOneSided => f32::from(stat).map(cdf_n01).map(T::from),
        TestKind::TwoSided => f32::from(stat.abs())
            .map(cdf_n01)
            .map(T::from)
            .map(|cdf_abs_stat| two * (one - cdf_abs_stat)),
    };
    Some(TestOutput {
        statistics: stat,
        pvalue,
    })
}
