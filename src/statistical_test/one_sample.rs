use num_traits::Float;

use crate::{cdf_n01, mean, variance};

use super::{TestOutput, TestTSide};

/// One sample Z-test (also called Gauss-test) with known population variance.
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
pub fn one_sample_ztest<S, T: Float>(
    sample: &S,
    mu0: T,
    sigma: T,
    test_type: TestTSide,
) -> Option<TestOutput<T, T>>
where
    for<'a> &'a S: IntoIterator<Item = &'a T>,
{
    let stat = if let Some((mean, n_sample)) = mean(sample) {
        if sigma > T::zero {
            (mean - mu0) * n_sample.sqrt() / sigma
        } else {
            return None;
        }
    } else {
        return None;
    };
    let one = T::one();
    let two = one + one;
    let pvalue = match test_type {
        TestTSide::UpperOneSided => f32::from(stat)
            .map(cdf_n01)
            .map(T::from)
            .map(|cdf_stat| one - cdf_stat),
        TestTSide::LowerOneSided => f32::from(stat).map(cdf_n01).map(T::from),
        TestTSide::TwoSided => f32::from(stat.abs())
            .map(cdf_n01)
            .map(T::from)
            .map(|cdf_abs_stat| two * (one - cdf_abs_stat)),
    };
    Some(TestOutput {
        statistics: stat,
        pvalue,
    })
}

// /// One sample T-test (also called Student-test) with unknown population
// /// variance.
// ///
// /// Assuming the samples X<sub>1</sub>, . . . , X<sub>n</sub> ∼ N(μ,
// /// σ<sup>2</sup>) are i.i.d., where N(μ, σ<sup>2</sup>) is the gaussian law
// /// with mean μ and **known variance** σ<sup>2</sup>. This function can test
// the /// following:
// ///
// /// (a) H<sub>0</sub> : μ ≤ μ<sub>0</sub> vs H<sub>1</sub> : μ >
// μ<sub>0</sub> /// (upper one-sided test).
// ///
// /// (b) H<sub>0</sub> : μ ≥ μ<sub>0</sub> vs H<sub>1</sub> : μ <
// μ<sub>0</sub> /// (lower one-sided test).
// ///
// /// (c) H<sub>0</sub> : μ = μ<sub>0</sub> vs H<sub>1</sub> : μ !=
// μ<sub>0</sub> /// (two-sided test).
// pub fn one_sample_ttest<S, T: Float>(
//     sample: &S,
//     mu0: T,
//     test_type: TestTSide,
// ) -> Option<TestOutput<T, T>>
// where
//     for<'a> &'a S: IntoIterator<Item = &'a T>,
// {
//     let stat = if let Some((var, mean, n_sample)) = variance(sample) {
//         if sigma > T::zero {
//             (mean - mu0) * n_sample.sqrt() / var.sqrt()
//         } else {
//             return None;
//         }
//     } else {
//         return None;
//     };
//     let one = T::one();
//     let two = one + one;
//     let pvalue = match test_type {
//         TestTSide::UpperOneSided => f32::from(stat)
//             .map(cdf_n01)
//             .map(T::from)
//             .map(|cdf_stat| one - cdf_stat),
//         TestTSide::LowerOneSided =>
// f32::from(stat).map(cdf_n01).map(T::from),         TestTSide::TwoSided =>
// f32::from(stat.abs())             .map(cdf_n01)
//             .map(T::from)
//             .map(|cdf_abs_stat| two * (one - cdf_abs_stat)),
//     };
//     Some(TestOutput {
//         statistics: stat,
//         pvalue,
//     })
// }
