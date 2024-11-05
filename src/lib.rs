pub mod binary_treatment;
pub mod conditional_binary_treatment;
pub mod data;
pub mod distribution;
pub mod statistical_test;
use num_traits::NumCast;

use distribution::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        binary_treatment::BinaryTreatment,
        cdf_t,
        conditional_binary_treatment::ConditionalBinaryTreatment,
        data::lalonde::*,
        statistical_test::{two_sample_homoscedastic_ttest, TestTSide},
    };
    use ndarray::*;
    #[test]
    fn bin_treat() {
        let income = Array::from_iter(RE78);
        let treatment = Array::from_iter(TREAT.iter().map(|v| *v as f32));
        let mut binary_treatment: BinaryTreatment<f32> = BinaryTreatment::new();
        binary_treatment.fit(&treatment, &income);
        println!("{:#?}", binary_treatment);
    }

    #[test]
    fn cond_bin_treat() {
        let mut binary_treatment: ConditionalBinaryTreatment<f32> =
            ConditionalBinaryTreatment::new();
        binary_treatment.fit(
            &[
                [0., 1., 0., 0., 1.],
                [0., 0., 0., 0., 1.],
                [1., 1., 1., 0., 0.],
            ],
            &[0., 1., 2.],
            &[
                [0., 10., 6., 0., 1.],
                [0., 20., 5., 8., 1.],
                [10., 1., 1., 0., 0.],
            ],
        );
        println!("{:#?}", binary_treatment);
    }

    #[test]
    fn student_cdf() {
        // let t = (25f64).sqrt() * (2800. - 3000.) / 600.;
        // println!("{t}");
        // println!("{:?}", 1. - cdf_t(t, 24.).unwrap());
        // println!("{:?}", cdf_t(t, 24.).unwrap());
        // println!("{:?}", 2. * (1. - cdf_t(t.abs(), 24.).unwrap()));
        // println!("{:?}", cdf_t(100., 10.));
        println!(
            "{:?}",
            two_sample_homoscedastic_ttest(
                0f64,
                (1180., 1353.75),
                (8., 8.),
                (26.32218f64.powi(2), 28.02661f64.powi(2)),
                TestTSide::TwoSided
            )
        );
        // println!("{}", lgamma(0.001));
        println!("{:?}", cdf_nt(1., 10., 1.));
    }
}
