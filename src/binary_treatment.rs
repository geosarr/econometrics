use std::ops::{AddAssign, Index};

use num_traits::Float;

#[derive(Debug)]
pub struct BinaryTreatment<F> {
    pub candidate_causal_effect: Option<F>,
    pub intercept: Option<F>,
}

impl<F> BinaryTreatment<F> {
    pub fn new() -> Self {
        Self {
            candidate_causal_effect: None,
            intercept: None,
        }
    }
    pub fn fit<D, Y>(&mut self, treatment: &D, output: &Y) -> bool
    where
        F: Float,
        for<'a> &'a D: IntoIterator<Item = &'a F>,
        Y: Index<usize, Output = F>,
        for<'a> &'a Y: IntoIterator<Item = &'a F>,
    {
        let mut sum_output_non_treat = F::zero();
        let mut n_non_treat = F::zero();
        let mut sum_output_treat = F::zero();
        let mut sum_treatment = F::zero();
        let mut n_treat = F::zero();
        let mut pos = 0;
        for treat in treatment {
            if *treat == F::zero() {
                sum_output_non_treat = sum_output_non_treat + output[pos];
                n_non_treat = n_non_treat + F::one();
            } else {
                sum_output_treat = sum_output_treat + output[pos];
                n_treat = n_treat + F::one();
                sum_treatment = sum_treatment + F::one();
            }
            pos += 1;
        }
        if (n_non_treat == F::zero()) | (n_treat == F::zero()) {
            return false;
        }
        let mean_output_treat = sum_output_treat / n_treat;
        let candidate_causal_effect = mean_output_treat - mean_output_non_treat;
        let n_sample = n_non_treat + n_treat;
        let mean_output = (sum_output_non_treat + sum_output_treat) / n_sample;
        let mean_treatment = sum_treatment / n_sample;
        self.intercept = Some(mean_output - candidate_causal_effect * mean_treatment);
        self.candidate_causal_effect = Some(candidate_causal_effect);
        true
    }
}
