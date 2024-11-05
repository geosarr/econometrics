use std::{
    fmt::Debug,
    ops::{AddAssign, DivAssign, Index},
};

use num_traits::{zero, Float};

use crate::binary_treatment::{self, BinaryTreatment};

#[derive(Debug)]
pub struct ConditionalBinaryTreatment<F> {
    pub candidate_causal_effect: Option<F>,
    pub intercept: Option<F>,
    pub sample_size: Option<F>,
}

impl<F> ConditionalBinaryTreatment<F> {
    pub fn new() -> Self {
        Self {
            candidate_causal_effect: None,
            intercept: None,
            sample_size: None,
        }
    }
    pub fn fit<D, Y>(&mut self, treatment: &[D], group: &[F], output: &[Y]) -> bool
    where
        F: Float + AddAssign + DivAssign + Debug,
        for<'a> &'a D: IntoIterator<Item = &'a F>,
        Y: Index<usize, Output = F>,
        for<'a> &'a Y: IntoIterator<Item = &'a F>,
    {
        let mut candidate = F::zero();
        let mut sample_size = F::zero();
        for (i, g) in group.iter().enumerate() {
            let mut binary_treatment = BinaryTreatment::new();
            if binary_treatment.fit(&treatment[i], &output[i]) {
                let (cce_i, n_sample) = (
                    binary_treatment.candidate_causal_effect.unwrap(),
                    binary_treatment.sample_size.unwrap(),
                );
                println!("{:?}, {:?}", n_sample, cce_i);
                candidate += cce_i * n_sample;
                sample_size += n_sample;
            } else {
                return false; // TODO: error handling
            }
        }
        if sample_size > F::zero() {
            self.sample_size = Some(sample_size);
            self.candidate_causal_effect = Some(candidate / sample_size) // ? weighted average of candidate causal effects.
        }
        true
    }
}
