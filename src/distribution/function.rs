use core::f64::consts::PI;

use num_traits::{Float, FloatConst};

/// Taken from https://en.wikipedia.org/wiki/Lanczos_approximation
const GAMMA_G: f64 = 7.;
const GAMMA_N: usize = 9;
const GAMMA_P: [f64; GAMMA_N] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];
const GAMMA_EPS: f64 = 1e-07;

pub(crate) fn lgamma<T: Float + FloatConst>(mut z: T) -> T {
    let one_half = T::from(0.5).unwrap();
    let one = T::from(1.).unwrap();
    if z < one_half {
        // Reflection formula
        T::PI().ln() - (T::PI() * z).sin().ln() - lgamma(one - z)
    } else {
        z = z - one;
        let mut x = T::from(GAMMA_P[0]).unwrap();
        for i in 1..GAMMA_N {
            x = x + T::from(GAMMA_P[i]).unwrap() / (z + T::from(i).unwrap());
        }
        let t = z + T::from(GAMMA_G).unwrap() + one_half;
        one_half * (T::LN_2() + T::PI().ln()) + (z + one_half) * t.ln() - t + x.ln()
    }
}
