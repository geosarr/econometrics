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

/// Computes the natural logarithm of the Gamma function using Lanczos approximation.
pub fn lngamma<T: Float + FloatConst>(mut z: T) -> T {
    let one_half = T::from(0.5).unwrap();
    let one = T::from(1.).unwrap();
    if z < one_half {
        // Reflection formula
        T::PI().ln() - (T::PI() * z).sin().ln() - lngamma(one - z)
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

/// Computes the complete Beta function.
/// - B(a, b) = ∫<sub>0</sub><sup>1</sup> t<sup>a-1</sup>(1-t)<sup>b-1</sup>dt
///
/// Adapted from C++ boost implementation.
pub fn beta(a: f64, b: f64) -> f64 {
    if (a <= 0.) | (b <= 0.) {
        panic!("beta inputs a and b should be > 0.")
    }
    let c = a + b;
    if (c == a) && (b < f64::EPSILON) {
        return 1. / b;
    } else if (c == b) && (a < f64::EPSILON) {
        return 1. / a;
    }
    if b == 1. {
        return 1. / a;
    } else if a == 1. {
        return 1. / b;
    } else if (c < f64::EPSILON) {
        return c / (b * a);
    }
    (lngamma(a) + lngamma(b) - lngamma(a + b)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lngamma_integers() {
        let x = (1..1000).collect::<Vec<usize>>();
        for a in x {
            let ln_facto = (1..a).fold(0., |acc, x| acc + (x as f64).ln());
            assert!((lngamma(a as f64) - ln_facto).abs() < 1e-10);
        }
    }

    #[test]
    fn beta_symmetry() {
        let x = (1..1000).collect::<Vec<usize>>();
        let y = (1..1000).collect::<Vec<usize>>();
        for a in &x {
            for b in &y {
                assert!((beta(*a as f64, *b as f64) - beta(*b as f64, *a as f64)).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn beta_pascal_identity() {
        let x = (1..1000).collect::<Vec<usize>>();
        let y = (1..1000).collect::<Vec<usize>>();
        for a in &x {
            for b in &y {
                assert!(
                    (beta(*a as f64, *b as f64)
                        - beta(*a as f64, *b as f64 + 1.)
                        - beta(*a as f64 + 1., *b as f64))
                    .abs()
                        < 1e-10
                );
            }
        }
    }
}
