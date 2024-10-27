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
const PI: f64 = 3.14159265358979323846264338327950288;

pub(crate) fn lgamma(mut z: f64) -> f64 {
    if z < 0.5 {
        // Reflection formula
        PI.ln() - (PI * z).sin().ln() - lgamma(1. - z)
    } else {
        z -= 1.;
        let mut x = GAMMA_P[0];
        for i in 1..GAMMA_N {
            x += GAMMA_P[i] / (z + (i as f64));
        }
        let t = z + GAMMA_G + 0.5;
        0.5 * (2. * PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
    }
}
