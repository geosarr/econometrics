use std::f64;

use super::lgamma;

/// Computes the Beta function B(a, b)
fn beta_function(a: f64, b: f64) -> f64 {
    (lgamma(a) + lgamma(b) - lgamma(a + b)).exp()
}

/// Computes the probability density function of Beta distribution
fn pdf_beta(x: f64, alpha: f64, beta: f64) -> f64 {
    if x <= 0.0 || x >= 1.0 {
        return 0.0;
    }
    let beta_const = 1.0 / beta_function(alpha, beta);
    beta_const * x.powf(alpha - 1.0) * (1.0 - x).powf(beta - 1.0)
}

/// Computes the cumulative distribution function of Beta distribution
/// using numerical integration (Simpson's rule)
pub(crate) fn cdf_beta(x: f64, alpha: f64, beta: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    // Simpson's rule implementation
    let n = 1000; // number of intervals (must be even)
    let h = x / n as f64;

    let mut sum = pdf_beta(0.0, alpha, beta) + pdf_beta(x, alpha, beta);

    // Even terms
    for i in (2..n).step_by(2) {
        let xi = i as f64 * h;
        sum += 2.0 * pdf_beta(xi, alpha, beta);
    }

    // Odd terms
    for i in (1..n).step_by(2) {
        let xi = i as f64 * h;
        sum += 4.0 * pdf_beta(xi, alpha, beta);
    }

    (h / 3.0) * sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cdf_beta_bounds() {
        let alpha = 2.0;
        let beta = 3.0;

        assert!((cdf_beta(0.0, alpha, beta) - 0.0).abs() < 1e-10);
        assert!((cdf_beta(1.0, alpha, beta) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cdf_beta_symmetry() {
        let x = 0.5;
        let alpha = 2.0;
        let beta = 2.0;

        // For alpha = beta, CDF(0.5) should be 0.5
        assert!((cdf_beta(x, alpha, beta) - 0.5).abs() < 1e-3);
    }
}
