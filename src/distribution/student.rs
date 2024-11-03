use core::f64::consts::PI;
use std::fmt::Debug;

use num_traits::{Float, FloatConst, NumCast};

use super::{cdf_beta, cdf_n01, lngamma};

/// Computes the CDF of the (central) Student's t-distribution with degree
/// of freedom `n`.
///
/// Adapted from this [crate][page]
///
/// page: https://docs.rs/distrs/latest/src/distrs/students_t.rs.html#28-116
pub(crate) fn cdf_t<T: num_traits::Float + num_traits::FloatConst>(x: T, n: T) -> Option<T> {
    let one = T::one();
    // TODO support n > 0.0
    if x.is_nan() || n.is_nan() || n < one {
        return Some(T::nan());
    }

    let zero = T::zero();
    if x == T::neg_infinity() {
        return Some(zero);
    }
    if x == T::infinity() {
        return Some(one);
    }

    if n == T::infinity() {
        return cdf_n01(x);
    }

    let (start, sign) = if x < zero { (zero, one) } else { (one, -one) };

    let mut z = one;
    let t = x * x;
    let mut y = t / n;
    let mut b = one + y;

    if n > n.floor() || (n >= T::from(20.0).unwrap() && t < n) || n > T::from(200.0).unwrap() {
        // asymptotic series for large or noninteger n
        if y > T::from(10e-6).unwrap() {
            y = b.ln();
        }
        let a = n - T::from(0.5).unwrap();
        b = T::from(48.0).unwrap() * a * a;
        y = y * a;
        y = (((((T::from(-0.4).unwrap() * y - T::from(3.3).unwrap()) * y
            - T::from(24.0).unwrap())
            * y
            - T::from(85.5).unwrap())
            / (T::from(0.8).unwrap() * y * y + T::from(100.0).unwrap() + b)
            + y
            + T::from(3.0).unwrap())
            / b
            + T::from(1.0).unwrap())
            * y.sqrt();
        return cdf_n01(-y).map(|cdf_y| one - (start + sign * cdf_y));
    }

    // make n mutable and int
    // n is int between 1 and 200 if made it here
    let mut n = <usize as NumCast>::from(n).unwrap();
    let pi = T::PI();
    let two = T::from(2.).unwrap();

    if n < 20 && t < T::from(4.0).unwrap() {
        // nested summation of cosine series
        y = y.sqrt();
        let mut a = y;
        if n == 1 {
            a = zero;
        }

        // loop
        if n > 1 {
            n -= 2;
            while n > 1 {
                a = T::from(n - 1).unwrap() / (b * T::from(n).unwrap()) * a + y;
                n -= 2;
            }
        }
        a = if n == 0 {
            a / b.sqrt()
        } else {
            (y.atan() + a / b) * (two / pi)
        };
        return Some(start + sign * (z - a) / two);
    }

    // tail series expanation for large t-values
    let mut a = b.sqrt();
    y = a * T::from(n).unwrap();
    let mut j = 0;
    while a != z {
        j += 2;
        z = a;
        y = y * T::from(j - 1).unwrap() / (b * T::from(j).unwrap());
        a = a + y / T::from(n + j).unwrap();
    }
    z = zero;
    y = zero;
    a = -a;

    // loop (without n + 2 and n - 2)
    while n > 1 {
        a = T::from(n - 1).unwrap() / (b * T::from(n).unwrap()) * a + y;
        n -= 2;
    }
    a = if n == 0 {
        a / b.sqrt()
    } else {
        (y.atan() + a / b) * (two / pi)
    };
    Some(start + sign * (z - a) / two)
}

/// Computes the probability mass function of Poisson distribution
fn poisson_pmf<T: Float + FloatConst>(lambda: T, k: T) -> T {
    let zero = T::zero();
    let one = T::one();
    if lambda <= zero {
        return if k == zero { one } else { zero };
    }
    let log_pmf = k * lambda.ln() - lambda - lngamma(k + one);
    log_pmf.exp()
}

/// Computes the CDF of the non-central Student's t-distribution
pub(crate) fn cdf_nt<T: Float + FloatConst + Debug>(t: T, df: T, ncp: T) -> Option<T> {
    let zero = T::from(0.0).unwrap();
    if !df.is_finite() || df <= zero {
        return None;
    }

    if ncp == zero {
        return cdf_t(t, df);
    }

    // Use series expansion for the non-central case
    let max_terms = T::from(100.).unwrap();
    let tolerance = T::from(1e-10).unwrap();

    let lambda = T::from(0.5).unwrap() * ncp.powi(2);
    let mut sum = Some(zero);
    let mut term = Some(zero);
    let mut k = zero;
    let one = T::from(1.).unwrap();
    let two = T::from(2.).unwrap();

    while k < max_terms {
        let poisson_term = poisson_pmf(lambda, k);
        let adjusted_t = t - ncp * (k + one).sqrt() / df.sqrt();
        let cdf_term = cdf_t(adjusted_t, df + two * k);
        println!("adj: {:?}, k: {:?}, sum: {:?}\n", adjusted_t, cdf_term, sum);

        term = cdf_term.map(|t| poisson_term * t);
        if term.is_none() {
            break;
        }
        if term.unwrap().abs() < tolerance {
            break;
        }
        sum = term.zip(sum).map(|(t, s)| t + s);

        k = k + one;
    }

    sum.map(|s| s.max(zero).min(one))
}

// const M_LN2: f64 = 0.693147180559945309417232121458;
// const M_SQRT_2DPI: f64 = 0.797884560802865355879892119869;
// const M_LN_SQRT_PI: f64 = 0.572364942924700087071713675677;

// pub(crate) fn cdf_nt(t: f64, df: f64, ncp: f64) -> Option<f64> {
//     let (mut albeta, mut a, mut b, mut del, mut errbd, mut lambda, mut rxb,
// mut tt, mut x) =         (0., 0., 0., 0., 0., 0., 0., 0., 0.);
//     let (mut geven, mut godd, mut p, mut q, mut s, mut tnc, mut xeven, mut
// xodd) =         (0., 0., 0., 0., 0., 0., 0., 0.);
//     let (mut it, mut negdel) = (0, true);

//     /* note - itrmax and errmax may be changed to suit one's needs. */
//     let itrmax: isize = 1000;
//     let errmax: f64 = 1e-12;

//     // if (df <= 0.0);
//     // if(ncp == 0.0) return pt(t, df, lower_tail, log_p);

//     if t.is_infinite() {
//         return if t < 0. { Some(0.) } else { Some(1.) };
//     }
//     if t >= 0. {
//         negdel = false;
//         tt = t;
//         del = ncp;
//     } else {
//         /* We deal quickly with left tail if extreme,
//         since pt(q, df, ncp) <= pt(0, df, ncp) = \Phi(-ncp) */
//         // if (ncp > 40) && (!log_p | !lower_tail) {
//         if ncp > 40. {
//             return Some(0.);
//         };
//         negdel = true;
//         tt = -t;
//         del = -ncp;
//     }

//     if (df > 4e5) | (del.powi(2) > 2. * M_LN2 * (-f64::MIN_EXP as f64)) {
//         /*-- 2nd part: if del > 37.62, then p=0 below
//         FIXME: test should depend on `df', `tt' AND `del' ! */
//         /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */
//         s = 1. / (4. * df);
//         let sigma = (1. + tt.powi(2) * 2. * s).sqrt();
//         let cdf = cdf_n01(tt * (1. - s) - del / sigma);
//         return if !negdel { cdf } else { cdf.map(|c| 1. - c) };
//     }

//     /* initialize twin series */
//     /* Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */
//     x = t * t;
//     rxb = df / (x + df); /* := (1 - x) {x below} -- but more accurately */
//     x = x / (x + df); /* in [0,1) */
//     if x > 0. {
//         /* <==>  t != 0 */
//         lambda = del.powi(2);
//         p = 0.5 * (-0.5 * lambda).exp();

//         if p == 0. {
//             /* underflow! */
//             /* ========== really use an other algorithm for this case !!! */
//             return Some(0.);
//         }
//         q = M_SQRT_2DPI * p * del;
//         s = 0.5 - p;
//         /* s = 0.5 - p = 0.5*(1 - exp(-.5 L)) =  -0.5*expm1(-.5 L)) */
//         if s < 1e-7 {
//             s = -0.5 * (-0.5 * lambda).exp_m1();
//         }
//         a = 0.5;
//         b = 0.5 * df;
//         /* rxb = (1 - x) ^ b   [ ~= 1 - b*x for tiny x --> see 'xeven' below]
//          * where '(1 - x)' =: rxb {accurately!} above */
//         rxb = rxb.powf(b);
//         albeta = M_LN_SQRT_PI + lgamma(b) - lgamma(0.5 + b);
//         xodd = cdf_beta(x, a, b);
//         godd = 2. * rxb * (a * x.ln() - albeta).exp();
//         tnc = b * x;
//         xeven = if tnc < f64::EPSILON { tnc } else { 1. - rxb };
//         geven = tnc * rxb;
//         tnc = p * xodd + q * xeven;

//         /* repeat until convergence or iteration limit */
//         it = 1;
//         while it <= itrmax {
//             a += 1.;
//             xodd -= godd;
//             xeven -= geven;
//             godd *= x * (a + b - 1.) / a;
//             geven *= x * (a + b - 0.5) / (a + 0.5);
//             let itf = it as f64;
//             p *= lambda / (2. * itf);
//             q *= lambda / (2. * itf + 1.);
//             tnc += p * xodd + q * xeven;
//             s -= p;
//             /* R 2.4.0 added test for rounding error here. */
//             if s < -1e-10 {
//                 /* happens e.g. for (t,df,ncp)=(40,10,38.5), after 799 it. */
//                 break;
//                 // goto finis;
//             }
//             if s <= 0. && it > 1 {
//                 break;
//                 // goto finis;
//             }
//             errbd = 2. * s * (xodd - godd);
//             if errbd.abs() < errmax {
//                 break;
//                 //  goto finis;/*convergence*/
//             }
//         }
//         /* non-convergence: */
//     } else {
//         /* x = t = 0 */
//         tnc = 0.;
//     }
//     //  finis:
//     let res = cdf_n01(-del).map(|cdf| cdf + tnc);

//     // lower_tail = lower_tail != negdel; /* xor */
//     // if(tnc > 1 - 1e-10 && lower_tail)
//     // ML_WARNING(ME_PRECISION, "pnt{final}");

//     return res.map(|r| r.min(1.)) /* Precaution */;
// }
