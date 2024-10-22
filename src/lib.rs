pub mod binary_treatment;
pub mod data;
pub mod statistical_test;

pub(crate) fn mean<S, T>(sample: &S) -> Option<(T, T)>
where
    T: num_traits::Float,
    for<'a> &'a S: IntoIterator<Item = &'a T>,
{
    let mut m = T::zero();
    let mut n = T::zero();
    for x in sample {
        m = m + *x;
        n = n + T::one();
    }
    if n > T::zero() {
        Some((m / n, n))
    } else {
        None
    }
}

pub(crate) fn variance<S, T>(sample: &S) -> Option<(T, T, T)>
where
    T: num_traits::Float,
    for<'a> &'a S: IntoIterator<Item = &'a T>,
{
    let (mean, mut n) = if let Some((m, s)) = mean(sample) {
        (m, s)
    } else {
        return None;
    };
    let mut var = T::zero();
    for x in sample {
        var = var + (*x - mean).powi(2);
        n = n + T::one();
    }
    if n > T::one() {
        let mean = mean / n;
        let var = (T::one() / n) * var;
        Some((var, mean, n))
    } else {
        None
    }
}

/// Computes the cumulative distribution function of the standard normal using
/// the formula of [Dia (2023)][lien].
///
/// [lien]: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4487559
pub(crate) fn cdf_n01(x: &f32) -> f32 {
    let x2 = x.powi(2);
    let abs_x = x.abs();
    let one_minus_cdf_abs_x = (0.39894228040143268 / (abs_x + 2.92678600515804815))
        * ((x2 + 8.42742300458043240 * abs_x + 18.38871225773938487)
            / (x2 + 5.81582518933527391 * abs_x + 8.97280659046817350))
        * ((x2 + 7.30756258553673541 * abs_x + 18.25323235347346525)
            / (x2 + 5.70347935898051437 * abs_x + 10.27157061171363079))
        * ((x2 + 5.66479518878470765 * abs_x + 18.61193318971775795)
            / (x2 + 5.51862483025707963 * abs_x + 12.72323261907760928))
        * ((x2 + 4.91396098895240075 * abs_x + 24.14804072812762821)
            / (x2 + 5.26184239579604207 * abs_x + 16.88639562007936908))
        * ((x2 + 3.83362947800146179 * abs_x + 11.61511226260603247)
            / (x2 + 4.92081346632882033 * abs_x + 24.12333774572479110))
        * ((-x2 / 2.).exp());
    if *x > 0. {
        1. - one_minus_cdf_abs_x
    } else {
        one_minus_cdf_abs_x
    }
}

/// Computes the cumulative distribution function of a Student with degree of
/// freedom `n`.
///
/// Adapted from this [crate][page]
///
/// page: https://docs.rs/distrs/latest/src/distrs/students_t.rs.html#28-116
pub(crate) fn cdf_student(x: f32, n: f32) -> f32 {
    use core::f32::consts::PI;

    // TODO support n > 0.0
    if x.is_nan() || n.is_nan() || n < 1.0 {
        return f32::NAN;
    }

    if x == f32::NEG_INFINITY {
        return 0.0;
    }

    if x == f32::INFINITY {
        return 1.0;
    }

    if n == f32::INFINITY {
        return cdf_n01(&(x as f32)) as f32;
    }

    let (start, sign) = if x < 0.0 { (0.0, 1.0) } else { (1.0, -1.0) };

    let mut z = 1.0;
    let t = x * x;
    let mut y = t / n;
    let mut b = 1.0 + y;

    if n > n.floor() || (n >= 20.0 && t < n) || n > 200.0 {
        // asymptotic series for large or noninteger n
        if y > 10e-6 {
            y = b.ln();
        }
        let a = n - 0.5;
        b = 48.0 * a * a;
        y *= a;
        y = (((((-0.4 * y - 3.3) * y - 24.0) * y - 85.5) / (0.8 * y * y + 100.0 + b) + y + 3.0)
            / b
            + 1.0)
            * y.sqrt();
        return start + sign * (cdf_n01(&(-y as f32)) as f32);
    }

    // make n mutable and int
    // n is int between 1 and 200 if made it here
    let mut n = n as u8;

    if n < 20 && t < 4.0 {
        // nested summation of cosine series
        y = y.sqrt();
        let mut a = y;
        if n == 1 {
            a = 0.0;
        }

        // loop
        if n > 1 {
            n -= 2;
            while n > 1 {
                a = (n - 1) as f32 / (b * n as f32) * a + y;
                n -= 2;
            }
        }
        a = if n == 0 {
            a / b.sqrt()
        } else {
            (y.atan() + a / b) * (2.0 / PI)
        };
        return start + sign * (z - a) / 2.0;
    }

    // tail series expanation for large t-values
    let mut a = b.sqrt();
    y = a * n as f32;
    let mut j = 0;
    while a != z {
        j += 2;
        z = a;
        y = y * (j - 1) as f32 / (b * j as f32);
        a += y / (n + j) as f32;
    }
    z = 0.0;
    y = 0.0;
    a = -a;

    // loop (without n + 2 and n - 2)
    while n > 1 {
        a = (n - 1) as f32 / (b * n as f32) * a + y;
        n -= 2;
    }
    a = if n == 0 {
        a / b.sqrt()
    } else {
        (y.atan() + a / b) * (2.0 / PI)
    };
    start + sign * (z - a) / 2.0
}

#[cfg(test)]
mod tests {
    use crate::{binary_treatment::BinaryTreatment, cdf_student, data::lalonde::*};
    use ndarray::*;
    #[test]
    fn binary_treatment() {
        let income = Array::from_iter(RE78);
        let treatment = Array::from_iter(TREAT.iter().map(|v| *v as f32));
        let mut binary_treatment: BinaryTreatment<f32> = BinaryTreatment::new();
        binary_treatment.fit(&treatment, &income);
        println!("{:#?}", binary_treatment);
    }

    #[test]
    fn student_cdf() {
        print!("{}", cdf_student(100., 10.))
    }
}
