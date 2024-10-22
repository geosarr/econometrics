use num_traits::NumCast;

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
/// the formula of [Dia (2023)][paper]. The number `x` should be convertible to
/// `f32`.
///
/// [paper]: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4487559
pub(crate) fn cdf_n01<T: num_traits::Float>(x: T) -> Option<T> {
    let abs_x = if let Some(x) = <f32 as NumCast>::from(x) {
        x
    } else {
        return None;
    };
    let x2 = abs_x.powi(2);
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
    if x > T::zero() {
        T::from(1. - one_minus_cdf_abs_x)
    } else {
        T::from(one_minus_cdf_abs_x)
    }
}

/// Computes the cumulative distribution function of a Student with degree of
/// freedom `n`.
///
/// Adapted from this [crate][page]
///
/// page: https://docs.rs/distrs/latest/src/distrs/students_t.rs.html#28-116
pub(crate) fn cdf_student<T: num_traits::Float + num_traits::FloatConst>(x: T, n: T) -> Option<T> {
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
        return cdf_n01(-y).map(|cdf_y| start + sign * cdf_y);
    }

    // make n mutable and int
    // n is int between 1 and 200 if made it here
    let mut n = <u8 as NumCast>::from(n).unwrap();
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
        print!("{:?}", cdf_student(100f32, 10.))
    }
}
