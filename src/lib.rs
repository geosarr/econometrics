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
        m = m + x;
        n = n + T::one();
    }
    if n > T::zero() {
        Some((m / n, n))
    } else {
        None
    }
}

/// Calcule la fonction de répartition
/// de la loi normale centrée réduite en utilisant
/// la formule proposée par [Dia (2023)][lien].
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
#[cfg(test)]
mod tests {
    use crate::{binary_treatment::BinaryTreatment, data::lalonde::*};
    use ndarray::*;
    #[test]
    fn binary_treatment() {
        let income = Array::from_iter(RE78);
        let treatment = Array::from_iter(TREAT.iter().map(|v| *v as f32));
        let mut binary_treatment: BinaryTreatment<f32> = BinaryTreatment::new();
        binary_treatment.fit(&treatment, &income);
        println!("{:#?}", binary_treatment);
    }
}
