use num_traits::NumCast;

/// Computes the cumulative distribution function of the standard normal using
/// the formula of [Dia (2023)][paper]. The number `x` should be convertible to
/// `f64`.
///
/// [paper]: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4487559
pub(crate) fn cdf_n01<T: num_traits::Float>(x: T) -> Option<T> {
    let abs_x = if let Some(x) = <f64 as NumCast>::from(x) {
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
