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
