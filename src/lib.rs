mod binary_treatment;
mod data;

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
