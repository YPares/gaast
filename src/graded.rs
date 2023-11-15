use super::grade_set::GradeSet;
use std::{collections::HashMap, rc::Rc};

/// The trait for all objects that are graded
pub trait Graded {
    /// Get the GradeSet of the object
    fn grade_set(&self) -> GradeSet;
    /// Get a slice to the components of the k-vector part, given k
    fn get_grade_slice(&self, k: usize) -> &[f64];
}

macro_rules! Graded_blanket_impls {
    ($($ref:tt),*) => {
        $(impl<T: Graded> Graded for $ref<T> {
            fn grade_set(&self) -> GradeSet {
                (**self).grade_set()
            }
            fn get_grade_slice(&self, k: usize) -> &[f64] {
                (**self).get_grade_slice(k)
            }
        })*
    };
}
Graded_blanket_impls!(Rc, Box);

impl Graded for f64 {
    fn grade_set(&self) -> GradeSet {
        GradeSet::g(0)
    }
    fn get_grade_slice(&self, k: usize) -> &[f64] {
        assert!(k == 0);
        std::slice::from_ref(self)
    }
}

/// A Vec-based multivector representation
pub struct DynSizedMV {
    /// Contains the components of each grade
    contents: HashMap<usize, Vec<f64>>,
    /// Gives the grades contained in `contents`
    grade_set: GradeSet,
}

impl DynSizedMV {
    /// Create a [`DynSizedMV`] from its components, sorted by grade
    pub fn from_map_of_comps(m: HashMap<usize, Vec<f64>>) -> Self {
        Self {
            grade_set: m
                .keys()
                .fold(GradeSet::g_any(), |acc, &k| acc + GradeSet::g(k)),
            contents: m,
        }
    }
}

impl Graded for DynSizedMV {
    fn grade_set(&self) -> GradeSet {
        self.grade_set.clone()
    }
    fn get_grade_slice(&self, k: usize) -> &[f64] {
        self.contents[&k].as_ref()
    }
}

/// Computes n! / (k! * (n-k)!)
pub(crate) const fn n_choose_k(n: usize, k: usize) -> usize {
    let mut res = n;
    let mut i = n - k + 1;
    loop {
        if i >= n {
            break;
        }
        res *= i;
        i += 1;
    }
    let mut k_fac = k;
    let mut i = 2;
    loop {
        if i >= k {
            break;
        }
        k_fac *= i;
        i += 1;
    }
    res / k_fac
}
