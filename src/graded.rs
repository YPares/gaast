//! The [`Graded`] trait, for types that contain or allow to read graded data
//! (raw multivectors)

use super::{algebra::n_choose_k, grade_set::*};
use std::{collections::HashMap, rc::Rc};

/// The trait for all objects that are graded, ie. from which we can extract an
/// array of components corresponding to some grade
pub trait Graded {
    /// Get the GradeSet of the object
    fn grade_set(&self) -> GradeSet;
    /// Get a slice to the components of the k-vector part, given k. The length
    /// of the slice must exactly correspond to what is expected for that grade
    fn grade_slice(&self, k: Grade) -> &[f64];
}

/// The trait for all objects that are graded and writeable
pub trait GradedMut: Graded {
    /// Create a multivector that contains all the components to hold data of
    /// given grade for a given vector space dimension. All components are
    /// initialized to zero
    fn init_null_mv(vec_space_dim: usize, gs: &GradeSet) -> Self;
    /// Get a mutable slice to the components of the k-vector part, given k. The
    /// length of the slice must exactly correspond to what is expected for that
    /// grade
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64];
    /// Multiply all the components of a given grade by -1
    fn negate_grade(&mut self, k: Grade) {
        for x in self.grade_slice_mut(k) {
            *x = -*x;
        }
    }
}

macro_rules! Graded_blanket_impls {
    ($($ref:tt),*) => {
        $(impl<T: Graded> Graded for $ref<T> {
            fn grade_set(&self) -> GradeSet {
                (**self).grade_set()
            }
            fn grade_slice(&self, k: Grade) -> &[f64] {
                (**self).grade_slice(k)
            }
        })*
    };
}
Graded_blanket_impls!(Rc, Box);

impl<T: GradedMut> GradedMut for Box<T> {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        (**self).grade_slice_mut(k)
    }
    fn init_null_mv(dim: usize, gs: &GradeSet) -> Self {
        Box::new(T::init_null_mv(dim, &gs))
    }
}

impl Graded for f64 {
    fn grade_set(&self) -> GradeSet {
        GradeSet::g(0)
    }
    fn grade_slice(&self, k: Grade) -> &[f64] {
        assert!(k == 0);
        std::slice::from_ref(self)
    }
}

impl GradedMut for f64 {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        assert!(k == 0);
        std::slice::from_mut(self)
    }
    fn init_null_mv(_dim: usize, gs: &GradeSet) -> Self {
        assert!(gs == &GradeSet::g(0));
        0.0
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
                .fold(GradeSet::g_none(), |acc, &k| acc + GradeSet::g(k)),
            contents: m,
        }
    }
}

impl Graded for DynSizedMV {
    fn grade_set(&self) -> GradeSet {
        self.grade_set.clone()
    }
    fn grade_slice(&self, k: Grade) -> &[f64] {
        self.contents[&k].as_ref()
    }
}

impl GradedMut for DynSizedMV {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        self.contents.get_mut(&k).unwrap()
    }
    fn init_null_mv(dim: usize, gs: &GradeSet) -> Self {
        let mut res = DynSizedMV {
            grade_set: gs.clone(),
            contents: HashMap::new(),
        };
        for k in res.grade_set.iter_grades() {
            res.contents.insert(k, vec![0.0; n_choose_k(dim, k)]);
        }
        res
    }
}
