//! The [`Graded`] trait, for types that contain graded data or act as
//! placeholders for graded data, and its subtypes to read and write grade
//! components (slices)

use super::{algebra::n_choose_k, grade_set::*};
use std::{collections::HashMap, rc::Rc, sync::Arc};

/// The trait for all objects from which we can query grades
pub trait Graded {
    /// Either directly an owned GradeSet (for Graded objects which do not
    /// directly hold on to a GradeSet and must reconstruct it on the fly) or a
    /// ref to a GradeSet (for Graded objects which do internally already have
    /// some GradeSet which can be reused read-only)
    type GradeSetOrRef<'a>: std::borrow::Borrow<GradeSet> + 'a
    where
        Self: 'a;
    /// Get the GradeSet of the object
    fn grade_set<'a>(&'a self) -> Self::GradeSetOrRef<'a>;
}

/// The identity implementation. A [`GradeSet`] just returns a reference to
/// itself
impl Graded for GradeSet {
    type GradeSetOrRef<'a> = &'a Self;
    fn grade_set<'a>(&'a self) -> Self::GradeSetOrRef<'a> {
        self
    }
}

/// The trait for all objects that are graded and contain readable data
/// (components) associated to each grade
pub trait GradedInput: Graded {
    /// Get a slice to the components of the k-vector part, given k. The length
    /// of the slice must exactly correspond to what is expected for that grade
    fn grade_slice(&self, k: Grade) -> &[f64];
}

/// The trait for all objects that are graded and contain writeable data
/// (components) associated to each grade
pub trait GradedOutput: Graded {
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
        $(
            impl<T: Graded> Graded for $ref<T> {
                type GradeSetOrRef<'a> = T::GradeSetOrRef<'a> where T: 'a;
                fn grade_set<'a>(&'a self) -> Self::GradeSetOrRef<'a> {
                    (**self).grade_set()
                }
            }
            impl<T: GradedInput> GradedInput for $ref<T> {
                fn grade_slice(&self, k: Grade) -> &[f64] {
                    (**self).grade_slice(k)
                }
            }
        )*
    };
}
Graded_blanket_impls!(Box, Rc, Arc);

impl<T: GradedOutput> GradedOutput for Box<T> {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        (**self).grade_slice_mut(k)
    }
    fn init_null_mv(dim: usize, gs: &GradeSet) -> Self {
        Box::new(T::init_null_mv(dim, &gs))
    }
    fn negate_grade(&mut self, k: Grade) {
        (**self).negate_grade(k);
    }
}

impl Graded for f64 {
    type GradeSetOrRef<'a> = GradeSet;
    fn grade_set(&self) -> Self::GradeSetOrRef<'_> {
        GradeSet::single(0)
    }
}
impl GradedInput for f64 {
    fn grade_slice(&self, k: Grade) -> &[f64] {
        assert!(k == 0);
        std::slice::from_ref(self)
    }
}
impl GradedOutput for f64 {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        assert!(k == 0);
        std::slice::from_mut(self)
    }
    fn init_null_mv(_dim: usize, gs: &GradeSet) -> Self {
        assert!(gs == &GradeSet::single(0));
        0.0
    }
}

/// A HashMap of Vecs raw multivector representation
///
/// Convenient for testing, but neither efficient nor compact in memory.
/// Notably, just storing a scalar part (ie. just a [`f64`]) requires 2
/// allocations (singleton hashmap and singleton vec)
#[derive(Debug, PartialEq)]
pub struct GradeMapMV(pub HashMap<Grade, Vec<f64>>);

impl Graded for GradeMapMV {
    type GradeSetOrRef<'a> = GradeSet;
    fn grade_set(&self) -> Self::GradeSetOrRef<'_> {
        self.0
            .keys()
            .fold(GradeSet::empty(), |acc, &k| acc.add_grade(k))
    }
}
impl GradedInput for GradeMapMV {
    fn grade_slice(&self, k: Grade) -> &[f64] {
        &self.0[&k]
    }
}
impl GradedOutput for GradeMapMV {
    fn grade_slice_mut(&mut self, k: Grade) -> &mut [f64] {
        self.0.get_mut(&k).unwrap()
    }
    fn init_null_mv(dim: usize, gs: &GradeSet) -> Self {
        let mut m = HashMap::new();
        for k in gs.iter() {
            m.insert(k, vec![0.0; n_choose_k(dim, k)]);
        }
        Self(m)
    }
}

#[macro_export]
/// A macro to create a [`GradeMapMV`], for instance in Cl(3):
///
/// `grade_map_mv!(0 => 1, 1 => 1 1 0, 3 => 1)`
///
/// this has a scalar part, a vector part, and a trivector part
macro_rules! grade_map_mv {
    ($($grade:expr => $($x:expr)+),+) => {{
        let mut m = std::collections::HashMap::new();
        $({
            let mut v = Vec::new();
            $(
                v.push($x as f64);
            )+
            m.insert($grade, v);
        })+
        crate::graded::GradeMapMV(m)
    }}
}
pub use grade_map_mv;

#[cfg(test)]
mod tests {
    use super::grade_map_mv;
    use crate::test_macros::*;

    simple_eqs! {
        hash_map_mv_eq: grade_map_mv!(1 => 1 2 3) => grade_map_mv!(1 => 1 2 3)
    }
}
