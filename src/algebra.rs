//! Describe a geometric algebra over which computations may be done
//!
//! An [`Algebra`] tells how many base vectors the underlying vector space has.
//! A [`MetricAlgebra`] additionally tells how these base vectors multiply under
//! the dot product

use crate::grade_set::Grade;

/// The description of a geometric algebra over some vector space (over the
/// field of real numbers, [`f64`] here)
pub trait Algebra {
    /// The dimensionality of the underlying vector space
    fn vec_space_dim(&self) -> usize;

    /// The number of components of a general multivector in this algebra
    fn algebraic_dim(&self) -> usize {
        (2 as usize).pow(self.vec_space_dim() as u32)
    }

    /// Number of grades in this algebra
    fn num_grades(&self) -> usize {
        self.vec_space_dim() + 1
    }

    /// The number of components of `k`-vectors in this algebra
    fn grade_dim(&self, k: Grade) -> usize {
        n_choose_k(self.vec_space_dim(), k)
    }
}

/// Computes n! / (k! * (n-k)!)
pub(crate) const fn n_choose_k(n: Grade, k: Grade) -> usize {
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

/// A metric geometric algebra over some vector space
pub trait MetricAlgebra: Algebra {
    /// Give the dot product of two base vectors (identified by index)
    fn base_vec_dot(&self, v1: usize, v2: usize) -> f64;
}

/// Representation of an [`Algebra`] as an array of orthogonal base vector squares
impl<const D: usize> Algebra for [f64; D] {
    fn vec_space_dim(&self) -> usize {
        D
    }
}

/// Representation of a [`MetricAlgebra`] as an array of orthogonal base vector
/// squares
impl<const D: usize> MetricAlgebra for [f64; D] {
    fn base_vec_dot(&self, v1: usize, v2: usize) -> f64 {
        if v1 == v2 {
            self[v1]
        } else {
            0.0
        }
    }
}

/// A [`MetricAlgebra`] over a euclidean vector space of some dimension N. Thus
/// we have N orthogonal base vectors, each one squaring to 1
///
/// Therefore, `OrthoEuclidN(X)` represents the exact same algebra as `[1;
/// X]`, only their representation (and therefore, size) in memory is different
pub struct OrthoEuclidN(
    /// Dimension N of the underlying vector space
    pub usize,
);

impl Algebra for OrthoEuclidN {
    fn vec_space_dim(&self) -> usize {
        self.0
    }
}

impl MetricAlgebra for OrthoEuclidN {
    fn base_vec_dot(&self, v1: usize, v2: usize) -> f64 {
        if v1 == v2 {
            1.0
        } else {
            0.0
        }
    }
}
