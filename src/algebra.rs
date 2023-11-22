//! Describe a geometric algebra over which computations may be done
//!
//! An [`Algebra`] tells how many base vectors the underlying vector space has.
//! A [`MetricAlgebra`] additionally tells how these base vectors multiply under
//! the dot product

use crate::grade_set::Grade;
use bitvec::prelude::*;
use num_bigint::BigUint;

// # TYPES & TRAITS //

/// The description of a geometric algebra over some vector space (over the
/// field of real numbers, [`f64`] here)
pub trait Algebra {
    /// The dimensionality of the underlying vector space
    fn vec_space_dim(&self) -> usize;

    /// Number of grades in this algebra
    fn num_grades(&self) -> usize {
        self.vec_space_dim() + 1
    }

    /// The number of basis blades in this algebra. Corresponds to the max
    /// number of f64 components of a general multivector in this algebra
    fn algebraic_dim(&self) -> BigUint {
        BigUint::from(2u32).pow(self.vec_space_dim() as u32)
    }

    /// The number of basis blades of grade `k` (ie. the number of components of
    /// `k`-vectors in this algebra)
    fn grade_dim(&self, k: Grade) -> usize {
        n_choose_k(self.vec_space_dim(), k)
    }
}

/// A metric geometric algebra over some vector space
pub trait MetricAlgebra: Algebra {
    /// Give the dot product of two base vectors (identified by index)
    fn base_vec_dot(&self, v1: usize, v2: usize) -> f64;
}

/// A basis blade in some algebra, of the form `1`, `eX`, `eX^eY`, `eX^eY^eZ`,
/// etc.
pub struct BasisBlade(BitVec<u8>);

impl std::fmt::Debug for BasisBlade {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for b in &self.0 {
            f.write_str(if *b {"1"} else {"0"})?;
        }
        Ok(())
    }
}

/// An algebra augmented with some information that we want to compute just once
/// and keep around
#[derive(Debug)]
pub struct ReadyAlgebra<A> {
    alg: A,
    /// Each basis blade used by this algebra associated to its index.
    /// Invariant: elements at position `n` in the outer vec should always have
    /// exactly `n` bits to 1. (So at pos 0, the only possibility in the
    /// singleton vec containing an empty bitfield)
    basis_blades: Vec<Vec<BasisBlade>>,
}

/// So the algebra methods can still be called through ReadyAlgebra
impl<A> std::ops::Deref for ReadyAlgebra<A> {
    type Target = A;
    fn deref(&self) -> &Self::Target {
        &self.alg
    }
}

impl<A: Algebra> ReadyAlgebra<A> {
    pub fn prepare(alg: A) -> Self {
        // Give a name (as a bitfield) to each basis blade of the algebra. This
        // follows the convention described in `Dorst, L., Fontijne, D., & Mann,
        // S. (2010). Geometric algebra for computer science: an object-oriented
        // approach to geometry. Elsevier`. Each bit represents a base vector of
        // the underlying vector space, so eg.:
        // - 0000000000000 is the scalar unique base blade `1`
        // - 1000000000000 is the first vector (e1)
        // - 0100000000000 is the second vector (e2)
        // - 1110000000000 is the trivector e1 ^ e2 ^ e3, etc.
        //
        // The length of each bitfield being the dim of the underlying vector
        // space. Note that bitvec considers the first bit to be the leftmost
        // (accesses are done like a vector, from left to right, regardless of
        // the underlying storage being MSB or LSB)
        let mut basis_blades: Vec<_> = (0..alg.num_grades())
            .map(|k| Vec::with_capacity(alg.grade_dim(k)))
            .collect();
        let mut cur_blade = BigUint::from(0u32);
        let max = alg.algebraic_dim();
        while cur_blade < max {
            let k = cur_blade.count_ones() as usize;
            let bv = BitVec::from_vec(cur_blade.to_bytes_le());
            basis_blades[k].push(BasisBlade(bv));
            cur_blade += 1u32;
        }
        ReadyAlgebra { alg, basis_blades }
    }

    pub fn basis_blade_from_coord(&self, k: Grade, pos_in_grade: usize) -> &BasisBlade {
        &self.basis_blades[k][pos_in_grade]
    }

    pub fn coord_from_basis_blade<'a>(&'a self, BasisBlade(b0): &'a BasisBlade) -> (Grade, usize) {
        let k = b0.count_ones() as usize;
        for (i, BasisBlade(b)) in self.basis_blades[k].iter().enumerate() {
            if b0 == b {
                return (k, i); // TODO. This is crap. To be improved
            }
        }
        panic!("coord_from_basis_blade: BasisBlade {b0} not found");
    }
}

fn canonical_reordering_sign(mut b1: BitVec<u8>, b2: &BitVec<u8>) -> f64 {
    b1.shift_left(1);
    let mut sum = 0;
    while b1.any() {
        sum += (b1.clone() & b2.clone()).count_ones();
        b1.shift_left(1);
    }
    if sum & 1 == 0 {
        1.0
    } else {
        -1.0
    }
}

impl<A: MetricAlgebra> ReadyAlgebra<A> {
    /// The geometric product of two basis blades restricted to cases where all
    /// basis vectors are orthogonal (diagonal metric). TODO: Diagonalize the
    /// Gram matrix when metric isn't diagonal.
    pub fn ortho_basis_blades_mul(
        &self,
        BasisBlade(b1): &BasisBlade,
        BasisBlade(b2): &BasisBlade,
    ) -> (BasisBlade, f64) {
        let mut coef = canonical_reordering_sign(b1.clone(), b2);
        for bit in (b1.clone() & b2).iter_ones() {
            coef *= self.alg.base_vec_dot(bit, bit);
        }
        (BasisBlade(b1.clone() ^ b2), coef)
    }
}

// # IMPLEMENTATIONS //

/// Representation of an [`Algebra`] as an array of orthogonal base vector squares
impl<const D: usize> Algebra for [f64; D] {
    fn vec_space_dim(&self) -> usize {
        D
    }
}

/// Representation of a [`MetricAlgebra`] as an array of orthogonal base vector
/// squares
impl<const D: usize> MetricAlgebra for [f64; D] {
    /// Gram Matrix
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
#[derive(Debug)]
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

// # UTILITY FUNCTIONS //

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
    let mut k_fac = if k == 0 { 1 } else { k };
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

// /// Computes `sum(for i = 0 to k)(n_choose_i)`
// ///
// /// ie. 1 + n + n(n-1)/2 + n(n-1)(n-2)/(2*3) + n(n-1)(n-3)/(2*3*4) + ...
// pub(crate) fn sum_n_choose_ks(n: Grade, k: Grade) -> BigUint {
//     let mut coef = BigUint::from(1u32);
//     let mut sum = BigUint::from(1u32);
//     for i in 1..=k {
//         coef *= (coef.clone() * (n - i + 1)) / i;
//         sum += coef.clone();
//     }
//     sum
// }

#[cfg(test)]
mod tests {
    use bitvec::prelude::*;
    use num_bigint::BigUint;

    #[test]
    fn bla() {
        assert!(BigUint::from(1u32) & BigUint::from(8u32) == BigUint::from(0u32));
    }
    #[test]
    fn plop() {
        let mut b = BitVec::<u8>::from_vec(BigUint::from(259u16).to_bytes_le());
        let c = bitvec![u8, Lsb0; 1; 8];
        b.shift_left(1);
        let x = b.get(0);
        println!("B: {x:?}");
        assert_eq!(b, c)
    }
}
