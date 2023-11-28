//! Describe a geometric algebra over which computations may be done
//!
//! An [`Algebra`] tells how many base vectors the underlying vector space has.
//! A [`MetricAlgebra`] additionally tells how these base vectors multiply under
//! the dot product

use super::{
    ast::{mv, GaExpr},
    grade_set::{Grade, GradeSet},
    graded::GradedDataMut,
};
use bitvec::prelude::*;
use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

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
            f.write_str(if *b { "1" } else { "0" })?;
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
    /// Prepare an Algebra so it is ready to be used for evaluating expressions
    pub fn from(alg: A) -> Self {
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
        let mut cur_blade = BigUint::zero();
        let max = alg.algebraic_dim();
        while cur_blade < max {
            let k = cur_blade.count_ones() as usize;
            let bv = BitVec::from_vec(cur_blade.to_bytes_le());
            basis_blades[k].push(BasisBlade(bv));
            cur_blade += 1u32;
        }
        ReadyAlgebra { alg, basis_blades }
    }

    pub fn basis_blade_from_indexes(&self, k: Grade, pos_in_grade: usize) -> &BasisBlade {
        &self.basis_blades[k][pos_in_grade]
    }

    pub fn indexes_from_basis_blade<'a>(
        &'a self,
        BasisBlade(b0): &'a BasisBlade,
    ) -> (Grade, usize) {
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
    pub fn ortho_basis_blades_gp(
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

impl<const D: usize> ReadyAlgebra<[f64; D]> {
    /// Return the base vectors of the algebra, whose number can be known
    /// statically here
    pub fn base_vec_exprs<T: GradedDataMut>(&self) -> [GaExpr<T>; D] {
        array_init::array_init(|i| {
            let mut v = T::init_null_mv(D, &GradeSet::single(1));
            v.grade_slice_mut(1)[i] = 1.0;
            mv(v)
        })
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

/// Given a word with exactly k bits at 1, generates the lexicographically next
/// bit permutation (the next smallest word with k bits at 1)
///
/// (we use the second version as it doesn't require a bit NOT operator, which
/// doesn't exist on variable-size BigInts)
pub fn next_bit_permutation(v: &BigInt) -> BigInt {
    if v == &BigInt::zero() {
        return v.clone();
    }
    let t: BigInt = (v.clone() | (v.clone() - 1)) + 1;
    t.clone() | ((((t.clone() & -t) / (v.clone() & -v)) >> 1) - 1)
}

/// Yield in ascending order all the `n`-bit words with exactly `k` bits equal
/// to 1
pub fn all_bit_permutations(n: usize, k: usize) -> impl Iterator<Item = BigUint> {
    let mut x = BigInt::zero();
    for i in 0..k {
        x.set_bit(i as u64, true)
    }
    (0..n_choose_k(n, k)).map(move |_| {
        let z = x.clone();
        x = next_bit_permutation(&z);
        z.to_biguint().unwrap()
    })
}

/// Directly computes the nth term of the lexicographic permutation suite of
/// `n`-bit words with exactly `k` bits equal to 1. `O(n)` time complexity.
///
/// `i` starts at `0`.
///
/// See first answer of https://math.stackexchange.com/questions/1304731/computing-the-n-textrmth-permutation-of-bits
/// 
/// See also
/// https://graphics.stanford.edu/%7Eseander/bithacks.html#NextBitPermutation
/// about a way to enumerate all bitfields in lexicographic order
pub fn idx_to_bitfield_permut(n: usize, mut k: usize, mut i: usize) -> BigUint {
    let mut res = BigUint::zero();
    for b in 1..=n {
        let z = n_choose_k(n - b, k);
        if i >= z {
            res.set_bit((n - b) as u64, true);
            i -= z;
            k -= 1; // We just set one bit to 1, so k-1 bits remain to be set to 1
        }
    }
    res
}

/// Does the reverse of [`idx_to_bitfield_permut`]: given a bitfield, finds its
/// index (starting at `0`) in the lexicographic permutation suite
pub fn bitfield_permut_to_idx(n: usize, v: &BigUint) -> usize {
    let mut res = 0;
    let mut k = v.count_ones() as usize;
    for b in 1..=n {
        let z = n_choose_k(n - b, k);
        if v.bit((n - b) as u64) {
            res += z;
            k -= 1;
        }
    }
    res
}

/// Binomial coefficient ("n choose k", or "k amongst n")
///
// #[memoize::memoize]
#[inline]
pub fn n_choose_k(n: usize, k: usize) -> usize {
    num_integer::binomial(n, k)
}

// /// Computes `sum(for i = 0 to k)(n_choose_i)`
// ///
// /// ie. 1 + n + n(n-1)/2 + n(n-1)(n-2)/(2*3) + n(n-1)(n-3)/(2*3*4) + ...
// pub(crate) fn sum_n_choose_ks(n: Grade, k: Grade) -> BigUint {
//     let mut coef = BigUint::one();
//     let mut sum = BigUint::one();
//     for i in 1..=k {
//         coef *= (coef.clone() * (n - i + 1)) / i;
//         sum += coef.clone();
//     }
//     sum
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_macros::*;

    simple_eqs! {
        n_choose_zero: n_choose_k(5, 0) => 1,
        zero_choose_zero: n_choose_k(0, 0) => 1,
        three_choose_two: n_choose_k(3, 2) => 3,
        bit_permuts_1: all_bit_permutations(4, 2).map(|i| i.to_u32_digits()[0]).collect::<Vec<_>>()
            => vec![0b11, 0b101, 0b110, 0b1001, 0b1010, 0b1100],
        bit_permuts_2: all_bit_permutations(10, 5).collect::<Vec<_>>()
            => (0..n_choose_k(10, 5)).map(|i| idx_to_bitfield_permut(10, 5, i)).collect::<Vec<_>>()
    }

    #[test]
    fn idx_bitfield_permut_roundtrip() {
        let indexes = (0..n_choose_k(10, 5)).collect::<Vec<_>>();
        let indexes2 = indexes.iter().map(|i| bitfield_permut_to_idx(10, &idx_to_bitfield_permut(10, 5, *i))).collect::<Vec<_>>();
        assert_eq!(indexes, indexes2);
    }

    #[test]
    fn bitfield_permut_idx_roundtrip() {
        let bitfields = (0..n_choose_k(9, 4)).map(|i| idx_to_bitfield_permut(9, 4, i)).collect::<Vec<_>>();
        let bitfields2 = bitfields.iter().map(|b| idx_to_bitfield_permut(9, 4, bitfield_permut_to_idx(9, b))).collect::<Vec<_>>();
        assert_eq!(bitfields, bitfields2);
    }
}
