//! Describe a geometric algebra over which computations may be done
//!
//! An [`Algebra`] tells how many base vectors the underlying vector space has.
//! A [`MetricAlgebra`] additionally tells how these base vectors multiply under
//! the dot product

use crate::GradeSet;

use super::grade_set::Grade;
use bitvec::prelude::*;

// # TYPES & TRAITS //

/// The description of a geometric algebra over some vector space (over the
/// field of real numbers, [`f64`] here)
pub trait Algebra {
    /// The dimensionality of the underlying vector space
    fn vec_space_dim(&self) -> usize;

    /// The maximal [`GradeSet`] available in this algebra
    fn full_grade_set(&self) -> GradeSet {
        (0..=self.vec_space_dim()).fold(GradeSet::empty(), |acc, k| acc.add_grade(k))
    }

    /// The number of basis blades of grade `k` (ie. the number of components of
    /// `k`-vectors in this algebra)
    fn grade_dim(&self, k: Grade) -> usize {
        n_choose_k(self.vec_space_dim(), k)
    }

    /// Given a grade and an index in a slice storing that grade's components,
    /// find the associated [`BasisBlade`]
    fn component_to_basis_blade(&self, Component { grade, index }: &Component) -> BasisBlade {
        BasisBlade(index_to_bitfield_permut(
            self.vec_space_dim(),
            *grade,
            *index,
        ))
    }

    /// Does the reverse: for some [`BasisBlade`], find its grade and its index
    /// in the slice of components for that grade
    fn basis_blade_to_component(&self, BasisBlade(b): &BasisBlade) -> Component {
        let grade = b.count_ones() as usize;
        let index = bitfield_permut_to_index(self.vec_space_dim(), grade, b);
        Component { grade, index }
    }
}

/// Given an algebra and a grade, yield all the basis blades of that grade in
/// this algebra
pub fn iter_basis_blades_of_grade(
    alg: &impl Algebra,
    grade: Grade,
) -> impl Iterator<Item = BasisBlade> + '_ {
    // TODO: reimplement that using
    // https://graphics.stanford.edu/%7Eseander/bithacks.html#NextBitPermutation
    (0..alg.grade_dim(grade))
        .map(move |index| alg.component_to_basis_blade(&Component { grade, index }))
}

/// A metric geometric algebra over some vector space
pub trait MetricAlgebra: Algebra {
    /// Give the dot product of two base vectors (identified by index). This
    /// function must always be symmetric (`base_vec_dot(u,v) ==
    /// base_vec_dot(v,u)`).
    ///
    /// Evaluated on each pair of base vectors, this gives the (symmetric) Gram
    /// matrix of the metric
    fn base_vec_dot(&self, v1: usize, v2: usize) -> f64;

    /// The geometric product of two basis blades, _restricted to cases where
    /// all basis vectors are orthogonal (diagonal Gram matrix)._ TODO:
    /// Diagonalize the Gram matrix in other cases
    fn ortho_basis_blades_gp(
        &self,
        BasisBlade(b1): &BasisBlade,
        BasisBlade(b2): &BasisBlade,
    ) -> (BasisBlade, f64) {
        let mut coef = canonical_reordering_sign(b1.clone(), b2);
        for bit in (b1.clone() & b2).iter_ones() {
            coef *= self.base_vec_dot(bit, bit);
        }
        (BasisBlade(b1.clone() ^ b2), coef)
    }
}

/// The position of a component (an item in a slice) in a graded object
#[derive(Debug)]
pub struct Component {
    pub grade: Grade,
    pub index: usize,
}

impl Component {
    pub fn from_tuple((grade, index): (Grade, usize)) -> Self {
        Self { grade, index }
    }
    pub fn to_tuple(self) -> (Grade, usize) {
        (self.grade, self.index)
    }
}

/// A basis blade in some algebra, of the form `1`, `eX`, `eX^eY`, `eX^eY^eZ`,
/// etc.
///
/// Internally, each basis blade is represented with a bitfield, and base
/// vectors are assumed to be orthogonal, so their outer and geometric products
/// are equal, and computable with a simple bitwise XOR. This follows the
/// recommendations of _Dorst, L., Fontijne, D., & Mann, S. (2010). Geometric
/// algebra for computer science: an object-oriented approach to geometry.
/// Elsevier_ (section 19.1, page 512 in the revised edition). Each bit
/// represents a basis _vector_ of the underlying vector space, so eg.:
///
/// - 0000000000000 : a blade with no basis vector. This is therefore the scalar
///   unique basis blade (`1`)
/// - 1000000000000 : first basis vector (`e1`)
/// - 0100000000000 : second basis vector (`e2`)
/// - 1110000000000 : "first" basis trivector (`e1 ^ e2 ^ e3`)
/// - 1011011000000 : the basis 5-vector `e1 ^ e3 ^ e4 ^ e6 ^ e7`
/// - etc.
///
/// The length of each bitfield is therefore just the dimension of the
/// underlying vector space.
///
/// The order of the vectors in each basis blade always follows the "natural"
/// order of the basis vectors it is composed of. That means eg. that the
/// trivector `e5 ^ e2 ^ e8` will be represented by the `e2 ^ e5 ^ e8` basis
/// trivector associated to a `-1` component.
///
/// There is just a variation compared to the book mentioned above: we start
/// from the leftmost bit (ie. the MSB represents `e1`), because gaast uses
/// variable-length bitvecs to represent basis blades, which behave more like
/// arrays than unsigned integers.
#[derive(Clone)]
pub struct BasisBlade(BitVec);

impl std::fmt::Debug for BasisBlade {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for b in &self.0 {
            f.write_str(if *b { "1" } else { "0" })?;
        }
        Ok(())
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
/// Therefore, `OrthoEuclidN(N)` represents the exact same algebra as `[1.0;
/// N]`, only their representation (and therefore, size) in memory is different
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

/// When two basis blades are multiplied with the geometric product, computes
/// whether we should add a minus sign or not (given anticommutativity of
/// orthogonal basis vectors)
pub(crate) fn canonical_reordering_sign(mut b1: BitVec, b2: &BitVec) -> f64 {
    let mut sum: i32 = 0;
    loop {
        b1.shift_left(1);
        sum += (b1.clone() & b2.clone()).count_ones() as i32;
        if b1.not_any() {
            break;
        }
    }
    (1 - (sum % 2) * 2) as f64 // -1 if sum is odd, +1 if even
}

/// Directly computes the nth term of the lexicographic permutation suite of
/// `n`-bit words with exactly `k` bits equal to 1. `O(n)` time complexity.
///
/// `i` starts at `0`.
///
/// See first answer of <https://math.stackexchange.com/questions/1304731/computing-the-n-textrmth-permutation-of-bits>
///
/// See also
/// <https://graphics.stanford.edu/%7Eseander/bithacks.html#NextBitPermutation>
/// about a way to enumerate all bitfields in lexicographic order
pub(crate) fn index_to_bitfield_permut(n: usize, mut k: usize, mut i: usize) -> BitVec {
    let mut res = bitvec![0; n];
    for b in 1..=n {
        let z = n_choose_k(n - b, k);
        if i >= z {
            res.set(n - b, true);
            i -= z;
            k -= 1; // We just set one bit to 1, so k-1 bits remain to be set to 1
        }
    }
    res
}

/// Does the reverse of [`idx_to_bitfield_permut`]: given a bitfield, finds its
/// index (starting at `0`) in the lexicographic permutation suite
pub(crate) fn bitfield_permut_to_index(n: usize, mut k: usize, v: &BitVec) -> usize {
    let mut res = 0;
    for b in 1..=n {
        let z = n_choose_k(n - b, k);
        if *v.get(n - b).as_deref().unwrap_or(&false) {
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
pub(crate) fn n_choose_k(n: usize, k: usize) -> usize {
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
        three_choose_two: n_choose_k(3, 2) => 3
    }

    #[test]
    fn idx_bitfield_permut_roundtrip() {
        let indexes = (0..n_choose_k(10, 5)).collect::<Vec<_>>();
        let indexes2 = indexes
            .iter()
            .map(|i| bitfield_permut_to_index(10, 5, &index_to_bitfield_permut(10, 5, *i)))
            .collect::<Vec<_>>();
        assert_eq!(indexes, indexes2);
    }

    #[test]
    fn bitfield_permut_idx_roundtrip() {
        let bitfields = (0..n_choose_k(9, 4))
            .map(|i| index_to_bitfield_permut(9, 4, i))
            .collect::<Vec<_>>();
        let bitfields2 = bitfields
            .iter()
            .map(|b| index_to_bitfield_permut(9, 4, bitfield_permut_to_index(9, 4, b)))
            .collect::<Vec<_>>();
        assert_eq!(bitfields, bitfields2);
    }
}
