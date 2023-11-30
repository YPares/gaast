//! Represent the effect of geometric algebra primitives over the grades of the
//! multivectors

use std::fmt::Debug;

use bitvec::prelude::*;

/// Represents the set of grades that can be contained in some multivector. You
/// can think of it as the "type" of a multivector, which supports the same
/// operation actual multivectors do, by mirroring the effects those operations
/// have on multivectors grades. Therefore, GradeSets can be added or multiplied
/// together, and this will yield the GradeSet of the multivector obtained by
/// addition/geometric multiplication of two multivectors. Other GA primitives
/// are provided (like exp & log), with some limitations indicated in the
/// methods' documentation. This allows to perform _grade inference_ on
/// multivector expressions _without_ having to actually compute them.
///
/// Implementation notes: This is implemented using a heap-allocated BitVec (not
/// a "big enough" statically-sized type like u64) in order to be truly agnostic
/// of the underlying vec space dimensionality N, and therefore of the final
/// number of different grades (N+1) of the GA it creates. The penalty it incurs
/// for vector spaces of lower dimensions (which will generate a GA with fewer than 64 grades)
/// is not evaluated yet
#[derive(Eq, Clone, Hash)]
pub struct GradeSet {
    bv: BitVec,
}

impl Debug for GradeSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.iter().collect::<Vec<_>>().fmt(f)
    }
}

impl PartialEq for GradeSet {
    /// Allows equality between bitvecs of different lengths: tests if they are
    /// equal up to some trailing zeroes
    fn eq(&self, other: &Self) -> bool {
        let (small, big) = sort_by_len(&self.bv, &other.bv);
        big[0..small.len()] == &small[..] && big[small.len()..].not_any()
    }
}

/// Grades are just regular positive integers. Alias introduced for clarity
pub type Grade = usize;

impl GradeSet {
    /// The `GradeSet` of _zero_. In GA, `0` is _polymorphic_: it's a scalar or
    /// any k-vector or linear combination of those at the same time. Whatever
    /// the grades contained in a multivector, it can _always_ be zero. And when
    /// a multivector has _no_ grades, then it can _only_ be zero.
    pub fn empty() -> Self {
        GradeSet { bv: BitVec::new() }
        // An empty bitvec is just treated as a bitvec full of zeroes
    }

    fn from_usize(k: Grade) -> Self {
        let mut bv = bitvec![0; k + 1];
        bv.set(k, true);
        GradeSet { bv }
    }

    // The grade of a k-vector. Can be negative, if so, the resulting GradeSet
    // is empty
    pub fn single(k: i64) -> Self {
        if k < 0 {
            Self::empty()
        } else {
            Self::from_usize(k as Grade)
        }
    }

    /// Grades ranging from x to y (incl)
    pub fn range(x: Grade, y: Grade) -> Self {
        let mut bv = bitvec![0; y + 1];
        for mut i in &mut bv.as_mut_bitslice()[x..=y] {
            *i = true;
        }
        GradeSet { bv }
    }

    /// GradeSet intersection: select the grades contained in both `self` and `rhs`.
    /// You can think of it as grade projection (extraction) performed on `self`,
    /// using `rhs` as the set of grades to keep
    pub fn intersection(self, rhs: Self) -> Self {
        GradeSet {
            bv: self.bv & rhs.bv,
        }
        // Note: the resulting bitvec will always have the same length as
        // `self.bv`
    }

    /// Iterate over each grade present in the GradeSet
    pub fn iter(&self) -> impl Iterator<Item = Grade> + Clone + '_ {
        self.bv.iter_ones()
    }

    /// Whether the GradeSet contains only even-graded parts
    pub fn all_even(&self) -> bool {
        self.iter().all(|x| x % 2 == 0)
    }

    /// Whether the GradeSet contains only odd-graded parts
    pub fn all_odd(&self) -> bool {
        self.iter().all(|x| x % 2 == 1)
    }

    /// Whether the GradeSet can represent a versor, ie. an object factorizable
    /// into a sequence of geometric products of `k` invertible vectors. An even
    /// (resp. odd) `k` will always result in a multivector with only even
    /// (resp. odd) grades, therefore versors can only be fully even-graded or
    /// fully odd-graded.
    ///
    /// We say _can_ because the parity condition above is necessary but not
    /// sufficient a condition so a multivector is a versor, factorizability is
    /// necessary too. This factorizability means that a versor will always
    /// square to a scalar.
    pub fn can_be_versor(&self) -> bool {
        self.all_even() || self.all_odd()
    }

    /// Whether the GradeSet contains no grades. If so, the expression it is
    /// attached to can only be equal to zero
    pub fn is_empty(&self) -> bool {
        self.bv.not_any()
    }

    /// Whether the GradeSet contains exactly one grade
    pub fn is_single(&self) -> bool {
        let mut iter = self.bv.iter_ones();
        if let None = iter.next() {
            return false;
        }
        for _ in iter {
            return false;
        }
        true
    }

    /// Whether the GradeSet contains the grade k
    pub fn contains(&self, k: Grade) -> bool {
        match self.bv.get(k) {
            None => false,
            Some(x) => *x,
        }
    }

    /// Whether this GradeSet fully contains another
    pub fn includes(self, other: GradeSet) -> bool {
        (self.bv.clone() | other.bv) == self.bv
    }

    /// Whether the GradeSet contains only the grade k
    pub fn is_just(&self, k: Grade) -> bool {
        self.contains(k) && self.is_single()
    }

    /// Add a grade to the set
    pub fn add_grade(mut self, k: Grade) -> Self {
        if k >= self.bv.len() {
            self.bv.resize(k + 1, false);
        }
        self.bv.set(k, true);
        self
    }

    /// Remove a grade from the set
    pub fn rm_grade(mut self, k: Grade) -> Self {
        if k < self.bv.len() {
            self.bv.set(k, false);
        }
        self
    }

    /// Exponential. Is defined only for **single**-graded k-vectors
    pub fn exp(self) -> Self {
        assert!(
            self.is_single(),
            "exp cannot be used on a multivector, only a k-vector"
        );
        Self::from_usize(0) + self
    }

    /// Logarithm. Is defined only for multivectors of the form \<A\>_0 + \<A\>_k
    pub fn log(self) -> Self {
        let other_grade = self.rm_grade(0);
        assert!(
            other_grade.is_single(),
            "log can only be used on multivectors of the form <A>_0 + <A>_k"
        );
        other_grade
    }

    /// Keeps only the lowest grade contained in the set
    pub fn min(&self) -> Option<Grade> {
        self.bv.first_one()
    }

    /// Keeps only the highest grade contained in the set
    pub fn max(&self) -> Option<Grade> {
        self.bv.last_one()
    }

    /// For each grade k in self and each grade g in other, yield (k,g)
    pub fn iter_cartesian_product<'a>(
        &'a self,
        other: &'a Self,
    ) -> impl Iterator<Item = (Grade, Grade)> + Clone + 'a {
        self.iter()
            .flat_map(move |kl| other.iter().map(move |kr| (kl, kr)))
    }

    // Apply a function to each pair of the cartesian product of 2 grade sets
    pub fn sum_map_cartesian_product(
        &self,
        other: &Self,
        mut f: impl FnMut(i64, i64) -> GradeSet,
    ) -> GradeSet {
        self.iter_cartesian_product(&other)
            .map(|(k1, k2)| f(k1 as i64, k2 as i64))
            .reduce(|a, b| a + b)
            .unwrap_or(GradeSet::empty())
    }

    /// Using `self` as a geometric product result, yields all the pairs of
    /// grades in `left` and `right` that will, when multiplied, contribute to
    /// at least one of the grades contained in self.
    ///
    /// Each element yielded is of the form (resulting_part_of_self,
    /// grade_in_left, grade_in_right)
    ///
    /// GREEDY O(N^2) IMPLEMENTATION FOR NOW
    pub fn iter_contribs_to_gp<'a>(
        &'a self,
        left: &'a Self,
        right: &'a Self,
    ) -> impl Iterator<Item = (GradeSet, Grade, Grade)> + Clone + 'a {
        left.iter_cartesian_product(right).filter_map(|(kl, kr)| {
            let contribs = self
                .clone()
                .intersection(Self::from_usize(kl) * Self::from_usize(kr));
            if !contribs.is_empty() {
                Some((contribs, kl, kr))
            } else {
                None
            }
        })
    }

    /// Uses [`Self::iter_contribs_to_gp`] to collect
    /// (contributing_grades_in_left, contributing_grades_in_right)
    pub fn parts_contributing_to_gp(&self, left: &Self, right: &Self) -> (Self, Self) {
        let mut filtered_left = GradeSet::empty();
        let mut filtered_right = GradeSet::empty();
        for (_, k_left, k_right) in self.iter_contribs_to_gp(left, right) {
            filtered_left = filtered_left.add_grade(k_left);
            filtered_right = filtered_right.add_grade(k_right);
        }
        (filtered_left, filtered_right)
    }

    #[inline]
    pub(crate) fn id(self) -> Self {
        self
    }
}

fn sort_by_len<T>(v1: T, v2: T) -> (T, T)
where
    T: std::borrow::Borrow<BitVec>,
{
    if v1.borrow().len() <= v2.borrow().len() {
        (v1, v2)
    } else {
        (v2, v1)
    }
}

impl std::ops::Add for GradeSet {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (small, big) = sort_by_len(self.bv, rhs.bv);
        GradeSet { bv: big | small }
    }
}

/// Outer product
impl std::ops::BitXor for GradeSet {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self::Output {
        self.sum_map_cartesian_product(&rhs, |k1, k2| GradeSet::single(k1 + k2))
    }
}

/// Inner product
impl std::ops::BitAnd for GradeSet {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self::Output {
        self.sum_map_cartesian_product(&rhs, |k1, k2| {
            if k1 == 0 || k2 == 0 {
                GradeSet::empty()
            } else {
                GradeSet::single((k1 - k2).abs())
            }
        })
    }
}

/// Left contraction
impl std::ops::Shl for GradeSet {
    type Output = Self;
    fn shl(self, rhs: Self) -> Self::Output {
        self.sum_map_cartesian_product(&rhs, |k1, k2| GradeSet::single(k2 - k1))
    }
}

/// Right contraction
impl std::ops::Shr for GradeSet {
    type Output = Self;
    fn shr(self, rhs: Self) -> Self::Output {
        self.sum_map_cartesian_product(&rhs, |k1, k2| GradeSet::single(k1 - k2))
    }
}

/// O(N^3) IMPLEMENTATION FOR NOW. Though its usage is limited to AST
/// construction (upwards grade inference), therefore this method is not used
/// when actually evaluating a GA expression
impl std::ops::Mul for GradeSet {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let (small, big) = sort_by_len(self.bv, rhs.bv);
        if small.len() == 0 {
            GradeSet { bv: small }
        } else {
            let mut res = bitvec![0; big.len() + small.len() - 1];
            for r in 0..res.len() {
                for i in 0..small.len() {
                    for j in 0..big.len() {
                        let m = (i as i32 - j as i32).abs();
                        if i + j >= r && m <= r as i32 && m % 2 == r as i32 % 2 {
                            let mut x = res.get_mut(r).unwrap();
                            *x = *x || (small[i] && big[j]);
                        }
                    }
                }
            }
            GradeSet { bv: res }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_macros::*;

    const E: fn() -> GradeSet = GradeSet::empty;
    const S: fn(i64) -> GradeSet = GradeSet::single;
    const R: fn(usize, usize) -> GradeSet = GradeSet::range;

    #[test]
    fn neq() {
        assert_ne!(S(3), S(4))
    }

    simple_eqs! {
        neg_grade_is_empty: S(-1) => E(),
        add_self_id: S(3) + S(3) => S(3),
        add_empty_id: S(3) + E() => S(3),
        mul_empty_absorb: S(3) * E() => E(),
        mul_scal_id: S(40) * S(0) => S(40),
        mul_vecs: S(1) * S(1) => S(0) + S(2),
        mul_bivec_quadvec: S(2) * S(4) => S(2) + S(4) + S(6),
        mul_trivec_quadvec: S(3) * S(4) => S(1) + S(3) + S(5) + S(7),
        mul_trivec_pentavec: S(3) * S(5) => S(2) + S(4) + S(6) + S(8),
        mul_vec_rotor: S(1) * (S(0) + S(2)) => S(1) + S(3),
        range: R(4,7) => S(4) + S(5) + S(6) + S(7),
        intersect: R(0,10).intersection(R(4,30)) => R(4,10),
        single_graded: (S(1) + S(1)).is_single() => true,
        not_single_graded: (S(1) + S(2)).is_single() => false,
        empty_not_single_graded: E().is_single() => false,
        empty_intersection_is_empty: S(0).intersection(S(1)).is_empty() => true,
        iter_grades: (S(1) + S(22) + S(10)).iter().collect::<Vec<_>>() => vec![1,10,22],
        parts_contributing_to_gp:
          S(0).parts_contributing_to_gp( &(S(1) + S(0) + S(2) + S(10))
                                        , &(S(0) + S(2) + S(6)) )
          => (S(0) + S(2), S(0) + S(2)),
        outer: (S(1) + S(4)) ^ (S(1) + S(5)) => S(2) + S(6) + S(5) + S(9),
        inner: (S(1) + S(4)) & (S(0) + S(1) + S(5)) => S(0) + S(4) + S(3) + S(1)
    }
}
