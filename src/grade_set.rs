//! Represent the effect of geometric algebra primitives over the grades of the
//! multivectors

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
#[derive(Debug, Eq, Clone)]
pub struct GradeSet(BitVec);

impl PartialEq for GradeSet {
    /// Allows equality between bitvecs of different lengths: tests if they are
    /// equal up to some trailing zeroes
    fn eq(&self, other: &Self) -> bool {
        let (small, big) = sort_by_len(&self.0, &other.0);
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
        GradeSet(BitVec::new())
        // An empty bitvec is just treated as a bitvec full of zeroes
    }

    /// The grade of a k-vector
    pub fn single(k: Grade) -> Self {
        let mut v = bitvec![0; k + 1];
        v.set(k, true);
        GradeSet(v)
    }

    /// Grades ranging from x to y (incl)
    pub fn range(x: Grade, y: Grade) -> Self {
        let mut v = bitvec![0; y + 1];
        for mut i in &mut v.as_mut_bitslice()[x..=y] {
            *i = true;
        }
        GradeSet(v)
    }

    /// Iterate over each grade present in the GradeSet
    pub fn iter(&self) -> impl Iterator<Item = Grade> + '_ {
        self.0.iter_ones()
    }

    /// Whether the GradeSet contains no grades. If so, the expression it is
    /// attached to can only be equal to zero
    pub fn is_empty(&self) -> bool {
        self.0.not_any()
    }

    /// Whether the GradeSet contains exactly one grade
    pub fn is_single(&self) -> bool {
        let mut iter = self.0.iter_ones();
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
        match self.0.get(k) {
            None => false,
            Some(x) => *x,
        }
    }

    /// Whether the GradeSet contains only the grade k
    pub fn is_just(&self, k: Grade) -> bool {
        self.contains(k) && self.is_single()
    }

    /// Add a grade to the set
    pub fn add_grade(mut self, k: Grade) -> Self {
        if k >= self.0.len() {
            self.0.resize(k+1, false);
        }
        self.0.set(k, true);
        self
    }

    /// Remove a grade from the set
    pub fn rm_grade(mut self, k: Grade) -> Self {
        if k < self.0.len() {
            self.0.set(k, false);
        }
        self
    }

    pub(crate) fn id(self) -> Self {
        self
    }

    /// Exponential. Is defined only for **single**-graded k-vectors
    pub fn exp(self) -> Self {
        assert!(
            self.is_single(),
            "exp cannot be used on a multivector, only a k-vector"
        );
        Self::single(0) + self
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

    /// Returns `a` and `b` restricted to their grades that will, when
    /// multiplied, affect those of `self`
    /// 
    /// NAIVE GREEDY IMPLEMENTATION FOR NOW. Though its usage is limited to the
    /// AST grade minimisation phase (downwards grade inference), therefore this
    /// method is not used when actually evaluating a GA expression
    pub fn grades_affecting_mul(&self, a: &Self, b: &Self) -> (Self, Self) {
        let mut ra = Self::empty();
        let mut rb = Self::empty();
        for ka in a.iter() {
            for kb in b.iter() {
                if (self.clone() & (Self::single(ka) * Self::single(kb))).0.any() {
                    // The product of ka and kb yields at least one grade that
                    // is in self
                    ra = ra + Self::single(ka);
                    rb = rb + Self::single(kb);
                }
            }
        }
        (ra, rb)
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
        let (small, big) = sort_by_len(self.0, rhs.0);
        GradeSet(big | small)
    }
}

/// GradeSet intersection: select the grades contained in both `self` and `rhs`.
/// You can think of it as grade projection (extraction) performed on `self`,
/// using `rhs` as the set of grades to keep
impl std::ops::BitAnd for GradeSet {
    type Output = GradeSet;
    fn bitand(self, rhs: Self) -> Self::Output {
        GradeSet(self.0 & rhs.0)
        // Note: the resulting bitvec will always have the same length as
        // `self.0`
    }
}

/// O(N^3) IMPLEMENTATION FOR NOW. Though its usage is limited to AST
/// construction (upwards grade inference), therefore this method is not used
/// when actually evaluating a GA expression
impl std::ops::Mul for GradeSet {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let (small, big) = sort_by_len(self.0, rhs.0);
        if small.len() == 0 {
            GradeSet(small)
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
            GradeSet(res)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const S: fn(Grade) -> GradeSet = GradeSet::single;
    const E: fn() -> GradeSet = GradeSet::empty;

    macro_rules! test_eqs {
        ($($test_name:ident : $a:expr => $b:expr),*) => {
            $(
                #[test]
                fn $test_name() {
                    assert_eq!($a, $b);
                }
            )*
        }
    }

    #[test]
    fn neq() {
        assert_ne!(S(3), S(4))
    }

    test_eqs!(
        add_self_id: S(3) + S(3) => S(3),
        add_empty_id: S(3) + E() => S(3),
        mul_empty_absorb: S(3) * E() => E(),
        mul_vecs: S(1) * S(1) => S(0) + S(2),
        mul_scal_id: S(40) * S(0) => S(40),
        mul_bivec_quadvec: S(2) * S(4) => S(2) + S(4) + S(6),
        mul_trivec_quadvec: S(3) * S(4) => S(1) + S(3) + S(5) + S(7),
        mul_trivec_pentavec: S(3) * S(5) => S(2) + S(4) + S(6) + S(8),
        mul_vec_rotor: S(1) * (S(0) + S(2)) => S(1) + S(3),
        range: GradeSet::range(4,6) => S(4) + S(5) + S(6),
        intersect: GradeSet::range(0,10) & GradeSet::range(4,6) => GradeSet::range(4,6),
        single_graded: (S(1) + S(1)).is_single() => true,
        not_single_graded: (S(1) + S(2)).is_single() => false,
        empty_not_single_graded: E().is_single() => false,
        iter_grades: (S(1) + S(22) + S(10)).iter().collect::<Vec<_>>() => vec![1,10,22],
        grades_affecting_mul:
          S(0).grades_affecting_mul(&(S(1) + S(0) + S(2) + S(10)), &(S(0) + S(2) + S(6)))
          => (S(0) + S(2), S(0) + S(2))
    );
}
