//! Represent the effect of geometric algebra primitives over the grades of the
//! multivectors

use bitvec::prelude::*;

/// Represents the set of grades that can be contained in some multivector. Can
/// be added or multiplied together, and will yield the grades of the
/// multivector obtained by addition/geometric multiplication of two
/// multivectors.
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

impl GradeSet {
    /// The grade of zero. (In GA, zero is polymorphic: it's a scalar, a vector,
    /// a bivector etc. at the same time)
    pub fn g_none() -> Self {
        GradeSet(BitVec::new())
        // A empty bitvec is just treated as a bitvec full of zeroes
    }

    /// The grade of a k-vector
    pub fn g(k: usize) -> Self {
        let mut v = bitvec![0; k + 1];
        v.set(k, true);
        GradeSet(v)
    }

    /// Grades ranging from x to y (incl)
    pub fn range(x: usize, y: usize) -> Self {
        let mut v = bitvec![0; y + 1];
        for mut i in &mut v.as_mut_bitslice()[x..=y] {
            *i = true;
        }
        GradeSet(v)
    }

    /// Grade projection: select part of the grades contained in self, using
    /// another GradeSet as a selector
    pub fn prj(self, grades: Self) -> Self {
        GradeSet(self.0 & grades.0)
    }

    /// Iterate over each grade present in the GradeSet
    pub fn iter_grades(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter_ones()
    }

    pub fn is_none(&self) -> bool {
        self.0.not_any()
    }

    /// Whether the GradeSet contains exactly one grade
    pub fn is_single_graded(&self) -> bool {
        let mut iter = self.0.iter_ones();
        if let None = iter.next() {
            return false;
        }
        for _ in iter {
            return false;
        }
        true
    }

    /// Remove a grade from the set
    pub fn rm_grade(mut self, k: usize) -> Self {
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
            self.is_single_graded(),
            "exp cannot be used on a multivector, only k-vector"
        );
        GradeSet::g(0) + self
    }

    /// Logarithm. Is defined only for multivectors of the form \<A\>_0 + \<A\>_k
    pub fn log(self) -> Self {
        let other_grade = self.rm_grade(0);
        assert!(
            other_grade.is_single_graded(),
            "log can only be used on multivectors of the form <A>_0 + <A>_k"
        );
        other_grade
    }

    /// Restricts `a` and `b` to the grades that will, when multiplied, affect
    /// those of `self`
    pub fn grades_affecting_mul(&self, a: &Self, b: &Self) -> (Self, Self) {
        let mut ra = GradeSet::g_none();
        let mut rb = GradeSet::g_none();
        for ka in a.iter_grades() {
            for kb in b.iter_grades() {
                if self.clone().prj(GradeSet::g(ka) * GradeSet::g(kb)).0.any() {
                    // The product of ka and kb yields at least one grade that
                    // is in self
                    ra = ra + GradeSet::g(ka);
                    rb = rb + GradeSet::g(kb);
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
    const G: fn(usize) -> GradeSet = GradeSet::g;
    const G_NONE: fn() -> GradeSet = GradeSet::g_none;

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
        assert_ne!(G(3), G(4))
    }

    test_eqs!(
        add_self_id: G(3) + G(3) => G(3),
        add_g_any_id: G(3) + G_NONE() => G(3),
        mul_g_any_absorb: G(3) * G_NONE() => G_NONE(),
        mul_vecs: G(1) * G(1) => G(0) + G(2),
        mul_scal_id: G(40) * G(0) => G(40),
        mul_bivec_quadvec: G(2) * G(4) => G(2) + G(4) + G(6),
        mul_trivec_quadvec: G(3) * G(4) => G(1) + G(3) + G(5) + G(7),
        mul_trivec_pentavec: G(3) * G(5) => G(2) + G(4) + G(6) + G(8),
        mul_vec_rotor: G(1) * (G(0) + G(2)) => G(1) + G(3),
        range: GradeSet::range(4,6) => G(4) + G(5) + G(6),
        project: GradeSet::range(0,10).prj(GradeSet::range(4,6)) => GradeSet::range(4,6),
        single_graded: (G(1) + G(1)).is_single_graded() => true,
        not_single_graded: (G(1) + G(2)).is_single_graded() => false,
        g_any_not_single_graded: G_NONE().is_single_graded() => false,
        iter_grades: (G(1) + G(22) + G(10)).iter_grades().collect::<Vec<_>>() => vec![1,10,22],
        grades_affecting_mul:
          G(0).grades_affecting_mul(&(G(1) + G(0) + G(2) + G(10)), &(G(0) + G(2) + G(6)))
          => (G(0) + G(2), G(0) + G(2))
    );
}
