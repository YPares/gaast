//! How to evaluate a GA expression into an actual multivector

use std::collections::HashMap;

use super::{algebra::*, ast::*, graded::*};
use AstNode as N;

type Cache<R> = HashMap<ExprId, R>;

impl<T: GradedInput> ReadyGaExpr<T> {
    /// Evaluates a [`GaExpr`]. The given [`MetricAlgebra`] must make sense with
    /// respect to the input values contained in the [`GaExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R>(&self, alg: &ReadyAlgebra<impl MetricAlgebra>) -> R
    where
        R: GradedInput + GradedOutput + Clone,
    {
        self.eval_with_cache(alg, &mut HashMap::new())
    }

    fn eval_with_cache<R>(
        &self,
        alg: &ReadyAlgebra<impl MetricAlgebra>,
        cache: &mut Cache<R>,
    ) -> R
    where
        R: GradedInput + GradedOutput + Clone,
    {
        if self.is_reused() {
            match cache.get(&self.identify()) {
                Some(r) => r.clone(),
                None => {
                    let mut res = R::init_null_mv(alg.vec_space_dim(), &self.grade_set());
                    self.add_to_res(alg, cache, &mut res);
                    cache.insert(self.identify(), res.clone());
                    res
                }
            }
        } else {
            let mut res = R::init_null_mv(alg.vec_space_dim(), &self.grade_set());
            self.add_to_res(alg, cache, &mut res);
            res
        }
    }

    fn add_to_res<R>(
        &self,
        alg: &ReadyAlgebra<impl MetricAlgebra>,
        cache: &mut Cache<R>,
        res: &mut R,
    ) where
        R: GradedInput + GradedOutput + Clone,
    {
        if self.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match self.ast_node() {
            N::GradedObj(input) => {
                res.add_grades_from(input, &self.grade_set());
            }
            N::Addition(e_left, e_right) => {
                e_left.add_to_res(alg, cache, res);
                e_right.add_to_res(alg, cache, res);
            }
            N::Negation(e) => {
                e.add_to_res(alg, cache, res);
                for k in self.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            N::GeometricProduct(e_left, e_right) => {
                self.eval_gp(e_left, e_right, alg, cache, res);
            }
            N::Reverse(e) => {
                e.add_to_res(alg, cache, res);
                for k in self.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::GradeInvolution(e) => {
                e.add_to_res(alg, cache, res);
                for k in self.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::ScalarInversion(e) => {
                e.add_to_res(alg, cache, res);
                let s = res.grade_slice_mut(0);
                s[0] = 1.0 / s[0];
            }
            N::GradeProjection(e, _) => {
                if self.grade_set() == e.grade_set() {
                    // Projection is a no-op: `res` is already what the
                    // underlying expr `e` expects
                    e.add_to_res(alg, cache, res);
                } else {
                    // Projection actually filters out stuff: `res` misses some
                    // grades to be readily used by the underlying expr
                    // evaluator. We need to allocate and copy
                    res.add_grades_from(&e.eval::<R>(alg), &self.grade_set());
                }
            }
            N::Exponential(_e) => todo!(),
            N::Logarithm(_e) => todo!(),
        }
    }

    fn eval_gp<R>(
        &self,
        e_left: &ReadyGaExpr<T>,
        e_right: &ReadyGaExpr<T>,
        alg: &ReadyAlgebra<impl MetricAlgebra>,
        cache: &mut Cache<R>,
        mv_res: &mut R,
    ) where
        R: GradedInput + GradedOutput + Clone,
    {
        let mv_left: R = e_left.eval_with_cache(alg, cache);
        let mv_right: R = e_right.eval_with_cache(alg, cache);
        for (_, k_left, k_right) in self
            .grade_set()
            .iter_contributions_to_gp(&e_left.grade_set(), &e_right.grade_set())
        {
            for (i_left, c_left) in mv_left.grade_slice(k_left).iter().enumerate() {
                for (i_right, c_right) in mv_right.grade_slice(k_right).iter().enumerate() {
                    let bb_left = alg.basis_blade_from_indexes(k_left, i_left);
                    let bb_right = alg.basis_blade_from_indexes(k_right, i_right);
                    let (bb_res, coef) = alg.ortho_basis_blades_gp(bb_left, bb_right);
                    let (k_res, i_res) = alg.indexes_from_basis_blade(&bb_res);
                    if self.grade_set().contains(k_res) {
                        mv_res.grade_slice_mut(k_res)[i_res] += c_left * c_right * coef;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{algebra::*, grade_map_mv, graded::GradeMapMV};
    use rstest::*;

    macro_rules! expr_eq {
        ($alg:ident, $a:expr, $b:expr) => {{
            let a = $a;
            let b = $b;
            assert_eq!(a.minimize_grades().eval::<GradeMapMV>(&$alg), b);
        }};
    }

    type Ega3 = ReadyAlgebra<[f64; 3]>;
    #[fixture]
    fn ega3() -> Ega3 {
        ReadyAlgebra::from([1.0, 1.0, 1.0])
    }

    type Pga2 = ReadyAlgebra<[f64; 3]>;
    #[fixture]
    fn pga2() -> Pga2 {
        ReadyAlgebra::from([0.0, 1.0, 1.0])
    }

    #[rstest]
    fn vecs_to_bivec(ega3: Ega3) {
        let [e1, e2, _] = ega3.base_vec_exprs::<GradeMapMV>();
        let ex = e2 ^ e1;
        expr_eq!(ega3, ex, grade_map_mv!(2 => -1 0 0));
    }

    #[rstest]
    fn vecs_to_trivec(ega3: Ega3) {
        let [e1, e2, e3] = ega3.base_vec_exprs::<GradeMapMV>();
        expr_eq!(ega3, e2 ^ e1 ^ e3, grade_map_mv!(3 => -1));
    }

    #[rstest]
    fn vec_norm(pga2: Pga2) {
        let [e0, e1, e2] = pga2.base_vec_exprs::<GradeMapMV>();
        expr_eq!(pga2, (e0 - 2 * e1 + e2).norm_sq(), grade_map_mv!(0 => 5));
    }

    #[rstest]
    fn projection(ega3: Ega3) {
        let [e1, e2, e3] = ega3.base_vec_exprs::<GradeMapMV>();
        let v = e1.clone() + e2.clone();
        let bv = 4 * e1 ^ e3;
        // project v onto bv:
        expr_eq!(ega3, (v & bv.clone()) & bv.inv(), grade_map_mv!(1 => 1 0 0));
    }
}
