//! How to evaluate a GA expression into an actual multivector

use super::{algebra::*, ast::*, graded::*};
use AstNode as N;

impl<T: GradedInput> GAExpr<T> {
    /// Evaluates a [`GAExpr`]. The given [`MetricAlgebra`] must make sense with
    /// respect to the input values contained in the [`GAExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R>(&self, alg: &ReadyAlgebra<impl MetricAlgebra>) -> R
    where
        R: GradedInput + GradedOutput,
    {
        let mut res = R::init_null_mv(alg.vec_space_dim(), &self.grade_set());
        self.add_to_res(alg, &mut res);
        res
    }

    fn add_to_res<R>(&self, alg: &ReadyAlgebra<impl MetricAlgebra>, res: &mut R)
    where
        R: GradedInput + GradedOutput,
    {
        if self.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match self.ast_node() {
            N::Val(input) => {
                res.add_grades_from(input, &self.grade_set());
            }
            N::Add(e_left, e_right) => {
                e_left.add_to_res(alg, res);
                e_right.add_to_res(alg, res);
            }
            N::Neg(e) => {
                e.add_to_res(alg, res);
                for k in self.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            N::Mul(e_left, e_right) => {
                self.perform_mul(e_left, e_right, alg, res);
            }
            N::Rev(e) => {
                e.add_to_res(alg, res);
                for k in self.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::GInvol(e) => {
                e.add_to_res(alg, res);
                for k in self.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::ScalarInv(e) => {
                e.add_to_res(alg, res);
                let s = res.grade_slice_mut(0);
                s[0] = 1.0 / s[0];
            }
            N::Prj(e, _) => {
                if *self.grade_set() == *e.grade_set() {
                    // Projection is a no-op: `res` is already what the
                    // underlying expr `e` expects
                    e.add_to_res(alg, res);
                } else {
                    // Projection actually filters out stuff: `res` misses some
                    // grades to be readily used by the underlying expr
                    // evaluator. We need to allocate and copy
                    res.add_grades_from(&e.eval::<R>(alg), &self.grade_set());
                }
            }
            N::Exp(_e) => todo!(),
            N::Log(_e) => todo!(),
        }
    }

    fn perform_mul<R>(
        &self,
        e_left: &GAExpr<T>,
        e_right: &GAExpr<T>,
        alg: &ReadyAlgebra<impl MetricAlgebra>,
        mv_res: &mut R,
    ) where
        R: GradedInput + GradedOutput,
    {
        let mv_left: R = e_left.eval(alg);
        let mv_right: R = e_right.eval(alg);
        for (_, k_left, k_right) in self
            .grade_set()
            .iter_contributions_to_mul(&e_left.grade_set(), &e_right.grade_set())
        {
            for (i_left, c_left) in mv_left.grade_slice(k_left).iter().enumerate() {
                for (i_right, c_right) in mv_right.grade_slice(k_right).iter().enumerate() {
                    let bb_left = alg.basis_blade_from_coord(k_left, i_left);
                    let bb_right = alg.basis_blade_from_coord(k_right, i_right);
                    let (bb_res, coef) = alg.ortho_basis_blades_mul(bb_left, bb_right);
                    let (k_res, i_res) = alg.coord_from_basis_blade(&bb_res);
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

    /// Tests that a GAExpr returns the correct result with AND without the use
    /// of the grade minimisation phase
    macro_rules! expr_eq {
        ($geom:ident, $a:expr, $b:expr) => {{
            let a = $a;
            let b = $b;
            assert_eq!(a.eval::<GradeMapMV>(&$geom), b);
            assert_eq!(a.minimize_grades().eval::<GradeMapMV>(&$geom), b);
        }};
    }

    //type Ex = GAExpr<GradeMapMV>;
    //const V: fn(GradeMapMV) -> Ex = GAExpr::val;

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
        expr_eq!(ega3, (e2 * e1).g(2), grade_map_mv!(2 => -1 0 0));
    }

    #[rstest]
    fn vecs_to_trivec(ega3: Ega3) {
        let [e1, e2, e3] = ega3.base_vec_exprs::<GradeMapMV>();
        expr_eq!(ega3, (e2 * e1 * e3).g(3), grade_map_mv!(3 => -1));
    }

    #[rstest]
    fn pga(pga2: Pga2) {
        let [e0, _, _] = pga2.base_vec_exprs::<GradeMapMV>();
        expr_eq!(pga2, (e0.clone() * e0).g(0), grade_map_mv!(0 => 0))
    }
}
