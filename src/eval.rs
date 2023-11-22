//! How to evaluate a GA expression into an actual multivector

use super::{algebra::*, ast::*, graded::*};
use std::borrow::Borrow as _;
use AstNode as N;

impl<T: GradedInput> GAExpr<T> {
    /// Evaluates a [`GAExpr`]. The given [`MetricAlgebra`] must make sense with
    /// respect to the input values contained in the [`GAExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R>(&self, ga: &ReadyAlgebra<impl MetricAlgebra>) -> R
    where
        R: GradedInput + GradedOutput,
    {
        let mut res = R::init_null_mv(ga.vec_space_dim(), &self.grade_set());
        self.add_to_res(ga, &mut res);
        res
    }

    fn add_to_res<R>(&self, ga: &ReadyAlgebra<impl MetricAlgebra>, res: &mut R)
    where
        R: GradedInput + GradedOutput,
    {
        if self.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match self.ast_node() {
            N::Val(input) => {
                let igs = input.grade_set();
                let rgs = self.grade_set().clone();
                for k in rgs.iter() {
                    if igs.borrow().contains(k) {
                        let input_slice = input.grade_slice(k);
                        let res_slice = res.grade_slice_mut(k);
                        for (r, i) in res_slice.iter_mut().zip(input_slice) {
                            *r = *r + i;
                        }
                    }
                }
            }
            N::Add(e_left, e_right) => {
                e_left.add_to_res(ga, res);
                e_right.add_to_res(ga, res);
            }
            N::Neg(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            N::Rev(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::GInvol(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::ScalarInv(e) => {
                e.add_to_res(ga, res);
                let s = res.grade_slice_mut(0);
                s[0] = 1.0 / s[0];
            }
            N::Prj(e, _) => {
                e.add_to_res(ga, res);
            }
            N::Mul(e_left, e_right) => {
                self.perform_mul(e_left, e_right, ga, res);
            }
            N::Exp(_e) => todo!(),
            N::Log(_e) => todo!(),
        }
    }

    fn perform_mul<R>(
        &self,
        e_left: &GAExpr<T>,
        e_right: &GAExpr<T>,
        ga: &ReadyAlgebra<impl MetricAlgebra>,
        res: &mut R,
    ) where
        R: GradedInput + GradedOutput,
    {
        let r_left: R = e_left.eval(ga);
        let r_right: R = e_right.eval(ga);
        for (_, k_left, k_right) in self
            .grade_set()
            .iter_contributions_to_mul(&e_left.grade_set(), &e_right.grade_set())
        {
            for (i_left, c_left) in r_left.grade_slice(k_left).iter().enumerate() {
                for (i_right, c_right) in r_right.grade_slice(k_right).iter().enumerate() {
                    let bb_left = ga.basis_blade_from_coord(k_left, i_left);
                    let bb_right = ga.basis_blade_from_coord(k_right, i_right);
                    let (bb_res, coef) = ga.ortho_basis_blades_mul(bb_left, bb_right);
                    let (k_res, i_res) = ga.coord_from_basis_blade(&bb_res);
                    if self.grade_set().contains(k_res) {
                        res.grade_slice_mut(k_res)[i_res] += c_left * c_right * coef;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{algebra::{ReadyAlgebra, OrthoEuclidN}, graded::HashMapMV, hash_map_mv, GAExpr};

    const V: fn(HashMapMV) -> GAExpr<HashMapMV> = GAExpr::val;

    #[test]
    fn vecs_to_bivec() {
        let a = hash_map_mv!(1 => 0 1 0);
        let b = hash_map_mv!(1 => 1 0 0);
        let alg = ReadyAlgebra::prepare(OrthoEuclidN(3));
        let e = (V(a) * V(b)).g(2);
        let r: HashMapMV = e.resolve_minimum_grades().eval(&alg);
        assert_eq!(r, hash_map_mv!(2 => -1 0 0));
    }
}
