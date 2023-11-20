//! How to evaluate a GA expression into an actual multivector

use super::{algebra::MetricAlgebra, ast::*, graded::*};
use std::borrow::Borrow as _;
use AstNode as N;

impl<T: GradedInput> GAExpr<T> {
    /// Evaluates a [`GAExpr`]. The given [`MetricAlgebra`] must make sense with
    /// respect to the input values contained in the [`GAExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R>(&self, ga: &impl MetricAlgebra) -> R
    where
        R: GradedInput + GradedOutput,
    {
        let mut res = R::init_null_mv(ga.vec_space_dim(), &self.grade_set());
        self.add_to_res(ga, &mut res);
        res
    }

    fn add_to_res<R>(&self, ga: &impl MetricAlgebra, res: &mut R)
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
                let r_left: R = e_left.eval(ga);
                let r_right: R = e_right.eval(ga);
                for (partial_gs, k_left, k_right) in self
                    .grade_set()
                    .iter_contributions_to_mul(&e_left.grade_set(), &e_right.grade_set())
                {
                    
                }
            }
            N::Exp(e) => todo!(),
            N::Log(e) => todo!(),
        }
    }
}
