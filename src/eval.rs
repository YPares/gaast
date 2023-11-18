//! How to evaluate a GA expression into an actual multivector

use super::{algebra::MetricAlgebra, ast::*, graded::*};
use AstNode::*;

impl<T: Graded> GAExpr<T> {
    /// Evaluates a [`GAExpr`]. The given [`MetricAlgebra`] must make sense with
    /// respect to the input values contained in the [`GAExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R: GradedMut>(&self, ga: &impl MetricAlgebra) -> R {
        let mut res = R::init_null_mv(ga.vec_space_dim(), &self.grade_set());
        self.add_to_res(ga, &mut res);
        res
    }

    fn add_to_res<R: GradedMut>(&self, ga: &impl MetricAlgebra, res: &mut R) {
        if self.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match self.ast_node() {
            Val(input) => {
                for k in res.grade_set().iter() {
                    if input.grade_set().contains(k) {
                        let input_slice = input.grade_slice(k);
                        let res_slice = res.grade_slice_mut(k);
                        for i in 0..res_slice.len() {
                            res_slice[i] += input_slice[i];
                        }
                    }
                }
            }
            Add(e1, e2) => {
                e1.add_to_res(ga, res);
                e2.add_to_res(ga, res);
            }
            Mul(e1, e2) => {
                let r1: R = e1.eval(ga);
                let r2: R = e2.eval(ga);
                for k in self.grade_set().iter() {}
            }
            Exp(e) => todo!(),
            Log(e) => todo!(),
            Neg(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            Rev(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            GInvol(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            ScalarInv(e) => {
                e.add_to_res(ga, res);
                let s = res.grade_slice_mut(0);
                s[0] = 1.0 / s[0];
            }
            Prj(e, _) => {
                e.add_to_res(ga, res);
            }
        }
    }
}
