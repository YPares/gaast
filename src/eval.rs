//! How to evaluate a GA expression into an actual multivector

use std::ops::Add;

use crate::grade_set::GradeSet;

use super::ast::*;
use super::graded::*;

/// A geometric algebra over some vector space
pub trait GA {
    /// The dimensionality of the underlying vector space
    fn vec_space_dim(&self) -> usize;

    /// The number of components of a multivector of this geometric algebra
    fn algebraic_dim(&self) -> usize {
        (2 as usize).pow(self.vec_space_dim() as u32)
    }

    /// Number of grades in this algebra
    fn num_grades(&self) -> usize {
        self.vec_space_dim() + 1
    }

    /// Give the number of components of k-vectors
    fn grade_dim(&self, k: usize) -> usize {
        n_choose_k(self.vec_space_dim(), k)
    }
}

/// A metric geometric algebra over some vector space
pub trait MetricGA: GA {
    /// Give the scalar product of two base vectors (identified by index)
    fn base_vec_scal_prod(&self, v1: usize, v2: usize) -> f64;
}

/// Representation of a [`GA`] as an array of orthogonal base vector squares
impl<const D: usize> GA for [f64; D] {
    fn vec_space_dim(&self) -> usize {
        D
    }
}

/// Representation of a [`MetricGA`] as an array of orthogonal base vector
/// squares
impl<const D: usize> MetricGA for [f64; D] {
    fn base_vec_scal_prod(&self, v1: usize, v2: usize) -> f64 {
        if v1 == v2 {
            self[v1]
        } else {
            0.0
        }
    }
}

impl<T: Graded> GAExpr<T> {
    /// Evaluates a [`GAExpr`]. The given [`MetricGA`] must make sense with
    /// respect to the input values contained in the [`GAExpr`], in terms of
    /// possible grades contained in those input values, and of number of
    /// components for each grade
    pub fn eval<R: GradedMut>(&self, ga: &impl MetricGA) -> R {
        let mut res = R::init_null_mv(ga.vec_space_dim(), &self.grade_set());
        self.add_to_res(ga, &mut res);
        res
    }

    fn add_to_res<R: GradedMut>(&self, ga: &impl MetricGA, res: &mut R) {
        if self.grade_set().has_no_grade() {
            return;
        }
        match self.ast_node() {
            AstNode::Val(input) => {
                for k in res.grade_set().iter_grades() {
                    let input_slice = input.grade_slice(k);
                    let res_slice = res.grade_slice_mut(k);
                    for i in 0..res_slice.len() {
                        res_slice[i] += input_slice[i];
                    }
                }
            }
            AstNode::Add(e1, e2) => {
                e1.add_to_res(ga, res);
                e2.add_to_res(ga, res);
            }
            AstNode::Mul(e1, e2) => {
                let r1: R = e1.eval(ga);
                let r2: R = e2.eval(ga);
                todo!()
            }
            AstNode::Exp(e) => todo!(),
            AstNode::Log(e) => todo!(),
            AstNode::Neg(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter_grades() {
                    res.negate_grade(k);
                }
            }
            AstNode::Rev(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter_grades() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            AstNode::GInvol(e) => {
                e.add_to_res(ga, res);
                for k in self.grade_set().iter_grades() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            AstNode::ScalarInv(e) => {
                e.add_to_res(ga, res);
                let s = res.grade_slice_mut(0);
                s[0] = 1.0 / s[0];
            }
            AstNode::Prj(e, _) => {
                e.add_to_res(ga, res);
            }
        }
    }
}
