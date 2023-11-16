//! How to evaluate a GA expression into an actual multivector

use crate::grade_set::GradeSet;

use super::ast::*;
use super::graded::*;

/// A geometric algebra over some vector space
pub trait GA {
    /// The dimensionality of the underlying vector space
    fn space_dim(&self) -> usize;

    /// The number of components of a multivector of this geometric algebra
    fn algebraic_dim(&self) -> usize {
        (2 as usize).pow(self.space_dim() as u32)
    }

    /// Number of grades in this algebra
    fn num_grades(&self) -> usize {
        self.space_dim() + 1
    }

    /// Give the number of components of k-vectors
    fn grade_dim(&self, k: usize) -> usize {
        n_choose_k(self.space_dim(), k)
    }
}

/// A metric geometric algebra over some vector space
pub trait MetricGA: GA {
    /// Give the scalar product of two base vectors (identified by index)
    fn base_vec_scal_prod(&self, v1: usize, v2: usize) -> f64;
}

/// Representation of a [`GA`] as an array of orthogonal base vector squares
impl<const D: usize> GA for [f64; D] {
    fn space_dim(&self) -> usize {
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
    /// Evaluates a [`GAExpr`]
    pub fn eval(&self, m: impl MetricGA) -> DynSizedMV {
        match self.ast() {
            AstNode::Val(x) => todo!(),
            AstNode::Add(_, _) => todo!(),
            AstNode::Mul(_, _) => todo!(),
            AstNode::Neg(_) => todo!(),
            AstNode::Exp(_) => todo!(),
            AstNode::Log(_) => todo!(),
            AstNode::Rev(_) => todo!(),
            AstNode::GInvol(_) => todo!(),
            AstNode::ScalarInv(_) => todo!(),
            AstNode::Prj(_, _) => todo!(),
        }
    }
}