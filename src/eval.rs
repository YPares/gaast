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

// impl<T> Ast<T> {
//     pub fn propagate_grades(&mut self, wanted: GradeSet) {
//         match self {
//             Ast::Val(_) => (),
//             Ast::Add(e1, e2) => {
//                 e1.grade_set.borrow_mut()
//             },
//             Ast::Mul(_, _) => todo!(),
//             Ast::Neg(_) => todo!(),
//             Ast::Exp(_) => todo!(),
//             Ast::Log(_) => todo!(),
//             Ast::Rev(_) => todo!(),
//             Ast::GInvol(_) => todo!(),
//             Ast::ScalarInv(_) => todo!(),
//         }
//     }
// }

// impl<T: Graded> GAExpr<T> {
//     /// Evaluates a [`GAExpr`]
//     pub fn eval(&self, m: impl MetricGA) -> DynSizedMV {
//         match self.ast.as_ref() {
//             Ast::Val(x) => todo!(),
//             Ast::Add(_, _) => todo!(),
//             Ast::Mul(_, _) => todo!(),
//             Ast::Neg(_) => todo!(),
//             Ast::Exp(_) => todo!(),
//             Ast::Log(_) => todo!(),
//             Ast::Rev(_) => todo!(),
//             Ast::GInvol(_) => todo!(),
//             Ast::ScalarInv(_) => todo!(),
//         }
//     }
// }
