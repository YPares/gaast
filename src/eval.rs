use super::ast::*;
use super::raw::RawMV;

/// Informs about:
/// 
/// - the dimensionality of the vector space
/// - how each base vector squares
pub struct Metric<const D: usize> {
    /// What each base vector squares to. Base vector are identified by position
    /// in the array. The length of the array tells the dimension of the vector
    /// space
    pub base_squares: [f64; D]
}

/// Evaluates a [`GAExpr`]
pub fn eval<T: Graded, const D: usize>(e: GAExpr<T>, m: Metric<D>) -> RawMV {
    match e.ast.as_ref() {
        Ast::Val(v) => todo!(),
        Ast::Add(_, _) => todo!(),
        Ast::Mul(_, _) => todo!(),
        Ast::Exp(_) => todo!(),
        Ast::Log(_) => todo!(),
        Ast::Rev(_) => todo!(),
        Ast::GInvol(_) => todo!(),
        Ast::ScalarInv(_) => todo!(),
    }
}
