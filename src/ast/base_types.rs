use crate::{algebra::Component, grade_set::*};
use std::{fmt::Debug, ops::Deref, rc::Rc};

/// The abstract syntax tree nodes representing geometric algebra primitive
/// operations. `T` is some raw multivector type, and `E` is a boxed type itself
/// containing an `AstNode`.
#[derive(Debug, Clone)]
pub enum AstNode<E, T> {
    /// Use of a raw multivector which exposes which grades it contains. To do
    /// so, most operations on `Ast` require `T: [Graded]`
    GradedObj(T),
    /// Multivector addition
    Addition(E, E),
    /// Some product
    Product(Product<E>),
    /// Multivector negation
    Negation(E),
    /// Multivector exponentiation
    Exponential(E),
    /// Multivector natural logarithm
    Logarithm(E),
    /// Grade projection (or "grade extraction")
    GradeProjection(E),
    /// Reverse (or "dagger")
    Reverse(E),
    /// Grade involution (or main involution)
    GradeInvolution(E),
    /// Operate only on the scalar part
    ScalarUnaryOp(ScalarUnaryOp, E),
}

/// Represents some product to perform. The actual individual
/// component-to-component multiplications to perform are initially empty, and
/// then updated during the AST specialization phase
#[derive(Debug, Clone)]
pub struct Product<E> {
    pub comp_muls_cell: Option<Vec<IndividualCompMul>>,
    pub(super) grades_to_produce: KVecsProductGradeSelection,
    pub left_expr: E,
    pub right_expr: E,
}

/// A component-to-component multiplication to perform on input data, and where
/// to store the result in output data
#[derive(Debug, Clone)]
pub struct IndividualCompMul {
    /// The component to read in the left operand
    pub left_comp: Component,
    /// The component to read in the right operand
    pub right_comp: Component,
    /// The component to update in the result
    pub result_comp: Component,
    /// The coefficient to apply to the product of the two input components
    pub coeff: f64,
}

/// When two k-vectors are multiplied together, select --given their grades--
/// which grades should be retained out of their geometric product. Is called
/// only during AST specialization phase
pub(super) struct KVecsProductGradeSelection(
    /// We use Rc and not Box so the closure can be directly reused by
    /// IsolatedGaExpr, as this closure cannot be cloned and is anyway immutable
    pub(super) Rc<dyn Fn((i64, i64)) -> GradeSet>,
);

impl Clone for KVecsProductGradeSelection {
    fn clone(&self) -> Self {
        KVecsProductGradeSelection(Rc::clone(&self.0))
    }
}

impl Debug for KVecsProductGradeSelection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("<KVecsProductGradeSelection>")
    }
}

impl KVecsProductGradeSelection {
    pub(super) fn get_fn(&self) -> impl Fn((i64, i64)) -> GradeSet + '_ {
        self.0.deref()
    }
}

#[derive(Hash, Debug, Clone)]
pub enum ScalarUnaryOp {
    Inversion,
    SquareRoot,
}
