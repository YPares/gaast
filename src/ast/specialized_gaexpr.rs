use super::{base_types::*, gaexpr::GradedNode};
use crate::{grade_set::*, Graded};
use std::{cell::Ref, fmt::Debug, rc::Rc};

/// A [`GaExpr`][crate::GaExpr] which is ready for compilation/evaluation
#[derive(Debug)]
#[repr(transparent)]
pub struct SpecializedGaExpr<T> {
    rc: Rc<GradedNode<Self, T>>,
}

/// Expressions are identifiable by reference, this allows to do caching when
/// evaluating an AST
#[derive(Hash, PartialEq, Eq)]
pub struct ExprId(
    /// A pointer that is just here for its hash/eq instances, it will never be dereferenced
    *const (),
);

impl<T> SpecializedGaExpr<T> {
    /// Whether that expression is used several times and could benefit from
    /// caching when evaluated
    pub fn is_reused(&self) -> bool {
        Rc::strong_count(&self.rc) >= 2
    }

    /// Get a hashable identifier for this expression. Expressions created via
    /// .clone() share the same identifier
    pub fn identify(&self) -> ExprId {
        ExprId(Rc::as_ptr(&self.rc).cast())
    }

    /// Get the node and the part of the AST it contains
    pub fn ast_node(&self) -> &AstNode<SpecializedGaExpr<T>, T> {
        &self.rc.ast_node
    }

    /// The underlying vector space of the algebra used by this expression
    pub fn vec_space_dim(&self) -> usize {
        *self
            .rc
            .vec_space_dim
            .get()
            .expect("vec_space_dim cell not set for this node")
    }
}

/// Get the minimal [`GradeSet`] inferred for this expression, and constrained
/// to what is available in the algebra given to [`GaExpr::specialize`][crate::GaExpr::specialize]
impl<T> Graded for SpecializedGaExpr<T> {
    type RefToGradeSet<'a> = Ref<'a, GradeSet> where T: 'a;
    fn grade_set(&self) -> Self::RefToGradeSet<'_> {
        // This corresponds to GaExpr::minimal_grade_set
        self.rc.minimal_grade_set.borrow()
    }
}
