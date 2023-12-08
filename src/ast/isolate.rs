use super::{base_types::*, gaexpr::*};
use crate::{GradeSet, Graded};
use std::{collections::HashMap, rc::Rc};
use AstNode as N;

/// Just an identifier. The container pointer is never dereferenced
#[derive(Hash, PartialEq, Eq, Clone, Copy)]
pub struct ExprId {
    pub(crate) ptr: *const (),
}

impl std::fmt::Debug for ExprId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.ptr.fmt(f)
    }
}

/// A node that is guaranteed not to be mutably shared between different
/// expressions, and can thus be safely mutated
#[derive(Debug)]
pub struct IsolatedNode<T> {
    /// The GradeSet inferred when constructing the AST
    pub(super) maximal_grade_set: GradeSet,
    /// The GradeSet inferred by specialization. Starts empty and receives
    /// updates as we go through the AST
    pub(super) minimal_grade_set: GradeSet,
    /// The dimension of the vec space of the algebra used by this node
    pub(super) vec_space_dim: Option<usize>,
    /// The operation performed by this node
    pub(super) ast_node: AstNode<ExprId, T>,
    /// How many times is the node referred to in the expression it belongs to
    pub(super) num_uses: u32,
}

/// Get the [`GradeSet`] inferred for this node by the specialization process
impl<T> Graded for IsolatedNode<T> {
    type RefToGradeSet<'a> = &'a GradeSet where T: 'a;
    /// The [`GradeSet`] inferred for this node by the specialization process
    fn grade_set(&self) -> Self::RefToGradeSet<'_> {
        &self.minimal_grade_set
    }
}
impl<T> IsolatedNode<T> {
    /// The dimension of the vector space used by the multivector this node
    /// evaluates to
    pub fn vec_space_dim(&self) -> usize {
        self.vec_space_dim.expect("vec_space_dim not set for node")
    }
    /// The operation that this node performs
    pub fn ast_node(&self) -> &AstNode<ExprId, T> {
        &self.ast_node
    }
    /// Tells if that node is used several times throughout the whole
    /// expression, and could therefore benefit from caching when evaluating
    pub fn is_used_several_times(&self) -> bool {
        self.num_uses >= 2
    }
}

pub(super) type NodeArena<T> = HashMap<ExprId, IsolatedNode<T>>;

/// A [`GaExpr`] where potentially mutable parts have been copied from their
/// original representation
#[derive(Debug)]
pub(super) struct IsolatedGaExpr<T> {
    pub(super) arena: NodeArena<T>,
    pub(super) root: ExprId,
}

impl<T> GaExpr<T> {
    pub(super) fn identify(&self) -> ExprId {
        ExprId {
            ptr: Rc::as_ptr(&self.rc).cast(),
        }
    }

    /// Turn the GaExpr into an arena-based storage, which can then be mutated
    pub(super) fn isolate(&self) -> IsolatedGaExpr<&T> {
        let mut arena = HashMap::new();
        self.rec_store_in_arena(&mut arena);
        IsolatedGaExpr {
            arena,
            root: self.identify(),
        }
    }

    fn rec_store_in_arena<'a>(&'a self, arena: &mut NodeArena<&'a T>) {
        if let Some(node) = arena.get_mut(&self.identify()) {
            // This node is already in the arena
            node.num_uses += 1;
            return;
        }
        let new_node = match &self.rc.ast_node {
            N::GradedObj(x) => N::GradedObj(x),
            N::Addition(left, right) => {
                left.rec_store_in_arena(arena);
                right.rec_store_in_arena(arena);
                N::Addition(left.identify(), right.identify())
            }
            N::Product(p) => {
                p.left_expr.rec_store_in_arena(arena);
                p.right_expr.rec_store_in_arena(arena);
                N::Product(Product {
                    comp_muls_cell: None,
                    grades_to_produce: p.grades_to_produce.clone(),
                    left_expr: p.left_expr.identify(),
                    right_expr: p.right_expr.identify(),
                })
            }
            N::Negation(e) => {
                e.rec_store_in_arena(arena);
                N::Negation(e.identify())
            }
            N::Exponential(e) => {
                e.rec_store_in_arena(arena);
                N::Exponential(e.identify())
            }
            N::Logarithm(e) => {
                e.rec_store_in_arena(arena);
                N::Logarithm(e.identify())
            }
            N::GradeProjection(e) => {
                e.rec_store_in_arena(arena);
                N::GradeProjection(e.identify())
            }
            N::Reverse(e) => {
                e.rec_store_in_arena(arena);
                N::Reverse(e.identify())
            }
            N::GradeInvolution(e) => {
                e.rec_store_in_arena(arena);
                N::GradeInvolution(e.identify())
            }
            N::ScalarUnaryOp(op, e) => {
                e.rec_store_in_arena(arena);
                N::ScalarUnaryOp(op.clone(), e.identify())
            }
        };
        arena.insert(
            self.identify(),
            IsolatedNode {
                maximal_grade_set: self.grade_set().clone(),
                minimal_grade_set: GradeSet::empty(),
                vec_space_dim: None,
                ast_node: new_node,
                num_uses: 1,
            },
        );
    }
}
