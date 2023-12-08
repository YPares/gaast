use crate::{algebra::Component, grade_set::*, Graded};
use std::{collections::HashMap, fmt::Debug, ops::Deref, rc::Rc};

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
    pub individual_comp_muls: Vec<IndividualCompMul>,
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
    /// IsolatedAst, as this closure cannot be cloned and is anyway immutable
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

/// Identifies some node of the AST, to allow for various nodes to point to the
/// same subexpression
#[derive(Hash, PartialEq, Eq, Clone, Copy)]
pub struct NodeId {
    /// Just used as an identifier. Is never dereferenced
    pub(crate) ptr: *const (),
}

impl std::fmt::Debug for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.ptr.fmt(f)
    }
}

/// An [`AstNode`] with associated metadata
#[derive(Debug)]
pub struct GradedNode<T> {
    /// The GradeSet inferred when constructing the AST
    pub(super) maximal_grade_set: GradeSet,
    /// The GradeSet inferred by specialization. Starts empty and receives
    /// updates as we go through the AST
    pub(super) minimal_grade_set: GradeSet,
    /// The dimension of the vec space of the algebra used by this node
    pub(super) vec_space_dim: usize,
    /// The operation performed by this node
    pub(super) ast_node: AstNode<NodeId, T>,
    /// How many times is the node referred to in the expression it belongs to
    pub(super) num_uses: u32,
    /// Whether this node has already been fully processed (minimal_grade_set
    /// fully infered and ready to use, and IndividualCompMuls resolved)
    pub(super) is_ready: bool,
}

/// Get the [`GradeSet`] inferred for this node by the specialization process
impl<T> Graded for GradedNode<T> {
    type RefToGradeSet<'a> = &'a GradeSet where T: 'a;
    /// The [`GradeSet`] inferred for this node by the specialization process
    fn grade_set(&self) -> Self::RefToGradeSet<'_> {
        &self.minimal_grade_set
    }
}
impl<T> GradedNode<T> {
    /// The dimension of the vector space used by the multivector this node
    /// evaluates to
    pub fn vec_space_dim(&self) -> usize {
        self.vec_space_dim
    }
    /// The operation that this node performs
    pub fn ast_node(&self) -> &AstNode<NodeId, T> {
        &self.ast_node
    }
    /// Tells if that node is used several times throughout the whole
    /// expression, and could therefore benefit from caching when evaluating
    pub fn is_used_several_times(&self) -> bool {
        self.num_uses >= 2
    }
}

pub(super) type NodeArena<T> = HashMap<NodeId, GradedNode<T>>;

/// An [`Expr`] that has been reified into an actual AST. This AST is guaranteed
/// not to be shared between different [`Expr`]s, it is therefore safe to mutate
/// it
#[derive(Debug)]
pub(super) struct ReifiedAst<T> {
    pub(super) arena: NodeArena<T>,
    pub(super) root_id: NodeId,
}
