//! Represent a GA expression using an abstract syntax tree

use crate::algebra::{iter_basis_blades_of_grade, MetricAlgebra};
use crate::graded::GradedDataMut;

use super::algebra::Component;
use super::{grade_set::*, graded::Graded};
use std::{
    cell::{OnceCell, Ref, RefCell},
    rc::Rc,
};
use AstNode as N;

/// The abstract syntax tree nodes representing geometric algebra primitive
/// operations. `T` is some raw multivector type, and `E` is a boxed type itself
/// containing an `AstNode`.
#[derive(Debug)]
pub enum AstNode<E, T> {
    /// Use of a raw multivector which exposes which grades it contains. To do
    /// so, most operations on `Ast` require `T: [Graded]`
    GradedObj(T),
    /// Multivector addition
    Addition(E, E),
    /// Geometric product. The actual individual component-to-component
    /// multiplications to perform are initially empty, and then updated during
    /// the AST preparation phase
    GeometricProduct(OnceCell<Vec<IndividualCompMul>>, E, E),
    /// Multivector negation
    Negation(E),
    /// Multivector exponentiation. See [`GaExpr::exp`] for limitations
    Exponential(E),
    /// Multivector natural logarithm. See [`GaExpr::log`] for limitations
    Logarithm(E),
    /// Grade projection (or "grade extraction"). The grade(s) to extract are
    /// stored here only for error-reporting reasons
    GradeProjection(E, GradeSet),
    /// Reverse (or "dagger")
    Reverse(E),
    /// Grade involution (or main involution)
    GradeInvolution(E),
    /// Operate only on the scalar part
    ScalarUnaryOp(ScalarUnaryOp, E),
}

/// A component-to-component multiplication to perform on input data, and where
/// to store the result in output data
#[derive(Debug)]
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

#[derive(Hash, Debug)]
pub enum ScalarUnaryOp {
    Inversion,
    SquareRoot,
}

/// Assign a [`GradeSet`] to some [`AstNode`]. This [`GradeSet`] can be modified
/// to be restricted further depending on where this [`AstNode`] is being used,
/// ie. which grades its use sites actually need
#[derive(Debug)]
struct GradedNode<E, T> {
    /// The grade set inferred at AST construction (inferred from the maximal
    /// grade sets of this node's subexpressions)
    maximal_grade_set: GradeSet,
    /// Starting at empty, is updated at AST preparation, to finally reflect on
    /// the minimal set of grades to allocate/compute at evaluation
    minimal_grade_set: RefCell<GradeSet>,
    /// Starting at empty, is set at AST preparation given the metric
    vec_space_dim: OnceCell<usize>,
    /// This node
    ast_node: AstNode<E, T>,
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [`GradeSet`] which can then
/// be further restrained depending on the use sites of this GaExpr, and then
/// used to optimize allocations while evaluating the expression
#[derive(Debug)]
#[repr(transparent)]
pub struct GaExpr<T> {
    rc: Rc<GradedNode<Self, T>>,
}

/// Create a GA expression that just returns some pre-evaluated [`Graded`]
/// value (multivector)
pub fn mv<T: Graded>(x: T) -> GaExpr<T> {
    let gs = x.grade_set().clone();
    GaExpr::new(gs, AstNode::GradedObj(x))
}

/// [`GaExpr`] uses [`Rc`] internally, so you can clone it extensively
impl<T> Clone for GaExpr<T> {
    fn clone(&self) -> Self {
        Self {
            rc: Rc::clone(&self.rc),
        }
    }
}

/// Multivector addition
impl<T, E: Into<GaExpr<T>>> std::ops::Add<E> for GaExpr<T> {
    type Output = Self;
    /// Multivector addition
    fn add(self, rhs: E) -> Self::Output {
        let e_rhs = rhs.into();
        let gs = self.grade_set().clone() + e_rhs.grade_set().clone();
        Self::new(gs, N::Addition(self, e_rhs))
    }
}

macro_rules! gaexpr_products {
    ($($trait:ident $method:ident $doc:literal),*) => {
        $(
        #[doc=$doc]
        impl<T, E: Into<GaExpr<T>>> std::ops::$trait<E> for GaExpr<T> {
            type Output = Self;
            #[doc=$doc]
            fn $method(self, rhs: E) -> Self::Output {
                let e_rhs = rhs.into();
                let gs = self.grade_set().clone().$method(e_rhs.grade_set().clone());
                Self::new(gs, N::GeometricProduct(OnceCell::new(), self, e_rhs))
            }
        }
        )*
    };
}

gaexpr_products! {
    Mul mul "Geometric product",
    BitXor bitxor "Outer product",
    BitAnd bitand "Inner product",
    Shl shl "Left contraction",
    Shr shr "Right contraction"
}

impl<T> std::ops::Neg for GaExpr<T> {
    type Output = Self;
    fn neg(self) -> Self {
        let gs = self.grade_set().clone();
        Self::new(gs, N::Negation(self))
    }
}

impl<T, E: Into<GaExpr<T>>> std::ops::Sub<E> for GaExpr<T> {
    type Output = Self;
    fn sub(self, rhs: E) -> Self::Output {
        self + -rhs.into()
    }
}

impl<T: GradedDataMut> From<f64> for GaExpr<T> {
    fn from(x: f64) -> GaExpr<T> {
        if x == 0.0 {
            return mv(T::init_null_mv(0, &GradeSet::empty()));
        }
        let mut scal_mv = T::init_null_mv(0, &GradeSet::single(0));
        scal_mv.grade_slice_mut(0)[0] = x;
        mv(scal_mv)
    }
}
impl<T: GradedDataMut> From<i64> for GaExpr<T> {
    #[inline]
    fn from(x: i64) -> GaExpr<T> {
        (x as f64).into()
    }
}

macro_rules! scalar_with_gaexpr_binary_ops {
    ($($t:ty),*) => {
        $(
        impl<T: GradedDataMut> std::ops::Add<GaExpr<T>> for $t {
            type Output = GaExpr<T>;
            #[inline]
            fn add(self, rhs: GaExpr<T>) -> Self::Output {
                Into::<GaExpr<T>>::into(self) + rhs
            }
        }
        impl<T: GradedDataMut> std::ops::Mul<GaExpr<T>> for $t {
            type Output = GaExpr<T>;
            #[inline]
            fn mul(self, rhs: GaExpr<T>) -> Self::Output {
                Into::<GaExpr<T>>::into(self) * rhs
            }
        }
        impl<T: GradedDataMut> std::ops::Div<$t> for GaExpr<T> {
            type Output = Self;
            fn div(self, rhs: $t) -> Self {
                self * (1.0/(rhs as f64))
            }
        }
        )*
    };
}
scalar_with_gaexpr_binary_ops!(f64, i64);

macro_rules! gaexpr_unary_ops {
    ($($fn_name:ident $ctor:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                let gs = self.grade_set().clone().$grade_op();
                Self::new(
                    gs,
                    N::$ctor(self),
                )
            }
        )*
    }
}

/// Get the [`GradeSet`] inferred for this expression. Note that it may contain
/// grades that won't exist in the final algebra, because that algebra isn't
/// known yet
impl<T> Graded for GaExpr<T> {
    type RefToGradeSet<'a> = &'a GradeSet where T: 'a;
    fn grade_set(&self) -> Self::RefToGradeSet<'_> {
        &self.rc.maximal_grade_set
    }
}

impl<T> GaExpr<T> {
    fn new(gs: GradeSet, ast_node: AstNode<Self, T>) -> Self {
        Self {
            rc: Rc::new(GradedNode {
                maximal_grade_set: gs,
                minimal_grade_set: RefCell::new(GradeSet::empty()),
                vec_space_dim: OnceCell::new(),
                ast_node,
            }),
        }
    }

    gaexpr_unary_ops! {
        rev Reverse id "Reverse (dagger)",
        ginvol GradeInvolution id "Grade involution (main involution)",
        exp Exponential exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Logarithm log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    }

    /// Return the expressions for the base vectors in a vector space of dim `D`
    pub fn basis_vectors<const D: usize>() -> [Self; D]
    where
        T: GradedDataMut,
    {
        array_init::array_init(|i| {
            let mut v = T::init_null_mv(D, &GradeSet::single(1));
            v.grade_slice_mut(1)[i] = 1.0;
            mv(v)
        })
    }

    /// Raise to some power. Shortcut for `exp(log(self) * p)`, usually with `p`
    /// evaluating to a scalar. Therefore, please refer to [`Self::log`] and
    /// [`Self::exp`] for limitations
    pub fn pow(self, p: impl Into<GaExpr<T>>) -> Self {
        GaExpr::exp(GaExpr::log(self) * p)
    }

    // Square root
    pub fn sqrt(self) -> Self
    where
        T: GradedDataMut,
    {
        if self.grade_set().is_just(0) {
            let gs = self.grade_set().clone();
            Self::new(gs, N::ScalarUnaryOp(ScalarUnaryOp::SquareRoot, self))
        } else {
            self.pow(0.5)
        }
    }

    /// Grade projection: a.g(k) = \<a\>_k
    pub fn g(self, k: i64) -> Self {
        self.gselect(|_| GradeSet::single(k))
    }

    /// Grade projection: select specific grades.
    pub fn gselect(self, f: impl FnOnce(&GradeSet) -> GradeSet) -> Self {
        let wanted = f(&self.grade_set());
        let gs = self.grade_set().clone().intersection(wanted);
        Self::new(gs.clone(), N::GradeProjection(self, gs))
    }

    /// Scalar product. Just a shortcut for `(self.rev() * rhs).g(0)`
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).g(0)
    }

    /// Clifford conjugate. Just a shortcut for `self.rev().ginvol()`
    pub fn conj(self) -> Self {
        self.rev().ginvol()
    }

    /// Inverse. In the case of a multivector, is just a shortcut for
    /// `self.clone().rev() * self.norm_sq().inv()`. If self is just a scalar,
    /// it will compute its regular inverse
    pub fn inv(self) -> Self {
        if self.grade_set().is_just(0) {
            // Regular scalar inversion
            let gs = self.grade_set().clone();
            Self::new(gs, N::ScalarUnaryOp(ScalarUnaryOp::Inversion, self))
        } else {
            self.clone().rev() * self.norm_sq().inv()
        }
    }

    /// Norm squared. Just a shortcut for `(self.clone().rev() * self).g(0)`
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }

    /// Tells which algebra this [`GaExpr`] is using, and recursively propagates
    /// wanted grades downwards so as to evaluate for each sub-expression only
    /// the grades that are necessary to compute the whole [`GaExpr`]
    ///
    /// The given [`MetricAlgebra`] must make sense with respect to the input
    /// values contained in the [`GaExpr`], in terms of possible grades
    /// contained in those input values, and of number of components for each
    /// grade
    pub fn prepare(self, alg: &impl MetricAlgebra) -> ReadyGaExpr<T> {
        // The process is split in two because sub-expressions may be used at
        // several places in the AST. First we collect all the requirements for
        // expressions throughout the whole tree, and add them to the minimal
        // grade sets of the corresponding nodes:
        self.rec_update_minimal_grade_sets(&self.grade_set());
        // Then for each node we apply and store at each node what is needed
        // from the algebra:
        self.rec_apply_algebra(alg);
        // GaExpr<T> and ReadyGaExpr<T> have the exact same memory
        // representation, therefore this is safe:
        unsafe { std::mem::transmute(self) }
    }

    fn rec_update_minimal_grade_sets(&self, wanted: &GradeSet) {
        // Here, we know the resulting grade of the operation represented by the
        // top node of the AST, and we recursively "undo" the AST it to find the
        // grades wanted in each subexpression (the minimal set of grades they
        // should evaluate to given the final result we want)
        self.rc
            .minimal_grade_set
            .replace_with(|old| old.clone() + wanted.clone());
        match &self.rc.ast_node {
            N::GradedObj(_) => {}
            N::GradeProjection(e, _)
            | N::Negation(e)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarUnaryOp(_, e) => {
                e.rec_update_minimal_grade_sets(wanted);
            }
            N::Addition(e1, e2) => {
                // <A + B>_k = <A>_k + <B>_k
                e1.rec_update_minimal_grade_sets(wanted);
                e2.rec_update_minimal_grade_sets(wanted);
            }
            N::GeometricProduct(_, e1, e2) => {
                // Find in e1 and e2 which grades, once multiplied, will affect
                // the grades in `wanted`
                let (e1_wanted, e2_wanted) =
                    wanted.parts_contributing_to_gp(&e1.grade_set(), &e2.grade_set());
                e1.rec_update_minimal_grade_sets(&e1_wanted);
                e2.rec_update_minimal_grade_sets(&e2_wanted);
            }
            N::Exponential(e) => e.rec_update_minimal_grade_sets(&wanted.clone().log()),
            N::Logarithm(e) => e.rec_update_minimal_grade_sets(&wanted.clone().exp()),
        }
    }

    fn minimal_grade_set(&self) -> Ref<'_, GradeSet> {
        self.rc.minimal_grade_set.borrow()
    }

    fn restrict_minimal_grade_set(&self, constraint: GradeSet) {
        self.rc
            .minimal_grade_set
            .replace_with(|orig| orig.clone().intersection(constraint));
    }

    fn rec_apply_algebra(&self, alg: &impl MetricAlgebra) {
        if let Err(_) = self.rc.vec_space_dim.set(alg.vec_space_dim()) {
            // Grade hints and algebra have already been applied for this node
            // (and its subnodes), because it is referred to in at least one
            // other part of the AST
            return;
        }
        self.restrict_minimal_grade_set(
            self.rc
                .maximal_grade_set
                .clone()
                .intersection(alg.full_grade_set()),
        );
        match &self.rc.ast_node {
            N::GradedObj(_) => {}
            N::Negation(e)
            | N::GradeProjection(e, _)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarUnaryOp(_, e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone());
            }
            N::Exponential(e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone().exp())
            }
            N::Logarithm(e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone().log())
            }
            N::Addition(e_left, e_right) => {
                e_left.rec_apply_algebra(alg);
                e_right.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(
                    e_left.minimal_grade_set().clone() + e_right.minimal_grade_set().clone(),
                );
            }
            N::GeometricProduct(comp_muls_cell, e_left, e_right) => {
                e_left.rec_apply_algebra(alg);
                e_right.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(
                    e_left.minimal_grade_set().clone() * e_right.minimal_grade_set().clone(),
                );
                // Now that the grades at play for this geometric product are
                // fully resolved, we can construct the set of
                // component-to-component multiplications that will be needed to
                // perform it:
                comp_muls_cell
                    .set(
                        iter_individual_comp_muls(
                            alg,
                            &self.minimal_grade_set(),
                            &e_left.minimal_grade_set(),
                            &e_right.minimal_grade_set(),
                        )
                        .collect(),
                    )
                    .expect("IndividualCompMul cell has already been set");
            }
        }
    }
}

fn iter_individual_comp_muls<'a>(
    alg: &'a impl MetricAlgebra,
    gs_product: &'a GradeSet,
    gs_left: &'a GradeSet,
    gs_right: &'a GradeSet,
) -> impl Iterator<Item = IndividualCompMul> + 'a {
    gs_product
        .iter_contribs_to_gp(gs_left, gs_right)
        .flat_map(move |(_, k_left, k_right)| {
            iter_basis_blades_of_grade(alg, k_left).flat_map(move |bb_left| {
                iter_basis_blades_of_grade(alg, k_right).filter_map(move |bb_right| {
                    let (bb_res, coeff) = alg.ortho_basis_blades_gp(&bb_left, &bb_right);
                    let result_comp = alg.basis_blade_to_component(&bb_res);
                    if gs_product.contains(result_comp.grade) {
                        Some(IndividualCompMul {
                            left_comp: alg.basis_blade_to_component(&bb_left),
                            right_comp: alg.basis_blade_to_component(&bb_right),
                            result_comp,
                            coeff,
                        })
                    } else {
                        None
                    }
                })
            })
        })
}

/// A [`GaExpr`] which is ready for compilation/evaluation
#[derive(Debug)]
#[repr(transparent)]
pub struct ReadyGaExpr<T> {
    rc: Rc<GradedNode<Self, T>>,
}

/// Expressions are identifiable by reference, this allows to do caching when
/// evaluating an AST
#[derive(Hash, PartialEq, Eq)]
pub struct ExprId(
    /// A pointer that is just here for its hash/eq instances, it will never be dereferenced
    *const (),
);

impl<T> ReadyGaExpr<T> {
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
    pub fn ast_node(&self) -> &AstNode<ReadyGaExpr<T>, T> {
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
/// to what is available in the algebra given to [`GaExpr::prepare`]
impl<T> Graded for ReadyGaExpr<T> {
    type RefToGradeSet<'a> = Ref<'a, GradeSet> where T: 'a;
    fn grade_set(&self) -> Self::RefToGradeSet<'_> {
        // This corresponds to GaExpr::minimal_grade_set
        self.rc.minimal_grade_set.borrow()
    }
}
