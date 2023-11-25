//! Represent a GA expression using an abstract syntax tree

use crate::graded::GradedOutput;

use super::{grade_set::*, graded::Graded};
use std::{borrow::Borrow as _, cell::RefCell, rc::Rc};
use AstNode as N;

/// The abstract syntax tree nodes representing geometric algebra primitive
/// operations. `T` is some raw multivector type, and `E` is a boxed type itself
/// containing an `AstNode`.
#[derive(Hash, Debug)]
pub enum AstNode<E, T> {
    /// Use of a raw multivector which exposes which grades it contains. To do
    /// so, most operations on `Ast` require `T: [Graded]`
    RawMultivector(T),
    /// Multivector addition
    Addition(E, E),
    /// Geometric product
    GeometricProduct(E, E),
    /// Multivector negation
    Negation(E),
    /// Multivector exponentiation. See [`GaExpr::exp`] for limitations
    Exponential(E),
    /// Multivector natural logarithm. See [`GaExpr::log`] for limitations
    Logarithm(E),
    /// Grade projection (or "grade extraction"). The grade to extract is stored
    /// here only for error-reporting reasons
    GradeProjection(E, GradeSet),
    /// Reverse (or "dagger")
    Reverse(E),
    /// Grade involution (or main involution)
    GradeInvolution(E),
    /// Invert the scalar part
    ScalarInversion(E),
}

/// Assign a [`GradeSet`] to some [`AstNode`]. This [`GradeSet`] can be modified
/// to be restricted further depending on where this [`AstNode`] is being used,
/// ie. which grades its use sites actually need
#[derive(Debug)]
pub struct GradedNode<E, T> {
    grade_set_cell: RefCell<GradeSet>,
    grade_set_hints: RefCell<Option<GradeSet>>,
    ast_node: AstNode<E, T>,
}

impl<E, T> GradedNode<E, T> {
    /// Get the [`GradeSet`] associated to this node
    pub fn grade_set(&self) -> impl std::ops::Deref<Target = GradeSet> + '_ {
        self.grade_set_cell.borrow()
    }
    /// Get the node and the part of the AST it contains
    pub fn ast_node(&self) -> &AstNode<E, T> {
        &self.ast_node
    }
    fn add_hint(&self, new: GradeSet) {
        self.grade_set_hints
            .replace_with(|mb_hints| match mb_hints {
                Some(old) => Some(old.clone() + new),
                None => Some(new),
            });
    }
    /// Important: apply_hints is idempotent
    fn apply_hints(&self) {
        match self.grade_set_hints.borrow().as_ref() {
            Some(hints) => {
                self.grade_set_cell
                    .replace_with(|g| g.clone().intersection(hints.clone()));
            }
            None => {}
        };
        self.grade_set_hints.replace(None);
    }
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [`GradeSet`] which can then
/// be further restrained depending on the use sites of this GaExpr, and then
/// used to optimize allocations while evaluating the expression
#[derive(Debug)]
pub struct GaExpr<T>(Rc<GradedNode<Self, T>>);

/// Create a GA expression that just returns some pre-evaluated [`Graded`]
/// value (multivector)
pub fn mv<T: Graded>(x: T) -> GaExpr<T> {
    let gs = x.grade_set().borrow().clone();
    GaExpr::wrap(gs, AstNode::RawMultivector(x))
}

/// [`GaExpr`] uses [`Rc`] internally, so you can safely clone it extensively
impl<T> Clone for GaExpr<T> {
    fn clone(&self) -> Self {
        GaExpr(Rc::clone(&self.0))
    }
}

impl<T> std::ops::Deref for GaExpr<T> {
    type Target = GradedNode<Self, T>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

macro_rules! gaexpr_binary_ops {
    ($($trait:ident $method:ident $ctor:ident $doc:literal),*) => {
        $(
        #[doc=$doc]
        impl<T, E: Into<GaExpr<T>>> std::ops::$trait<E> for GaExpr<T> {
            type Output = Self;
            #[doc=$doc]
            fn $method(self, rhs: E) -> Self::Output {
                let e_rhs = rhs.into();
                let gs = self.grade_set().clone().$method(e_rhs.grade_set().clone());
                Self::wrap(gs, N::$ctor(self, e_rhs))
            }
        }
        )*
    };
}
gaexpr_binary_ops! {
    Add add Addition "Multivector addition",
    Mul mul GeometricProduct "Geometric product",
    BitXor bitxor GeometricProduct "Outer product",
    BitAnd bitand GeometricProduct "Inner product",
    Shl shl GeometricProduct "Left contraction",
    Shr shr GeometricProduct "Right contraction"
}

impl<T> std::ops::Neg for GaExpr<T> {
    type Output = Self;
    fn neg(self) -> Self {
        let gs = self.grade_set().clone();
        Self::wrap(gs, N::Negation(self))
    }
}

impl<T, E: Into<GaExpr<T>>> std::ops::Sub<E> for GaExpr<T> {
    type Output = Self;
    fn sub(self, rhs: E) -> Self::Output {
        self + -rhs.into()
    }
}

impl<T: GradedOutput> From<f64> for GaExpr<T> {
    fn from(x: f64) -> GaExpr<T> {
        if x == 0.0 {
            return mv(T::init_null_mv(0, &GradeSet::empty()));
        }
        let mut scal_mv = T::init_null_mv(0, &GradeSet::single(0));
        scal_mv.grade_slice_mut(0)[0] = x;
        mv(scal_mv)
    }
}
impl<T: GradedOutput> From<i64> for GaExpr<T> {
    #[inline]
    fn from(x: i64) -> GaExpr<T> {
        (x as f64).into()
    }
}

macro_rules! scalar_with_gaexpr_binary_ops {
    ($($t:ty),*) => {
        $(
        impl<T: GradedOutput> std::ops::Add<GaExpr<T>> for $t {
            type Output = GaExpr<T>;
            #[inline]
            fn add(self, rhs: GaExpr<T>) -> Self::Output {
                Into::<GaExpr<T>>::into(self) + rhs
            }
        }
        impl<T: GradedOutput> std::ops::Mul<GaExpr<T>> for $t {
            type Output = GaExpr<T>;
            #[inline]
            fn mul(self, rhs: GaExpr<T>) -> Self::Output {
                Into::<GaExpr<T>>::into(self) * rhs
            }
        }
        impl<T: GradedOutput> std::ops::Div<$t> for GaExpr<T> {
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
                Self::wrap(
                    gs,
                    N::$ctor(self),
                )
            }
        )*
    }
}

impl<T> GaExpr<T> {
    fn wrap(gs: GradeSet, ast: AstNode<Self, T>) -> Self {
        Self(Rc::new(GradedNode {
            grade_set_cell: RefCell::new(gs),
            grade_set_hints: RefCell::new(None),
            ast_node: ast,
        }))
    }

    gaexpr_unary_ops! {
        rev Reverse id "Reverse (dagger)",
        ginvol GradeInvolution id "Grade involution (main involution)",
        exp Exponential exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Logarithm log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    }

    /// Raise to some power. Shortcut for `exp(log(self) * p)`, usually with `p`
    /// evaluating to a scalar. Therefore, please refer to [`Self::log`] and
    /// [`Self::exp`] for limitations
    pub fn pow(self, p: GaExpr<T>) -> Self {
        GaExpr::exp(GaExpr::log(self) * p)
    }

    /// Grade projection: a.g(k) = \<a\>_k
    pub fn g(self, k: i64) -> Self {
        self.gselect(|_| GradeSet::single(k))
    }

    /// Grade projection: select specific grades.
    pub fn gselect(self, f: impl FnOnce(&GradeSet) -> GradeSet) -> Self {
        let wanted = f(&self.grade_set());
        let gs = self.grade_set().clone().intersection(wanted);
        Self::wrap(gs.clone(), N::GradeProjection(self, gs))
    }

    /// Scalar product. Just a shortcut for `(self.rev() * rhs).prj(0)`
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).g(0)
    }

    /// Clifford conjugate. Just a shortcut for `self.rev().ginvol()`
    pub fn conj(self) -> Self {
        self.rev().ginvol()
    }

    /// Inverse. Just a shortcut for `self.clone().rev() * self.norm_sq().inv()`
    pub fn inv(self) -> Self {
        if self.grade_set().is_just(0) {
            // Regular scalar inversion
            let gs = self.grade_set().clone();
            Self::wrap(gs, N::ScalarInversion(self))
        } else {
            self.clone().rev() * self.norm_sq().inv()
        }
    }

    /// Norm squared. Just a shortcut for `(self.clone().rev() * self).g(0)`
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }

    /// Recursively propagate wanted grades downwards so as to evaluate for each
    /// sub-expression only the grades that are necessary to compute the whole
    /// [`GaExpr`]
    pub fn minimize_grades(self) -> Self {
        // The process is split in two because sub-expressions may be used at
        // several places in the AST. First we collect all the requirements for
        // expressions throughout the whole tree, as hints to be applied later:
        self.propagate_grade_hints(&self.grade_set());
        // Then for each node we apply the hints collected previously:
        self.apply_grade_hints();
        self
    }

    fn propagate_grade_hints(&self, wanted: &GradeSet) {
        // Here, we know the resulting grade of the operation represented by the
        // top node of the AST, and we "undo" it to find the grades wanted in
        // its operands
        self.add_hint(wanted.clone());
        match self.ast_node() {
            N::RawMultivector(_) => {}
            N::GradeProjection(e, _)
            | N::Negation(e)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarInversion(e) => {
                e.propagate_grade_hints(wanted);
            }
            N::Addition(e1, e2) => {
                // <A + B>_k = <A>_k + <B>_k
                e1.propagate_grade_hints(wanted);
                e2.propagate_grade_hints(wanted);
            }
            N::GeometricProduct(e1, e2) => {
                // Find in e1 and e2 which grades, once multiplied, will affect
                // the grades in `wanted`
                let (e1_wanted, e2_wanted) =
                    wanted.parts_contributing_to_gp(&e1.grade_set(), &e2.grade_set());
                e1.propagate_grade_hints(&e1_wanted);
                e2.propagate_grade_hints(&e2_wanted);
            }
            N::Exponential(e) => e.propagate_grade_hints(&wanted.clone().log()),
            N::Logarithm(e) => e.propagate_grade_hints(&wanted.clone().exp()),
        }
    }

    fn apply_grade_hints(&self) {
        self.apply_hints();
        match self.ast_node() {
            N::RawMultivector(_) => {}
            N::Negation(e)
            | N::GradeProjection(e, _)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarInversion(e)
            | N::Exponential(e)
            | N::Logarithm(e) => {
                e.apply_grade_hints();
            }
            N::Addition(e1, e2) | N::GeometricProduct(e1, e2) => {
                e1.apply_grade_hints();
                e2.apply_grade_hints();
            }
        }
    }
}
