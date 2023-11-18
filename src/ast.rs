//! Represent a GA expression using an abstract syntax tree

use super::{grade_set::*, graded::Graded};
use std::cell::RefCell;
use std::rc::Rc;
use AstNode::*;

/// The abstract syntax tree nodes representing geometric algebra primitive
/// operations. `T` is some raw multivector type, and `E` is a boxed type itself
/// containing an `AstNode`.
#[derive(Hash, Clone, Debug)]
pub enum AstNode<E, T> {
    /// Use of a raw multivector which exposes which grades it contains. To do
    /// so, most operations on `Ast` require `T: [Graded]`
    Val(T),
    /// Multivector addition
    Add(E, E),
    /// Geometric product
    Mul(E, E),
    /// Multivector negation
    Neg(E),
    /// Multivector exponentiation. See [`GAExpr::exp`] for limitations
    Exp(E),
    /// Multivector natural logarithm. See [`GAExpr::log`] for limitations
    Log(E),
    /// Grade projection (or "grade extraction"). The grade to extract is stored
    /// here only for error-reporting reasons
    Prj(E, Grade),
    /// Reverse (or "dagger")
    Rev(E),
    /// Grade involution (or main involution)
    GInvol(E),
    /// Invert the scalar part
    ScalarInv(E),
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
                    .replace_with(|g| g.clone().prj(hints.clone()));
            }
            None => {}
        };
        self.grade_set_hints.replace(None);
    }
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [`GradeSet`] which can then
/// be further restrained depending on the use sites of this GAExpr, and then
/// used to optimize allocations while evaluating the expression
///
/// `GAExpr` dereferences to a [`GradedNode`], so you can call [`GradedNode`]
/// methods on a `GAExpr`
#[derive(Clone, Debug)]
pub struct GAExpr<T>(Rc<GradedNode<Self, T>>);

impl<T> std::ops::Deref for GAExpr<T> {
    type Target = GradedNode<Self, T>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> std::ops::Add for GAExpr<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let gs = self.grade_set().clone() + rhs.grade_set().clone();
        Self::wrap(gs, Add(self, rhs))
    }
}

impl<T> std::ops::Mul for GAExpr<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let gs = self.grade_set().clone() * rhs.grade_set().clone();
        Self::wrap(gs, Mul(self, rhs))
    }
}

macro_rules! unary_ops {
    ($($fn_name:ident $ctor:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                let gs = self.grade_set().clone().$grade_op();
                Self::wrap(
                    gs,
                    $ctor(self),
                )
            }
        )*
    }
}

impl<T: Graded> GAExpr<T> {
    /// Create a GA expression from a raw input multivector value
    pub fn val(x: T) -> Self {
        Self::wrap(x.grade_set(), AstNode::Val(x))
    }

    /// To some floating-point power. `p` must therefore evaluate to a scalar.
    /// Refer to [`Self::log`] for limitations
    pub fn pow(self, p: GAExpr<T>) -> Self {
        assert!(
            p.grade_set().is_just(0),
            "Pow can only be applied if the expression in exponent evaluates to a scalar"
        );
        GAExpr::exp(GAExpr::log(self) * p)
    }
}

impl<T> GAExpr<T> {
    fn wrap(gs: GradeSet, ast: AstNode<Self, T>) -> Self {
        Self(Rc::new(GradedNode {
            grade_set_cell: RefCell::new(gs),
            grade_set_hints: RefCell::new(None),
            ast_node: ast,
        }))
    }

    unary_ops!(
        rev Rev id "Reverse (dagger)",
        ginvol GInvol id "Grade involution (main involution)",
        exp Exp exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Log log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    );

    /// Grade projection: a.prj(k) = \<a\>_k
    pub fn prj(self, k: Grade) -> Self {
        let gs = self.grade_set().clone().prj(GradeSet::single(k));
        Self::wrap(gs, Prj(self, k))
    }

    /// Scalar product
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).prj(0)
    }

    /// Clifford conjugate (combination of reverse & grade involution)
    pub fn conj(self) -> Self {
        self.rev().ginvol()
    }

    /// Recursively propagate wanted grades downwards so as to evaluate for each
    /// sub-expression only the grades that are necessary to compute the whole
    /// [`GAExpr`]
    pub fn resolve_minimum_grades(self) -> Self {
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
            Val(_) => {}
            Prj(e, _) | Neg(e) | Rev(e) | GInvol(e) | ScalarInv(e) => {
                e.propagate_grade_hints(wanted);
            }
            Add(e1, e2) => {
                // <A + B>_k = <A>_k + <B>_k
                e1.propagate_grade_hints(wanted);
                e2.propagate_grade_hints(wanted);
            }
            Mul(e1, e2) => {
                // Find in e1 and e2 which grades, once multiplied, will affect
                // the grades in `wanted`
                let (e1_wanted, e2_wanted) =
                    wanted.grades_affecting_mul(&e1.grade_set(), &e2.grade_set());
                e1.propagate_grade_hints(&e1_wanted);
                e2.propagate_grade_hints(&e2_wanted);
            }
            Exp(e) => e.propagate_grade_hints(&wanted.clone().log()),
            Log(e) => e.propagate_grade_hints(&wanted.clone().exp()),
        }
    }

    fn apply_grade_hints(&self) {
        self.apply_hints();
        match self.ast_node() {
            Val(_) => {}
            Neg(e) | Prj(e, _) | Rev(e) | GInvol(e) | ScalarInv(e) | Exp(e) | Log(e) => {
                e.apply_grade_hints();
            }
            Add(e1, e2) | Mul(e1, e2) => {
                e1.apply_grade_hints();
                e2.apply_grade_hints();
            }
        }
    }
}

impl<T: Clone> GAExpr<T> {
    /// Inverse
    pub fn inv(self) -> Self {
        if self.grade_set().is_just(0) {
            // Regular scalar inversion
            let gs = self.grade_set().clone();
            Self::wrap(gs, ScalarInv(self))
        } else {
            self.clone().rev() * self.norm_sq().inv()
        }
    }

    /// Norm squared
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }
}

impl<T: Graded + Clone> std::ops::Neg for GAExpr<T> {
    type Output = Self;
    fn neg(self) -> Self {
        let gs = self.grade_set().clone();
        Self::wrap(gs, Neg(self))
    }
}

impl<T: Graded + Clone> std::ops::Sub for GAExpr<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}
