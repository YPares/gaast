//! Represent a GA expression using an abstract syntax tree

use super::grade_set::*;
use super::graded::Graded;
use std::cell::RefCell;
use std::rc::Rc;
use Ast::*;

#[derive(Hash, Clone)]
pub(crate) enum Ast<E, T> {
    Val(T),
    Add(E, E),
    /// The geometric product
    Mul(E, E),
    /// Negate
    Neg(E),
    Exp(E),
    Log(E),
    Prj(E, usize),
    /// The reverse (or "dagger")
    Rev(E),
    /// The grade involution (or main involution)
    GInvol(E),
    /// Regular scalar inversion _on a expression that evaluates to a scalar_
    ScalarInv(E),
}

/// Assign a [`GradeSet`] to some [`Ast`]
pub(crate) struct GradedAst<E, T> {
    grade_set_cell: RefCell<GradeSet>,
    grade_set_hints: RefCell<Option<GradeSet>>,
    ast: Ast<E, T>,
}

impl<E, T> GradedAst<E, T> {
    pub(crate) fn grade_set(&self) -> GradeSet {
        self.grade_set_cell.borrow().clone()
    }
    pub(crate) fn ast(&self) -> &Ast<E, T> {
        &self.ast
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
#[derive(Clone)]
pub struct GAExpr<T>(pub(crate) Rc<GradedAst<Self, T>>);

impl<T> std::ops::Add for GAExpr<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::wrap(self.0.grade_set() + rhs.0.grade_set(), Add(self, rhs))
    }
}

impl<T> std::ops::Mul for GAExpr<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::wrap(self.0.grade_set() * rhs.0.grade_set(), Mul(self, rhs))
    }
}

macro_rules! unary_ops {
    ($($fn_name:ident $ctor:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                Self::wrap(
                    self.0.grade_set().$grade_op(),
                    $ctor(self),
                )
            }
        )*
    }
}

impl<T> GAExpr<T> {
    fn wrap(gs: GradeSet, ast: Ast<Self, T>) -> Self {
        Self(Rc::new(GradedAst {
            grade_set_cell: RefCell::new(gs),
            grade_set_hints: RefCell::new(None),
            ast,
        }))
    }

    unary_ops!(
        rev Rev id "Reverse (dagger)",
        ginvol GInvol id "Grade involution (main involution)",
        exp Exp exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Log log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    );

    /// Grade projection: a.prj(k) = \<a\>_k
    pub fn prj(self, k: usize) -> Self {
        Self::wrap(self.0.grade_set().prj(GradeSet::g(k)), Prj(self, k))
    }

    /// Scalar product
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).prj(0)
    }

    /// Recursively propagate wanted grades downwards so as to evaluate for each
    /// sub-expression only the grades that are necessary to compute the whole
    /// [`GAExpr`]
    pub fn resolve_minimum_grades(self) -> Self {
        // The process is split in two because sub-expressions may be used at
        // several places in the AST. First we collect all the requirements for
        // expressions throughout the whole tree, as hints to be applied later:
        self.propagate_grade_hints(self.0.grade_set());
        // Then for each node we apply the hints collected previously:
        self.apply_grade_hints();
        self
    }

    fn propagate_grade_hints(&self, wanted: GradeSet) {
        self.0.add_hint(wanted.clone());
        match self.0.ast() {
            Val(_) => {}
            Neg(e) | Prj(e, _) | Rev(e) | GInvol(e) | ScalarInv(e) => {
                e.propagate_grade_hints(wanted);
            }
            Add(e1, e2) => {
                // <A + B>_k = <A>_k + <B>_k
                e1.propagate_grade_hints(wanted.clone());
                e2.propagate_grade_hints(wanted);
            }
            Mul(e1, e2) => {
                // Find in e1 and e2 which grades, once multiplied, will affect
                // the grades in `wanted`
                let (e1_wanted, e2_wanted) =
                    wanted.grades_affecting_mul(&e1.0.grade_set(), &e2.0.grade_set());
                e1.propagate_grade_hints(e1_wanted);
                e2.propagate_grade_hints(e2_wanted);
            }
            Exp(_) => todo!(),
            Log(_) => todo!(),
        }
    }

    fn apply_grade_hints(&self) {
        self.0.apply_hints();
        match self.0.ast() {
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

impl<T: Graded> GAExpr<T> {
    /// Create a GA expression from a raw value
    pub fn val(x: T) -> Self {
        Self::wrap(x.grade_set(), Ast::Val(x))
    }

    /// To some floating-point power. `p` must therefore evaluate to a scalar.
    /// Refer to [`Self::log`] for limitations
    pub fn pow(self, p: GAExpr<T>) -> Self {
        assert!(p.0.grade_set() == GradeSet::g(0));
        GAExpr::exp(GAExpr::log(self) * p)
    }
}

impl<T: Clone> GAExpr<T> {
    /// Inverse
    pub fn inv(self) -> Self {
        if self.0.grade_set() == GradeSet::g(0) {
            // Regular scalar inversion
            Self::wrap(self.0.grade_set(), ScalarInv(self))
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
        Self::wrap(self.0.grade_set(), Neg(self))
    }
}

impl<T: Graded + Clone> std::ops::Sub for GAExpr<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}
