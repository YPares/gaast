//! Represent a GA expression using an abstract syntax tree

use super::grade_set::*;
use super::graded::Graded;
use std::rc::Rc;
use Ast::*;

#[derive(Hash)]
pub(crate) enum Ast<E, T> {
    Val(T),
    Add(E, E),
    /// The geometric product
    Mul(E, E),
    /// Negate
    Neg(E),
    Exp(E),
    Log(E),
    /// The reverse (or "dagger")
    Rev(E),
    /// The grade involution (or main involution)
    GInvol(E),
    /// Regular scalar inversion _on a expression that evaluates to a scalar_
    ScalarInv(E),
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [GradeSet] which can then be
/// further restrained depending on the use sites of this GAExpr, and then used
/// to optimize allocations while evaluating the expression
#[derive(Clone)]
pub struct GAExpr<T> {
    pub(crate) grade_set: GradeSet,
    pub(crate) ast: Rc<Ast<GAExpr<T>, T>>,
}

impl<T> std::ops::Add for GAExpr<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            grade_set: self.grade_set.clone() + rhs.grade_set.clone(),
            ast: Rc::new(Add(self, rhs)),
        }
    }
}

impl<T> std::ops::Mul for GAExpr<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            grade_set: self.grade_set.clone() * rhs.grade_set.clone(),
            ast: Rc::new(Mul(self, rhs)),
        }
    }
}

macro_rules! unary_ops {
    ($($fn_name:ident $ctor:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                Self {
                    grade_set: self.grade_set.clone().$grade_op(),
                    ast: Rc::new($ctor(self))
                }
            }
        )*
    }
}

impl<T> GAExpr<T> {
    unary_ops!(
        rev Rev id "Reverse (dagger)",
        ginvol GInvol id "Grade involution (main involution)",
        exp Exp exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Log log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    );

    /// Grade projection: a.prj(k) = \<a\>_k
    pub fn prj(mut self, k: usize) -> Self {
        self.grade_set = GradeSet::g(k);
        self
    }

    /// Scalar product
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).prj(0)
    }
}

impl<T: Graded> GAExpr<T> {
    /// Create a GA expression from a raw value
    pub fn val(x: T) -> Self {
        Self {
            grade_set: x.grade_set(),
            ast: Rc::new(Ast::Val(x)),
        }
    }

    /// To some floating-point power. `p` must therefore evaluate to a scalar.
    /// Refer to [`Self::log`] for limitations
    pub fn pow(self, p: GAExpr<T>) -> Self {
        assert!(p.grade_set == GradeSet::g(0));
        GAExpr::exp(GAExpr::log(self) * p)
    }
}

impl<T: Clone> GAExpr<T> {
    /// Inverse
    pub fn inv(self) -> Self {
        if self.grade_set == GradeSet::g(0) {
            // Regular scalar inversion
            Self {
                grade_set: self.grade_set.clone(),
                ast: Rc::new(ScalarInv(self)),
            }
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
        Self {
            grade_set: self.grade_set.clone(),
            ast: Rc::new(Neg(self)),
        }
    }
}

impl<T: Graded + Clone> std::ops::Sub for GAExpr<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}
