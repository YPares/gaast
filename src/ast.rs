//! Represent a GA expression using an abstract syntax tree

use super::grade_set::*;
use std::rc::Rc;
use Ast::*;

/// The trait for all objects that are graded
pub trait Graded {
    /// Get the GradeSet of the object
    fn grade_set(&self) -> GradeSet;
    /// Create an instance of the graded object with just a grade 0 part
    fn from_scalar(x: f64) -> Self;
}

macro_rules! Graded_blanket_impls {
    ($($ref:tt),*) => {
        $(impl<T: Graded> Graded for $ref<T> {
            fn grade_set(&self) -> GradeSet {
                (**self).grade_set()
            }
            fn from_scalar(x: f64) -> Self {
                $ref::new(T::from_scalar(x))
            }
        })*
    };
}
Graded_blanket_impls!(Rc, Box);

impl Graded for f64 {
    fn grade_set(&self) -> GradeSet {
        GradeSet::g(0)
    }
    fn from_scalar(x: f64) -> Self {
        x
    }
}

pub(crate) enum Ast<T> {
    Val(T),
    Add(GAExpr<T>, GAExpr<T>),
    /// The geometric product
    Mul(GAExpr<T>, GAExpr<T>),
    Exp(GAExpr<T>),
    Log(GAExpr<T>),
    /// The reverse (or "dagger")
    Rev(GAExpr<T>),
    /// The grade involution (or main involution)
    GInvol(GAExpr<T>),
    /// Regular scalar inversion _on a expression that evaluates to a scalar_
    ScalarInv(GAExpr<T>),
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [GradeSet] which can then be
/// used to optimize allocations while evaluating the expression
#[derive(Clone)]
pub struct GAExpr<T> {
    pub(crate) grade_set: GradeSet,
    pub(crate) ast: Rc<Ast<T>>,
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

    /// To some floating-point power. Refer to [`Self::log`] for limitations
    pub fn powf(self, p: f64) -> Self {
        GAExpr::exp(GAExpr::log(self) * GAExpr::val(T::from_scalar(p)))
    }

    /// Square root. Refer to [`Self::log`] for limitations
    pub fn sqrt(self) -> Self {
        self.powf(0.5)
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
    fn neg(self) -> Self::Output {
        self * GAExpr::val(T::from_scalar(-1.0))
    }
}

impl<T: Graded + Clone> std::ops::Sub for GAExpr<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}
