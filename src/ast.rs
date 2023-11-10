use super::grade_set::*;
use std::rc::Rc;
use Ast::*;

enum Ast {
    Val(Rc<dyn Graded>),
    Add(AstRef, AstRef),
    /// The geometric product
    Mul(AstRef, AstRef),
    Exp(AstRef),
    Log(AstRef),
    /// The reverse (or "dagger")
    Rev(AstRef),
    /// The grade involution (or main involution)
    GInvol(AstRef),
    /// Invert the scalar part of an expression
    ScalarInv(AstRef),
}

type AstRef = Rc<Ast>;

#[derive(Clone)]
/// An AST representing primitive GA operations
pub struct GAExpr {
    /// The grade set is recursively inferred from the AST when constructing the
    /// GAExpr
    grade_set: GradeSet,
    ast: AstRef,
}

impl std::ops::Add for GAExpr {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        GAExpr {
            grade_set: self.grade_set + rhs.grade_set,
            ast: Rc::new(Add(self.ast, rhs.ast)),
        }
    }
}

impl std::ops::Mul for GAExpr {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        GAExpr {
            grade_set: self.grade_set * rhs.grade_set,
            ast: Rc::new(Mul(self.ast, rhs.ast)),
        }
    }
}

impl std::ops::Neg for GAExpr {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self * GAExpr::val(-1.0)
    }
}

impl std::ops::Sub for GAExpr {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

macro_rules! grade_op {
    ($_:ident, id) => {};
    ($x:ident, $fn_name:ident) => {
        $x.grade_set = $x.grade_set.$fn_name();
    };
}
macro_rules! unary_ops {
    ($($fn_name:ident $ctor:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(mut self) -> Self {
                grade_op!(self, $grade_op);
                self.ast = Rc::new($ctor(self.ast));
                self
            }
        )*
    }
}

impl GAExpr {
    /// Create a GA expression from a raw value
    pub fn val(x: impl Graded + 'static) -> Self {
        GAExpr {
            grade_set: x.grade_set(),
            ast: Rc::new(Val(Rc::new(x))),
        }
    }

    unary_ops!(
        rev Rev id "Reverse (dagger)",
        ginvol GInvol id "Grade involution (main involution)",
        exp Exp exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors that square to a scalar",
        log Log log "Natural logarithm"
    );

    /// Inverse
    pub fn inv(mut self) -> Self {
        if self.grade_set == GradeSet::g(0) {
            // Regular scalar inversion
            self.ast = Rc::new(ScalarInv(self.ast));
            self
        } else {
            self.clone().rev() * self.norm_sq().inv()
        }
    }

    /// Grade projection
    pub fn prj(mut self, k: usize) -> Self {
        self.grade_set = GradeSet::g(k);
        self
    }

    /// Scalar product
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).prj(0)
    }

    /// Norm squared
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }

    /// To some floating-point power
    pub fn powf(self, p: f64) -> Self {
        GAExpr::exp(GAExpr::log(self) * GAExpr::val(p))
    }

    /// Square root
    pub fn sqrt(self) -> Self {
        self.powf(0.5)
    }
}
