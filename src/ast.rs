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

/// Grade preserving operations
macro_rules! pub_gpos {
    ($($fn_name:ident : $ctor:ident)*) => {
        $(
            pub fn $fn_name(mut self) -> Self {
                self.ast = Rc::new($ctor(self.ast));
                self
            }
        )*
    };
}

impl GAExpr {
    /// Create a GA expression from a raw value
    pub fn val(x: impl Graded + 'static) -> Self {
        GAExpr {
            grade_set: x.grade_set(),
            ast: Rc::new(Val(Rc::new(x))),
        }
    }

    pub_gpos!(exp:Exp log:Log rev:Rev ginvol:GInvol);

    /// Inverse
    pub fn inv(mut self) -> Self {
        if self.grade_set == GradeSet::g(0) {
            // Regular scalar inversion
            self.ast = Rc::new(ScalarInv(self.ast));
            self
        } else {
            self.clone() * self.norm_sq().inv()
        }
    }

    /// Grade projection
    pub fn prj(mut self, k: usize) -> Self {
        self.grade_set = GradeSet::g(k);
        self
    }

    /// Norm squared
    pub fn norm_sq(self) -> Self {
        (self.clone().rev() * self).prj(0)
    }
}
