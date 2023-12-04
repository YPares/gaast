use super::base_types::{AstNode as N, *};
use crate::{grade_set::*, graded::GradedDataMut, Graded};
use std::{
    cell::{OnceCell, RefCell},
    fmt::Debug,
    rc::Rc,
};

/// Assign a [`GradeSet`] to some [`AstNode`]. This [`GradeSet`] can be modified
/// to be restricted further depending on where this [`AstNode`] is being used,
/// ie. which grades its use sites actually need
#[derive(Debug)]
pub(super) struct GradedNode<E, T> {
    /// The grade set inferred at AST construction (inferred from the maximal
    /// grade sets of this node's subexpressions)
    pub(super) maximal_grade_set: GradeSet,
    /// Starting at empty, is updated during AST specialization, to finally
    /// reflect on the minimal set of grades to allocate/compute at evaluation
    pub(super) minimal_grade_set: RefCell<GradeSet>,
    /// Starting at empty, is set during AST specialization given the algebra
    pub(super) vec_space_dim: OnceCell<usize>,
    /// This node
    pub(super) ast_node: AstNode<E, T>,
}

/// An AST representing primitive GA operations, which allows to infer the
/// grades contained in the multivector returned by the expression without
/// having to evaluate it. This is represented by a [`GradeSet`] which can then
/// be further restrained depending on the use sites of this GaExpr, and then
/// used to optimize allocations while evaluating the expression
#[derive(Debug)]
#[repr(transparent)]
pub struct GaExpr<T> {
    pub(super) rc: Rc<GradedNode<Self, T>>,
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
    ($($doc:literal $trait:ident $method:ident ($fn:expr)),*) => {
        $(
        #[doc=$doc]
        impl<T, E: Into<GaExpr<T>>> std::ops::$trait<E> for GaExpr<T> {
            type Output = Self;
            #[doc=$doc]
            fn $method(self, rhs: E) -> Self::Output {
                self.product(rhs.into(), $fn)
            }
        }
        )*
    };
}

gaexpr_products! {
    "Geometric product" Mul mul
        (|(k1, k2)| GradeSet::single(k1) * GradeSet::single(k2)),
    "Outer product" BitXor bitxor
        (|(k1, k2)| GradeSet::single(k1 + k2)),
    "Inner product" BitAnd bitand
        (|(k1, k2)| {
            if k1 == 0 || k2 == 0 {
                GradeSet::empty()
            } else {
                GradeSet::single((k1 - k2).abs())
            }
        }),
    "Left contraction" Shl shl
        (|(k1, k2)| GradeSet::single(k2 - k1)),
    "Right contraction" Shr shr
        (|(k1, k2)| GradeSet::single(k1 - k2))
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

    /// Compute a product (based on the geometric product) of two operands,
    /// given a function that will select grades from the result of each
    /// individual k-vector to g-vector product
    pub fn product(
        self,
        rhs: Self,
        grades_to_produce: impl Fn((i64, i64)) -> GradeSet + 'static,
    ) -> Self {
        let gs = iter_grade_sets_cp(self.grade_set(), rhs.grade_set())
            .map(&grades_to_produce)
            .collect();
        Self::new(
            gs,
            N::Product(Product {
                comp_muls_cell: OnceCell::new(),
                grades_to_produce: KVecsProductGradeSelection(Box::new(grades_to_produce)),
                left_expr: self,
                right_expr: rhs,
            }),
        )
    }

    /// Raise to some power. Shortcut for `exp(log(self) * p)`. Therefore,
    /// please refer to [`Self::log`] and [`Self::exp`] for limitations
    pub fn pow(self, p: impl Into<GaExpr<T>>) -> Self {
        GaExpr::exp(GaExpr::log(self) * p)
    }

    /// Square root
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
        Self::new(gs, N::GradeProjection(self))
    }

    /// Scalar product. Just a shortcut for `(self.rev() * rhs).g(0)`
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).g(0)
    }

    /// Clifford conjugate. Just a shortcut for `self.rev().ginvol()`
    pub fn conj(self) -> Self {
        self.rev().ginvol()
    }

    /// Versor inverse. Just a shortcut for `self.clone().rev() *
    /// self.norm_sq().sinv()`. Will default to [`Self::sinv`] when self is just a
    /// scalar
    pub fn vinv(self) -> Self {
        if self.grade_set().is_just(0) {
            self.sinv()
        } else {
            self.clone().rev() * self.norm_sq().sinv()
        }
    }

    /// Inverts only the scalar part
    pub fn sinv(self) -> Self {
        let gs = self.grade_set().clone();
        Self::new(gs, N::ScalarUnaryOp(ScalarUnaryOp::Inversion, self))
    }

    /// Norm squared. Just a shortcut for `(self.clone().rev() * self).g(0)`
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }
}
