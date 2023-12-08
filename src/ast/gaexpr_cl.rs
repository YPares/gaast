use super::{
    base_types::*,
    isolate::{IsolatedGaExpr, NodeArena},
    ExprId, IsolatedNode,
};
use crate::{algebra::MetricAlgebra, grade_set::*, graded::GradedDataMut, Graded};
use std::{collections::HashMap, rc::Rc};

struct Builder<'a, T> {
    algebra: &'a dyn MetricAlgebra,
    arena: &'a mut NodeArena<T>,
}

impl<'a, T> Builder<'a, T> {
    fn add_node(&mut self, node_id: ExprId, (ast_node, node_gs): (AstNode<ExprId, T>, GradeSet)) {
        self.arena.insert(
            node_id,
            IsolatedNode {
                maximal_grade_set: node_gs.intersection(self.algebra.full_grade_set()),
                minimal_grade_set: GradeSet::empty(),
                vec_space_dim: Some(self.algebra.vec_space_dim()),
                ast_node,
                num_uses: 1,
            },
        );
    }
}

pub struct GaExpr<'a, T> {
    build: Rc<dyn Fn(ExprId, &mut Builder<'_, T>) -> () + 'a>,
}

impl<T> Clone for GaExpr<'_, T> {
    fn clone(&self) -> Self {
        Self {
            build: Rc::clone(&self.build),
        }
    }
}

enum Wrapper<'a, T> {
    Node((AstNode<ExprId, T>, GradeSet)),
    Expr(GaExpr<'a, T>),
}

impl<'a, T: 'a> GaExpr<'a, T> {
    pub(super) fn isolate(self, alg: impl MetricAlgebra) -> IsolatedGaExpr<T> {
        let mut arena = HashMap::new();
        let (root, _) = self.build_or_reuse(&mut Builder {
            algebra: &alg,
            arena: &mut arena,
        });
        IsolatedGaExpr { arena, root }
    }

    /// Adds a node to the arena (or does nothing if it already exists) and
    /// returns its identifier and GradeSet
    fn build_or_reuse(&self, b: &mut Builder<'_, T>) -> (ExprId, GradeSet) {
        let id = ExprId {
            ptr: Rc::as_ptr(&self.build).cast(),
        };
        match b.arena.get_mut(&id) {
            None => (self.build)(id, b),
            Some(node) => {
                node.num_uses += 1;
            }
        }
        (id, b.arena[&id].maximal_grade_set.clone())
    }

    fn new(f: impl Fn(&mut Builder<'_, T>) -> (AstNode<ExprId, T>, GradeSet) + 'a) -> Self {
        Self {
            build: Rc::new(move |this_id, b| {
                let node_and_gs = f(b);
                b.add_node(this_id, node_and_gs)
            }),
        }
    }

    /// Enables to wrap a subexpression either in a new node, or by reusing
    /// pre-existing functions returning GaExprs
    fn wrap(self, f: impl Fn(Self, ExprId, GradeSet) -> Wrapper<'a, T> + 'a) -> Self {
        Self {
            build: Rc::new(move |wrapper_id, b| {
                let (self_id, self_gs) = self.build_or_reuse(b);
                match f(self.clone(), self_id, self_gs) {
                    Wrapper::Node(wrapper_node_and_gs) => {
                        b.add_node(wrapper_id, wrapper_node_and_gs)
                    }
                    Wrapper::Expr(wrapper_expr) => {
                        (wrapper_expr.build)(wrapper_id, b);
                        // We correct for the extra use we added when calling
                        // self.build_or_reuse above:
                        b.arena.get_mut(&self_id).unwrap().num_uses -= 1;
                    }
                }
            }),
        }
    }

    pub fn product(
        self,
        rhs: Self,
        grades_to_produce: impl Fn((i64, i64)) -> GradeSet + 'static,
    ) -> Self {
        let grades_to_produce = KVecsProductGradeSelection(Rc::new(grades_to_produce));
        Self::new(move |b| {
            let (left_id, left_gs) = self.build_or_reuse(b);
            let (right_id, right_gs) = rhs.build_or_reuse(b);
            (
                AstNode::Product(Product {
                    comp_muls_cell: None,
                    grades_to_produce: grades_to_produce.clone(),
                    left_expr: left_id,
                    right_expr: right_id,
                }),
                iter_grade_sets_cp(&left_gs, &right_gs)
                    .map(grades_to_produce.0.as_ref())
                    .collect(),
            )
        })
    }
}

pub fn mv<'a, T: Graded + Clone + 'a>(x: T) -> GaExpr<'a, T> {
    GaExpr::new(move |_| (AstNode::GradedObj(x.clone()), x.grade_set().clone()))
}

macro_rules! gaexpr_products {
    ($($doc:literal $trait:ident $method:ident ($fn:expr)),*) => {
        $(
        #[doc=$doc]
        impl<'a, T: 'a, E: Into<GaExpr<'a, T>>> std::ops::$trait<E> for GaExpr<'a, T> {
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

impl<'a, T: 'a, E: Into<GaExpr<'a, T>>> std::ops::Add<E> for GaExpr<'a, T> {
    type Output = Self;
    fn add(self, rhs: E) -> Self::Output {
        let rhs = rhs.into();
        Self::new(move |b| {
            let (left_id, left_gs) = self.build_or_reuse(b);
            let (right_id, right_gs) = rhs.build_or_reuse(b);
            (AstNode::Addition(left_id, right_id), left_gs + right_gs)
        })
    }
}

impl<'a, T: 'a> std::ops::Neg for GaExpr<'a, T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(move |b| {
            let (id, gs) = self.build_or_reuse(b);
            (AstNode::Negation(id), gs)
        })
    }
}

impl<'a, T: 'a, E: Into<GaExpr<'a, T>>> std::ops::Sub<E> for GaExpr<'a, T> {
    type Output = Self;
    fn sub(self, rhs: E) -> Self::Output {
        self + -rhs.into()
    }
}

impl<'a, T: GradedDataMut + Clone + 'a> From<f64> for GaExpr<'a, T> {
    fn from(x: f64) -> Self {
        if x == 0.0 {
            return mv(T::init_null_mv(0, &GradeSet::empty()));
        }
        let mut scal_mv = T::init_null_mv(0, &GradeSet::single(0));
        scal_mv.grade_slice_mut(0)[0] = x;
        mv(scal_mv)
    }
}
impl<'a, T: GradedDataMut + Clone + 'a> From<i64> for GaExpr<'a, T> {
    #[inline]
    fn from(x: i64) -> Self {
        (x as f64).into()
    }
}

macro_rules! scalar_with_gaexpr_binary_ops {
    ($($t:ty),*) => {
        $(
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Add<GaExpr<'a, T>> for $t {
            type Output = GaExpr<'a, T>;
            #[inline]
            fn add(self, rhs: GaExpr<'a, T>) -> Self::Output {
                Into::<GaExpr<'a, T>>::into(self) + rhs
            }
        }
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Mul<GaExpr<'a, T>> for $t {
            type Output = GaExpr<'a, T>;
            #[inline]
            fn mul(self, rhs: GaExpr<'a, T>) -> Self::Output {
                Into::<GaExpr<T>>::into(self) * rhs
            }
        }
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Div<$t> for GaExpr<'a, T> {
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
    ($($fn_name:ident $variant:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                Self::new(move |b| {
                    let (id, gs) = self.build_or_reuse(b);
                    (AstNode::$variant(id), gs.$grade_op())
                })
            }
        )*
    }
}

impl<'a, T: 'a> GaExpr<'a, T> {
    gaexpr_unary_ops! {
        rev Reverse id "Reverse (dagger)",
        ginvol GradeInvolution id "Grade involution (main involution)",
        exp Exponential exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Logarithm log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    }

    pub fn pow(self, p: impl Into<Self>) -> Self {
        GaExpr::exp(GaExpr::log(self) * p)
    }

    pub fn sqrt(self) -> Self
    where
        T: GradedDataMut + Clone,
    {
        self.wrap(|this, this_id, this_gs| {
            if this_gs.is_just(0) {
                Wrapper::Node((
                    AstNode::ScalarUnaryOp(ScalarUnaryOp::SquareRoot, this_id),
                    this_gs.clone(),
                ))
            } else {
                Wrapper::Expr(this.pow(0.5))
            }
        })
    }

    /// Grade projection: a.g(k) = \<a\>_k
    pub fn g(self, k: i64) -> Self {
        self.gselect(move |_| GradeSet::single(k))
    }

    /// Grade projection: select specific grades.
    pub fn gselect(self, get_wanted_grades: impl Fn(&GradeSet) -> GradeSet + 'a) -> Self {
        Self::new(move |b| {
            let (id, gs) = self.build_or_reuse(b);
            (
                AstNode::GradeProjection(id),
                get_wanted_grades(&gs).intersection(gs),
            )
        })
    }

    /// Clifford conjugate. Just a shortcut for `self.rev().ginvol()`
    pub fn conj(self) -> Self {
        self.rev().ginvol()
    }

    /// Scalar product. Just a shortcut for `(self.rev() * rhs).g(0)`
    pub fn scal(self, rhs: Self) -> Self {
        (self.rev() * rhs).g(0)
    }

    /// Norm squared. Just a shortcut for `(self.clone().rev() * self).g(0)`
    pub fn norm_sq(self) -> Self {
        self.clone().scal(self)
    }

    /// Inverts only the scalar part
    pub fn sinv(self) -> Self {
        Self::new(move |b| {
            let (id, gs) = self.build_or_reuse(b);
            (AstNode::ScalarUnaryOp(ScalarUnaryOp::Inversion, id), gs)
        })
    }

    /// Versor inverse. Just a shortcut for `self.clone().rev() *
    /// self.norm_sq().sinv()`. Will default to [`Self::sinv`] when self is just a
    /// scalar
    pub fn vinv(self) -> Self {
        self.wrap(|this, _, this_gs| {
            Wrapper::Expr(if this_gs.is_just(0) {
                this.sinv()
            } else {
                this.clone().rev() * this.norm_sq().sinv()
            })
        })
    }
}
