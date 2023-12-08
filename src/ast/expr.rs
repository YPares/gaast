use super::{base_types::*, GradedNode};
use crate::{algebra::MetricAlgebra, grade_set::*, graded::GradedDataMut, Graded};
use std::{collections::HashMap, rc::Rc};

/// Just a type grouping references to the metric and to the node arena together
struct Builder<'a, T> {
    algebra: &'a dyn MetricAlgebra,
    arena: &'a mut NodeArena<T>,
}

impl<'a, T> Builder<'a, T> {
    /// Add a node to the arena. Called by the rec_reify closures of the Exprs
    fn add_node(&mut self, node_id: NodeId, (ast_node, node_gs): (AstNode<NodeId, T>, GradeSet)) {
        self.arena.insert(
            node_id,
            GradedNode {
                maximal_grade_set: node_gs.intersection(self.algebra.full_grade_set()),
                minimal_grade_set: GradeSet::empty(),
                vec_space_dim: self.algebra.vec_space_dim(),
                ast_node,
                num_uses: 1,
                is_ready: false,
            },
        );
    }
}

/// The base type to use to build GA expressions
pub struct Expr<'a, T> {
    /// A closure that adds to the arena the nodes needed to process this
    /// expression.
    ///
    /// We use a closure-based representation for two reasons:
    /// - follow-up grade inference operations mutate the AST, and a flat
    ///   storage (a hashmap of nodes) is more convenient to do so. But this
    ///   storage should be abstracted from the user. It's the job of the
    ///   closure stored here to actually put its nodes into that type of
    ///   storage
    /// - The AST cannot be reified until the algebra is known, because some
    ///   nodes may depend on which algebra is used. By using closures, we can
    ///   defer actual construction to when the algebra is known, and pass it
    ///   along as we add the actual nodes to the arena storage
    run: Rc<dyn Fn(NodeId, &mut Builder<'_, T>) -> () + 'a>,
}

/// [`Expr`] uses [`Rc`] internally so clones are actually cheap
impl<T> Clone for Expr<'_, T> {
    fn clone(&self) -> Self {
        Self {
            run: Rc::clone(&self.run),
        }
    }
}

enum Wrapper<'a, T> {
    Node((AstNode<NodeId, T>, GradeSet)),
    Expr(Expr<'a, T>),
}

impl<'a, T: 'a> Expr<'a, T> {
    /// Given the algebra, call the Expr closures and get a mutable AST
    pub(super) fn reify(self, alg: &impl MetricAlgebra) -> ReifiedAst<T> {
        let mut arena = HashMap::new();
        let (root_id, _) = self.reify_or_reuse(&mut Builder {
            algebra: alg,
            arena: &mut arena,
        });
        ReifiedAst { arena, root_id }
    }

    /// Adds a node to the arena (or does nothing if it already exists),
    /// increments its number of uses, and returns its identifier and GradeSet
    fn reify_or_reuse(&self, b: &mut Builder<'_, T>) -> (NodeId, GradeSet) {
        let id = NodeId {
            ptr: Rc::as_ptr(&self.run).cast(),
        };
        match b.arena.get_mut(&id) {
            None => (self.run)(id, b),
            Some(node) => {
                node.num_uses += 1;
            }
        }
        (id, b.arena[&id].maximal_grade_set.clone())
    }

    fn new(f: impl Fn(&mut Builder<'_, T>) -> (AstNode<NodeId, T>, GradeSet) + 'a) -> Self {
        Self {
            run: Rc::new(move |this_id, b| {
                let node_and_gs = f(b);
                b.add_node(this_id, node_and_gs)
            }),
        }
    }

    /// Enables to wrap a subexpression either in a new node, or by reusing
    /// pre-existing functions returning Exprs
    fn wrap(self, f: impl Fn(Self, NodeId, GradeSet) -> Wrapper<'a, T> + 'a) -> Self {
        Self {
            run: Rc::new(move |wrapper_id, b| {
                let (self_id, self_gs) = self.reify_or_reuse(b);
                match f(self.clone(), self_id, self_gs) {
                    Wrapper::Node(wrapper_node_and_gs) => {
                        b.add_node(wrapper_id, wrapper_node_and_gs)
                    }
                    Wrapper::Expr(wrapper_expr) => {
                        (wrapper_expr.run)(wrapper_id, b);
                        // We correct for the extra use we counted when calling
                        // self.build_or_reuse above (as wrapper_expr.run will
                        // count it):
                        b.arena.get_mut(&self_id).unwrap().num_uses -= 1;
                    }
                }
            }),
        }
    }

    /// Compute some product of two multivector operands, from the geometric
    /// products of their individual k-vector parts. A grade projection will be
    /// applied to each one of these individual products. The exact grade
    /// projection to perform is given by the `grades_to_produce` argument,
    /// which is given the grades of the 2 k-vector parts, and outputs the set
    /// of grades to project out of their product
    pub fn product(
        self,
        rhs: Self,
        grades_to_produce: impl Fn((i64, i64)) -> GradeSet + 'static,
    ) -> Self {
        let grades_to_produce = KVecsProductGradeSelection(Rc::new(grades_to_produce));
        Self::new(move |b| {
            let (left_id, left_gs) = self.reify_or_reuse(b);
            let (right_id, right_gs) = rhs.reify_or_reuse(b);
            (
                AstNode::Product(Product {
                    individual_comp_muls: vec![],
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

    /// Returns the expressions for the basis vectors in a vector space of dim
    /// `D`
    pub fn basis_vectors<const D: usize>() -> [Self; D]
    where
        T: GradedDataMut + Clone,
    {
        array_init::array_init(|i| {
            let mut v = T::init_null_mv(D, &GradeSet::single(1));
            v.grade_slice_mut(1)[i] = 1.0;
            mv(v)
        })
    }
}

/// Create an [`Expr`] that just evaluates to some pre-existing [`Graded`] value
/// (ie. a multivector)
pub fn mv<'a, T: Graded + Clone + 'a>(x: T) -> Expr<'a, T> {
    Expr::new(move |_| (AstNode::GradedObj(x.clone()), x.grade_set().clone()))
}

macro_rules! expr_products {
    ($($doc:literal $trait:ident $method:ident ($fn:expr)),*) => {
        $(
        #[doc=$doc]
        impl<'a, T: 'a, E: Into<Expr<'a, T>>> std::ops::$trait<E> for Expr<'a, T> {
            type Output = Self;
            #[doc=$doc]
            fn $method(self, rhs: E) -> Self::Output {
                self.product(rhs.into(), $fn)
            }
        }
        )*
    };
}
expr_products! {
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

/// Multivector addition
impl<'a, T: 'a, E: Into<Expr<'a, T>>> std::ops::Add<E> for Expr<'a, T> {
    type Output = Self;
    fn add(self, rhs: E) -> Self::Output {
        let rhs = rhs.into();
        Self::new(move |b| {
            let (left_id, left_gs) = self.reify_or_reuse(b);
            let (right_id, right_gs) = rhs.reify_or_reuse(b);
            (AstNode::Addition(left_id, right_id), left_gs + right_gs)
        })
    }
}

/// Multivector negation
impl<'a, T: 'a> std::ops::Neg for Expr<'a, T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(move |b| {
            let (id, gs) = self.reify_or_reuse(b);
            (AstNode::Negation(id), gs)
        })
    }
}

/// Multivector subtraction
impl<'a, T: 'a, E: Into<Expr<'a, T>>> std::ops::Sub<E> for Expr<'a, T> {
    type Output = Self;
    fn sub(self, rhs: E) -> Self::Output {
        self + -rhs.into()
    }
}

impl<'a, T: GradedDataMut + Clone + 'a> From<f64> for Expr<'a, T> {
    fn from(x: f64) -> Self {
        if x == 0.0 {
            return mv(T::init_null_mv(0, &GradeSet::empty()));
        }
        let mut scal_mv = T::init_null_mv(0, &GradeSet::single(0));
        scal_mv.grade_slice_mut(0)[0] = x;
        mv(scal_mv)
    }
}
impl<'a, T: GradedDataMut + Clone + 'a> From<i64> for Expr<'a, T> {
    #[inline]
    fn from(x: i64) -> Self {
        (x as f64).into()
    }
}

macro_rules! scalar_with_expr_binary_ops {
    ($($t:ty),*) => {
        $(
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Add<Expr<'a, T>> for $t {
            type Output = Expr<'a, T>;
            #[inline]
            fn add(self, rhs: Expr<'a, T>) -> Self::Output {
                Into::<Expr<'a, T>>::into(self) + rhs
            }
        }
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Mul<Expr<'a, T>> for $t {
            type Output = Expr<'a, T>;
            #[inline]
            fn mul(self, rhs: Expr<'a, T>) -> Self::Output {
                Into::<Expr<T>>::into(self) * rhs
            }
        }
        impl<'a, T: GradedDataMut + Clone + 'a> std::ops::Div<$t> for Expr<'a, T> {
            type Output = Self;
            fn div(self, rhs: $t) -> Self {
                self * (1.0/(rhs as f64))
            }
        }
        )*
    };
}
scalar_with_expr_binary_ops!(f64, i64);

macro_rules! expr_unary_ops {
    ($($fn_name:ident $variant:ident $grade_op:ident $doc:literal),*) => {
        $(
            #[doc=$doc]
            pub fn $fn_name(self) -> Self {
                Self::new(move |b| {
                    let (id, gs) = self.reify_or_reuse(b);
                    (AstNode::$variant(id), gs.$grade_op())
                })
            }
        )*
    }
}

impl<'a, T: 'a> Expr<'a, T> {
    expr_unary_ops! {
        rev Reverse id "Reverse (dagger)",
        ginvol GradeInvolution id "Grade involution (main involution)",
        exp Exponential exp "Exponential. IMPORTANT: Is defined only for **single**-graded k-vectors",
        log Logarithm log "Natural logarithm. IMPORTANT: Is defined only for multivectors of the form \\<A\\>_0 + \\<A\\>_k"
    }

    /// Raises to some power. Shortcut for `exp(log(self) * p)`. Therefore,
    /// please refer to [`Self::log`] and [`Self::exp`] for limitations
    pub fn pow(self, p: impl Into<Self>) -> Self {
        Expr::exp(Expr::log(self) * p)
    }

    /// Square root. See [`Self::pow`]
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
            let (id, gs) = self.reify_or_reuse(b);
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
            let (id, gs) = self.reify_or_reuse(b);
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
