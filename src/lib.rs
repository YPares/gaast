//! # Geometric Algebra Abstract Syntax Tree
//!
//! Contruct geometric algebra expressions and perform grade inference on them,
//! in order to limit allocations and computations to what is needed to get the
//! wanted result. The main type to build and combine these expressions is the
//! [`Expr`] type, which exposes methods for the geometric algebra primitives.
//!
//! [`Expr`] is agnostic over which types are actually used to store input and
//! result multivectors. This is achived via the [`Graded`] trait and its
//! subtraits. It is also agnostic over which actual vector space is used (both
//! in terms of dimension and metric), via the [`algebra::MetricAlgebra`] trait.
//! The only hard requirement is that this vector space must be over the field
//! of reals ([`f64`] here).
//!
//! Computing GA expressions with `gaast` is done fully at runtime, in 4 phases:
//!
//! - 1: **Expression construction**. Combine raw multivectors into a geometric
//!   algebra expression using the operators and methods provided by [`Expr`].
//!   The actual algebra is not known at this point, therefore internally, at
//!   this step, the expressions are represented by closures, expecting to
//!   receive an object representing an algebra to create the actual AST.
//! - 2: **AST reification**. We feed the algebra and evaluate the
//!   aforementioned closures. This step also does _upwards grade inference_:
//!   the algebra and the grades of the input data are known, so we can, as we
//!   construct the AST, infer for each subexpression the grades contained in
//!   the multivector to which it will evaluate to.
//! - 3: **AST specialization**. This step does _downwards grade inference_.
//!   It's the same process than during phase 2, but in reverse: knowing the
//!   grade(s) that our entire expression evaluates to (notably if this
//!   expression uses some grade projections), we propagate the wanted grades
//!   downwards so as to update each AST node with only the grades that will
//!   need to be read or computed in the end. Individual component-to-component
//!   multiplications to perform for each multivector product are fully resolved
//!   at this phase, and therefore the metric is no longer needed after it.
//! - 4: **AST evaluation**. Given the annotated AST resulting from the previous
//!   phases, allocates, reads and computes only what is needed to get the final
//!   multivector to which the full AST evaluates to.
//!
//! To mark these various stages, different types are used to reflect the result
//! of each phase. Phases 2 & 3 are done by the [`Expr::specialize`] function,
//! which transforms the [`Expr`] into a [`SpecializedAst`], so trying to
//! evaluate an expression that hasn't been specialized yet will result in a
//! compile-time type error.
//!
//! The main point of this approach is that phases 1, 2 & 3 can be cached. They
//! don't need to access the actual data: they just need to know the grades and
//! the metric in use. So as long as the final data remains compliant with that,
//! the same AST can be "pre-compiled" and reused several times. In the future,
//! codegen could be done to compile an AST after phase 3 to LLVM or to a GPU
//! kernel or shader. Enough of the internal representation of the AST is made
//! public so that external crates can implement this.
//!
//! The second point of this approach is to limit the amount of primitives that
//! must be special-cased. Basically, we can just implement a general notion of
//! product (a combination of 2 multivectors that results in a specific set of
//! component-to-component multiplications, the most complete case being the
//! geometric product for which every possible component-to-component
//! multiplication is performed). [`Expr`] provides a generic [`Expr::product`]
//! method to which we just give a function that will tell, for a given pair of
//! k-vectors in the input, which grades of their geometric product are of
//! interest to us. Operators `*` (geometric), `^` (outer), `&` (inner), `<<`
//! (left contraction), and `>>` (right contraction) are just provided as
//! shortcuts for the most common products. But they don't need any special
//! implementation, because we just specify each time to `gaast` what we want,
//! and then let it optimize away the parts we do not need.
//!
//! `gaast` will also provide exponential and logarithm of multivectors as
//! primitives, though with some restrictions. However, these restrictions are
//! checked before evaluation actually starts, to permit early failures.

pub mod algebra;
pub mod ast;
pub mod grade_set;
pub mod graded;

pub use ast::{mv, Expr, SpecializedAst};
pub use grade_set::{Grade, GradeSet};
pub use graded::Graded;

#[cfg(feature = "eval")]
mod eval;

#[cfg(test)]
pub(crate) mod test_macros {
    macro_rules! simple_eqs {
        {$($test_name:ident : $a:expr => $b:expr),+} => {
          mod simple_eqs {
            use super::*;
            $(
                #[test]
                fn $test_name() {
                    assert_eq!($a, $b);
                }
            )+
          }
        }
    }
    pub(crate) use simple_eqs;
}
