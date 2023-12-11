//! # Geometric Algebra Abstract Syntax Tree
//!
//! Construct geometric algebra expressions and perform grade inference on them,
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
//! - 1: **Expression construction**. We combine raw multivectors into a
//!   geometric algebra expression, using the operators and methods provided by
//!   [`Expr`]. The actual algebra and metric are not known yet. To account for
//!   this, expressions at this point are represented by closures ("lambdas"),
//!   that defer the actual AST construction to when the algebra will be given.
//! - 2: **AST reification**. We feed the algebra (metric) and call the
//!   aforementioned closures. This step does _upwards grade inference_ at the
//!   same time: since the algebra and the grades of the input data are known,
//!   we can infer for each subexpression the grades that will be contained in
//!   the multivector it should evaluate to.
//! - 3: **AST specialization**. This step does _downwards grade inference_. It
//!   is very similar to the process just above, only in reverse. Knowing the
//!   grade(s) that our entire expression should evaluate to (notably if this
//!   expression uses some grade projections), we propagate the wanted grades
//!   downwards. We update each AST node, annotating it with only the grades
//!   that will need to be read or computed in the end. We also resolve the
//!   actual individual component-to-component multiplications to perform for
//!   each multivector product along the way. Therefore the algebra is no longer
//!   needed after this phase.
//! - 4: **AST evaluation**. Given the annotated AST resulting from the above,
//!   we now read, allocate and compute only what is needed to get the final
//!   multivector to which the full AST evaluates to.
//!
//! Different types are used to represent the result of each phase. Phases 2 & 3
//! are both done by the [`Expr::specialize`] function, which transforms the
//! [`Expr`] into a [`SpecializedAst`]. Thus trying to evaluate an expression
//! that has not been specialized yet will result in a compile-time error.
//!
//! The main point of this approach is that phases 1, 2 & 3 can be done just
//! once. They don't need to access the actual data: they just need to know the
//! grades and the metric in use. So as long as the final data remains compliant
//! with that, the same AST can be "pre-compiled" and reused several times. In
//! the future, codegen could be done to compile a specialized AST to LLVM or to
//! a GPU kernel or shader. Enough of the internal representation of the AST is
//! made public so that external crates can implement this.
//!
//! The second point of this approach is to limit the amount of primitives that
//! must be special-cased. Basically, we can just implement a general notion of
//! product: a combination of 2 multivectors that results in a specific set of
//! component-to-component multiplications, the most complete case being the
//! geometric product for which every possible component-to-component
//! multiplication is performed. [`Expr`] provides a generic [`Expr::product`]
//! method to which we give simple a function that will tell, for a given pair
//! of k-vectors in the input, which grades of their geometric product are of
//! interest to us. Then, operators like `*` (geometric product), `^` (outer
//! product), `&` (inner product), `<<` (left contraction), and `>>` (right
//! contraction) can simply be one-liners based on that common `product` method.
//! `gaast` provides these operators because they are very common in
//! applications, but they could just as easily be implemented in client code.
//! They don't need any special implementation because we can just specify for
//! each one what we want out of it, and let `gaast` optimize away the parts we
//! do not need.
//!
//! This approach is also what lets `gaast` be agnostic of the actual
//! representation of both input and result multivectors. For instance, if your
//! application already has a preferred vector representation, then you don't
//! need to convert back and forth. The input multivectors that the AST refers
//! to can be either owned data (ie. be stored in the AST directly) or
//! references to data stored elsewhere. Also, input multivectors can be of
//! different types (e.g. so they can be tailored to which grades they contain)
//! by using something like `Box<dyn GradedData>` or some custom enum, and still
//! interoperate in the same AST. They can even be external resources read in a
//! lazy fashion.
//! 
//! `gaast` will also provide exponential and logarithm of multivectors as
//! primitives, though with some restrictions. However, the nature of `gaast`
//! enables it to check for these conditions before evaluation actually starts,
//! which enables early failures.

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
