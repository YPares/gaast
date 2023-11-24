/*!
# Geometric Algebra Abstract Syntax Tree

Contruct geometric algebra expressions and perform grade inference on them,
in order to limit computations and RAM usage to only what is needed to get
the wanted result. The main type to build and combine these expressions is
[`GAExpr`], which exposes methods for the geometric algebra primitives.

[`GAExpr`] is agnostic over which types are actually used to store input and
result multivectors. This is achived via the [`Graded`] trait and its
subtraits. It is also agnostic over which actual vector space is used (both
in terms of dimension and metric), via the [`algebra::MetricAlgebra`] trait.
The only hard requirement is that this vector space must be over the field
of reals ([`f64`] here).

Computing GA expressions with `gaast` is done fully at runtime, in 3 phases:

- 1: **AST construction**. This step does _upwards grade inference_: the
  grades contained in input data are queried (_without_ actually reading
  that data), and then we go up the AST, inferring the grades of the result
  of each AST node (ie. each intermediate operation), and therefore of the
  full expression. We also issue errors and warnings during that phase: for
  instance if a grade projection operator is guaranteed to always return 0
  (because it projects to a grade that doesn't exist in the expression it is
  applied to).
- 2: **Grade minimisation**. This step does _downwards grade inference_.
  It's the same process, but in reverse: knowing the grade(s) that our
  entire expression evaluates to (notably if this expression uses some grade
  projections), we propagate these wanted grades downwards so as to update
  each AST node with only the grades that will need to be read or computed
  in the end.
- 3: **Evaluation**. Given the annotated AST resulting from the two previous
  phases, allocates, reads and computes only what is needed to get the final
  multivector to which the full GAExpr evaluates to.

The main point of this approach is that phases 1 & 2 can be cached. They
don't need to access the actual data: they just need to know the grades and
the metric in use. So as long as the final data remains compliant with that,
the same AST can be "pre-compiled" and reused several times. In the future,
codegen could be done to compile an AST after phase 2 to LLVM or to a GPU
kernel or shader. Enough of the internal representation of the AST is made
public so that external crates can implement this.

The second point of this approach is to limit the amount of primitives that
must be implemented directly. Basically, we can just implement the geometric
product in its most generic form, and let `gaast` optimize away the parts we
do not need each time we use it. For instance, there is no need to provide a
specialized implementation for the outer product and left & right
contractions: all of these correspond to a geometric product followed by a
grade projection, and the grade inference will make it so that only the part
of the geometric product that is actually needed is computed.

`gaast` also provides exponential and logarithm of multivectors as
primitives, though with some restrictions. However, these restrictions are
checked during phase 1 to permit early failures.
*/

pub mod algebra;
pub mod ast;
pub mod grade_set;
pub mod graded;

pub use ast::{mv, GaExpr};
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
