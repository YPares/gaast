# Geometric Algebra Abstract Syntax Tree

Contruct geometric algebra expressions and perform grade inference on them,
before direct evaluation or code generation, in order to limit allocations and
computations to what is needed to get the wanted result. For instance, let's
take the following expression:

$D = \langle AB + C \rangle_{2}$

where $A$, $B$ and $C$ are arbitrary multivector-valued sub-expressions or
literals. To evaluate $D$, we will only read, allocate or compute:

- the grade 2 part of $C$
- the parts of $A$ and $B$ which have a grade that will contribute to the grade
  2 part of the geometric product $AB$

To do so, GA expressions are processed in three phases:

- 1: AST construction (declare input multivectors and combine them into an
  expression with the provided operators)
- 2: AST specialization (apply the metric, perform grade inference and resolve
  individual component-to-component operations)
- 3: AST evaluation (actually read the input data and run the operations
  resolved at the previous steps)

`gaast` is still pretty experimental. It aims at applications dealing with
metric & finite-dimensional vector spaces, but which must cope with a broad
range of possible dimensionalities. For instance: data mining, machine learning
and multi-dimensional signal processing, or just general GA
teaching/visualization. It therefore does not aim at being the best fit for
fixed and low-dimension applications like computer graphics or physics, for
which an implementation tailored to a specific algebra and dimension (like 3D
Plane-based or Conformal GA) would be more efficient.

Please refer to the documentation in [`lib.rs`](src/lib.rs) for more
information.

## Implemented

_When applicable, "**MV**" means that the feature is implemented for expressions
evaluating to arbitrary multivectors, and "**VSR**" means it is only implemented
for expressions evaluating to versors._

- Linear combinations of expressions, raw multivectors and scalar literals
  (**MV**)
- Common GA products (**MV**)
- Reverse and grade involution (**MV**)
- Inverse (**VSR**)
- Arbitrary diagonal metrics
- AST construction, specialization and evaluation
- Sub-expression identification and caching: if the same expression is reused
  twice in the AST, its result is computed only once and reused. Note that for
  now this is limited to the usage of these subexpressions in products (where
  intermediary allocations are needed), because the advantage of such caching
  for operations that don't require intermediary allocations (sums, reverse,
  grade involutions...) isn't clear.
- Unit testing of the grade inference logic and of the evaluation of some basic
  expressions in EGA3

Quite a few implementation ideas are taken from the book _"Geometric Algebra for
Computer Science"_, by Leo Dorst, Daniel Fontijne and Stephen Mann.

## Roadmap and limitations

`gaast` is still very much a work in progress, so there is some more work to be
done:

### Immediate roadmap

- Versor exponentiation & logarithm
- Change the inputs used by a SpecializedGaExpr (ie. re-using a "precompiled"
  expression with different inputs that respect the same "schema"). It's doable
  right know by using a custom mutable datatype that implements the `Graded`
  trait, but a more convenient API should be built
- Gram matrix diagonalization (to work on algebras with a non-diagonal metric,
  ie. with a basis of non-orthogonal vectors)
- More extensive test suite
- Better handling of sparsity (see next section)

### Caveats

- This is Rust. Expect to use `.clone()` a lot on subexpressions that you want
  to refer to several times in your final GA expression. Note that no memory is
  actually copied, so it's mostly an API concern. `clone` also ensures the two
  subexpressions will share the same identifier, which allows caching. Aside
  from this, the API should be pretty concise, and thanks to the genericity of
  Rust operators, the need for explicit casts should be limited.
- For now, cloned subexpressions should ONLY be reused as subexpressions of the
  SAME final expression, they should NOT be used in different expressions that
  will end up being specialized and evaluated separately. This is because, for
  now, the `GaExpr` and `SpecializedGaExpr` types simply share memory (so there
  is aliasing).
- Grade sets and basis blades are represented by dynamically-sized bitvectors
  (to be agnostic of vector space dimension), which are of course much slower
  than stack-allocated bitfields like `u64`. However, this bitvector
  manipulation is limited to phases 1 & 2, so it should not impact the speed of
  the actual computations.
- Multivectors are read and written in a "one array per grade" fashion. This
  imposes a dense storage whatever the grade. Enabling a sparse storage for some
  grades would be necessary for higher-dimension spaces.

### Potential future work

- JIT compilation of `SpecializedGaExpr`s with LLVM thanks to
  https://github.com/TheDan64/inkwell (see e.g.
  https://createlang.rs/01_calculator/jit_intro.html)
- Investigate opportunities for parallelization (eg. with
  https://github.com/rayon-rs/rayon)
- Bindings to other languages (eg. Python via https://pyo3.rs)

## Why Rust?

Because of its very good trade-offs between expressiveness, performance and
safety (ie., guarantees about your code not running straight into a segfault).
Rust traits and general approach to polymorphism offer very good abstraction
powers, while incurring only small to non-existent runtime performance
penalties, which is wonderful for a tool like `gaast` which wants to make as
little assumptions as possible. Rust is also good when it comes to
interoperability, as it can produce shared libraries with a standard C ABI with
no extra dependencies. Tooling and ecosystem in general are also great.
