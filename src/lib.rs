//! # Geometric Algebra Abstract Syntax Tree
//!
//! Allows to contruct geometric algebra expressions and perform grade inference
//! on them, in order to limit computations and RAM usage to only what is needed
//! to get the wanted result. The main type to build and combine these
//! expressions is [`GAExpr`], which exposes methods for the geometric algebra
//! primitives.
//!
//! [`GAExpr`] is agnostic over which type is actually used to store input
//! multivectors. This is achived via the [`Graded`] trait.

pub mod algebra;
pub mod ast;
pub mod eval;
pub mod grade_set;
pub mod graded;

pub use ast::GAExpr;
pub use grade_set::{Grade, GradeSet};
pub use graded::{DynSizedMV, Graded};
