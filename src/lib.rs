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

pub mod ast;
pub mod eval;
pub mod grade_set;
pub mod graded;

pub use ast::GAExpr;
pub use graded::{DynSizedMV, Graded};

macro_rules! pub_all {
    (
        $(#[$meta:meta])*
        struct $name:ident {
            $(
                $field:ident : $typ:ty
            ),*
        }
    ) => {
        $(#[$meta])*
        pub struct $name {
            $(
                pub $field : $typ
            ),*
        }
    };
}
macro_rules! pub_crate_all {
    (
        $(#[$meta:meta])*
        struct $name:ident {
            $(
                $field:ident : $typ:ty
            ),*
        }
    ) => {
        $(#[$meta])*
        pub(crate) struct $name {
            $(
                pub(crate) $field : $typ
            ),*
        }
    };
}
pub(crate) use pub_all;
pub(crate) use pub_crate_all;
