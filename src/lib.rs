//! GAAST
pub mod ast;
pub mod eval;
pub mod grade_set;
pub mod graded;

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
