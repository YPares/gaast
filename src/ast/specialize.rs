use super::base_types::*;
use crate::{
    algebra::{iter_basis_blades_of_grade, MetricAlgebra},
    grade_set::*,
    Expr,
};
use AstNode as N;

/// An [`Expr`][crate::Expr] which is ready for compilation/evaluation
#[derive(Debug)]
pub struct SpecializedAst<T> {
    ast: ReifiedAst<T>,
}

impl<T> SpecializedAst<T> {
    /// The identifier of the root node of the expression
    pub fn root_id(&self) -> NodeId {
        self.ast.root_id
    }

    /// Get a node from its [`NodeId`]
    pub fn get_node(&self, node_id: NodeId) -> &GradedNode<T> {
        &self.ast.arena[&node_id]
    }
}

impl<'a, T: Clone> Expr<'a, T> {
    /// Tells which algebra this [`Expr`] is using, and recursively propagates
    /// wanted grades downwards so as to evaluate for each sub-expression only
    /// the grades that are necessary to compute the whole [`Expr`].
    ///
    /// The given [`MetricAlgebra`] must make sense with respect to the input
    /// values contained in the [`Expr`], in terms of possible grades
    /// contained in those input values, and of number of components for each
    /// grade
    pub fn specialize(self, alg: &impl MetricAlgebra) -> SpecializedAst<T> {
        // First, we reify the expression, to go from a closure-based repr to an
        // actual tree stored in an arena
        let mut e = self.reify(alg);
        let root_gs = e.arena[&e.root_id].maximal_grade_set.clone();
        // Then the rest of the process is split in two because sub-expressions
        // may be used at several places in the AST. First we collect all the
        // requirements for expressions throughout the whole tree, and add them
        // to the minimal grade sets of the corresponding nodes:
        rec_update_minimal_grade_sets(&mut e.arena, e.root_id, root_gs);
        // Then for each node we apply and store at each node what is needed
        // from the algebra:
        rec_apply_algebra(&mut e.arena, e.root_id, alg);
        SpecializedAst { ast: e }
    }
}

fn rec_update_minimal_grade_sets<T: Clone>(
    arena: &mut NodeArena<T>,
    this_id: NodeId,
    wanted: GradeSet,
) {
    // Here, we know the resulting grade of the operation represented by the top
    // node of the AST, and we recursively "undo" the AST it to find the grades
    // wanted in each subexpression (the minimal set of grades they should
    // evaluate to given the final result we want).
    //
    // If a node is refered to in several parts of the AST, all the grade
    // requirements will be combined. So that node can be evaluated just once,
    // and it will return everything needed throughout the AST
    arena.get_mut(&this_id).unwrap().minimal_grade_set += wanted.clone();

    match arena[&this_id].ast_node.clone() {
        N::GradedObj(_) => {}
        N::GradeProjection(e)
        | N::Negation(e)
        | N::Reverse(e)
        | N::GradeInvolution(e)
        | N::ScalarUnaryOp(_, e) => rec_update_minimal_grade_sets(arena, e, wanted),
        N::Addition(left_expr, right_expr) => {
            // <A + B>_k = <A>_k + <B>_k
            rec_update_minimal_grade_sets(arena, left_expr, wanted.clone());
            rec_update_minimal_grade_sets(arena, right_expr, wanted);
        }
        N::Product(p) => {
            // Find in left_expr and right_expr which grades, once
            // multiplied, will affect the grades in `wanted`
            let (left_wanted_gs, right_wanted_gs) = wanted.parts_contributing_to_product(
                p.grades_to_produce.get_fn(),
                &arena[&p.left_expr].maximal_grade_set,
                &arena[&p.right_expr].maximal_grade_set,
            );
            rec_update_minimal_grade_sets(arena, p.left_expr, left_wanted_gs);
            rec_update_minimal_grade_sets(arena, p.right_expr, right_wanted_gs);
        }
        N::Exponential(e) => rec_update_minimal_grade_sets(arena, e, wanted.log()),
        N::Logarithm(e) => rec_update_minimal_grade_sets(arena, e, wanted.exp()),
    };
}

fn rec_apply_algebra<T: Clone>(
    arena: &mut NodeArena<T>,
    this_id: NodeId,
    alg: &impl MetricAlgebra,
) {
    if arena[&this_id].is_ready {
        // Algebra has already been applied for this node (and its subnodes),
        // because it is referred to in at least one other part of the AST
        assert!(
            arena[&this_id].is_used_several_times(),
            "Algebra was already applied to a node that is referred to only once"
        );
        return;
    }
    {
        let this = arena.get_mut(&this_id).unwrap();
        this.is_ready = true;
        assert!(
            this.maximal_grade_set
                .includes(this.minimal_grade_set.clone()),
            "Inferred minimal grade set contains grades not available in maximal grade set"
        );
    }
    match arena[&this_id].ast_node.clone() {
        N::GradedObj(_) => {}
        N::Negation(e)
        | N::GradeProjection(e)
        | N::Reverse(e)
        | N::GradeInvolution(e)
        | N::ScalarUnaryOp(_, e) => {
            rec_apply_algebra(arena, e, alg);
        }
        N::Addition(left_expr, right_expr) => {
            rec_apply_algebra(arena, left_expr, alg);
            rec_apply_algebra(arena, right_expr, alg);
        }
        N::Product(p) => {
            rec_apply_algebra(arena, p.left_expr, alg);
            rec_apply_algebra(arena, p.right_expr, alg);
            let gs_left = arena[&p.left_expr].minimal_grade_set.clone();
            let gs_right = arena[&p.right_expr].minimal_grade_set.clone();

            let comp_muls = arena[&this_id]
                .minimal_grade_set
                .iter_contribs_to_product(p.grades_to_produce.get_fn(), &gs_left, &gs_right)
                .flat_map(|individual_grades_and_contribs| {
                    iter_comp_muls_for_kvectors_prod(alg, individual_grades_and_contribs)
                })
                .collect();
            let this = arena.get_mut(&this_id).unwrap();
            match &mut this.ast_node {
                N::Product(p) => {
                    p.individual_comp_muls = comp_muls;
                }
                _ => panic!("Should be a Product node"),
            };
        }
        N::Exponential(e) => {
            rec_apply_algebra(arena, e, alg);
        }
        N::Logarithm(e) => {
            rec_apply_algebra(arena, e, alg);
        }
    }
}

fn iter_comp_muls_for_kvectors_prod(
    alg: &impl MetricAlgebra,
    (k_left, k_right, contribs): (usize, usize, GradeSet),
) -> impl Iterator<Item = IndividualCompMul> + '_ {
    iter_basis_blades_of_grade(alg, k_left).flat_map(move |bb_left| {
        let contribs = contribs.clone();
        iter_basis_blades_of_grade(alg, k_right).filter_map(move |bb_right| {
            let (bb_res, coeff) = alg.ortho_basis_blades_gp(&bb_left, &bb_right);
            let result_comp = alg.basis_blade_to_component(&bb_res);
            if contribs.contains(result_comp.grade) {
                Some(IndividualCompMul {
                    left_comp: alg.basis_blade_to_component(&bb_left),
                    right_comp: alg.basis_blade_to_component(&bb_right),
                    result_comp,
                    coeff,
                })
            } else {
                None
            }
        })
    })
}
