use super::{base_types::*, gaexpr::*, specialized_gaexpr::*};
use crate::{
    algebra::{iter_basis_blades_of_grade, MetricAlgebra},
    grade_set::*,
    Graded,
};
use std::cell::Ref;
use AstNode as N;

impl<T> GaExpr<T> {
    /// Tells which algebra this [`GaExpr`] is using, and recursively propagates
    /// wanted grades downwards so as to evaluate for each sub-expression only
    /// the grades that are necessary to compute the whole [`GaExpr`].
    ///
    /// The given [`MetricAlgebra`] must make sense with respect to the input
    /// values contained in the [`GaExpr`], in terms of possible grades
    /// contained in those input values, and of number of components for each
    /// grade
    pub fn specialize(self, alg: &impl MetricAlgebra) -> SpecializedGaExpr<T> {
        // The process is split in two because sub-expressions may be used at
        // several places in the AST. First we collect all the requirements for
        // expressions throughout the whole tree, and add them to the minimal
        // grade sets of the corresponding nodes:
        self.rec_update_minimal_grade_sets(&self.grade_set());
        // Then for each node we apply and store at each node what is needed
        // from the algebra:
        self.rec_apply_algebra(alg);
        // GaExpr<T> and SpecializedGaExpr<T> have the exact same memory
        // representation, therefore this is safe:
        unsafe { std::mem::transmute(self) }
    }

    fn rec_update_minimal_grade_sets(&self, wanted: &GradeSet) {
        // Here, we know the resulting grade of the operation represented by the
        // top node of the AST, and we recursively "undo" the AST it to find the
        // grades wanted in each subexpression (the minimal set of grades they
        // should evaluate to given the final result we want)
        self.rc
            .minimal_grade_set
            .replace_with(|gs| gs.clone() + wanted.clone());
        // If a node is refered to in several parts of the AST, all the
        // grade requirements will be combined. So that node can be
        // evaluated just once, and it will return everything needed
        // throughout the AST
        match &self.rc.ast_node {
            N::GradedObj(_) => {}
            N::GradeProjection(e)
            | N::Negation(e)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarUnaryOp(_, e) => {
                e.rec_update_minimal_grade_sets(wanted);
            }
            N::Addition(left_expr, right_expr) => {
                // <A + B>_k = <A>_k + <B>_k
                left_expr.rec_update_minimal_grade_sets(wanted);
                right_expr.rec_update_minimal_grade_sets(wanted);
            }
            N::Product(Product {
                left_expr,
                right_expr,
                grades_to_produce,
                ..
            }) => {
                // Find in left_expr and right_expr which grades, once
                // multiplied, will affect the grades in `wanted`
                let (left_wanted_gs, right_wanted_gs) = wanted.parts_contributing_to_product(
                    &grades_to_produce.0,
                    &left_expr.grade_set(),
                    &right_expr.grade_set(),
                );
                left_expr.rec_update_minimal_grade_sets(&left_wanted_gs);
                right_expr.rec_update_minimal_grade_sets(&right_wanted_gs);
            }
            N::Exponential(e) => e.rec_update_minimal_grade_sets(&wanted.clone().log()),
            N::Logarithm(e) => e.rec_update_minimal_grade_sets(&wanted.clone().exp()),
        }
    }

    fn minimal_grade_set(&self) -> Ref<'_, GradeSet> {
        self.rc.minimal_grade_set.borrow()
    }

    fn restrict_minimal_grade_set(&self, constraint: GradeSet) {
        self.rc
            .minimal_grade_set
            .replace_with(|orig| orig.clone().intersection(constraint));
    }

    fn rec_apply_algebra(&self, alg: &impl MetricAlgebra) {
        if let Err(_) = self.rc.vec_space_dim.set(alg.vec_space_dim()) {
            // Grade hints and algebra have already been applied for this node
            // (and its subnodes), because it is referred to in at least one
            // other part of the AST
            return;
        }
        self.restrict_minimal_grade_set(alg.full_grade_set());
        assert!(
            self.rc
                .maximal_grade_set
                .includes(self.minimal_grade_set().clone()),
            "Inferred minimal grade set contains grades not available in maximal grade set"
        );
        match &self.rc.ast_node {
            N::GradedObj(_) => {}
            N::Negation(e)
            | N::GradeProjection(e)
            | N::Reverse(e)
            | N::GradeInvolution(e)
            | N::ScalarUnaryOp(_, e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone());
            }
            N::Addition(left_expr, right_expr) => {
                left_expr.rec_apply_algebra(alg);
                right_expr.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(
                    left_expr.minimal_grade_set().clone() + right_expr.minimal_grade_set().clone(),
                );
            }
            N::Product(Product {
                comp_muls_cell,
                left_expr,
                right_expr,
                grades_to_produce,
            }) => {
                left_expr.rec_apply_algebra(alg);
                right_expr.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(
                    iter_grade_sets_cp(
                        &left_expr.minimal_grade_set(),
                        &right_expr.minimal_grade_set(),
                    )
                    .map(&grades_to_produce.0)
                    .collect(),
                );
                // Now that the grades at play for this product are fully
                // resolved, we can construct the set of component-to-component
                // multiplications that will be needed to perform it:
                comp_muls_cell
                    .set(
                        self.minimal_grade_set()
                            .iter_contribs_to_product(
                                &grades_to_produce.0,
                                &left_expr.minimal_grade_set(),
                                &right_expr.minimal_grade_set(),
                            )
                            .flat_map(|individual_grades_and_contribs| {
                                iter_comp_muls_for_kvectors_prod(
                                    alg,
                                    individual_grades_and_contribs,
                                )
                            })
                            .collect(),
                    )
                    .expect("IndividualCompMul cell has already been set");
            }
            N::Exponential(e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone().exp())
            }
            N::Logarithm(e) => {
                e.rec_apply_algebra(alg);
                self.restrict_minimal_grade_set(e.minimal_grade_set().clone().log())
            }
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
