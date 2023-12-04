//! How to evaluate a GA expression into an actual multivector

use std::collections::HashMap;

use crate::{ast::*, graded::*};
use AstNode as N;

type Cache<R> = HashMap<ExprId, R>;

impl<T: GradedData + std::fmt::Debug> SpecializedGaExpr<T> {
    /// Evaluates a [`GaExpr`]
    pub fn eval<R>(&self) -> R
    where
        R: GradedDataMut + Clone,
    {
        self.eval_with_cache(self.root(), &mut HashMap::new())
    }

    fn eval_with_cache<R>(&self, this_id: ExprId, cache: &mut Cache<R>) -> R
    where
        R: GradedDataMut + Clone,
    {
        let this = self.get_node(this_id);
        if this.is_used_several_times() {
            match cache.get(&this_id) {
                Some(r) => r.clone(),
                None => {
                    let mut res = R::init_null_mv(this.vec_space_dim(), this.grade_set());
                    self.add_to_res(this_id, cache, &mut res);
                    cache.insert(this_id, res.clone());
                    res
                }
            }
        } else {
            let mut res = R::init_null_mv(this.vec_space_dim(), this.grade_set());
            self.add_to_res(this_id, cache, &mut res);
            res
        }
    }

    fn add_to_res<R>(&self, this_id: ExprId, cache: &mut Cache<R>, res: &mut R)
    where
        R: GradedDataMut + Clone,
    {
        let this = self.get_node(this_id);
        if this.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match this.ast_node() {
            N::GradedObj(input) => {
                res.add_grades_from(input, this.grade_set());
            }
            N::Addition(left_expr, right_expr) => {
                self.add_to_res(*left_expr, cache, res);
                self.add_to_res(*right_expr, cache, res);
            }
            N::Negation(e) => {
                self.add_to_res(*e, cache, res);
                for k in this.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            N::Product(Product {
                comp_muls_cell,
                left_expr,
                right_expr,
                ..
            }) => {
                let mv_left: R = self.eval_with_cache(*left_expr, cache);
                let mv_right: R = self.eval_with_cache(*right_expr, cache);
                for mul in comp_muls_cell
                    .as_ref()
                    .expect("IndividualCompMul cell has not been set")
                {
                    let val_left = mv_left.grade_slice(mul.left_comp.grade)[mul.left_comp.index];
                    let val_right =
                        mv_right.grade_slice(mul.right_comp.grade)[mul.right_comp.index];
                    let val_result =
                        &mut res.grade_slice_mut(mul.result_comp.grade)[mul.result_comp.index];
                    *val_result += val_left * val_right * mul.coeff;
                }
            }
            N::Reverse(e) => {
                self.add_to_res(*e, cache, res);
                for k in this.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::GradeInvolution(e) => {
                self.add_to_res(*e, cache, res);
                for k in this.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::ScalarUnaryOp(op, e) => {
                self.add_to_res(*e, cache, res);
                let s = res.grade_slice_mut(0);
                s[0] = match op {
                    ScalarUnaryOp::Inversion => 1.0 / s[0],
                    ScalarUnaryOp::SquareRoot => s[0].sqrt(),
                }
            }
            N::GradeProjection(e) => self.add_to_res(*e, cache, res),
            N::Exponential(_e) => todo!(),
            N::Logarithm(_e) => todo!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{algebra::OrthoEuclidN, grade_map_mv, graded::GradeMapMV, GaExpr};

    macro_rules! expr_eq {
        ($alg:ident, $a:expr, $b:expr) => {
            let expr = $a;
            let specialized = expr.specialize(&$alg);
            assert_eq!(specialized.eval::<GradeMapMV>(), $b);
        };
    }

    type E = GaExpr<GradeMapMV>;
    const EGA3: OrthoEuclidN = OrthoEuclidN(3);
    const PGA2: [f64; 3] = [0.0, 1.0, 1.0];

    #[test]
    fn vecs_to_bivec() {
        let [e1, e2, _] = E::basis_vectors();
        expr_eq!(EGA3, e1 ^ e2, grade_map_mv!(2 => 1 0 0));
    }

    #[test]
    fn vecs_to_trivec() {
        let [e1, e2, e3] = E::basis_vectors();
        expr_eq!(EGA3, e2 ^ e1 ^ e3, grade_map_mv!(3 => -1));
    }

    #[test]
    fn vec_norm() {
        let [e0, e1, e2] = E::basis_vectors();
        expr_eq!(PGA2, (e0 - 2 * e1 + e2).norm_sq(), grade_map_mv!(0 => 5));
    }

    #[test]
    fn projection() {
        let [e1, e2, e3] = E::basis_vectors();
        let v = e1.clone() + e2.clone();
        let bv = 4 * e1 ^ e3;
        // project v onto bv:
        expr_eq!(
            EGA3,
            (v & bv.clone()) & bv.vinv(),
            grade_map_mv!(1 => 1 0 0)
        );
    }
}
