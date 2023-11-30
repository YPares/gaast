//! How to evaluate a GA expression into an actual multivector

use std::collections::HashMap;

use super::{ast::*, graded::*};
use AstNode as N;

type Cache<R> = HashMap<ExprId, R>;

impl<T: GradedData> ReadyGaExpr<T> {
    /// Evaluates a [`GaExpr`]
    pub fn eval<R>(&self) -> R
    where
        R: GradedDataMut + Clone,
    {
        self.eval_with_cache(&mut HashMap::new())
    }

    fn eval_with_cache<R>(&self, cache: &mut Cache<R>) -> R
    where
        R: GradedDataMut + Clone,
    {
        if self.is_reused() {
            match cache.get(&self.identify()) {
                Some(r) => r.clone(),
                None => {
                    let mut res = R::init_null_mv(self.vec_space_dim(), &self.grade_set());
                    self.add_to_res(cache, &mut res);
                    cache.insert(self.identify(), res.clone());
                    res
                }
            }
        } else {
            let mut res = R::init_null_mv(self.vec_space_dim(), &self.grade_set());
            self.add_to_res(cache, &mut res);
            res
        }
    }

    fn add_to_res<R>(&self, cache: &mut Cache<R>, res: &mut R)
    where
        R: GradedDataMut + Clone,
    {
        if self.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match self.ast_node() {
            N::GradedObj(input) => {
                res.add_grades_from(input, &self.grade_set());
            }
            N::Addition(e_left, e_right) => {
                e_left.add_to_res(cache, res);
                e_right.add_to_res(cache, res);
            }
            N::Negation(e) => {
                e.add_to_res(cache, res);
                for k in self.grade_set().iter() {
                    res.negate_grade(k);
                }
            }
            N::GeometricProduct(individual_muls_cell, e_left, e_right) => {
                let mv_left: R = e_left.eval_with_cache(cache);
                let mv_right: R = e_right.eval_with_cache(cache);
                for mul in individual_muls_cell
                    .get()
                    .expect("IndividualCoordMul cell has not been set")
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
                e.add_to_res(cache, res);
                for k in self.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::GradeInvolution(e) => {
                e.add_to_res(cache, res);
                for k in self.grade_set().iter() {
                    if k % 2 == 1 {
                        res.negate_grade(k);
                    }
                }
            }
            N::ScalarUnaryOp(op, e) => {
                e.add_to_res(cache, res);
                let s = res.grade_slice_mut(0);
                s[0] = match op {
                    ScalarUnaryOp::Inversion => 1.0 / s[0],
                    ScalarUnaryOp::SquareRoot => s[0].sqrt(),
                }
            }
            N::GradeProjection(e, _) => e.add_to_res(cache, res),
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
            assert_eq!($a.prepare(&$alg).eval::<GradeMapMV>(), $b);
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
        expr_eq!(EGA3, (v & bv.clone()) & bv.inv(), grade_map_mv!(1 => 1 0 0));
    }
}
