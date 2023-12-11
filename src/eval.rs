//! How to evaluate a GA expression into an actual multivector

use std::collections::HashMap;

use crate::{ast::*, graded::*, GradeSet};
use AstNode as N;

type Cache<R> = HashMap<NodeId, R>;

impl<T: GradedData + std::fmt::Debug> SpecializedAst<T> {
    /// Evaluates a [`SpecializedAst`]
    pub fn eval<R>(&self) -> R
    where
        R: GradedDataMut,
    {
        let mut cache = HashMap::new();
        self.store_in_cache(self.root_id(), &mut cache);
        return cache.remove(&self.root_id()).unwrap();
    }

    fn store_in_cache<R>(&self, this_id: NodeId, cache: &mut Cache<R>)
    where
        R: GradedDataMut,
    {
        let this = self.get_node(this_id);
        if let None = cache.get(&this_id) {
            cache.insert(
                this_id,
                R::init_null_mv(this.vec_space_dim(), this.grade_set()),
            );
            self.add_to_res(this_id, this_id, cache);
        }
    }

    fn add_to_res<R>(&self, res_id: NodeId, this_id: NodeId, cache: &mut Cache<R>)
    where
        R: GradedDataMut,
    {
        let this = self.get_node(this_id);
        if this.grade_set().is_empty() {
            // self necessarily evaluates to zero, no need to go further
            return;
        }
        match this.ast_node() {
            N::GradedObj(input) => {
                cache
                    .get_mut(&res_id)
                    .unwrap()
                    .add_grades_from(input, this.grade_set());
            }
            N::Addition(left_expr, right_expr) => {
                self.add_to_res(res_id, *left_expr, cache);
                self.add_to_res(res_id, *right_expr, cache);
            }
            N::Negation(e) => {
                self.add_to_res(res_id, *e, cache);
                for k in this.grade_set().iter() {
                    cache.get_mut(&res_id).unwrap().negate_grade(k);
                }
            }
            N::Product(Product {
                individual_comp_muls,
                left_expr,
                right_expr,
                ..
            }) => {
                self.store_in_cache(*left_expr, cache);
                self.store_in_cache(*right_expr, cache);
                // We move res out of the cache so we can mutate it:
                let mut res = std::mem::replace(
                    cache.get_mut(&res_id).unwrap(),
                    R::init_null_mv(0, &GradeSet::empty()),
                );
                // Cache is thus only borrowed immutably:
                let left = &cache[left_expr];
                let right = &cache[right_expr];
                for mul in individual_comp_muls {
                    let val_left = left.grade_slice(mul.left_comp.grade)[mul.left_comp.index];
                    let val_right = right.grade_slice(mul.right_comp.grade)[mul.right_comp.index];
                    let val_result =
                        &mut res.grade_slice_mut(mul.result_comp.grade)[mul.result_comp.index];
                    *val_result += val_left * val_right * mul.coeff;
                }
                // We move res back into the cache:
                std::mem::swap(cache.get_mut(&res_id).unwrap(), &mut res);
            }
            N::Reverse(e) => {
                self.add_to_res(res_id, *e, cache);
                for k in this.grade_set().iter() {
                    if (k * (k - 1) / 2) % 2 == 1 {
                        cache.get_mut(&res_id).unwrap().negate_grade(k);
                    }
                }
            }
            N::GradeInvolution(e) => {
                self.add_to_res(res_id, *e, cache);
                for k in this.grade_set().iter() {
                    if k % 2 == 1 {
                        cache.get_mut(&res_id).unwrap().negate_grade(k);
                    }
                }
            }
            N::ScalarUnaryOp(op, e) => {
                self.add_to_res(res_id, *e, cache);
                let s = cache.get_mut(&res_id).unwrap().grade_slice_mut(0);
                s[0] = match op {
                    ScalarUnaryOp::Inversion => 1.0 / s[0],
                    ScalarUnaryOp::SquareRoot => s[0].sqrt(),
                }
            }
            N::GradeProjection(e) => self.add_to_res(res_id, *e, cache),
            N::Exponential(_e) => todo!(),
            N::Logarithm(_e) => todo!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{algebra::OrthoEuclidN, grade_map_mv, graded::GradeMapMV, Expr};

    macro_rules! expr_eq {
        ($alg:ident, $a:expr, $b:expr) => {
            let expr = $a;
            let specialized = expr.specialize(&$alg);
            assert_eq!(specialized.eval::<GradeMapMV>(), $b);
        };
    }

    type E = Expr<'static, GradeMapMV>;
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
