#[derive(Debug)]
pub struct GradeSet(Vec<bool>);

impl GradeSet {
    pub fn g_null() -> Self {
        GradeSet(Vec::new())
    }
    pub fn g(x: usize) -> Self {
        let mut v = vec![false; x + 1];
        v[x] = true;
        GradeSet(v)
    }
}

fn sort_vecs(GradeSet(v1): GradeSet, GradeSet(v2): GradeSet) -> (Vec<bool>, Vec<bool>) {
    if v1.len() <= v2.len() {
        (v1, v2)
    } else {
        (v2, v1)
    }
}

impl std::ops::Add for GradeSet {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (small, mut big) = sort_vecs(self, rhs);
        for i in 0..small.len() {
            big[i] = big[i] || small[i];
        }
        GradeSet(big)
    }
}

impl std::ops::Mul for GradeSet {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let (small, big) = sort_vecs(self, rhs);
        if small.len() == 0 { return GradeSet(small) };
        let mut res = vec![false; big.len()];
        for r in 0..res.len() {
            for i in 0..small.len() {
                for j in 0..big.len() {
                    let m = (i as i32 - j as i32).abs();
                    if i+j >= r && m <= r as i32 && m % 2 == r as i32 % 2 {
                        res[r] = res[r] || (small[i] && big[j]);
                    }
                }
            }
        }
        GradeSet(res)
    }
}
