#[derive(Debug)]
pub struct Translator {
    pub scalar: f32,
    pub v_bivector: [f32; 3],
}

impl Translator {
    pub fn new([t1, t2, t3]: &[f32; 3]) -> Self {
        Self {
            scalar: 1.,
            v_bivector: [t1 * 0.5, t2 * 0.5, t3 * 0.5],
        }
    }
    pub fn sqrt(&self) -> Self {
        let ts = self.scalar;
        let tv = self.v_bivector;
        let fac = 2.0f32.sqrt() / (2. * (ts + 1.).sqrt());
        Self {
            scalar: 0.5 * (2. * ts + 2.).sqrt(),
            v_bivector: [tv[0] * fac, tv[1] * fac, tv[2] * fac],
        }
    }

    pub fn apply_to_point(&self, p: &super::Point) -> super::Point {
        let v = self.v_bivector;
        let s = self.scalar;
        let p = p.trivector;
        super::Point {
            trivector: [
                p[0] * s * s,
                s * (-2. * p[0] * v[0] + p[1] * s),
                s * (-2. * p[0] * v[1] + p[2] * s),
                s * (-2. * p[0] * v[2] + p[3] * s),
            ],
        }
    }

    pub fn mul_rotor(&self, r: &super::Rotor) -> super::Motor {
        let v = self.v_bivector;
        let e = r.e_bivector;
        let ts = self.scalar;
        let rs = r.scalar;
        super::Motor {
            scalar: rs * ts,
            pseudo: e[0] * v[0] + e[1] * v[1] + e[2] * v[2],
            e_bivector: [e[0] * ts, e[1] * ts, e[2] * ts],
            v_bivector: [
                e[1] * v[2] - e[2] * v[1] + rs * v[0],
                -e[0] * v[2] + e[2] * v[0] + rs * v[1],
                e[0] * v[1] - e[1] * v[0] + rs * v[2],
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn new() {
        let t = Translator::new(&[1., 2., -3.]);
        let p = Point::new(&[0., 0., 0.]);
        assert_eq!(t.apply_to_point(&p), Point::new(&[1., 2., -3.]));
    }
}
