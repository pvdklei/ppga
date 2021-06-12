use crate::error;

#[derive(Debug)]
pub struct Motor {
    pub scalar: f32,
    pub v_bivector: [f32; 3],
    pub e_bivector: [f32; 3],
    pub pseudo: f32,
}

impl Motor {
    /// Creates a motor that moves every param (e.g., a) to its destination (e.g., a_).
    /// If this transformation cannot be represented by a motor (i.e., is not orthogonal),
    /// then a invalid (possibly filled with NaN or Inf values) motor is returned.
    ///
    /// See PGA4CS page 62 for more info.
    pub fn from_point_correspondences(
        a: &super::Point,
        a_: &super::Point,
        b: &super::Point,
        b_: &super::Point,
        c: &super::Point,
        c_: &super::Point,
    ) -> Self {
        let v_a = a_.div(&a).sqrt();
        let b_a = v_a.apply_to_point(&b);

        let from = super::join::points(&a_, &b_a);
        let to = super::join::points(&a_, &b_);

        let v_ba = to.div(&from).sqrt().mul_translator(&v_a);
        let c_ba = v_ba.apply_to_point(&c);

        let from = super::join::three_points(&a_, &b_, &c_ba);
        let to = super::join::three_points(&a_, &b_, &c_);

        to.div(&from).sqrt().mul(&v_ba)
    }

    pub fn into_rotor_checked(&self) -> Result<super::Rotor, error::CastError<Self>> {
        if self.v_bivector.iter().any(|e| e.abs() > 0.1) || self.pseudo.abs() > 0.1 {
            return Err(error::CastError::new(self, "Rotor"));
        }
        Ok(self.into_rotor_unchecked())
    }

    pub fn into_rotor_unchecked(&self) -> super::Rotor {
        super::Rotor {
            scalar: self.scalar,
            e_bivector: self.e_bivector,
        }
    }

    pub fn norm(&self) -> f32 {
        (self.e_bivector[0] * self.e_bivector[0]
            + self.e_bivector[1] * self.e_bivector[1]
            + self.e_bivector[2] * self.e_bivector[2]
            + self.scalar * self.scalar)
            .sqrt()
    }

    pub fn normalize(&self) -> Self {
        let fac = 1. / self.norm();
        Self {
            scalar: self.scalar * fac,
            pseudo: self.pseudo * fac,
            v_bivector: {
                let v = na::Vector3::from(self.v_bivector);
                (v * fac).into()
            },
            e_bivector: {
                let e = na::Vector3::from(self.e_bivector);
                (e * fac).into()
            },
        }
    }

    pub fn sqrt(&self) -> Self {
        let s = self.scalar;
        let ps = self.pseudo;
        let e = self.e_bivector;
        let v = self.v_bivector;
        let everywhere = (-ps + 2. * s + 2.) * 2.0f32.sqrt();
        let fac = 1. / (4. * (s + 1.).powf(1.5));
        let everywhere_fac = everywhere * fac;
        Self {
            scalar: everywhere / (4. * (s + 1.).sqrt()),
            pseudo: ps * everywhere_fac,
            v_bivector: [
                v[0] * everywhere_fac,
                v[1] * everywhere_fac,
                v[2] * everywhere_fac,
            ],
            e_bivector: [
                e[0] * everywhere_fac,
                e[1] * everywhere_fac,
                e[2] * everywhere_fac,
            ],
        }
    }

    /// PGA4CS page 69
    pub fn ln(&self) -> super::Line {
        let w = super::Line::from(self).div_scalar(self.scalar);
        let wrev = w.reverse();
        let wdotwrev = super::inner::lines(&w, &wrev);
        let wmeetwrev = super::meet::lines(&w, &wrev);
        let sqrt_wdotwrev = wdotwrev.sqrt();

        let a = sqrt_wdotwrev.atan();
        let we = w.e_bivector;
        let wv = w.v_bivector;
        super::Line {
            e_bivector: [
                a * we[0] / sqrt_wdotwrev,
                a * we[1] / sqrt_wdotwrev,
                a * we[2] / sqrt_wdotwrev,
            ],
            v_bivector: [
                a * wv[0] / sqrt_wdotwrev
                    + 0.5 * a * we[0] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev)
                    - 0.5 * we[0] * wmeetwrev.0 / wdotwrev,
                a * wv[1] / sqrt_wdotwrev
                    + 0.5 * a * we[1] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev)
                    - 0.5 * we[1] * wmeetwrev.0 / wdotwrev,
                a * wv[2] / sqrt_wdotwrev
                    + 0.5 * a * we[2] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev)
                    - 0.5 * we[2] * wmeetwrev.0 / wdotwrev,
            ],
        }

        // let bi = wmeetwrev.div_scalar(2. * sqrt_wdotwrev);
        // let we = w.e_bivector;
        // let wv = w.v_bivector;
        // let l = super::Line {
        //     e_bivector: [
        //         we[0] / sqrt_wdotwrev,
        //         we[1] / sqrt_wdotwrev,
        //         we[2] / sqrt_wdotwrev,
        //     ],
        //     v_bivector: [
        //         wv[0] / sqrt_wdotwrev + 0.5 * we[0] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev),
        //         wv[1] / sqrt_wdotwrev + 0.5 * we[1] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev),
        //         wv[2] / sqrt_wdotwrev + 0.5 * we[2] * wmeetwrev.0 / (sqrt_wdotwrev * wdotwrev),
        //     ],
        // };
        // let lv = l.v_bivector;
        // let le = l.e_bivector;
        // super::Line {
        //     e_bivector: [a * le[0], a * le[1], a * le[2]],
        //     v_bivector: [
        //         a * lv[0] - bi.0 * le[0],
        //         a * lv[1] - bi.0 * le[1],
        //         a * lv[2] - bi.0 * le[2],
        //     ],
        // }
    }

    pub fn mul(&self, other: &Self) -> Self {
        let s1 = self.scalar;
        let s2 = other.scalar;
        let ps1 = self.pseudo;
        let ps2 = other.pseudo;
        let v1 = self.v_bivector;
        let v2 = other.v_bivector;
        let e1 = self.e_bivector;
        let e2 = other.e_bivector;
        Self {
            scalar: -e1[0] * e2[0] - e1[1] * e2[1] - e1[2] * e2[2] + s1 * s2,
            pseudo: e1[0] * v2[0]
                + e1[1] * v2[1]
                + e1[2] * v2[2]
                + e2[0] * v1[0]
                + e2[1] * v1[1]
                + e2[2] * v1[2]
                + ps1 * s2
                + ps2 * s1,
            v_bivector: [
                -e1[0] * ps2 - e1[1] * v2[2] + e1[2] * v2[1] - e2[0] * ps1 + e2[1] * v1[2]
                    - e2[2] * v1[1]
                    + s1 * v2[0]
                    + s2 * v1[0],
                e1[0] * v2[2] - e1[1] * ps2 - e1[2] * v2[0] - e2[0] * v1[2] - e2[1] * ps1
                    + e2[2] * v1[0]
                    + s1 * v2[1]
                    + s2 * v1[1],
                -e1[0] * v2[1] + e1[1] * v2[0] - e1[2] * ps2 + e2[0] * v1[1]
                    - e2[1] * v1[0]
                    - e2[2] * ps1
                    + s1 * v2[2]
                    + s2 * v1[2],
            ],
            e_bivector: [
                e1[0] * s2 - e1[1] * e2[2] + e1[2] * e2[1] + e2[0] * s1,
                e1[0] * e2[2] + e1[1] * s2 - e1[2] * e2[0] + e2[1] * s1,
                -e1[0] * e2[1] + e1[1] * e2[0] + e1[2] * s2 + e2[2] * s1,
            ],
        }
    }

    pub fn mul_translator(&self, t: &super::Translator) -> Self {
        let ts = t.scalar;
        let tv = t.v_bivector;
        let ms = self.scalar;
        let mps = self.pseudo;
        let mv = self.v_bivector;
        let me = self.e_bivector;
        Self {
            scalar: ms * ts,
            pseudo: me[0] * tv[0] + me[1] * tv[1] + me[2] * tv[2] + mps * ts,
            v_bivector: [
                -me[1] * tv[2] + me[2] * tv[1] + ms * tv[0] + mv[0] * ts,
                me[0] * tv[2] - me[2] * tv[0] + ms * tv[1] + mv[1] * ts,
                -me[0] * tv[1] + me[1] * tv[0] + ms * tv[2] + mv[2] * ts,
            ],
            e_bivector: [me[0] * ts, me[1] * ts, me[2] * ts],
        }
    }

    pub fn apply_to_point(&self, p: &super::Point) -> super::Point {
        let ms = self.scalar;
        let mps = self.pseudo;
        let mv = self.v_bivector;
        let me = self.e_bivector;
        let p = p.trivector;
        super::Point {
            trivector: [
                p[0] * (me[0] * me[0] + me[1] * me[1] + me[2] * me[2] + ms * ms),
                me[0] * me[0] * p[1] + 2. * me[0] * me[1] * p[2] + 2. * me[0] * me[2] * p[3]
                    - 2. * me[0] * mps * p[0]
                    - me[1] * me[1] * p[1]
                    - 2. * me[1] * ms * p[3]
                    + 2. * me[1] * mv[2] * p[0]
                    - me[2] * me[2] * p[1]
                    + 2. * me[2] * ms * p[2]
                    - 2. * me[2] * mv[1] * p[0]
                    + ms * ms * p[1]
                    - 2. * ms * mv[0] * p[0],
                -me[0] * me[0] * p[2] + 2. * me[0] * me[1] * p[1] + 2. * me[0] * ms * p[3]
                    - 2. * me[0] * mv[2] * p[0]
                    + me[1] * me[1] * p[2]
                    + 2. * me[1] * me[2] * p[3]
                    - 2. * me[1] * mps * p[0]
                    - me[2] * me[2] * p[2]
                    - 2. * me[2] * ms * p[1]
                    + 2. * me[2] * mv[0] * p[0]
                    + ms * ms * p[2]
                    - 2. * ms * mv[1] * p[0],
                -me[0] * me[0] * p[3] + 2. * me[0] * me[2] * p[1] - 2. * me[0] * ms * p[2]
                    + 2. * me[0] * mv[1] * p[0]
                    - me[1] * me[1] * p[3]
                    + 2. * me[1] * me[2] * p[2]
                    + 2. * me[1] * ms * p[1]
                    - 2. * me[1] * mv[0] * p[0]
                    + me[2] * me[2] * p[3]
                    - 2. * me[2] * mps * p[0]
                    + ms * ms * p[3]
                    - 2. * ms * mv[2] * p[0],
            ],
        }
    }
    pub fn into_klein(&self) -> [[f32; 4]; 2] {
        let m = self;
        [
            [m.scalar, m.e_bivector[0], m.e_bivector[1], m.e_bivector[2]],
            [m.pseudo, m.v_bivector[0], m.v_bivector[1], m.v_bivector[2]],
        ]
    }
}

impl PartialEq for Motor {
    fn eq(&self, other: &Self) -> bool {
        self.e_bivector
            .iter()
            .zip(other.e_bivector.iter())
            .all(|(a, b)| (a - b).abs() < 0.01)
            && self
                .v_bivector
                .iter()
                .zip(other.v_bivector.iter())
                .all(|(a, b)| (a - b).abs() < 0.01)
            && (self.scalar - other.scalar).abs() < 0.01
            && (self.pseudo - other.pseudo).abs() < 0.01
    }
}

impl From<&super::Translator> for Motor {
    fn from(t: &super::Translator) -> Self {
        Self {
            scalar: t.scalar,
            e_bivector: [0.0f32; 3],
            v_bivector: t.v_bivector,
            pseudo: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[allow(dead_code)]
    fn test_motor() -> super::Motor {
        let p1 = Plane::new(3., &na::Vector3::new(1., 3., -2.).normalize().into());
        let p2 = Plane::new(-2., &na::Vector3::new(2., -3., -2.).normalize().into());
        p1.div(&p2).sqrt()
    }
    #[allow(dead_code)]
    fn rotating_test_motor() -> super::Motor {
        let p1 = Point::new(&[0., 0., 0.]);
        let p2 = Point::new(&[-6.4, 9.1, 0.4]);
        let p3 = Point::new(&[0.0, -10., 3.]);
        let l1 = join::points(&p1, &p2);
        let l2 = join::points(&p1, &p3);
        l1.div(&l2).sqrt()
    }
    #[allow(dead_code)]
    fn translating_test_motor() -> super::Motor {
        let p1 = Point::new(&[0., 0., 0.]);
        let p2 = Point::new(&[-6.4, 9.1, 0.4]);
        Motor::from(&p1.div(&p2).sqrt())
    }

    #[test]
    fn point_corr() {
        let a = Point::new(&[2., 3., 5.]);
        let b = Point::new(&[2., 8., 7.]);
        let c = Point::new(&[3., -2., 1.]);

        let rot = Rotor::new(2.6, &na::Vector3::from([1.2, 1., 0.]).normalize().into());
        let trans = Point::new(&[6., 4., 1.])
            .div(&Point::new(&[2., 0., 9.]))
            .sqrt();
        let rot = rot.mul_translator(&trans);

        let a_ = rot.apply_to_point(&a);
        let b_ = rot.apply_to_point(&b);
        let c_ = rot.apply_to_point(&c);

        let v_cba = Motor::from_point_correspondences(&a, &a_, &b, &b_, &c, &c_);
        assert_eq!(v_cba, rot);
    }

    #[test]
    fn logarithm1() {
        let m = test_motor().normalize();
        let m_ = m.ln().exp();
        println!("{:?}", m);
        println!("{:?}", m_);

        let p = Point::new(&[0.1, -3.9, 0.2]);
        let res = m.apply_to_point(&p);
        let res_ = m_.apply_to_point(&p);

        assert_eq!(res, res_);
    }
    #[test]
    fn logarithm2() {
        let l = Line {
            v_bivector: [2., -3., 1.],
            e_bivector: [-9., 4., 3.5],
        };
        let l_ = l.exp().ln();
        assert_eq!(l, l_);
    }
}
