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

    pub fn neg(&self) -> Self {
        Self {
            scalar: -self.scalar,
            pseudo: -self.pseudo,
            v_bivector: (-na::Vec3::from(self.v_bivector)).into(),
            e_bivector: (-na::Vec3::from(self.e_bivector)).into(),
        }
    }

    pub fn norm(&self) -> f32 {
        (self.e_bivector[0] * self.e_bivector[0]
            + self.e_bivector[1] * self.e_bivector[1]
            + self.e_bivector[2] * self.e_bivector[2]
            + self.scalar * self.scalar)
            .sqrt()
    }
    pub fn norm_squared(&self) -> f32 {
        self.e_bivector[0] * self.e_bivector[0]
            + self.e_bivector[1] * self.e_bivector[1]
            + self.e_bivector[2] * self.e_bivector[2]
            + self.scalar * self.scalar
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

    /// Check whether the motor has no grade 4 part, and therefore
    /// called 'simple'. E.g., rotation that is not around the origin.
    pub fn is_simple(&self) -> bool {
        self.pseudo.abs() < 0.001
    }

    pub fn sqrt(&self) -> Self {
        let s = self.scalar;
        let ps = self.pseudo;
        let e = self.e_bivector;
        let v = self.v_bivector;
        fn sqrt(x: f32) -> f32 {
            x.sqrt()
        }
        // TODO: Optimize this
        Self {
            scalar: sqrt(2. * s + 2.) / 2.,
            pseudo: sqrt(2.) * ps / (4. * sqrt(s + 1.)),
            e_bivector: [
                sqrt(2.) * e[0] / (2. * sqrt(s + 1.)),
                sqrt(2.) * e[1] / (2. * sqrt(s + 1.)),
                sqrt(2.) * e[2] / (2. * sqrt(s + 1.)),
            ],
            v_bivector: [
                sqrt(2.) * e[0] * ps / (4. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * s * v[0] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * v[0] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.))),
                sqrt(2.) * e[1] * ps / (4. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * s * v[1] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * v[1] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.))),
                sqrt(2.) * e[2] * ps / (4. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * s * v[2] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.)))
                    + sqrt(2.) * v[2] / (2. * (s * sqrt(s + 1.) + sqrt(s + 1.))),
            ],
        }
    }

    /// Square root for a simple motor (no grade 4 part)
    pub fn ssqrt(&self) -> Self {
        self.add_scalar(1.).normalize()
    }

    /// PGA4CS page 69
    pub fn ln_4cs(&self) -> super::Line {
        let w = super::Line::from(self).div_scalar(self.scalar);
        let wrev = w.reverse();
        let wdotwrev = super::inner::lines(&w, &wrev);
        let wmeetwrev = super::meet::lines(&w, &wrev);
        let sqrt_wdotwrev = wdotwrev.sqrt();
        let a = sqrt_wdotwrev.atan();

        println!("Dot {:?}", wdotwrev);
        println!("SqrtDot {:?}", sqrt_wdotwrev);
        println!("Meet {:?}", wmeetwrev);
        println!("Alpha {:?}", a);

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
    }

    /// SIGGRAPH Course Notes 8.1.6.
    pub fn ln(&self) -> super::Line {
        let s1 = self.scalar;
        let p1 = self.pseudo;
        let b = super::Line::from(self);
        let bdb = -super::inner::lines(&b, &b);
        // When self is a translator nothing has to be done
        if bdb.abs() < 0.001 {
            return super::Line::from(self);
        }
        let s2 = bdb.sqrt();
        let p2 = super::meet::lines(&b, &b).mul_scalar(-2. * s2);

        let (u, v) = if s1.abs() > 0.001 {
            let u = s2.atan2(s1);
            let v = p2.div_scalar(s1);
            (u, v)
        } else {
            let u = (-p1).atan2(p2.0);
            let v = super::PseudoScalar(-p1 / s2);
            (u, v)
        };
        let be = b.e_bivector;
        let bv = b.v_bivector;
        let s2sq = s2 * s2;
        super::Line {
            e_bivector: [be[0] * u / s2, be[1] * u / s2, be[2] * u / s2],
            v_bivector: [
                -be[0] * p2.0 * u / s2sq - be[0] * v.0 / s2 + bv[0] * u / s2,
                -be[1] * p2.0 * u / s2sq - be[1] * v.0 / s2 + bv[1] * u / s2,
                -be[2] * p2.0 * u / s2sq - be[2] * v.0 / s2 + bv[2] * u / s2,
            ],
        }
    }

    /// Just like a log, maps a motor to a bivector and can be reversed.
    /// Tingelstad 2018 (https://link.springer.com/article/10.1007/s00006-018-0850-2)
    pub fn cayley(&self) -> super::Line {
        super::Line::from(&self.neg().add_scalar(1.0).div(&self.add_scalar(1.0)))
    }

    /// Again a map to the bivectors (Tingelstad, 2018)
    pub fn outer_ln(&self) -> super::Line {
        unimplemented!();
    }
    /// A simplified version on Tingelstad (2018) Steven de
    /// Keninck sent to me in an image.
    pub fn outer_ln_steven(&self) -> super::Line {
        super::Line::from(self).div_scalar(self.scalar)
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

    pub fn div(&self, other: &Self) -> Self {
        self.mul(&other.inverse())
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

    pub fn apply_to<T: Applicable>(&self, g: &T) -> T {
        g.apply(self)
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

    pub fn apply_to_plane(&self, p: &super::Plane) -> super::Plane {
        let pvec = p.vector;
        let ms = self.scalar;
        let mps = self.pseudo;
        let mv = self.v_bivector;
        let me = self.e_bivector;
        super::Plane {
            vector: [
                me[0] * me[0] * pvec[0] + 2. * me[0] * mps * pvec[1] + 2. * me[0] * mv[1] * pvec[3]
                    - 2. * me[0] * mv[2] * pvec[2]
                    + me[1] * me[1] * pvec[0]
                    + 2. * me[1] * mps * pvec[2]
                    - 2. * me[1] * mv[0] * pvec[3]
                    + 2. * me[1] * mv[2] * pvec[1]
                    + me[2] * me[2] * pvec[0]
                    + 2. * me[2] * mps * pvec[3]
                    + 2. * me[2] * mv[0] * pvec[2]
                    - 2. * me[2] * mv[1] * pvec[1]
                    + ms * ms * pvec[0]
                    + 2. * ms * mv[0] * pvec[1]
                    + 2. * ms * mv[1] * pvec[2]
                    + 2. * ms * mv[2] * pvec[3],
                me[0] * me[0] * pvec[1]
                    + 2. * me[0] * me[1] * pvec[2]
                    + 2. * me[0] * me[2] * pvec[3]
                    - me[1] * me[1] * pvec[1]
                    - 2. * me[1] * ms * pvec[3]
                    - me[2] * me[2] * pvec[1]
                    + 2. * me[2] * ms * pvec[2]
                    + ms * ms * pvec[1],
                -me[0] * me[0] * pvec[2]
                    + 2. * me[0] * me[1] * pvec[1]
                    + 2. * me[0] * ms * pvec[3]
                    + me[1] * me[1] * pvec[2]
                    + 2. * me[1] * me[2] * pvec[3]
                    - me[2] * me[2] * pvec[2]
                    - 2. * me[2] * ms * pvec[1]
                    + ms * ms * pvec[2],
                -me[0] * me[0] * pvec[3] + 2. * me[0] * me[2] * pvec[1]
                    - 2. * me[0] * ms * pvec[2]
                    - me[1] * me[1] * pvec[3]
                    + 2. * me[1] * me[2] * pvec[2]
                    + 2. * me[1] * ms * pvec[1]
                    + me[2] * me[2] * pvec[3]
                    + ms * ms * pvec[3],
            ],
        }
    }

    pub fn apply_to_line(&self, l: &super::Line) -> super::Line {
        let ms = self.scalar;
        let mps = self.pseudo;
        let mv = self.v_bivector;
        let me = self.e_bivector;
        let le = l.e_bivector;
        let lv = l.v_bivector;
        super::Line {
            e_bivector: [
                le[0] * me[0] * me[0] - le[0] * me[1] * me[1] - le[0] * me[2] * me[2]
                    + le[0] * ms * ms
                    + 2. * le[1] * me[0] * me[1]
                    + 2. * le[1] * me[2] * ms
                    + 2. * le[2] * me[0] * me[2]
                    - 2. * le[2] * me[1] * ms,
                2. * le[0] * me[0] * me[1] - 2. * le[0] * me[2] * ms - le[1] * me[0] * me[0]
                    + le[1] * me[1] * me[1]
                    - le[1] * me[2] * me[2]
                    + le[1] * ms * ms
                    + 2. * le[2] * me[0] * ms
                    + 2. * le[2] * me[1] * me[2],
                2. * le[0] * me[0] * me[2] + 2. * le[0] * me[1] * ms - 2. * le[1] * me[0] * ms
                    + 2. * le[1] * me[1] * me[2]
                    - le[2] * me[0] * me[0]
                    - le[2] * me[1] * me[1]
                    + le[2] * me[2] * me[2]
                    + le[2] * ms * ms,
            ],
            v_bivector: [
                2. * le[0] * me[0] * mv[0]
                    - 2. * le[0] * me[1] * mv[1]
                    - 2. * le[0] * me[2] * mv[2]
                    - 2. * le[0] * mps * ms
                    + 2. * le[1] * me[0] * mv[1]
                    + 2. * le[1] * me[1] * mv[0]
                    - 2. * le[1] * me[2] * mps
                    + 2. * le[1] * ms * mv[2]
                    + 2. * le[2] * me[0] * mv[2]
                    + 2. * le[2] * me[1] * mps
                    + 2. * le[2] * me[2] * mv[0]
                    - 2. * le[2] * ms * mv[1]
                    + lv[0] * me[0] * me[0]
                    - lv[0] * me[1] * me[1]
                    - lv[0] * me[2] * me[2]
                    + lv[0] * ms * ms
                    + 2. * lv[1] * me[0] * me[1]
                    + 2. * lv[1] * me[2] * ms
                    + 2. * lv[2] * me[0] * me[2]
                    - 2. * lv[2] * me[1] * ms,
                2. * le[0] * me[0] * mv[1] + 2. * le[0] * me[1] * mv[0] + 2. * le[0] * me[2] * mps
                    - 2. * le[0] * ms * mv[2]
                    - 2. * le[1] * me[0] * mv[0]
                    + 2. * le[1] * me[1] * mv[1]
                    - 2. * le[1] * me[2] * mv[2]
                    - 2. * le[1] * mps * ms
                    - 2. * le[2] * me[0] * mps
                    + 2. * le[2] * me[1] * mv[2]
                    + 2. * le[2] * me[2] * mv[1]
                    + 2. * le[2] * ms * mv[0]
                    + 2. * lv[0] * me[0] * me[1]
                    - 2. * lv[0] * me[2] * ms
                    - lv[1] * me[0] * me[0]
                    + lv[1] * me[1] * me[1]
                    - lv[1] * me[2] * me[2]
                    + lv[1] * ms * ms
                    + 2. * lv[2] * me[0] * ms
                    + 2. * lv[2] * me[1] * me[2],
                2. * le[0] * me[0] * mv[2] - 2. * le[0] * me[1] * mps
                    + 2. * le[0] * me[2] * mv[0]
                    + 2. * le[0] * ms * mv[1]
                    + 2. * le[1] * me[0] * mps
                    + 2. * le[1] * me[1] * mv[2]
                    + 2. * le[1] * me[2] * mv[1]
                    - 2. * le[1] * ms * mv[0]
                    - 2. * le[2] * me[0] * mv[0]
                    - 2. * le[2] * me[1] * mv[1]
                    + 2. * le[2] * me[2] * mv[2]
                    - 2. * le[2] * mps * ms
                    + 2. * lv[0] * me[0] * me[2]
                    + 2. * lv[0] * me[1] * ms
                    - 2. * lv[1] * me[0] * ms
                    + 2. * lv[1] * me[1] * me[2]
                    - lv[2] * me[0] * me[0]
                    - lv[2] * me[1] * me[1]
                    + lv[2] * me[2] * me[2]
                    + lv[2] * ms * ms,
            ],
        }
    }

    pub fn inverse(&self) -> Self {
        self.reverse().div_scalar(self.norm_squared())
    }

    pub fn mul_scalar(&self, s: f32) -> Self {
        Self {
            scalar: self.scalar * s,
            pseudo: self.pseudo * s,
            e_bivector: (na::Vector3::from(self.e_bivector) * s).into(),
            v_bivector: (na::Vector3::from(self.v_bivector) * s).into(),
        }
    }
    pub fn div_scalar(&self, s: f32) -> Self {
        self.mul_scalar(1. / s)
    }
    pub fn into_klein(&self) -> [[f32; 4]; 2] {
        let m = self;
        [
            [m.scalar, m.e_bivector[0], m.e_bivector[1], m.e_bivector[2]],
            [m.pseudo, m.v_bivector[0], m.v_bivector[1], m.v_bivector[2]],
        ]
    }
    pub fn is_similar_to(&self, d: f32, other: &Self) -> bool {
        let p = super::Point::random().normalize();
        self.apply_to_point(&p)
            .is_similar_to(d, &other.apply_to_point(&p))
    }
    pub fn random() -> Self {
        Self {
            scalar: rand::random(),
            pseudo: rand::random(), // or should this be 0??
            e_bivector: na::Vector3::new_random().into(),
            v_bivector: na::Vector3::new_random().into(),
        }
    }

    pub fn squared(&self) -> Self {
        self.mul(&self)
    }
    pub fn add_scalar(&self, s: f32) -> Self {
        Self {
            scalar: self.scalar + s,
            ..*self
        }
    }
    pub fn reverse(&self) -> Self {
        Self {
            e_bivector: (-na::Vec3::from(self.e_bivector)).into(),
            v_bivector: (-na::Vec3::from(self.v_bivector)).into(),
            ..*self
        }
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

impl From<&super::Line> for Motor {
    fn from(l: &super::Line) -> Self {
        Self {
            scalar: 0.,
            pseudo: 0.,
            v_bivector: l.v_bivector,
            e_bivector: l.e_bivector,
        }
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

impl From<&super::Rotor> for Motor {
    fn from(r: &super::Rotor) -> Self {
        Self {
            scalar: r.scalar,
            e_bivector: r.e_bivector,
            v_bivector: [0.; 3],
            pseudo: 0.,
        }
    }
}

pub trait Applicable {
    fn apply(&self, m: &Motor) -> Self;
}

impl Applicable for super::Point {
    fn apply(&self, m: &Motor) -> Self {
        m.apply_to_point(self)
    }
}
impl Applicable for super::Line {
    fn apply(&self, m: &Motor) -> Self {
        m.apply_to_line(self)
    }
}
impl Applicable for super::Plane {
    fn apply(&self, m: &Motor) -> Self {
        m.apply_to_plane(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[allow(dead_code)]
    fn test_motor() -> super::Motor {
        let p1 = Plane::new(3., &na::Vector3::new_random().normalize().into());
        let p2 = Plane::new(-2., &na::Vector3::new_random().normalize().into());
        let p3 = Plane::new(-5., &na::Vector3::new_random().normalize().into());
        let p4 = Plane::new(20., &na::Vector3::new_random().normalize().into());
        let p5 = Plane::random().nnormalize();
        let m = p1.move_to(&p2).mul(&p3.move_to(&p4)).mul(&p3.move_to(&p5));
        m
    }
    fn test_motor2() -> super::Motor {
        // let p1 = Plane::new(1., &na::Vec3::new(1., 5., -1.2).normalize().into());
        // let p2 = Plane::new(6., &na::Vec3::new(0., 1., -1.).normalize().into());
        // let p3 = Plane::new(-4., &na::Vec3::new(1., 41., 0.).normalize().into());
        let p1 = Plane::new(1., &na::Vec3::new_random().normalize().into());
        let p2 = Plane::new(6., &na::Vec3::new_random().normalize().into());
        let p3 = Plane::new(-4., &na::Vec3::new_random().normalize().into());
        p1.move_to(&p2).mul(&p2.move_to(&p3))
    }
    #[allow(dead_code)]
    fn rotating_test_motor() -> super::Motor {
        let p1 = Point::new(&[0., 0., 0.]);
        let p2 = Point::new(&[-6.4, 9.1, 0.4]);
        let p3 = Point::new(&[0.0, -10., 3.]);
        let l1 = join::points(&p1, &p2).normalize();
        let l2 = join::points(&p1, &p3).normalize();
        l1.div(&l2).sqrt()
    }
    #[allow(dead_code)]
    fn translating_test_motor() -> super::Motor {
        let p1 = Point::new(&[0., 0., 0.]);
        let p2 = Point::new(&[-6.4, 9.1, 0.4]);
        Motor::from(&p1.div(&p2).sqrt())
    }

    #[test]
    fn sqrt1() {
        // Taking the square root and then squaring
        // should leave it unchanged
        let m = test_motor();
        let m_ = m.sqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn sqrt2() {
        let m = test_motor2();
        let m_ = m.sqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn sqrt3() {
        let m = rotating_test_motor();
        let m_ = m.sqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
        let m = translating_test_motor();
        let m_ = m.sqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn sqrt4() {
        // For a translator and rotor the simple sqrt can be used.
        let m = rotating_test_motor();
        let m_ = m.ssqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
        let m = translating_test_motor();
        let m_ = m.ssqrt().squared();
        assert_eq!(m, m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn sqrt5() {
        // Taking the square root should only perform half the transformation.
        // so a full turn becomes a half turn. Twice the square root should do
        // a quarter of the rotation, etc
        let r = Rotor::new(std::f32::consts::PI * 2., &[0., 1., 0.]);
        let m = Motor::from(&r);
        let p = m.sqrt().apply_to_point(&Point::new(&[1., 0., 0.]));
        assert_eq!(p, Point::new(&[-1., 0., 0.]));
        let p = m.sqrt().sqrt().apply_to_point(&Point::new(&[1., 0., 0.]));
        assert_eq!(p, Point::new(&[0., 0., -1.]));
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
        let m = test_motor2();
        let m_ = m.ln().exp();
        println!("{:?}", m);
        println!("{:?}", m_);
        // assert_eq!(m, m_);
        let p = Point::new(&[2., 3., 4.]);
        println!("{:?}", p);
        println!("{:?}", m.apply_to_point(&p));
        println!("{:?}", m_.apply_to_point(&p));
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn logarithm2() {
        let l = join::points(&Point::random(), &Point::random());
        let l_ = l.exp().ln();
        println!("{:?}", l.exp());
        println!("{:?}", l_.exp());
        assert_eq!(l, l_);
    }

    #[test]
    fn logrithm3() {
        let m = rotating_test_motor();
        assert_eq!(m, m.ln().exp());
    }
    #[test]
    fn logrithm4() {
        let m = translating_test_motor();
        assert_eq!(m, m.ln().exp());
    }

    #[test]
    fn cayley() {
        let m = test_motor();
        assert_eq!(m, m.cayley().cayley())
    }

    #[test]
    fn apply_plane() {
        let t = Translator::new(&[1., 0., 0.]);
        let r = Rotor::new(std::f32::consts::PI, &[0., 1., 0.]);
        let m = t.mul_rotor(&r);
        let p = Plane::new(0., &[1., 0., 0.]);
        assert_eq!(m.apply_to_plane(&p), Plane::new(-1., &[-1., 0., 0.]))
    }

    #[test]
    fn outer_exp_steven() {
        let m = test_motor();
        let m_ = m.outer_ln_steven().outer_exp_steven().normalize();
        println!("{:?}", m);
        println!("{:?}", m_);
        assert!(m.is_similar_to(0.1, &m_));
        let m = rotating_test_motor();
        let m_ = m.outer_ln_steven().outer_exp_steven().normalize();
        println!("{:?}", m);
        println!("{:?}", m_);
        assert!(m.is_similar_to(0.1, &m_));
        let m = translating_test_motor();
        let m_ = m.outer_ln_steven().outer_exp_steven().normalize();
        println!("{:?}", m);
        println!("{:?}", m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
}
