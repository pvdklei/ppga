#[derive(Debug, Copy, Clone)]
pub struct Line {
    pub e_bivector: [f32; 3],
    pub v_bivector: [f32; 3],
}

impl Line {
    /// PGA4CS page 29
    pub fn new(point: &[f32; 3], dir: &[f32; 3]) -> Self {
        Self {
            v_bivector: [
                point[2] * dir[1] - point[1] * dir[2],
                point[0] * dir[2] - point[2] * dir[0],
                point[1] * dir[0] - point[0] * dir[1],
            ],
            e_bivector: *dir,
        }
    }

    pub fn random() -> Self {
        Self {
            v_bivector: na::Vec3::new_random().into(),
            e_bivector: na::Vec3::new_random().into(),
        }
    }

    pub fn zero() -> Self {
        Self {
            v_bivector: [0.0; 3],
            e_bivector: [0.0; 3],
        }
    }

    pub fn dual(&self) -> Self {
        Self {
            v_bivector: self.e_bivector,
            e_bivector: self.v_bivector,
        }
    }

    pub fn reverse(&self) -> Self {
        self.neg()
    }

    pub fn neg(&self) -> Self {
        Self {
            e_bivector: (-na::Vector3::from(self.e_bivector)).into(),
            v_bivector: (-na::Vector3::from(self.v_bivector)).into(),
        }
    }

    pub fn inverse(&self) -> Self {
        let e = self.e_bivector;
        let v = self.v_bivector;
        let fac = 1. / (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        Self {
            e_bivector: [-e[0] * fac, -e[1] * fac, -e[2] * fac],
            v_bivector: [-v[0] * fac, -v[1] * fac, -v[2] * fac],
        }
    }

    pub fn norm(&self) -> f32 {
        na::Vec3::from(self.e_bivector).norm()
    }
    pub fn inorm(&self) -> f32 {
        na::Vec3::from(self.v_bivector).norm()
    }
    pub fn normalize(&self) -> Self {
        Self {
            v_bivector: self.v_bivector,
            e_bivector: na::Vec3::from(self.e_bivector).normalize().into(),
        }
    }

    pub fn is_zero(&self) -> bool {
        na::Vec3::from(self.e_bivector).is_empty() && na::Vec3::from(self.v_bivector).is_empty()
    }

    /// PGA4CS chapter 5.6 and 7
    /// Has some numerical stability issues, don't know why.
    /// Also no exceptions handled, f.i. when e == 0
    /// (vanishing line / translation).
    pub fn exp_unstable(&self) -> super::Motor {
        let (e, v) = self.decompose();
        if e.is_zero() {
            return super::Motor::from(self).add_scalar(1.);
        }
        let half_phi = -e.norm();
        let e_hat = e.normalize();
        let cos_half_phi = half_phi.cos();
        let sin_half_phi = half_phi.sin();

        let eucl_e = e_hat.e_bivector;
        let eucl_v = e_hat.v_bivector;
        let van_v = v.v_bivector;

        super::Motor {
            scalar: cos_half_phi,
            pseudo: -eucl_e[0] * sin_half_phi * van_v[0]
                - eucl_e[1] * sin_half_phi * van_v[1]
                - eucl_e[2] * sin_half_phi * van_v[2],
            v_bivector: [
                cos_half_phi * van_v[0] - eucl_e[1] * sin_half_phi * van_v[2]
                    + eucl_e[2] * sin_half_phi * van_v[1]
                    - eucl_v[0] * sin_half_phi,
                cos_half_phi * van_v[1] + eucl_e[0] * sin_half_phi * van_v[2]
                    - eucl_e[2] * sin_half_phi * van_v[0]
                    - eucl_v[1] * sin_half_phi,
                cos_half_phi * van_v[2] - eucl_e[0] * sin_half_phi * van_v[1]
                    + eucl_e[1] * sin_half_phi * van_v[0]
                    - eucl_v[2] * sin_half_phi,
            ],
            e_bivector: [
                -eucl_e[0] * sin_half_phi,
                -eucl_e[1] * sin_half_phi,
                -eucl_e[2] * sin_half_phi,
            ],
        }
    }

    /// Just like the exponent, maps to a motor space and back (inverse).
    pub fn cayley(&self) -> super::Motor {
        let x = super::Motor::from(self);
        x.neg().add_scalar(1.0).inverse().mul(&x.add_scalar(1.0))
    }

    /// SIGGRAPH Course Notes 8.1.3 & 8.1.4.
    /// Implementation taken from ganja.js. More stable.
    pub fn exp(&self) -> super::Motor {
        let bdb = -super::inner::lines(&self, &self);
        let u = bdb.sqrt();
        if u < 0.001 {
            return super::Motor::from(self).add_scalar(1.);
        }
        let v = super::meet::lines(&self, &self).mul_scalar(-2. * u);
        let cu = u.cos();
        let su = u.sin();
        let be = self.e_bivector;
        let bv = self.v_bivector;
        let u2 = u * u;
        super::Motor {
            scalar: cu,
            pseudo: -su * v.0,
            v_bivector: [
                -be[0] * cu * v.0 / u - be[0] * su * v.0 / u2 + bv[0] * su / u,
                -be[1] * cu * v.0 / u - be[1] * su * v.0 / u2 + bv[1] * su / u,
                -be[2] * cu * v.0 / u - be[2] * su * v.0 / u2 + bv[2] * su / u,
            ],
            e_bivector: [be[0] * su / u, be[1] * su / u, be[2] * su / u],
        }
    }

    /// PGA4CS chapter 5.6
    /// Decomposes a line into a vanishing and euclidian line,
    /// so that
    ///     self = vanishing + euclidian.
    /// These parts commute, so that
    ///     vanishing * euclidian = euclidian * vanishing
    pub fn decompose(&self) -> (Line, Line) {
        let rev = self.reverse();
        let bdb = super::inner::lines(&self, &rev);
        if bdb.abs() < 0.001 {
            return (Line::zero(), *self);
        }
        let bmb = super::meet::lines(&self, &rev);
        let be = self.e_bivector;
        let bv = self.v_bivector;
        let van = Line {
            v_bivector: [
                -0.5 * be[0] * bmb.0 / bdb,
                -0.5 * be[1] * bmb.0 / bdb,
                -0.5 * be[2] * bmb.0 / bdb,
            ],
            e_bivector: [0.0; 3],
        };
        let eucl = Line {
            v_bivector: [
                bv[0] + 0.5 * be[0] * bmb.0 / bdb,
                bv[1] + 0.5 * be[1] * bmb.0 / bdb,
                bv[2] + 0.5 * be[2] * bmb.0 / bdb,
            ],
            e_bivector: [be[0], be[1], be[2]],
        };
        (eucl, van)
    }

    pub fn mul(&self, other: &Self) -> super::Motor {
        let e1 = self.e_bivector;
        let v1 = self.v_bivector;
        let e2 = other.e_bivector;
        let v2 = other.v_bivector;
        super::Motor {
            scalar: -e1[0] * e2[0] - e1[1] * e2[1] - e1[2] * e2[2],
            v_bivector: [
                -e1[1] * v2[2] + e1[2] * v2[1] + e2[1] * v1[2] - e2[2] * v1[1],
                e1[0] * v2[2] - e1[2] * v2[0] - e2[0] * v1[2] + e2[2] * v1[0],
                -e1[0] * v2[1] + e1[1] * v2[0] + e2[0] * v1[1] - e2[1] * v1[0],
            ],
            e_bivector: [
                -e1[1] * e2[2] + e1[2] * e2[1],
                e1[0] * e2[2] - e1[2] * e2[0],
                -e1[0] * e2[1] + e1[1] * e2[0],
            ],
            pseudo: e1[0] * v2[0]
                + e1[1] * v2[1]
                + e1[2] * v2[2]
                + e2[0] * v1[0]
                + e2[1] * v1[1]
                + e2[2] * v1[2],
        }
    }

    pub fn div(&self, other: &Self) -> super::Motor {
        self.mul(&other.inverse())
    }

    pub fn mul_scalar(&self, s: f32) -> Self {
        let we = na::Vector3::from(self.e_bivector) * s;
        let wv = na::Vector3::from(self.v_bivector) * s;
        Self {
            e_bivector: we.into(),
            v_bivector: wv.into(),
        }
    }

    pub fn div_scalar(&self, s: f32) -> Self {
        self.mul_scalar(1. / s)
    }

    pub fn add(&self, other: &Self) -> Self {
        Self {
            e_bivector: (na::Vector3::from(self.e_bivector) + na::Vector3::from(other.e_bivector))
                .into(),
            v_bivector: (na::Vector3::from(self.v_bivector) + na::Vector3::from(other.v_bivector))
                .into(),
        }
    }

    pub fn move_to(&self, dest: &Self) -> super::Motor {
        dest.div(&self).ssqrt()
    }
}

impl From<&super::Motor> for Line {
    fn from(m: &super::Motor) -> Self {
        Self {
            e_bivector: m.e_bivector,
            v_bivector: m.v_bivector,
        }
    }
}

impl PartialEq for Line {
    fn eq(&self, other: &Line) -> bool {
        self.exp() == other.exp()
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn new() {
        let l = super::Line::new(&[4., 3., 2.], &[2., 7., 6.]);
        let l_ = super::Line {
            v_bivector: [-4., 20., -22.],
            e_bivector: [2., 7., 6.],
        };
        assert_eq!(l, l_);

        let l = super::Line::new(&[5., 3., 4.], &[3., 2., 1.]);
        let l_ = super::Line {
            v_bivector: [5., -7., -1.],
            e_bivector: [3., 2., 1.],
        };
        assert_eq!(l, l_);
    }

    #[test]
    fn dual() {
        let d = super::Line {
            v_bivector: [2., 3., 4.],
            e_bivector: [8., 5., 2.],
        }
        .dual();

        let d_ = super::Line {
            v_bivector: [8., 5., 2.],
            e_bivector: [2., 3., 4.],
        };
        assert_eq!(d, d_);

        let l = super::Line::new(&[4., -2., 9.], &[-4., 6., 3.]);
        assert_eq!(l, l.dual().dual());
    }

    #[test]
    fn decompose() {
        // decomposition should add up to original result,
        // and to parts should commute.
        let l1 = Line::random();
        let (e, v) = l1.decompose();
        assert_eq!(l1, e.add(&v));
        assert_eq!(e.mul(&v), v.mul(&e));
    }

    #[test]
    fn move_to() {
        let l1 = Line::random().mul_scalar(6.5).normalize();
        let l2 = Line::random().mul_scalar(-10.7).normalize();
        let l3 = Line::random().mul_scalar(-3.).normalize();
        let m1 = l1.move_to(&l2);
        let m2 = l2.move_to(&l3);
        let p = Point::random();
        let r1 = l2.exp();
        let r2 = m1.apply_to_line(&l1).exp();
        println!("{:?}", r1);
        println!("{:?}", r2);
        assert_eq!(l2, m1.apply_to_line(&l1));
        assert_eq!(l3, m2.apply_to_line(&l2));
    }

    #[test]
    fn neg() {
        let l = super::Line {
            v_bivector: [1., 1., 1.],
            e_bivector: [1., 1., 1.],
        };
        let ln = super::Line {
            v_bivector: [-1., -1., -1.],
            e_bivector: [-1., -1., -1.],
        };
        assert_eq!(l.neg(), ln);
    }

    #[test]
    fn exp_sanity1() {
        // rotated 360 degrees around normalized line,
        // should do nothing.
        let l = super::Line {
            v_bivector: [0.0f32; 3],
            e_bivector: na::Vector3::new(9., -2.4, 1.).normalize().into(),
        };
        let angle = 2. * std::f32::consts::PI;
        let v = l.mul_scalar(-angle * 0.5).exp();
        let p = Point::new(&[3., -9., 2.8]);
        assert_eq!(p, v.apply_to_point(&p));
    }

    #[test]
    fn exp_sanity2() {
        // 1, 0, 0 rotated 180 deg around y axis should become -1, 0, 0
        // and translating it 0, 1, 1, should make -1, 1, 1
        let t = super::Line {
            v_bivector: [0., 1., 1.],
            e_bivector: [0.0; 3],
        };
        let l = super::Line {
            v_bivector: [0.0f32; 3],
            e_bivector: [0., 1., 0.],
        };
        let angle = std::f32::consts::PI;
        let v = t
            .mul_scalar(0.5)
            .exp()
            .mul(&l.mul_scalar(-angle * 0.5).exp());
        let p = Point::new(&[1., 0., 0.]);
        assert_eq!(Point::new(&[-1., 1., 1.]), v.apply_to_point(&p));
    }

    #[test]
    fn exp_sanity3() {
        // Should translate by 0, 1, 0
        let l = super::Line {
            e_bivector: [0.0f32; 3],
            v_bivector: [0., 1., 0.],
        };
        let v = l.mul_scalar(0.5).exp();
        let p = Point::new(&[1., 0., 0.]);
        assert_eq!(Point::new(&[1., 1., 0.]), v.apply_to_point(&p));
    }
}
