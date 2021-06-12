#[derive(Debug, PartialEq)]
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

    pub fn exp(&self) -> super::Motor {
        let mut l = na::Vector3::from(self.e_bivector);
        let ao2 = -l.norm();
        l = l.normalize();

        let exp_v_v = self.v_bivector; // and its scalar component is 1
        let exp_e_s = ao2.cos();
        let exp_e_e = -ao2.sin() * l;

        super::Motor {
            scalar: exp_e_s,
            pseudo: exp_e_e[0] * exp_v_v[0] + exp_e_e[1] * exp_v_v[1] + exp_e_e[2] * exp_v_v[2],
            v_bivector: [
                exp_e_e[1] * exp_v_v[2] - exp_e_e[2] * exp_v_v[1] + exp_e_s * exp_v_v[0],
                -exp_e_e[0] * exp_v_v[2] + exp_e_e[2] * exp_v_v[0] + exp_e_s * exp_v_v[1],
                exp_e_e[0] * exp_v_v[1] - exp_e_e[1] * exp_v_v[0] + exp_e_s * exp_v_v[2],
            ],
            e_bivector: exp_e_e.into(),
        }
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
}

impl From<&super::Motor> for Line {
    fn from(m: &super::Motor) -> Self {
        Self {
            e_bivector: m.e_bivector,
            v_bivector: m.v_bivector,
        }
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
        let l = super::Line {
            v_bivector: [0.0f32; 3],
            e_bivector: [0., 1., 0.],
        };
        let angle = std::f32::consts::PI;
        let v = l.mul_scalar(-angle * 0.5).exp();
        let p = Point::new(&[1., 0., 0.]);
        assert_eq!(Point::new(&[-1., 0., 0.]), v.apply_to_point(&p));
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
