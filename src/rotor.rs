#[derive(Debug, Copy, Clone)]
pub struct Rotor {
    pub scalar: f32,
    pub e_bivector: [f32; 3],
}

impl Rotor {
    pub fn new(a: f32, axis: &[f32; 3]) -> Self {
        let ha = 0.5 * a;
        let sha = ha.sin();
        Self {
            scalar: ha.cos(),
            e_bivector: [-sha * axis[0], -sha * axis[1], -sha * axis[2]],
        }
    }

    pub fn random() -> Self {
        Self::new(
            rand::random::<f32>() * std::f32::consts::PI * 2.0,
            &na::Vec3::new_random().into(),
        )
    }

    /// Creates a rotor out of a base transformation (e.g., matrix columns).
    /// Note that the base vectors must be normalized and orthogonal to each other.
    /// If not, this method will not panic, but returns an invalid rotor.
    pub fn from_base(e1: &[f32; 3], e2: &[f32; 3], e3: &[f32; 3]) -> Self {
        // let e1_ = super::Point::new(e1);
        // let e2_ = super::Point::new(e2);
        // let e3_ = super::Point::new(e3);
        // let e1 = super::Point::x();
        // let e2 = super::Point::y();
        // let e3 = super::Point::z();
        let e1_ = super::Plane::new(0.0, e1).normalize();
        let e2_ = super::Plane::new(0.0, e2).normalize();
        let e3_ = super::Plane::new(0.0, e3).normalize();
        let e1 = super::Plane::yz().normalize();
        let e2 = super::Plane::zx().normalize();
        let e3 = super::Plane::xy().normalize();

        // TODO: Make something simpler than three point corrs, like two line corrs
        //let from = super::join::three_points(&e1, &e2, &e3);
        //let to = super::join::three_points(&e1_, &e2_, &e3_);
        //let r = to.div(&from).sqrt().into_rotor_unchecked();

        let m = super::Motor::from_plane_correspondences(&e1, &e1_, &e2, &e2_, &e3, &e3_);
        let r = m.into_rotor_unchecked();
        r
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
            e_bivector: {
                let e = na::Vector3::from(self.e_bivector);
                (e * fac).into()
            },
        }
    }

    pub fn sqrt(&self) -> Self {
        let s = self.scalar;
        let e = self.e_bivector;
        let fac = 2.0f32.sqrt() / (2.0 * (s + 1.).sqrt());
        Self {
            scalar: 0.5 * (2. * s + 2.).sqrt(),
            e_bivector: [e[0] * fac, e[1] * fac, e[2] * fac],
        }
    }

    pub fn neg(&self) -> Self {
        Self {
            scalar: -self.scalar,
            e_bivector: (-na::Vec3::from(self.e_bivector)).into(),
        }
    }

    // Course notes chapter 8
    pub fn ln(&self) -> super::Line {
        let e = na::Vec3::from(self.e_bivector);
        let s2 = e.norm_squared().sqrt();
        let u = s2.atan2(self.scalar);
        super::Line {
            v_bivector: [0.0; 3],
            e_bivector: (e * u / s2).into(),
        }
    }
    pub fn outer_ln(&self) -> super::Line {
        super::Line {
            e_bivector: (na::Vec3::from(self.e_bivector) / self.scalar).into(),
            v_bivector: [0.; 3],
        }
    }
    pub fn qtangent_ln(&self) -> [f32; 3] {
        if self.scalar.is_sign_negative() {
            (-na::Vec3::from(self.e_bivector)).into()
        } else {
            self.e_bivector
        }
    }
    pub fn qtangent_exp(q: &[f32; 3]) -> Self {
        let q = na::Vec3::from(*q);
        Self {
            e_bivector: q.into(),
            scalar: (1. - q.dot(&q)).sqrt(),
        }
    }

    pub fn apply_to_point(&self, p: &super::Point) -> super::Point {
        let p = p.trivector;
        let e = self.e_bivector;
        let s = self.scalar;

        super::Point {
            trivector: [
                e[0] * e[0] * p[0] + e[1] * e[1] * p[0] + e[2] * e[2] * p[0] + p[0] * s * s,
                e[0] * (e[0] * p[1] + e[1] * p[2] + e[2] * p[3])
                    - e[1] * (-e[0] * p[2] + e[1] * p[1] + p[3] * s)
                    + e[2] * (e[0] * p[3] - e[2] * p[1] + p[2] * s)
                    + s * (-e[1] * p[3] + e[2] * p[2] + p[1] * s),
                e[0] * (-e[0] * p[2] + e[1] * p[1] + p[3] * s)
                    + e[1] * (e[0] * p[1] + e[1] * p[2] + e[2] * p[3])
                    - e[2] * (-e[1] * p[3] + e[2] * p[2] + p[1] * s)
                    + s * (e[0] * p[3] - e[2] * p[1] + p[2] * s),
                -e[0] * (e[0] * p[3] - e[2] * p[1] + p[2] * s)
                    + e[1] * (-e[1] * p[3] + e[2] * p[2] + p[1] * s)
                    + e[2] * (e[0] * p[1] + e[1] * p[2] + e[2] * p[3])
                    + s * (-e[0] * p[2] + e[1] * p[1] + p[3] * s),
            ],
        }
    }

    pub fn mul_translator(&self, t: &super::Translator) -> super::Motor {
        let ts = t.scalar;
        let tv = t.v_bivector;
        let rs = self.scalar;
        let re = self.e_bivector;
        super::Motor {
            scalar: rs * ts,
            pseudo: re[0] * tv[0] + re[1] * tv[1] + re[2] * tv[2],
            e_bivector: [re[0] * ts, re[1] * ts, re[2] * ts],
            v_bivector: [
                -re[1] * tv[2] + re[2] * tv[1] + rs * tv[0],
                re[0] * tv[2] - re[2] * tv[0] + rs * tv[1],
                -re[0] * tv[1] + re[1] * tv[0] + rs * tv[2],
            ],
        }
    }
}

impl From<&Rotor> for [f32; 4] {
    fn from(r: &Rotor) -> [f32; 4] {
        [r.scalar, r.e_bivector[0], r.e_bivector[1], r.e_bivector[2]]
    }
}
impl From<Rotor> for [f32; 4] {
    fn from(r: Rotor) -> [f32; 4] {
        [r.scalar, r.e_bivector[0], r.e_bivector[1], r.e_bivector[2]]
    }
}

impl PartialEq for Rotor {
    fn eq(&self, other: &Self) -> bool {
        (self.scalar - other.scalar).abs() < 0.03
            && self
                .e_bivector
                .iter()
                .zip(other.e_bivector.iter())
                .all(|(a, b)| (a - b).abs() < 0.3)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn sanity1() {
        let p = Point::new(&[3., 4., 5.]);
        let rot = Rotor::new(2. * std::f32::consts::PI, &[4., -3., 1.3]);
        assert_eq!(rot.apply_to_point(&p), p);
    }

    #[test]
    fn base() {
        let r = Rotor::from_base(&[0., 1., 0.], &[1., 0., 0.], &[0., 0., 1.]);
        let p = Point::random();
        println!("{:?}", r);
        println!("{:?}", p.eucl());
        let p_ = r.apply_to_point(&p);
        let mut p = p.eucl();
        p.swap(0, 1);
        assert_eq!(p, p_.eucl())
    }

    #[test]
    fn for_steven() {
        let r = Rotor::from_base(&[-0.95, 0.0, 0.31], &[0.0, 1.0, 0.0], &[-0.31, 0.0, -0.95])
            .normalize();
        let pos = [0.7, 3.0, -2.6];
        let t = Translator::new(&pos);
        let m = t.mul_rotor(&r).normalize();

        let m_ = m.ln().exp().normalize();
        println!("{:?}", m);
        println!("{:?}", m.ln());
        println!("{:?}", m_);
        let p = Point::random();
        println!("True Point After Motor {:?}", m.apply_to(&p));
        println!("Point After Log Motor {:?}", m_.apply_to(&p));
        assert_eq!(m, m_);
    }

    #[test]
    fn ln() {
        let r = Rotor::random().normalize();
        let m = Motor::from(&r);
        let m_ = r.ln().exp();
        println!("{:?}", m);
        println!("{:?}", m_);
        assert!(m.is_similar_to(0.01, &m_));
    }
    #[test]
    fn outer_ln() {
        let r = Rotor::random().normalize();
        let m = Motor::from(&r);
        let m_ = r.outer_ln().outer_exp();
        println!("{:?}", m);
        println!("{:?}", m_);
        assert!(m.is_similar_to(0.01, &m_));
    }

    #[test]
    fn qtangent() {
        let r = Rotor::random().normalize();
        let r_ = Rotor::qtangent_exp(&r.qtangent_ln());
        assert!(r == r_ || r == r_.neg());
    }
}
