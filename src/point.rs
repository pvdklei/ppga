#[derive(Debug)]
pub struct Point {
    pub trivector: [f32; 4],
}

impl Point {
    /// Simply (x, y, z) but on basis { e023, e031, e012 },
    /// no neg since our base is diffent, and one times e123.
    /// Based on PGA4CS.
    pub fn new([x, y, z]: &[f32; 3]) -> Self {
        Self {
            trivector: [1.0, -x, -y, -z],
        }
    }

    /// Directions, or points at infinity, have zero for e123.
    pub fn inf([x, y, z]: &[f32; 3]) -> Self {
        Self {
            trivector: [0.0, -x, -y, -z],
        }
    }

    pub fn x() -> Self {
        Self::new(&[1., 0., 0.])
    }
    pub fn y() -> Self {
        Self::new(&[0., 1., 0.])
    }
    pub fn z() -> Self {
        Self::new(&[0., 0., 1.])
    }

    pub fn r3(&self) -> [f32; 3] {
        [-self.trivector[1], -self.trivector[2], -self.trivector[3]]
    }

    pub fn random() -> Self {
        use rand::prelude::*;
        let mut rng = rand::thread_rng();
        Self {
            trivector: [1., rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()],
        }
    }

    pub fn random_dir() -> Self {
        let mut x = Self::random();
        x.trivector[0] = 0.0;
        x
    }

    pub fn dual(&self) -> super::Plane {
        super::Plane {
            vector: self.trivector,
        }
    }

    pub fn neg(&self) -> Self {
        Self {
            trivector: (-na::Vector4::from_row_slice(&self.trivector)).into(),
        }
    }

    pub fn norm(&self) -> f32 {
        self.trivector[0]
    }
    pub fn normalize(&self) -> Self {
        Self {
            trivector: (na::Vec4::from(self.trivector) / self.norm()).into(),
        }
    }

    pub fn is_inf(&self) -> bool {
        self.trivector[0].abs() < 0.01
    }

    pub fn inverse(&self) -> Self {
        let p = self.trivector;
        let fac = 1. / (p[0] * p[0]);
        Self {
            trivector: [-1. / p[0], -p[1] * fac, -p[2] * fac, -p[3] * fac],
        }
    }

    pub fn mul(&self, other: &Self) -> super::Translator {
        let p1 = self.trivector;
        let p2 = other.trivector;
        super::Translator {
            scalar: -p1[0] * p2[0],
            v_bivector: [
                -p1[0] * p2[1] + p1[1] * p2[0],
                -p1[0] * p2[2] + p1[2] * p2[0],
                -p1[0] * p2[3] + p1[3] * p2[0],
            ],
        }
    }

    pub fn div(&self, other: &Self) -> super::Translator {
        self.mul(&other.inverse())
    }

    pub fn move_to(&self, dest: &Self) -> super::Translator {
        dest.div(&self).sqrt()
    }
    pub fn is_similar_to(&self, d: f32, other: &Self) -> bool {
        !self
            .trivector
            .iter()
            .zip(other.trivector.iter())
            .any(|(a, b)| (a - b).abs() > d)
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.is_similar_to(0.01, other)
    }
}

impl From<[f32; 3]> for Point {
    fn from(x: [f32; 3]) -> Self {
        Self::new(&x)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn new() {
        let p = super::Point::new(&[4., 3., 9.]);
        let p_ = super::Point {
            trivector: [1., -4., -3., -9.],
        };
        assert_eq!(p, p_);
    }

    #[test]
    fn dual() {
        let p = super::Point {
            trivector: [5., 2., 3., 1.],
        }
        .dual();
        let p_ = Plane {
            vector: [5., 2., 3., 1.],
        };
        assert_eq!(p, p_);

        let p = super::Point::new(&[4., -2., 7.]);
        assert_eq!(p, p.dual().dual());
    }

    #[test]
    fn move_to() {
        let p1 = Point::random().normalize();
        let p2 = Point::random().normalize();
        let p3 = Point::random().normalize();
        let m1 = p1.move_to(&p2);
        let m2 = p2.move_to(&p3);
        assert_eq!(p2, m1.apply_to_point(&p1));
        assert_eq!(p3, m2.apply_to_point(&p2));
    }
}
