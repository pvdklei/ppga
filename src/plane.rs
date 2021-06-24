#[derive(Debug)]
pub struct Plane {
    pub vector: [f32; 4],
}

impl Plane {
    /// So that d*n is on the plane. Convention made by PGA4CS.
    /// Standard form: n_1x + n_2y + n_3z = d.
    pub fn new(d: f32, n: &[f32; 3]) -> Self {
        Self {
            vector: [d, n[0], n[1], n[2]],
        }
    }

    pub fn random() -> Self {
        let n = na::Vec3::new_random().normalize().into();
        Self::new(rand::random(), &n)
    }
    pub fn yz() -> Self {
        Self::new(0.0, &[1., 0., 0.])
    }
    pub fn zx() -> Self {
        Self::new(0.0, &[0., 1., 0.])
    }
    pub fn xy() -> Self {
        Self::new(0.0, &[0., 0., 1.])
    }

    pub fn dual(&self) -> super::Point {
        super::Point {
            trivector: self.vector,
        }
    }

    pub fn normalize(&self) -> Self {
        let n = self.norm();
        Self {
            vector: (na::Vec4::from(self.vector) / n).into(),
        }
    }
    /// Normalizes only the normal vector
    pub fn nnormalize(&self) -> Self {
        let n = na::Vec3::from_row_slice(&self.vector[1..=3]).normalize();
        Self {
            vector: [self.vector[0], n[0], n[1], n[2]],
        }
    }
    pub fn norm(&self) -> f32 {
        na::Vec3::from_row_slice(&self.vector[1..=3]).norm()
    }

    pub fn neg(&self) -> Self {
        Self {
            vector: (-na::Vector4::from_row_slice(&self.vector)).into(),
        }
    }

    pub fn inverse(&self) -> Self {
        let v = self.vector;
        let fac = 1. / (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
        Self {
            vector: [v[0] * fac, v[1] * fac, v[2] * fac, v[3] * fac],
        }
    }

    pub fn mul_scalar(&self, s: f32) -> Self {
        Self {
            vector: (s * na::Vec4::from(self.vector)).into(),
        }
    }

    pub fn div_scalar(&self, s: f32) -> Self {
        self.mul_scalar(1. / s)
    }

    pub fn mul(&self, other: &Self) -> super::Motor {
        let v1 = self.vector;
        let v2 = other.vector;
        super::Motor {
            scalar: v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3],
            v_bivector: [
                v1[0] * v2[1] - v1[1] * v2[0],
                v1[0] * v2[2] - v1[2] * v2[0],
                v1[0] * v2[3] - v1[3] * v2[0],
            ],
            e_bivector: [
                v1[2] * v2[3] - v1[3] * v2[2],
                -v1[1] * v2[3] + v1[3] * v2[1],
                v1[1] * v2[2] - v1[2] * v2[1],
            ],
            pseudo: 0.,
        }
    }

    pub fn div(&self, other: &Self) -> super::Motor {
        self.mul(&other.inverse())
    }

    pub fn move_to(&self, dest: &Self) -> super::Motor {
        dest.div(&self).sqrt()
    }
}

impl PartialEq for Plane {
    fn eq(&self, other: &Self) -> bool {
        !self
            .vector
            .iter()
            .zip(other.vector.iter())
            .any(|(a, b)| (a - b).abs() > 0.01)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn new() {
        let p = super::Plane::new(1., &[3., 2., 4.]);
        let p_ = super::Plane {
            vector: [1., 3., 2., 4.],
        };
        assert_eq!(p, p_);
    }

    #[test]
    fn dual() {
        let p = super::Plane {
            vector: [2., 4., 5., 6.],
        }
        .dual();
        let p_ = Point {
            trivector: [2., 4., 5., 6.],
        };
        assert_eq!(p, p_);

        let p = super::Plane::new(-5., &[4., -3., 7.]);
        assert_eq!(p, p.dual().dual());
    }

    #[test]
    fn move_to() {
        let p1 = Plane::random().mul_scalar(-4.5).nnormalize();
        let p2 = Plane::random().mul_scalar(18.).nnormalize();
        let p3 = Plane::random().mul_scalar(-3.9).nnormalize();

        let m1 = p1.move_to(&p2);
        let m2 = p2.move_to(&p3);
        assert_eq!(p2, m1.apply_to_plane(&p1));
        assert_eq!(p3, m2.apply_to_plane(&p2));
    }
}
