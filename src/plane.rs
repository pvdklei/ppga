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

    pub fn dual(&self) -> super::Point {
        super::Point {
            trivector: self.vector,
        }
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
    use crate::Point;
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
}
