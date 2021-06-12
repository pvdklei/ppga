/// p1 ^ p1 = (ae0 + be1 + ce2 + de3) ^ (xe0 + ye1 + ze2 + we3)
///         = e01(ay - bx)
///           + e02(az - cx)
///           + e03(aw - dx)
///           + e23(cw - dz)
///           + e31(dy - bw)
///           + e12(bz - cy)
pub fn planes(p1: &super::Plane, p2: &super::Plane) -> super::Line {
    super::Line {
        v_bivector: [
            p1.vector[0] * p2.vector[1] - p1.vector[1] * p2.vector[0],
            p1.vector[0] * p2.vector[2] - p1.vector[2] * p2.vector[0],
            p1.vector[0] * p2.vector[3] - p1.vector[3] * p2.vector[0],
        ],
        e_bivector: [
            p1.vector[2] * p2.vector[3] - p1.vector[3] * p2.vector[2],
            p1.vector[3] * p2.vector[1] - p1.vector[1] * p2.vector[3],
            p1.vector[1] * p2.vector[2] - p1.vector[2] * p2.vector[1],
        ],
    }
}

/// p ^ L = L ^ p
///       = (b0e01 + b1e02 + b2e03 + b3e23 + b4e31 + b5e12)
///         ^ (v0e0 + v1e1 + v2e2 + v3e3)
///       = 0123(v1b3 + v2b4 + v3b5)
///         + e032(-v0b3 + v2b2 - v3b1)
///         + e013(-v0b4 - v1b2 + v3b0)
///         + e021(-v0b3 + v1b1 - v2b0)
pub fn plane_with_line(p: &super::Plane, l: &super::Line) -> super::Point {
    super::Point {
        trivector: [
            p.vector[1] * l.e_bivector[0]
                + p.vector[2] * l.e_bivector[1]
                + p.vector[3] * l.e_bivector[2],
            -p.vector[0] * l.e_bivector[0] + p.vector[2] * l.v_bivector[2]
                - p.vector[3] * l.v_bivector[1],
            -p.vector[0] * l.e_bivector[1] - p.vector[1] * l.v_bivector[2]
                + p.vector[3] * l.v_bivector[0],
            -p.vector[0] * l.e_bivector[2] + p.vector[1] * l.v_bivector[1]
                - p.vector[2] * l.v_bivector[0],
        ],
    }
}

pub fn lines(b1: &super::Line, b2: &super::Line) -> super::PseudoScalar {
    let e1 = b1.e_bivector;
    let e2 = b2.e_bivector;
    let v1 = b1.v_bivector;
    let v2 = b2.v_bivector;
    super::PseudoScalar(
        e1[0] * v2[0]
            + e1[1] * v2[1]
            + e1[2] * v2[2]
            + e2[0] * v1[0]
            + e2[1] * v1[1]
            + e2[2] * v1[2],
    )
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn sanity1() {
        let pxy = Plane::new(1., &[0., 0., 1.]);
        let pxz = Plane::new(1., &[0., 1., 0.]);
        let l = Line::new(&[0., 1., 1.], &[-1., 0., 0.]);
        assert_eq!(meet::planes(&pxy, &pxz), l);
    }

    #[test]
    fn sanity2() {
        let p1 = Plane::new(1., &[1., 0., 0.]);
        let p2 = Plane::new(1., &[0., 1., 0.]);
        let p3 = Plane::new(1., &[0., 0., 1.]);
        assert_eq!(
            meet::plane_with_line(&p1, &meet::planes(&p2, &p3)),
            Point::new(&[1., 1., 1.])
        );
    }

    #[test]
    fn sanity3() {
        let p = Plane::new(4., &[0., 1., 0.]);
        let l = Line::new(&[4., 0., 1.], &[0., 1., 0.]);
        assert_eq!(meet::plane_with_line(&p, &l), Point::new(&[4., 4., 1.]));
    }
}
