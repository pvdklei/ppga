pub fn points(p1: &super::Point, p2: &super::Point) -> super::Line {
    let t1 = p1.trivector;
    let t2 = p2.trivector;
    super::Line {
        v_bivector: [
            -t1[2] * t2[3] + t1[3] * t2[2],
            t1[1] * t2[3] - t1[3] * t2[1],
            -t1[1] * t2[2] + t1[2] * t2[1],
        ],
        e_bivector: [
            -t1[0] * t2[1] + t1[1] * t2[0],
            -t1[0] * t2[2] + t1[2] * t2[0],
            -t1[0] * t2[3] + t1[3] * t2[0],
        ],
    }
}

pub fn line_to_point(l: &super::Line, p: &super::Point) -> super::Plane {
    let t1 = p.trivector;
    let vb = l.v_bivector;
    let eb = l.e_bivector;
    super::Plane {
        vector: [
            -t1[1] * vb[0] - t1[2] * vb[1] - t1[3] * vb[2],
            eb[1] * t1[3] - eb[2] * t1[2] + t1[0] * vb[0],
            -eb[0] * t1[3] + eb[2] * t1[1] + t1[0] * vb[1],
            eb[0] * t1[2] - eb[1] * t1[1] + t1[0] * vb[2],
        ],
    }
}

pub fn three_points(p1: &super::Point, p2: &super::Point, p3: &super::Point) -> super::Plane {
    line_to_point(&points(p1, p2), p3)
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn sanity1() {
        let p1 = Point::new(&[3., 4., 5.]);
        let p2 = Point::new(&[8., 3., 2.]);
        let l = join::points(&p1, &p2);
        assert_eq!(l, Line::new(&[8., 3., 2.], &[8. - 3., 3. - 4., 2. - 5.]));
    }
    #[test]
    fn sanity2() {
        let l = Line::new(&[0., 0., 0.], &[1., 0., 0.]);
        let p = Point::new(&[1., 0., -1.]);
        let p = join::line_to_point(&l, &p);
        assert_eq!(p, Plane::new(0., &[0., -1., 0.]));
    }

    // The meet operation is defined as a & b = (a* ^ b*)-*.
    // So meet(&p1.dual(), &p2.dual()).dual().neg() should also work.

    #[test]
    fn dual_meet_points() {
        let p1 = Point::new(&[3., -2.3, 1.7]);
        let p2 = Point::new(&[8., -7.3, -1.7]);
        assert_eq!(
            join::points(&p1, &p2),
            meet::planes(&p1.dual(), &p2.dual()).dual().neg()
        )
    }
    #[test]
    fn dual_meet_point_line() {
        let p1 = Point::new(&[3., -2.3, 1.7]);
        let l1 = Line::new(&[3., -1.7, 3.4], &[-10., 2., 6.]);
        assert_eq!(
            join::line_to_point(&l1, &p1),
            meet::plane_with_line(&p1.dual(), &l1.dual()).dual().neg()
        )
    }
}
