pub fn lines(l1: &super::Line, l2: &super::Line) -> f32 {
    let e1 = l1.e_bivector;
    let e2 = l2.e_bivector;
    -e1[0] * e2[0] - e1[1] * e2[1] - e1[2] * e2[2]
}
