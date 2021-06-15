use ppga::*;

fn main() {
    let p1 = Point::new(&[1., -10., 3.]);
    let p2 = Point::new(&[5., -2., -8.]);
    let p3 = Point::new(&[2., -6., 0.]);
    let p4 = Point::new(&[-5., 2., -7.]);
    let l1 = join::points(&p1, &p2);
    let l2 = join::points(&p3, &p4);
    let m = l2.div(&l1).sqrt();
    let m_ = m.cayley_ln().cayley_exp();
    let p = Point::random();
    println!("{:?}", m);
    println!("{:?}", m.cayley_ln());
    println!("{:?}", m_);
    println!("{:?}", m.is_similar_to(0.3, &m_));

    println!("{:?}", m.apply_to(&p));
    println!("{:?}", m_.apply_to(&p));
}
