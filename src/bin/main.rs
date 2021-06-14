use ppga::*;

type f32 = f64;

fn main() {
    let p1 = Plane::new(1., &[1., 0., 0.]);
    let p2 = Plane::new(6., &na::Vec3::new(0., 0., 1.).into());
    let p3 = Plane::new(-4., &na::Vec3::new(0., 1., 0.).into());
    let m = p2.div(&p1).ssqrt().mul(&p3.div(&p2).ssqrt());
    let p = Point::new(&[2., 3., 4.]);
    println!("{:?}", m);
    println!("{:?}", m.ln().exp());
    // println!("{:?}", m.apply_to_point(&p));
    // println!("{:?}", m_.apply_to_point(&p));
    // println!("{:?}", m_.ln().exp().apply_to_point(&p));
}
