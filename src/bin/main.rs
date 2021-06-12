use ppga::*;

fn main() {
    //let a = Point::new(&[2., 3., 5.]);
    //let b = Point::new(&[2., 8., 7.]);
    //let c = Point::new(&[3., -2., 1.]);
    //let a = [3., 4., 1.];
    //let b = [2., 3., 1.];
    //let c = [-4., 7., -8.];
    let a = [1., 0., 0.];
    let b = [0., 0., 1.];
    let c = [0., 1., 0.];

    //let t = a.div(&a).sqrt();
    //println!("{:#?}", t.apply_to_point(&a));

    let r = Rotor::from_base(&a, &b, &c);
    println!("{:#?}", r);

    println!("{:#?}", r.apply_to_point(&Point::x()));
}
