//! This crate contains an implementation of plane-based (projective)
//! geometric algebra, i.e., G(3,0,1). It uses the memory model and
//! conventions described in Dorst 2020 (https://bivector.net/PGA4CS.html).
//!
//! Most calculations are generated using the python version of enki's
//! ganja.js (https://github.com/enkimute/ganja.js/) and sympy. See the
//! generate folder for details.
//!
//! Memory model:
//!     scalar = 1
//!     pseudo = e0123
//!     vector = { e0, e1, e2, e3 }
//!     v_bivector = { e01, e02, e03 }
//!     e_bivector = { e23, e31, e12 }
//!     trivector = { e123, e032, e013, e021 }

mod line;
mod motor;
mod plane;
mod point;
mod rotor;
mod translator;

mod error;

pub mod inner;
pub mod join;
pub mod meet;

pub use line::Line;
pub use motor::Motor;
pub use plane::Plane;
pub use point::Point;
pub use rotor::Rotor;
pub use translator::Translator;

#[allow(non_upper_case_globals)]
pub const vector: [&'static str; 4] = ["e0", "e1", "e2", "e3"];
#[allow(non_upper_case_globals)]
pub const v_bivector: [&'static str; 3] = ["e01", "e02", "e03"];
#[allow(non_upper_case_globals)]
pub const e_bivector: [&'static str; 3] = ["e23", "e31", "e12"];
#[allow(non_upper_case_globals)]
pub const trivector: [&'static str; 4] = ["e123", "e032", "e013", "e021"];

#[derive(Debug, Copy, Clone)]
pub struct PseudoScalar(pub f32);

impl PseudoScalar {
    pub fn mul_scalar(&self, s: f32) -> Self {
        Self(self.0 * s)
    }
    pub fn div_scalar(&self, s: f32) -> Self {
        self.mul_scalar(1. / s)
    }
}
