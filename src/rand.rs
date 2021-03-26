use glam::{vec3, vec4, Vec3, Vec4};
use rand::rngs::ThreadRng;
use rand::Rng;

use crate::AABB;

/// Returns a random vector of length <= 1 using an rng
pub fn random_unit_sphere(rng: &mut ThreadRng) -> Vec3 {
    random_unit_vector(rng) * rng.gen_range(0. ..1.)
}

/// Returns a random vector of length 1 using an rng
pub fn random_unit_vector(rng: &mut ThreadRng) -> Vec3 {
    let x = rng.gen_range(f32::MIN..f32::MAX);
    let y = rng.gen_range(f32::MIN..f32::MAX);
    let z = rng.gen_range(f32::MIN..f32::MAX);
    vec3(x, y, z).normalize()
}

pub fn random_in_box(rng: &mut ThreadRng, bounds: AABB) -> Vec4 {
    let AABB { min, max } = bounds;
    let x = rng.gen_range(min.x..max.x);
    let y = rng.gen_range(min.y..max.y);
    let z = rng.gen_range(min.z..max.z);
    vec4(x, y, z, 1.)
}
