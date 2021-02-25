use glam::f32::vec3;
use glam::{Mat3, Vec3};
use rand::rngs::ThreadRng;
use rand::Rng;

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    origin: Vec3,
    direction: Vec3,
}

/// Represents a ray tracing ray
impl Ray {
    pub fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray {
            origin,
            direction: direction.normalize(),
        }
    }
    pub fn origin(&self) -> Vec3 {
        self.origin
    }
    pub fn direction(&self) -> Vec3 {
        self.direction
    }
    /// Computes the position `P(t)` of this ray at time `t`, according to the
    /// function `P(t) = origin + direction * t`.
    /// ```rust
    ///use glam::Vec3;
    ///use ray_tracer::Ray;
    ///
    ///let origin = Vec3::zero();
    ///let direction = Vec3::new(1., 0., 0.);
    ///
    ///let r = Ray::new(origin, direction);
    ///
    ///let p_1 = r.at(0.5);
    ///
    ///assert_eq!(Vec3::new(0.5, 0., 0.), p_1);
    /// ```
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin + t * self.direction
    }
}

pub trait Origin {
    fn ray_through(&self, through: Vec3) -> Ray;
}

impl Origin for Vec3 {
    fn ray_through(&self, through: Vec3) -> Ray {
        Ray::new(*self, through - *self)
    }
}

pub struct Sphere {
    pub center: Vec3,
    pub radius: f32,
}

impl Surface for Sphere {
    /// Checks if ray intersects this sphere
    fn hit(&self, ray: &Ray) -> Option<HitRecord> {
        let oc = ray.origin() - self.center;
        let a = ray.direction().dot(ray.direction());
        let b = 2. * oc.dot(ray.direction());
        let c = oc.dot(oc) - self.radius.powi(2);
        let discriminant = b * b - 4. * a * c;

        if discriminant < 0. {
            None
        } else {
            let t = (-b - discriminant.sqrt()) / (2. * a);
            Some(HitRecord {
                t,
                p: ray.at(t),
                n: SurfaceNormal::new(ray.direction, ray.at(t) - self.center),
            })
        }
    }
    fn color(&self) -> Vec3 {
        vec3(1., 0., 0.)
    }
}

/// A camera which uses perspective projection
pub struct PerspCamera {
    /// Location of the center of the viewplane
    pub location: Vec3,

    //    /// Rotation around x, y, and z axis' represented by `rotation`'s x, y, and z
    //    /// corrdinates respectively
    //    pub rotation: Vec3,
    /// Distance from eye to viewplane
    pub focal_length: f32,

    /// width of the viewplane divided by height of the viewplane
    pub aspect_ratio: f32,

    pub zoom: f32,
}

impl Camera for PerspCamera {
    fn ray_through(&self, u: f32, v: f32) -> Ray {
        let origin =
            self.location + self.zoom * self.focal_length * Vec3::unit_z();
        origin.ray_through(
            self.zoom * vec3(u * self.aspect_ratio, v, 0.) + self.location,
        )
    }
}

pub trait Camera {
    fn ray_through(&self, u: f32, v: f32) -> Ray;
}

/// A camera which uses orthographic projection
pub struct OrthoCamera {
    /// Location of the center of the viewplane
    pub location: Vec3,

    //    /// Rotation around x, y, and z axis' represented by `rotation`'s x, y, and z
    //    /// corrdinates respectively
    //    pub rotation: Vec3,
    /// width of the viewplane divided by height of the viewplane
    pub aspect_ratio: f32,

    pub zoom: f32,
}

impl Camera for OrthoCamera {
    fn ray_through(&self, u: f32, v: f32) -> Ray {
        let origin =
            self.zoom * vec3(u * self.aspect_ratio, v, 0.) + self.location;
        let direction = -Vec3::unit_z();
        Ray::new(origin, direction)
    }
}

pub struct Plane {
    /// A point on this plane
    pub point: Vec3,
    /// A vector orthoganal to this plane
    pub normal: Vec3,
}

use SurfaceNormal::*;

impl Surface for Plane {
    fn hit(&self, ray: &Ray) -> Option<HitRecord> {
        let dn = ray.direction().dot(self.normal);
        if dn == 0. {
            None
        } else {
            let t = (self.point - ray.origin).dot(self.normal) / dn;

            Some(HitRecord {
                t,
                p: ray.at(t),
                n: SurfaceNormal::new(ray.direction(), self.normal),
            })
        }
    }
    fn color(&self) -> Vec3 {
        vec3(0.4, 0.4, 0.4)
    }
}

/// A normal vector which includes info about which side of the surface it is pointing
pub enum SurfaceNormal {
    Outer(Vec3),
    Inner(Vec3),
}

impl SurfaceNormal {
    fn new(direction: Vec3, normal: Vec3) -> SurfaceNormal {
        if normal.dot(direction) < 0. {
            Outer(normal)
        } else {
            Inner(normal)
        }
    }
}

pub struct HitRecord {
    /// The time at which the surface was hit by the ray
    pub t: f32,
    /// The point at which the surface was hit by the ray
    pub p: Vec3,

    /// The normal vector to the surface at p
    pub n: SurfaceNormal,
}

pub trait Surface {
    fn hit(&self, ray: &Ray) -> Option<HitRecord>;
    fn color(&self) -> Vec3;
}

type Triangle = (Vec3, Vec3, Vec3);

impl Surface for Triangle {
    fn hit(&self, ray: &Ray) -> Option<HitRecord> {
        let &(p0, p1, p2) = self;
        let d = ray.direction();
        let o = ray.origin();

        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let n = e1.cross(e2);

        ray_tri_intersect_mt(o, d, p0, p1, p2).map(|t| HitRecord {
            t,
            p: ray.at(t),
            n: Outer(n),
        })
    }

    fn color(&self) -> Vec3 {
        Vec3::new(0., 1., 0.)
    }
}

fn ray_tri_intersect_mt(
    o: Vec3,
    d: Vec3,
    p0: Vec3,
    p1: Vec3,
    p2: Vec3,
) -> Option<f32> {
    let e1 = p1 - p0;
    let e2 = p2 - p0;

    let q = d.cross(e2);
    let a = e1.dot(q);
    let eps = 10f32.powi(-5);
    if -eps < a && a < eps {
        None
    } else {
        let f = 1. / a;
        let s = o - p0;
        let u = f * s.dot(q);
        let r = s.cross(e1);
        let v = f * d.dot(r);
        let t = f * e2.dot(r);
        if u < 0. || v < 0. || u + v > 1. {
            None
        } else {
            Some(t)
        }
    }
}

fn ray_tri_intersect(
    o: Vec3,
    d: Vec3,
    p0: Vec3,
    p1: Vec3,
    p2: Vec3,
) -> Option<f32> {
    let e1 = p1 - p0;
    let e2 = p2 - p0;
    let s = o - p0;

    let m = Mat3::from_cols(-d, e1, e2);
    let eps = 10f32.powi(-5);
    let a = m.determinant();

    if -eps < a && a < eps {
        None
    } else {
        let t = Mat3::from_cols(s, e1, e2).determinant();
        let u = Mat3::from_cols(-d, s, e2).determinant();
        let v = Mat3::from_cols(-d, e1, s).determinant();
        if u < 0. || v < 0. || u + v > 1. {
            None
        } else {
            Some(t)
        }
    }
}

/// Returns a random vector of length <= 1 using an rng
pub fn random_unit_sphere(rng: &mut ThreadRng) -> Vec3 {
    random_unit_vector(rng) * rng.gen_range(0. ..1.)
}

/// Returns a random vector of length 1 using an rng
pub fn random_unit_vector(rng: &mut ThreadRng) -> Vec3 {
    let resolution = 100.;
    let x = rng.gen_range(-resolution..resolution);
    let y = rng.gen_range(-resolution..resolution);
    let z = rng.gen_range(-resolution..resolution);
    Vec3::new(x, y, z).normalize()
}
