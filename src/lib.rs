use glam::f32::vec3;
use glam::Vec3;

pub struct Ray {
    origin: Vec3,
    direction: Vec3,
}

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

pub struct PerspCamera {
    /// Location of the center of the viewport
    pub location: Vec3,

    //    /// Rotation around x, y, and z axis' represented by `rotation`'s x, y, and z
    //    /// corrdinates respectively
    //    pub rotation: Vec3,
    /// Distance from eye to viewport
    pub focal_length: f32,

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

pub struct OrthoCamera {
    /// Location of the center of the viewport
    pub location: Vec3,

    //    /// Rotation around x, y, and z axis' represented by `rotation`'s x, y, and z
    //    /// corrdinates respectively
    //    pub rotation: Vec3,
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
                t: t,
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
