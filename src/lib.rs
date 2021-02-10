use glam::f32::vec3;
use glam::Vec3;

pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
}

impl Ray {
    /// Computes the position `P(t)` of this ray at time `t`, according to the
    /// function `P(t) = origin + direction * t`.
    /// ```rust
    ///use glam::Vec3;
    ///use ray_tracer::Ray;
    ///
    ///let origin = Vec3::zero();
    ///let direction = Vec3::new(1., 2., 3.);
    ///
    ///let r = Ray { origin, direction };
    ///
    ///let p_1 = r.at(0.5);
    ///
    ///assert_eq!(Vec3::new(0.5, 1., 1.5), p_1);
    /// ```
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin + t * self.direction
    }
}

pub struct Sphere {
    pub center: Vec3,
    pub radius: f32,
}

impl Sphere {
    /// Checks if ray intersects this sphere
    pub fn intersects(&self, ray: &Ray) -> f32 {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(ray.direction);
        let b = 2. * oc.dot(ray.direction);
        let c = oc.dot(oc) - self.radius.powi(2);
        let discriminant = b * b - 4. * a * c;
        if discriminant < 0. {
            -1.
        } else {
            (-b - discriminant.sqrt()) / (2. * a)
        }
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
        let direction = self.zoom * vec3(u * self.aspect_ratio, v, 0.)
            + self.location
            - origin;
        Ray { origin, direction }
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
        Ray { origin, direction }
    }
}
