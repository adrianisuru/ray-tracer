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
