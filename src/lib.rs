#![feature(total_cmp)]
use glam::f32::vec2;
use glam::f32::vec3;
use glam::{Mat3, Mat4, Vec2, Vec3, Vec4};
use rand::rngs::ThreadRng;
use rand::Rng;
use std::f32;
use std::iter::FromIterator;
use std::rc::Rc;

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    origin: Vec4,
    direction: Vec4,
}

/// Represents a ray tracing ray
impl Ray {
    pub fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray {
            origin: (origin, 1.).into(),
            direction: (direction, 0.).into(),
        }
    }
    pub fn origin(&self) -> Vec4 {
        self.origin
    }
    pub fn direction(&self) -> Vec4 {
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
    pub fn at(&self, t: f32) -> Vec4 {
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
    pub center: Vec4,
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
                c: vec3(1., 0., 0.),
            })
        }
    }

    fn aabb(&self) -> Option<AABB> {
        Some(AABB {
            min: self.center - (Vec3::splat(self.radius), 0.).into(),
            max: self.center + (Vec3::splat(self.radius), 0.).into(),
        })
    }
}

pub struct Viewplane {
    pub hres: u32,
    pub vres: u32,
    pub s: f32,
}

impl Viewplane {
    /// Returns xy coordinates in camera space from row column and pixel offsets
    pub fn view(&self, c: u32, r: u32, px: f32, py: f32) -> Vec2 {
        let s = self.s as f32;
        let hres = self.hres as f32;
        let vres = self.vres as f32;
        let i = vec2(c as f32, vres - (r as f32));
        let res = vec2(hres, vres);
        let p = vec2(px, py);

        let c = c as f32;
        let r = r as f32;

        vec2(
            s * (c - (hres / 2.) + px),
            s * (vres - r - (vres / 2.) + py),
        )
    }
}

/// A camera which uses perspective projection
pub struct PerspCamera {
    /// Origin of all rays generated with this camera, origin of this camera's
    /// local euclidian frame
    eye: Vec3,
    /// Distance between eye and the viewplane
    distance: f32,
    /// Orthonormal basis of this camera's local euclidian frame
    onb: Mat3,
}

impl PerspCamera {
    pub fn new(eye: Vec3, lookat: Vec3, up: Vec3) -> PerspCamera {
        let w = (eye - lookat).normalize();
        let u = up.cross(w).normalize();
        let v = w.cross(u).normalize();

        PerspCamera {
            eye,
            distance: (eye - lookat).length(),
            onb: Mat3::from_cols(u, v, w),
        }
    }
}

impl Camera for PerspCamera {
    fn ray_through(&self, v: Vec2) -> Ray {
        let eye = self.eye;
        let onb = self.onb;
        let d = self.distance;
        let c: Vec3 = (v, -d).into();

        Ray::new(eye, onb * c + eye)
    }
}

pub trait Camera {
    fn ray_through(&self, v: Vec2) -> Ray;
}

// TODO this camera is broken
/// A camera which uses orthographic projection
pub struct OrthoCamera {
    /// Location of the center of the viewplane
    pub eye: Vec3,
    onb: Mat3,
}

impl OrthoCamera {
    pub fn new(eye: Vec3, lookat: Vec3, up: Vec3) -> OrthoCamera {
        let w = (eye - lookat).normalize();
        let u = up.cross(w).normalize();
        let v = w.cross(u).normalize();

        OrthoCamera {
            eye,
            onb: Mat3::from_cols(u, v, w),
        }
    }
}

impl Camera for OrthoCamera {
    fn ray_through(&self, v: Vec2) -> Ray {
        let onb = self.onb;
        let o: Vec3 = (v, 0.).into();
        Ray::new(onb * o, onb * (Vec3::Z))
    }
}

pub struct Plane {
    /// A point on this plane
    pub point: Vec4,
    /// A vector orthoganal to this plane
    pub normal: Vec4,
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
                c: vec3(0.3, 0.3, 0.3),
            })
        }
    }

    fn aabb(&self) -> Option<AABB> {
        None
    }
}

/// A normal vector which includes info about which side of the surface it is pointing
#[derive(Debug)]
pub enum SurfaceNormal {
    Outer(Vec4),
    Inner(Vec4),
}

impl SurfaceNormal {
    fn new(direction: Vec4, normal: Vec4) -> SurfaceNormal {
        if normal.dot(direction) < 0. {
            Outer(normal)
        } else {
            Inner(normal)
        }
    }
}

#[derive(Debug)]
pub struct HitRecord {
    /// The time at which the surface was hit by the ray
    pub t: f32,
    /// The point at which the surface was hit by the ray
    pub p: Vec4,

    /// The normal vector to the surface at p
    pub n: SurfaceNormal,

    /// The color of the surface at p
    pub c: Vec3,
}

pub trait Surface {
    /// Returns information about where ray intersects this surface (or None)
    fn hit(&self, ray: &Ray) -> Option<HitRecord>;

    /// Computes an axis aligned bounding box for this surface
    fn aabb(&self) -> Option<AABB>;
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

        ray_tri_intersect_mt(o.into(), d.into(), p0, p1, p2).map(|t| {
            HitRecord {
                t,
                p: ray.at(t),
                n: Outer((n, 0.).into()),
                c: vec3(0., 0., 1.),
            }
        })
    }

    fn aabb(&self) -> Option<AABB> {
        todo!()
    }
}

/// Möller-Trumbore ray triangle intersection
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

/// Axis aligned bounding box
#[derive(Copy, Clone, Debug)]
pub struct AABB {
    /// Lower left back of this bounding box
    min: Vec4,

    /// Upper right front of this bounding box
    max: Vec4,
}

impl AABB {
    fn hit(&self, ray: &Ray) -> bool {
        let d = Mat3::from_diagonal(ray.direction().into());
        let inv = d.inverse();
        let b_min: Vec3 = (self.min - ray.origin()).into();
        let b_max: Vec3 = (self.max - ray.origin()).into();
        let t_min: Vec3 = d.inverse() * b_min;
        let t_max: Vec3 = d.inverse() * b_max;

        let t0 = vec3(
            t_min.x.min(t_max.x),
            t_min.y.min(t_max.y),
            t_min.z.min(t_max.z),
        );

        let t = vec3(
            t_min.x.max(t_max.x),
            t_min.y.max(t_max.y),
            t_min.z.max(t_max.z),
        );

        t0.max_element() < t.min_element()
    }

    fn center(self) -> Vec4 {
        (self.max - self.min) * 0.5 + self.min
    }
}

/// Node of a Bounding Volume Hierarchy.
#[derive(Clone)]
pub enum BVH {
    Node(AABB, Box<BVH>, Box<BVH>),
    Leaf(Rc<dyn Surface>),
}

impl Surface for BVH {
    fn hit(&self, ray: &Ray) -> Option<HitRecord> {
        match self {
            Self::Leaf(surface) => surface.hit(&ray),
            Self::Node(bv, left, right) => {
                if bv.hit(&ray) {
                    let hit = [left, right]
                        .iter()
                        .filter_map(|s| {
                            s.aabb().and_then(|bb| {
                                if bb.hit(&ray) {
                                    s.hit(&ray)
                                } else {
                                    None
                                }
                            })
                        })
                        .min_by(|hit1, hit2| hit1.t.total_cmp(&hit2.t));
                    println!("{:?}", hit);
                    hit
                } else {
                    None
                }
            }
        }
    }

    fn aabb(&self) -> Option<AABB> {
        match self {
            Self::Node(bv, _, _) => Some(*bv),
            Self::Leaf(s) => s.aabb(),
        }
    }
}

use BVH::*;
impl BVH {
    /// Construct a new BVH using the midpoint method and cutting on widest axis
    pub fn new(surfaces: Vec<Rc<dyn Surface>>) -> BVH {
        if surfaces.is_empty() {
            panic!()
        };

        let leaves: Vec<BVH> =
            surfaces.iter().map(|s| Leaf(s.clone())).collect();

        recursive_build_bvh(leaves)
    }

    pub fn center(&self) -> Vec4 {
        match self {
            Self::Node(bv, _, _) => bv.center(),
            Self::Leaf(surface) => surface.aabb().unwrap().center(),
        }
    }
}

/// Builds a bvh recursively
fn recursive_build_bvh(nodes: Vec<BVH>) -> BVH {
    if (nodes.len() == 1) {
        let mut nodes = nodes.to_vec();
        return nodes.pop().unwrap();
    }
    // TODO clean this up
    // get bounds of all nodes
    let min_x = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .min_by(|x1, x2| x1.total_cmp(x2))
        .unwrap();
    let min_y = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .min_by(|y1, y2| y1.total_cmp(y2))
        .unwrap();
    let min_z = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .min_by(|z1, z2| z1.total_cmp(z2))
        .unwrap();
    let max_x = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .max_by(|x1, x2| x1.total_cmp(x2))
        .unwrap();
    let max_y = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .max_by(|y1, y2| y1.total_cmp(y2))
        .unwrap();
    let max_z = nodes
        .iter()
        .map(|node| match node {
            Leaf(surface) => surface.aabb().unwrap(),
            Node(bv, _, _) => *bv,
        })
        .map(AABB::center)
        .map(|v| v.x)
        .max_by(|z1, z2| z1.total_cmp(z2))
        .unwrap();

    let min = vec3(min_x, min_y, min_z);
    let max = vec3(max_x, max_y, max_z);
    let mid = min.lerp(max, 0.5);
    let spread = max - min;

    let mut m = spread.x.abs();
    let mut v = Vec3::X;
    if m < spread.y.abs() {
        m = spread.y.abs();
        v = Vec3::Y;
    }
    if m < spread.z.abs() {
        m = spread.z.abs();
        v = Vec3::Z;
    }

    let v: Vec4 = (v, 0.).into();

    let mut nodes = nodes.to_vec();
    nodes.sort_by(|bvh1, bvh2| {
        bvh1.center().dot(v).total_cmp(&bvh2.center().dot(v))
    });

    let i = nodes.len() / 2;
    let s = Sphere {
        center: Vec4::W,
        radius: 1.,
    };

    let left = recursive_build_bvh(nodes[..i].to_vec());
    let right = recursive_build_bvh(nodes[i..].to_vec());

    let lbv = left.aabb().unwrap();
    let rbv = right.aabb().unwrap();

    let bv = AABB {
        min: (
            f32::min(lbv.min.x, rbv.min.x),
            f32::min(lbv.min.y, rbv.min.y),
            f32::min(lbv.min.z, rbv.min.z),
            1.,
        )
            .into(),
        max: (
            f32::max(lbv.max.x, rbv.max.x),
            f32::max(lbv.max.y, rbv.max.y),
            f32::max(lbv.max.z, rbv.max.z),
            1.,
        )
            .into(),
    };

    Node(bv, Box::new(left), Box::new(right))
}

//impl FromIterator<Rc<dyn Surface>> for BVH {
//    fn from_iter<T: IntoIterator<Item = Rc<dyn Surface>>>(iter: T) -> Self {}
//}
