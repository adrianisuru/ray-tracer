use crate::Ray;
use glam::f32::{vec2, vec3, Mat3, Mat4, Vec2, Vec3, Vec4};

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

pub struct Viewplane {
    pub hres: u32,
    pub vres: u32,
    pub s: f32,
}
