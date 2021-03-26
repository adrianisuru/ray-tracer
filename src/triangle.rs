use glam::f32::{vec3, Vec3};

/// A triangle based on triangle::lib32::Triangle but designed to work well with glam.
#[derive(Debug)]
pub struct Triangle {
    pub a: Vec3,
    pub b: Vec3,
    pub c: Vec3,
}

impl From<triangle::lib32::Triangle> for Triangle {
    fn from(t: triangle::lib32::Triangle) -> Self {
        Triangle {
            a: point_to_vec3(t.a),
            b: point_to_vec3(t.b),
            c: point_to_vec3(t.c),
        }
    }
}

impl From<Triangle> for triangle::lib32::Triangle {
    fn from(t: Triangle) -> Self {
        triangle::lib32::Triangle {
            a: vec3_to_point(t.a),
            b: vec3_to_point(t.b),
            c: vec3_to_point(t.c),
        }
    }
}

/// Converts the given point to a Vec3
fn point_to_vec3(p: triangle::lib32::Point) -> Vec3 {
    vec3(p.x, p.y, p.z)
}

/// Converts the given Vec3 to a point
fn vec3_to_point(v: Vec3) -> triangle::lib32::Point {
    triangle::lib32::Point {
        x: v.x,
        y: v.y,
        z: v.z,
    }
}

/// Below code taken from https://github.com/p4ymak/triangle
impl Triangle {
    ///Returns two opposite points of axis-aligned bounding box.
    pub fn bounding_box(&self) -> [Vec3; 2] {
        let mut c_x = [self.a.x, self.b.x, self.c.x];
        let mut c_y = [self.a.y, self.b.y, self.c.y];
        let mut c_z = [self.a.z, self.b.z, self.c.z];
        c_x.sort_by(|i, j| i.partial_cmp(j).unwrap());
        c_y.sort_by(|i, j| i.partial_cmp(j).unwrap());
        c_z.sort_by(|i, j| i.partial_cmp(j).unwrap());
        return [vec3(c_x[0], c_y[0], c_z[0]), vec3(c_x[2], c_y[2], c_z[2])];
    }

    ///Gets angles of the triangle.
    pub fn angles(&self) -> Option<[f32; 3]> {
        if self.is_collinear() {
            return None;
        }
        let [la, lb, lc] = self.sides();
        let alpha =
            ((lb.powi(2) + lc.powi(2) - la.powi(2)) / (2.0 * lb * lc)).acos();
        let beta =
            ((la.powi(2) + lc.powi(2) - lb.powi(2)) / (2.0 * la * lc)).acos();
        let gamma = std::f32::consts::PI - alpha - beta;
        return Some([alpha, beta, gamma]);
    }

    ///Gets area of the triangle.
    pub fn area(&self) -> f32 {
        let s = self.semiperimeter();
        let [la, lb, lc] = self.sides();
        return (s * (s - la) * (s - lb) * (s - lc)).sqrt();
    }

    ///Converts barycentric coordinates of given point to cartesian coordinate system.
    pub fn barycentric_to_cartesian(&self, pt: &Vec3) -> Vec3 {
        let x = pt.x * self.a.x + pt.y * self.b.x + pt.z * self.c.x;
        let y = pt.x * self.a.y + pt.y * self.b.y + pt.z * self.c.y;
        let z = pt.x * self.a.z + pt.y * self.b.z + pt.z * self.c.z;
        return vec3(x, y, z);
    }

    ///Converts cartesian coordinates of given point to barycentric coordinate system.
    pub fn cartesian_to_barycentric(&self, pt: &Vec3) -> Vec3 {
        let v0 = self.b - self.a;
        let v1 = self.c - self.a;
        let v2 = *pt - self.a;
        let den = 1.0 / (v0.x * v1.y - v1.x * v0.y);
        let v = (v2.x * v1.y - v1.x * v2.y) * den;
        let w = (v0.x * v2.y - v2.x * v0.y) * den;
        let u = 1.0 - v - w;
        return vec3(u, v, w);
    }

    ///Gets centroid of the triangle.
    pub fn centroid(&self) -> Vec3 {
        return vec3(
            (self.a.x + self.b.x + self.c.x) / 3.0,
            (self.a.y + self.b.y + self.c.y) / 3.0,
            (self.a.z + self.b.z + self.c.z) / 3.0,
        );
    }

    ///Gets radius of a circle that passes through all of the triangle's vertices, so called
    ///circumradius.
    pub fn circumradius(&self) -> Option<f32> {
        if self.is_collinear() {
            return None;
        }
        return Some(
            self.sides().iter().product::<f32>() / (4.0 * self.area()),
        );
    }

    ///Checks whether a given point lies inside the triangle.
    pub fn has_point(&self, pt: Vec3) -> bool {
        fn sign(a: &Vec3, b: &Vec3, c: &Vec3) -> f32 {
            ((a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y)) as f32
        }
        let d1 = sign(&pt, &self.a, &self.b);
        let d2 = sign(&pt, &self.b, &self.c);
        let d3 = sign(&pt, &self.c, &self.a);
        let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
        let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
        return !(has_neg && has_pos);
    }

    ///Gets the heights of the triangle.
    pub fn heights(&self) -> Option<[f32; 3]> {
        if self.is_collinear() {
            return None;
        }
        let double_area = 2.0 * self.area();
        let [la, lb, lc] = self.sides();
        return Some([double_area / la, double_area / lb, double_area / lc]);
    }

    ///Gets radius of a circle which is tangent to each side of the triangle, so called inradius.
    pub fn inradius(&self) -> Option<f32> {
        if self.is_collinear() {
            return None;
        }
        return Some(self.area() / self.semiperimeter());
    }

    ///Checks if points of triangle are collinear.
    pub fn is_collinear(&self) -> bool {
        return self.area() == 0.0;
    }

    ///Checks if the triangle is equilateral.
    pub fn is_equilateral(&self) -> bool {
        let sides = self.sides();
        return sides[0] == sides[1] && sides[1] == sides[2];
    }

    ///Checks if the triangle is golden or sublime.
    pub fn is_golden(&self) -> bool {
        if !self.is_isosceles() {
            return false;
        }
        let mut sides = self.sides();
        sides.sort_by(|a, b| a.partial_cmp(&b).unwrap());
        let min = sides[0];
        let max = sides[2];
        return max / min == (1.0 + 5.0_f32.sqrt()) / 2.0;
    }

    ///Checks if the triangle is isosceles.
    pub fn is_isosceles(&self) -> bool {
        let sides = self.sides();
        return sides[0] == sides[1]
            || sides[1] == sides[2]
            || sides[2] == sides[0];
    }

    ///Checks if the triangle is right-angled.
    pub fn is_right(&self) -> bool {
        if self.is_collinear() {
            return false;
        }
        let angles = self.angles().unwrap();
        let half_pi = std::f32::consts::PI / 2.0;
        return angles[0] == half_pi
            || angles[1] == half_pi
            || angles[2] == half_pi;
    }

    ///Gets medians of the triangle.
    pub fn medians(&self) -> [f32; 3] {
        let [la, lb, lc] = self.sides();
        let ma =
            (2.0 * lb.powi(2) + 2.0 * lc.powi(2) - la.powi(2)).sqrt() / 2.0;
        let mb =
            (2.0 * lc.powi(2) + 2.0 * la.powi(2) - lb.powi(2)).sqrt() / 2.0;
        let mc =
            (2.0 * la.powi(2) + 2.0 * lb.powi(2) - lc.powi(2)).sqrt() / 2.0;
        return [ma, mb, mc];
    }

    ///Gets normal of the triangle, depending on vertices order.
    pub fn normal(&self) -> Option<Vec3> {
        if self.is_collinear() {
            return None;
        }
        let u = self.b - self.a;
        let v = self.c - self.a;
        let n = vec3(
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x,
        );
        return Some(n.normalize());
    }

    ///Gets perimeter of the triangle.
    pub fn perimeter(&self) -> f32 {
        return self.sides().iter().sum();
    }

    ///Gets distance from ray origin to intersection with triangle. Möller & Trumbore algorithm.
    pub fn ray_intersection(
        &self,
        ray_orig: &Vec3,
        ray_dir: &Vec3,
    ) -> Option<f32> {
        if self.is_collinear() {
            return None;
        }

        let e1 = self.b - self.a;
        let e2 = self.c - self.a;
        let pvec = ray_dir.cross(e2);
        let det = e1.dot(pvec);
        if det.abs() < f32::MIN {
            return None;
        }

        let inv_det = 1.0 / det;
        let tvec = *ray_orig - self.a;
        let u = tvec.dot(pvec) * inv_det;
        if u < 0.0 || u > 1.0 {
            return None;
        }

        let qvec = tvec.cross(e1);
        let v = ray_dir.dot(qvec) * inv_det;
        if v < 0.0 || (u + v) > 1.0 {
            return None;
        }

        return Some(e2.dot(qvec) * inv_det);
    }

    ///Gets semiperimeter of the triangle.
    pub fn semiperimeter(&self) -> f32 {
        return self.perimeter() / 2.0;
    }

    ///Gets lengths of sides opposite to points.
    pub fn sides(&self) -> [f32; 3] {
        return [
            (self.b - self.c).length(),
            (self.c - self.a).length(),
            (self.a - self.b).length(),
        ];
    }
}
