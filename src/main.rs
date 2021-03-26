#![feature(total_cmp, or_patterns)]

use clap::{App, Arg};
use glam::f32::vec3;
use glam::f32::vec4;
use glam::f32::Vec4;

use glam::f32::Mat3;
use glam::f32::Vec3;
use glam::Vec4Swizzles;
use image::{ImageBuffer, Rgb};
use rand::Rng;
use ray_tracer::camera::{Camera, OrthoCamera, PerspCamera, Viewplane};
use ray_tracer::rand::random_in_box;
use ray_tracer::triangle::Triangle;
use ray_tracer::{
    HitRecord, MeshTriangle, Origin, Plane, Ray, Sphere, Surface,
    SurfaceNormal, AABB, BVH,
};

use std::rc::Rc;
use tobj::load_obj;

fn main() {
    // This is used to validate that the arguments given for width and height
    // are indeed integers
    let u32validator = |x: std::string::String| {
        x.parse::<u32>().map(|_| {}).map_err(|e| e.to_string())
    };

    // This uses the library clap to parse command line arguments
    let matches = App::new(env!("CARGO_PKG_NAME"))
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about(env!("CARGO_PKG_DESCRIPTION"))
        .arg({
            let value_name = "file";
            Arg::with_name("output")
                .short("o")
                .value_name(value_name)
                .takes_value(true)
                //.validator(|x| std::fs::File::open(x).map_err())
                .help(&format!("Write output to <{}>", value_name))
        })
        .arg(
            Arg::with_name("dimensions")
                .short("d")
                .takes_value(true)
                .value_names(&["width", "height"])
                .validator(u32validator)
                .help("Set output dimensions to <width>,<height>"),
        )
        .arg(
            Arg::with_name("orthographic")
                .long("orthographic")
                .help("Use orthographic projection not default (perspective)"),
        )
        .arg(
            Arg::with_name("jittering")
                .long("jittering")
                .help("Render each pixel by jittering <jitters> times")
                .takes_value(true)
                .validator(u32validator)
                .value_name("jitters"),
        )
        .get_matches();

    // Get the values we need from command line args
    let output = "a.png";
    let output = matches.value_of("output").unwrap_or(output);
    let values = matches.values_of("dimensions").unwrap();
    let orthographic = matches.is_present("orthographic");
    let jittering = matches.value_of("jittering");
    let jittering = jittering.map(str::parse::<u32>).map(Result::ok).flatten();

    let width = 500;
    let height = 500;

    // Handle empty values for width and height and convert width and height strings to u32
    let (width, height) = match values
        .take(2)
        .map(str::parse::<u32>)
        .filter_map(Result::ok)
        .collect::<Vec<_>>()
        .as_slice()
    {
        &[w, h, ..] => (w, h),
        &[w, ..] => (w, w),
        _ => (width, height),
    };

    let mut surfaces: Vec<Rc<dyn Surface>> = Vec::new();

    //let mesh_bvh = gen_mesh("teapot.obj");
    //surfaces.push(Rc::new(mesh_bvh));

    let min = vec4(-1., -1., -1., 1.);
    let max = vec4(1., 1., 1., 1.);

    let num_spheres = 100;
    let spheres = gen_spheres(
        AABB { min, max },
        num_spheres,
        1. / (num_spheres as f32).sqrt(),
    );
    let spheres = BVH::new(spheres);

    surfaces.push(Rc::new(spheres));

    let persp_camera = PerspCamera::new(
        vec3(0., 0., 10.),
        vec3(0., 0., -90000.),
        vec3(0., 1., 0.),
    );

    let ortho_camera = OrthoCamera::new(
        vec3(0., 0., 0.),
        vec3(0., 0., -10.),
        vec3(0., 1., 0.),
    );

    let view_plane = Viewplane {
        hres: width,
        vres: height,
        s: 25.,
    };

    // A single point light
    // TODO: make this work for multiple lights
    let light = vec3(0., 1000., 0.);

    // ImageBuffer can generate images by iteration over pixels
    // (x, y) is the pixel and the return of the closure is the Rgb color value
    let imgbuf = ImageBuffer::from_fn(width, height, |row, col| {
        // we generate one ray and one color for every sample
        let samples = 1;
        let colors: Vec<Vec3> = (0..samples)
            .map(|_| {
                let view = view_plane.view(row, col, 0.5, 0.5);
                let r = {
                    if orthographic {
                        ortho_camera.ray_through(view)
                    } else {
                        persp_camera.ray_through(view)
                    }
                };

                // iterate over all surfaces and find the surface which was intersected by the ray
                // first (but not at a negative time value). Then use this surface to generate a
                // color. Or if we dont intersect any surface, return a neutral color
                surfaces
                    .iter()
                    .filter_map(|surface| {
                        // filter out negative t values
                        surface.hit(&r).filter(|hit| hit.t >= 0.)
                    })
                    // get the closest intersected surface by taking the min of the t values
                    .min_by(|hita, hitb| hita.t.total_cmp(&hitb.t))
                    // if there was really a hit of some surface then find out what the color was
                    .map(|hit| phong_shade(hit, r, light))
                    // if we didnt hit any surface give it a neutral color
                    .unwrap_or({
                        // Sky
                        let t = 0.5 * (r.direction().normalize().y + 1.);
                        vec3(1., 1., 1.).lerp(vec3(0.5, 0.7, 1.), t)
                    })
            })
            .collect();

        // Average of all the colors calculated for each jitter
        let c: Vec3 = colors.iter().sum::<Vec3>() * (1. / samples as f32);

        // RgbImage generates the image from u8 slices so we have to convert c
        let r = (c.x * u8::MAX as f32) as u8;
        let g = (c.y * u8::MAX as f32) as u8;
        let b = (c.z * u8::MAX as f32) as u8;
        Rgb([r, g, b])
    });

    imgbuf.save(output).unwrap();
}

/// Calculates color based on hit information, ray, light position with phong reflection model
fn phong_shade(hit: HitRecord, r: Ray, light: Vec3) -> Vec3 {
    // V, L, N, and R for Phong reflection model
    let v = hit.p.xyz().ray_through(r.origin());
    let v = v.direction().normalize();

    // send a ray though the point light
    let l = hit.p.xyz().ray_through(light);
    let l = l.direction().normalize();

    let SurfaceNormal::Outer(n) | SurfaceNormal::Inner(n) = hit.n;

    let n = n.xyz().normalize();

    // reflection of l over n
    let r = 2. * l.dot(n) * n - l;
    let r = r.normalize();

    let light_color = vec3(1.2, 1.2, 1.2);

    // ambient diffuse and specular components
    let ambient = 0.5 * light_color;
    let diffuse = 0.0f32.max(n.dot(l)) * light_color;
    let specular = 0.0f32.max(v.dot(r)).powi(8) * light_color;

    let phong = ambient + diffuse + specular;

    Mat3::from_diagonal(phong) * hit.c
}

/// Generates a mesh BVH from an obj file
fn gen_mesh(filename: &str) -> BVH {
    // parse the obj file
    let objs = load_obj(filename, true).unwrap();

    // assuming there is at least one model in the obj file we can get the first model
    let models = objs.0;
    let model = models.iter().next().unwrap();

    let mut mesh = model.mesh.clone();
    let normals = &mut mesh.normals;
    let positions = &mesh.positions;
    let indices = &mesh.indices;
    assert!(indices.len() != 0);

    normals.resize(positions.len(), 0.);

    // Compute area weighted normals
    for idxs in indices.chunks(3) {
        let astart = 3 * idxs[0] as usize;
        let bstart = 3 * idxs[1] as usize;
        let cstart = 3 * idxs[2] as usize;

        let aend = astart + 3;
        let bend = bstart + 3;
        let cend = cstart + 3;

        let arange = astart..aend;
        let brange = bstart..bend;
        let crange = cstart..cend;

        let a = &positions[arange.clone()];
        let b = &positions[brange.clone()];
        let c = &positions[crange.clone()];

        let a = Vec3::from_slice_unaligned(a);
        let b = Vec3::from_slice_unaligned(b);
        let c = Vec3::from_slice_unaligned(c);

        let triangle = Triangle { a, b, c };
        let normal = triangle.normal().unwrap_or(Vec3::ZERO) * triangle.area();

        let asum =
            Vec3::from_slice_unaligned(&normals[arange.clone()]) + normal;
        let bsum =
            Vec3::from_slice_unaligned(&normals[brange.clone()]) + normal;
        let csum =
            Vec3::from_slice_unaligned(&normals[crange.clone()]) + normal;

        asum.write_to_slice_unaligned(&mut normals[arange]);
        bsum.write_to_slice_unaligned(&mut normals[brange]);
        csum.write_to_slice_unaligned(&mut normals[crange]);
    }

    // Normalize area weighted normals
    for normal in &mut normals.chunks_mut(3) {
        let n = Vec3::from_slice_unaligned(normal);
        n.normalize().write_to_slice_unaligned(normal);
    }

    let indices = &mesh.indices;
    let mesh = &mesh;

    let mut triangles: Vec<Rc<dyn Surface>> = indices
        .chunks(3)
        .enumerate()
        .map(|(i, _)| Rc::new(MeshTriangle::new(mesh, i)) as Rc<dyn Surface>)
        .collect();

    BVH::new(triangles)
}

/// Generates a sphere bvh
fn gen_spheres(
    bounds: AABB,
    num_spheres: u32,
    radius: f32,
) -> Vec<Rc<dyn Surface>> {
    let mut rng = rand::thread_rng();

    (0..num_spheres)
        .map(|_| {
            let center = random_in_box(&mut rng, bounds);
            let color = random_in_box(
                &mut rng,
                AABB {
                    min: Vec4::ZERO,
                    max: Vec4::ONE,
                },
            )
            .xyz();

            let sphere = Sphere {
                center,
                radius,
                color,
            };
            Rc::new(sphere) as Rc<dyn Surface>
        })
        .collect()
}
