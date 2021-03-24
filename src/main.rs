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
use ray_tracer::{
    random_unit_vector, HitRecord, Origin, Plane, Ray, Sphere, Surface,
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

    // A few surfaces
    let sphere0 = Sphere {
        center: vec4(0., 0., -1., 1.),
        radius: 0.5,
        color: vec3(1., 0., 1.),
    };
    let sphere1 = Sphere {
        center: vec4(-1., 0., -0.5, 1.),
        radius: 0.1,
        color: vec3(1., 0., 1.),
    };
    let sphere2 = Sphere {
        center: vec4(1., 0., -0.5, 1.),
        radius: 0.1,
        color: vec3(1., 0., 1.),
    };
    let sphere3 = Sphere {
        center: vec4(0., 0.7, -0.5, 1.),
        radius: 0.1,
        color: vec3(1., 0., 1.),
    };

    let floor = Plane {
        point: vec4(0., -0.1, 0., 1.),
        normal: Vec4::Y,
    };
    let wall = Plane {
        point: vec4(0., 0., -3., 1.),
        normal: Vec4::Z,
    };

    let triangle = (
        vec3(0.7, -0.1, -0.5),
        vec3(0.75, 0.9, -1.),
        vec3(0.68, -0.1, -1.5),
    );

    let objs = load_obj("tepot.obj", false);

    let aspect_ratio = width as f32 / height as f32;

    let persp_camera = PerspCamera::new(
        vec3(0., 0., 1.),
        vec3(0., 0., -200.),
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
        s: 1.,
    };

    // Since we may extend this to models unknown at compile time, use a reference
    // counted pointer to each of our surfaces
    let mut surfaces: Vec<Rc<dyn Surface>> = Vec::new();
    //surfaces.push(Rc::new(sphere0));
    //surfaces.push(Rc::new(sphere1));
    //surfaces.push(Rc::new(sphere2));
    //surfaces.push(Rc::new(sphere3));
    //surfaces.push(Rc::new(floor));
    //surfaces.push(Rc::new(triangle));
    //surfaces.push(Rc::new(wall));

    let min = vec4(-1., -1., -1., 1.);
    let max = vec4(1., 1., -2., 1.);

    let num_spheres = 10;
    let bounds = AABB { min, max };

    let (w, h, d) = (max - min).xyz().into();
    let w = w.abs();
    let h = h.abs();
    let d = d.abs();

    let mut rng = rand::thread_rng();

    let radius = 1. / 10.;

    println!("radius: {}", radius);

    let surfaces: Vec<Rc<dyn Surface>> = (0..num_spheres)
        .map(|_| {
            let x = min.x + rng.gen_range(0.0..w);
            let y = min.y + rng.gen_range(0.0..h);
            let z = min.z + rng.gen_range(0.0..d);

            let center = vec4(x, y, z, 1.);
            let color = random_unit_vector(&mut rng);

            let sphere = Sphere {
                center,
                radius,
                color,
            };
            Rc::new(sphere) as Rc<dyn Surface>
        })
        .collect();

    let bvh = BVH::new(surfaces);

    let bvh = Rc::new(bvh);

    // A single point light
    // TODO: make this work for multiple lights
    let light = vec3(5., 5., 0.);

    // ImageBuffer can generation images by iteration over pixels
    // (x, y) is the pixel and the return of the closure is the Rgb color value
    let imgbuf = ImageBuffer::from_fn(width, height, |row, col| {
        // If we have jittering then store the nummber in jitters
        let jitters = match jittering {
            Some(jitters) => jitters,
            None => 1,
        };
        // we generate one ray and one color for every sample
        let colors: Vec<Vec3> = (0..jitters)
            .map(|_| {
                let view = view_plane.view(row, col, 0.5, 0.5);
                let r = {
                    if orthographic {
                        todo!();
                        ortho_camera.ray_through(view)
                    } else {
                        persp_camera.ray_through(view)
                    }
                };

                let surfaces = vec![bvh.clone()];

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
        let c: Vec3 = colors.iter().sum::<Vec3>() * (1. / jitters as f32);

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

    let n = match hit.n {
        SurfaceNormal::Outer(n) => -n,
        SurfaceNormal::Inner(n) => n,
    };

    let n = n.xyz().normalize();

    // reflection of l over n
    let r = l - 2. * l.dot(n) * n;
    let r = r.normalize();

    let light_color = vec3(1.2, 1.2, 1.2);

    // ambient diffuse and specular components
    let ambient = 0.5 * light_color;
    let diffuse = 0.0f32.max(n.dot(l)) * light_color;
    let specular = 0.0f32.max(v.dot(r)).powi(8) * light_color;

    let phong = ambient + diffuse + specular;

    Mat3::from_diagonal(phong) * hit.c
}
