#![feature(total_cmp, or_patterns)]

use clap::{App, Arg};
use glam::f32::vec3;
use glam::Vec3;
use image::{ImageBuffer, Rgb};
use rand::Rng;
use ray_tracer::{
    random_unit_sphere, Camera, Origin, OrthoCamera, PerspCamera, Plane, Ray,
    Sphere, Surface, SurfaceNormal,
};
use std::rc::Rc;

fn main() {
    let u32validator = |x: std::string::String| {
        x.parse::<u32>().map(|_| {}).map_err(|e| e.to_string())
    };
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

    let output = "a.png";
    let output = matches.value_of("output").unwrap_or(output);
    let values = matches.values_of("dimensions").unwrap();
    let orthographic = matches.is_present("orthographic");
    let jittering = matches.value_of("jittering");

    let width = 500;
    let height = 500;
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

    let sphere0 = Sphere {
        center: vec3(0., 0., -1.),
        radius: 0.5,
    };
    let sphere1 = Sphere {
        center: vec3(-1., 0., -0.5),
        radius: 0.1,
    };
    let sphere2 = Sphere {
        center: vec3(1., 0., -0.5),
        radius: 0.1,
    };
    let sphere3 = Sphere {
        center: vec3(0., 0.7, -0.5),
        radius: 0.1,
    };

    let floor = Plane {
        point: vec3(0., -0.2, 0.),
        normal: Vec3::unit_y(),
    };
    let wall = Plane {
        point: vec3(0., 0., -3.),
        normal: Vec3::unit_z(),
    };

    let triangle = (
        vec3(0., 0.5, -1.),
        vec3(-0.5, -0.5, -1.),
        vec3(0.5, -0.5, -1.),
    );

    let aspect_ratio = width as f32 / height as f32;
    let zoom = 1.;

    let persp_camera = PerspCamera {
        location: vec3(0., 0., 0.),
        //rotation: vec3(0., 0., 0.),
        focal_length: 1.,
        aspect_ratio,
        zoom,
    };
    let ortho_camera = OrthoCamera {
        location: vec3(0., 0., 0.),
        //rotation: vec3(0., 0., 0.),
        aspect_ratio,
        zoom,
    };

    let mut surfaces: Vec<Rc<dyn Surface>> = Vec::new();
    surfaces.push(Rc::new(sphere0));
    surfaces.push(Rc::new(sphere1));
    surfaces.push(Rc::new(sphere2));
    surfaces.push(Rc::new(sphere3));
    surfaces.push(Rc::new(floor));
    //surfaces.push(Rc::new(triangle));
    //surfaces.push(Rc::new(wall));

    let light = vec3(10., 10., 0.);

    let jittering = jittering.map(str::parse::<u32>).map(Result::ok).flatten();
    let mut rng = rand::thread_rng();
    let imgbuf = ImageBuffer::from_fn(width, height, |x, y| {
        let jitters = match jittering {
            Some(jitters) => jitters,
            None => 1,
        };
        let c: Vec<Vec3> = (0..jitters)
            .map(|_| {
                let x = x as f32;
                let y = y as f32;

                let dx: f32;
                let dy: f32;

                match jittering {
                    Some(_) => {
                        dx = rng.gen_range(0. ..1.);
                        dy = rng.gen_range(0. ..1.);
                    }
                    None => {
                        dx = 0.;
                        dy = 0.;
                    }
                };
                let u = ((x + dx) / width as f32) - 0.5;
                let v = (1. - (y + dy) / height as f32) - 0.5;

                let r = if orthographic {
                    ortho_camera.ray_through(u, v)
                } else {
                    persp_camera.ray_through(u, v)
                };

                surfaces
                    .iter()
                    .filter_map(|surface| {
                        surface
                            .hit(&r)
                            .filter(|hit| hit.t >= 0.)
                            .map(|hit| (surface, hit))
                    })
                    .min_by(|(_a, hita), (_b, hitb)| hita.t.total_cmp(&hitb.t))
                    .and_then(|(surface, hit)| {
                        // V, L, N, and R for Phong reflection model
                        let v = hit.p.ray_through(r.origin());
                        let v = v.direction();
                        let l = hit.p.ray_through(light);
                        let SurfaceNormal::Inner(n) | SurfaceNormal::Outer(n) =
                            hit.n;
                        let n = n.normalize();
                        let r = r.direction() - 2. * r.direction().dot(n) * n;
                        let r = r.normalize();

                        let lt = (light - l.origin()).length();

                        let c = match surfaces
                            .iter()
                            .filter(|s| !Rc::ptr_eq(s, surface))
                            .filter_map(|s| {
                                s.hit(&l).filter(|h| h.t >= 0.).map(|h| (s, h))
                            })
                            .min_by(|(a, ha), (b, hb)| ha.t.total_cmp(&hb.t))
                            .filter(|(s, h)| h.t < lt)
                        {
                            Some((s, hit)) => {
                                let shadow_percent = 0.4;

                                surface
                                    .color()
                                    .lerp(Vec3::zero(), shadow_percent)
                            }
                            None => {
                                // Phong reflection model
                                surface.color()
                                    + surface.color() * (l.direction().dot(n))
                                    + surface.color() * r.dot(v)
                            }
                        };

                        Some(c)
                    })
                    .or_else(|| {
                        // Sky
                        let t = 0.5 * (r.direction().normalize().y + 1.);
                        Some(vec3(1., 1., 1.).lerp(vec3(0.5, 0.7, 1.), t))
                    })
            })
            .filter_map(|a| a)
            .collect();

        let c: Vec3 = c.iter().sum::<Vec3>() * (1. / jitters as f32);

        let r = (c.x * u8::MAX as f32) as u8;
        let g = (c.y * u8::MAX as f32) as u8;
        let b = (c.z * u8::MAX as f32) as u8;
        Rgb([r, g, b])
    });

    imgbuf.save(output).unwrap();
}
