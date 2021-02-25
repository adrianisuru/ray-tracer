#![feature(total_cmp, or_patterns)]

use clap::{App, Arg};
use glam::f32::vec3;
use glam::Vec3;
use image::{ImageBuffer, Rgb};
use rand::Rng;
use ray_tracer::{
    Camera, HitRecord, Origin, OrthoCamera, PerspCamera, Plane, Sphere,
    Surface, SurfaceNormal,
};
use std::rc::Rc;

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
        point: vec3(0., -0.1, 0.),
        normal: Vec3::unit_y(),
    };
    let wall = Plane {
        point: vec3(0., 0., -3.),
        normal: Vec3::unit_z(),
    };

    let triangle = (
        vec3(0.7, -0.1, -0.5),
        vec3(0.75, 0.9, -1.),
        vec3(0.68, -0.1, -1.5),
    );

    let aspect_ratio = width as f32 / height as f32;

    let persp_camera = PerspCamera {
        location: vec3(0., 0., 0.),
        //rotation: vec3(0., 0., 0.),
        focal_length: 1.,
        aspect_ratio,
        zoom: 1.,
    };
    let ortho_camera = OrthoCamera {
        location: vec3(0., 0., 0.),
        //rotation: vec3(0., 0., 0.),
        aspect_ratio,
        zoom: 2.,
    };

    // Since we may extend this to models unknown at compile time, use a reference
    // counted pointer to each of our surfaces
    let mut surfaces: Vec<Rc<dyn Surface>> = Vec::new();
    surfaces.push(Rc::new(sphere0));
    surfaces.push(Rc::new(sphere1));
    surfaces.push(Rc::new(sphere2));
    surfaces.push(Rc::new(sphere3));
    surfaces.push(Rc::new(floor));
    surfaces.push(Rc::new(triangle));
    //surfaces.push(Rc::new(wall));

    // A single point light
    // TODO: make this work for multiple lights
    let light = vec3(10., 10., -0.5);

    // Rng for jittering
    let mut rng = rand::thread_rng();

    // ImageBuffer can generation images by iteration over pixels
    // (x, y) is the pixel and the return of the closure is the Rgb color value
    let imgbuf = ImageBuffer::from_fn(width, height, |x, y| {
        // If we have jittering then store the nummber in jitters
        let jitters = match jittering {
            Some(jitters) => jitters,
            None => 1,
        };
        // we generate one ray and one color for every sample
        let colors: Vec<Vec3> = (0..jitters)
            .map(|_| {
                let x = x as f32;
                let y = y as f32;

                // if we use jittering, the ray is not sent directly through the center of the
                // pixel
                let dx: f32;
                let dy: f32;

                match jittering {
                    Some(_) => {
                        dx = rng.gen_range(0. ..1.);
                        dy = rng.gen_range(0. ..1.);
                    }
                    // if no jittering then we shoot a ray through the point exactly
                    None => {
                        dx = 0.;
                        dy = 0.;
                    }
                };

                // convert to camera coordinates
                let u = ((x + dx) / width as f32) - 0.5;
                let v = (1. - (y + dy) / height as f32) - 0.5;

                // our ray r depends on which camera we use
                let r = if orthographic {
                    ortho_camera.ray_through(u, v)
                } else {
                    persp_camera.ray_through(u, v)
                };

                // iterate over all surfaces and find the surface which was intersected by the ray
                // first (but not at a negative time value). Then use this surface to generate a
                // color. Or if we dont intersect any surface, return a neutral color
                surfaces
                    .iter()
                    .filter_map(|surface| {
                        // filter out negative t values
                        surface
                            .hit(&r)
                            .filter(|hit| hit.t >= 0.)
                            .map(|hit| (surface, hit)) // we want the surface and the hitrecord
                    })
                    // get the closest intersected surface by taking the min of the t values
                    .min_by(|(_a, hita), (_b, hitb)| hita.t.total_cmp(&hitb.t))
                    // if there was really a hit of some surface then find out what the color was
                    .map(|(surface, hit)| {
                        // V, L, N, and R for Phong reflection model
                        let v = hit.p.ray_through(r.origin());
                        let v = v.direction();

                        // send a ray though the point light
                        let l = hit.p.ray_through(light);

                        // TODO: handle inner and outer differently
                        let SurfaceNormal::Inner(n) | SurfaceNormal::Outer(n) =
                            hit.n;
                        let n = n.normalize();

                        // reflection of l over n
                        let r = r.direction() - 2. * r.direction().dot(n) * n;
                        let r = r.normalize();

                        // Send a ray to the point light for shadow generation
                        let lt = (light - l.origin()).length();

                        let shadow_percent = 0.5;
                        // Loop over all the surfaces except the current one
                        surfaces
                            .iter()
                            .filter(|s| !Rc::ptr_eq(s, surface)) // ignore the current surface
                            .filter_map(|s| {
                                s.hit(&l).filter(|h| h.t >= 0.).map(|h| (s, h)) // filter negative t values
                            })
                            .min_by(|(_, ha), (_, hb)| ha.t.total_cmp(&hb.t)) // closest intersected surface
                            .filter(|(_, h)| h.t < lt) // if we hit the light first then filter out
                            .map(|_| {
                                // calculate shadows

                                surface
                                    .color()
                                    .lerp(Vec3::zero(), shadow_percent)
                            })
                            // if we do not hit a suface than shade via phong reflection model
                            .unwrap_or({
                                let s = shadow_percent;
                                // Phong reflection model
                                surface.color() * Vec3::one()
                                    + surface.color()
                                        * (l.direction().dot(n))
                                        * (Vec3::one() * s)
                                    + (surface.color() * r.dot(v))
                                        * (Vec3::one() * s)
                            })
                    })
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
