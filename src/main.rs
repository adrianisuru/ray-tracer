use clap::{App, Arg};
use glam::f32::vec3;
use glam::Vec3;
use image::{ImageBuffer, Rgb};
use ray_tracer::{Camera, OrthoCamera, PerspCamera, Sphere};

fn main() {
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
                .validator(|x| {
                    x.parse::<u32>().map(|_| {}).map_err(|e| e.to_string())
                })
                .help("Set output dimensions to <width>,<height>"),
        )
        .arg(
            Arg::with_name("orthographic")
                .long("orthographic")
                .help("Use orthographic projection not default (perspective)"),
        )
        .get_matches();

    let output = "a.png";
    let output = matches.value_of("output").unwrap_or(output);
    let values = matches.values_of("dimensions").unwrap();
    let orthographic = matches.is_present("orthographic");

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

    let sphere = Sphere {
        center: vec3(0., 0., -1.),
        radius: 0.5,
    };

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

    let imgbuf = ImageBuffer::from_fn(width, height, |x, y| {
        let u = (x as f32 / width as f32) - 0.5;
        let v = (1. - y as f32 / height as f32) - 0.5;

        let r = if orthographic {
            ortho_camera.ray_through(u, v)
        } else {
            persp_camera.ray_through(u, v)
        };

        let t = sphere.intersects(&r);

        let c = if t > 0. {
            let n = (r.at(t) + Vec3::unit_z()).normalize();
            0.5 * (n + Vec3::one())
        } else {
            let t = 0.5 * (r.direction.normalize().y + 1.);
            (1. - t) * Vec3::one() + t * vec3(0.5, 0.7, 1.0)
        };

        let r = (c.x * u8::MAX as f32) as u8;
        let g = (c.y * u8::MAX as f32) as u8;
        let b = (c.z * u8::MAX as f32) as u8;
        Rgb([r, g, b])
    });

    imgbuf.save(output).unwrap();
}
