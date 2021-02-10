use clap::{App, Arg};
use glam::f32::vec3;
use glam::Vec3;
use image::{ImageBuffer, Rgb, RgbImage};
use ray_tracer::Ray;
use ray_tracer::Sphere;

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
        .get_matches();

    let output = "a.png";
    let output = matches.value_of("output").unwrap_or(output);
    let values = matches.values_of("dimensions").unwrap();

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

    let mut imgbuf = ImageBuffer::new(width, height);

    // Camera
    let aspect_ratio = width as f32 / height as f32;
    let viewort_height = 2.;
    let viewport_width = aspect_ratio * viewort_height;
    let focal_length = 1.;

    let origin = Vec3::zero();
    let horizontal = viewport_width * Vec3::unit_x();
    let vertical = viewort_height * Vec3::unit_y();
    let diagonal = focal_length * Vec3::unit_z();
    let upper_left = origin + horizontal / 2. + vertical / 2. - diagonal;

    let sphere = Sphere {
        center: vec3(0., 0., -1.),
        radius: 0.5,
    };
    imgbuf.enumerate_pixels_mut().for_each(|(x, y, pixel)| {
        // Coordinates
        let u = x as f32 / width as f32;
        let v = y as f32 / height as f32;

        let direction = upper_left - u * horizontal - v * vertical - origin;
        let r = Ray { origin, direction };

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
        *pixel = Rgb([r, g, b]);
    });

    imgbuf.save(output).unwrap();
}
