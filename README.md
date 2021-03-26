```toml

[dependencies]
image = "0.23.12" # for png generation
clap = "2.33.3" # for argument parsing
glam = "0.13.0" # for linear algebra
rand = "0.8.3" # for random number generation
tobj = "2.0.4" # for obj loading
triangle = "0.1.351" # for barycentric interpolation
```
```
$ cargo run -- -h
ray-tracer 0.1.0
Adrian Herath <adrianisuru@gmail.com>


USAGE:
    ray-tracer [FLAGS] [OPTIONS]

FLAGS:
    -h, --help            Prints help information
        --orthographic    Use orthographic projection not default (perspective)
    -V, --version         Prints version information

OPTIONS:
    -d <width> <height>          Set output dimensions to <width>,<height>
        --jittering <jitters>    Render each pixel by jittering <jitters> times
    -o <file>                    Write output to <file>

```

Documentation is easier to read if generated with `cargo doc --open`
