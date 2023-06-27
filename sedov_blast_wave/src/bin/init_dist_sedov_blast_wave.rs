// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/sedov_blast_wave.csv";
    let n:u32 = 100; // Number of Particles
    let w:f64 = 1.; // width
    let l:f64 = 1.; // large
    let x0:f64 = -0.5; // initial x
    let y0:f64 = -0.5; // initial y
    let rho:f64 = 1.0; // density
    let h = w*l / n as f64; // Smoothing length
    if let Err(err) = sphfunctions::init_square(path, n, rho, h, w, l, x0, y0){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}