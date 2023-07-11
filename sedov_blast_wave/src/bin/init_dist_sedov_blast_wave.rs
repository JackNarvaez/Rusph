// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/sedov_blast_wave.csv";
    let n:u32 = 2500; // Number of Particles
    let w:f64 = 1.; // Width
    let l:f64 = 1.; // Large
    let x0:f64 = 0.; // Initial x
    let y0:f64 = 0.; // Initial y
    let rho:f64 = 1.0; // Density
    let h = 2.*(w*l / n as f64).sqrt(); // Smoothing length
    if let Err(err) = sphfunctions::init_square(path, n, rho, h, w, l, x0, y0){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}