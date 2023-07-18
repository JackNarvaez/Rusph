// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/sedov_blast_wave.csv";
    let n:u32 = 32*32; // Number of Particles
    let w:f64 = 1.; // Width
    let l:f64 = 1.; // Large
    let x0:f64 = 0.; // Initial x
    let y0:f64 = 0.; // Initial y
    let rho:f64 = 1.0; // Density
    let eta: f64 = 1.2;
    let d: i32 = 2; // Dimensions
    let h: f64 = eta*(w*l / n as f64).powf(1./d as f64); // Smoothing length
    let dist: usize = 1;
    if let Err(err) = sphfunctions::init_square(path, n, rho, h, w, l, x0, y0, dist){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}