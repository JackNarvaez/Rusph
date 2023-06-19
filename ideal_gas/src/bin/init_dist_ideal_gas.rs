// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/ideal_gas.csv";
    let n:u32 = 100; // Number of Particles
    let w:f64 = 1.; // width
    let l:f64 = 1.; // large
    let rho:f64 = 1.0; // density
    let h = 0.1; // Smoothing length
    if let Err(err) = sphfunctions::init_square(path, n, rho, h, w, l){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}