// Initialize a set of particle to simulate a Toy Star system in 2D.

use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/accretion_disk.csv";
    let n:u32 = 100; // Number of Particles
    let r:f64 = 0.75; // Star's radius
    let rho:f64 = 1.0; // Density
    let h = 0.04 /(n as f64 /1000.).sqrt(); // Smoothing length
    let (x0, y0) = (0.0, 0.0); // Circle's center
    if let Err(err) = sphfunctions::init_random_circle(path, n, r, rho, h, x0, y0) {
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}