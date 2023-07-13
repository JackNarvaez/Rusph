// Initialize a set of particle to simulate a Toy Star system in 2D.

use std::{
    error::Error,
    process,
};

use sphfunctions;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/toy_star_2D.csv";
    let d: i32 = 2; // Dimension
    let eta: f64 = 1.2;
    let n: u32 = 1000; // Number of Particles
    let r: f64 = 0.75; // Star's radius
    let m: f64 = 2.0; // Star's mass
    let rho: f64 = m/(PI*r*r); // Density
    let h: f64 = 0.1*eta*(m/(n as f64 * rho)).powf(1./d as f64); // Smoothing length
    let (x0, y0) = (0.0, 0.0); // Circle's center
    if let Err(err) = sphfunctions::init_random_circle(path, n, r, rho, h, x0, y0) {
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}