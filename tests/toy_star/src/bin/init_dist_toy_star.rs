// Initial setting for the Toy Star Problem in 3D.

use std::{
    error::Error,
    process,
};

use structures::Particle;
use datafunctions;
use partdistribution;
use sphfunctions::h_by_density;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path = "./Toystar/Ini_00.csv";
    let input_file = "./tests/toy_star/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta:f64  = input[0];         // Dimensionless constant specifying the smoothing length
    let _d: i32  = input[2] as i32;  // Dimensions
    let m: f64   = input[4];         // Star's mass
    let r: f64   = input[5];         // Star's radius
    
    let x0:f64   = input[6];         // Star's center (x-coordinate)
    let y0:f64   = input[7];         // Star's center (y-coordinate)
    let z0:f64   = input[8];         // Star's center (z-coordinate)
    let vx0:f64  = input[9];         // Star's velocity (x-coordinate)
    let vy0:f64  = input[10];        // Star's velocity (y-coordinate)
    let vz0:f64  = input[11];        // Star's velocity (z-coordinate)
    let u0:f64   = input[12];        // Initial energy
    
    let mut n:u32= input[17] as u32; // Resolution

    let mut particles :Vec<Particle> = Vec::new();
    let x_i: f64 = x0-r;
    let y_i: f64 = y0-r;
    let z_i: f64 = z0-r;
    let dm : f64 = 2.*r;
    let rsq: f64 = r*r;

    let rho:f64  = 3. * m/(4.*PI*r*r*r); // Density
    partdistribution::init_dist_hcp(&mut particles, n, n, n, rho, eta, dm, dm, dm, x_i, y_i, z_i, 1.0);
    particles.retain(|particle| distance_center(&particle, x0, y0, z0) < rsq);

    n = particles.len() as u32;
    let dm:f64 = m/n as f64;
    let h: f64 = h_by_density(dm, rho, eta);
    for ii in 0..n as usize {
        particles[ii].h = h;
        particles[ii].vx = vx0;
        particles[ii].vy = vy0;
        particles[ii].vz = vz0;
        particles[ii].u = u0;
    }

    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn distance_center(particle: &Particle, x0: f64, y0: f64, z0: f64) -> f64 {
    let xtem: f64 = particle.x - x0;
    let ytem: f64 = particle.y - y0;
    let ztem: f64 = particle.z - z0;
    return xtem*xtem + ytem*ytem + ztem*ztem;
}