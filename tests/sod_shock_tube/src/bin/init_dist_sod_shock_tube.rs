// Initial setting for the Sod Shock Tube Problem in 3D

use std::{
    f64,
    error::Error,
    process,
};

use structures::Particle;
use datafunctions;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path: &str      = "./Sodtube/Ini_00.csv";
    let input_file: &str= "./tests/sod_shock_tube/input";
    
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length
    let gamma: f64  = input[1];         // Heat capacity ratio
    
    let x0: f64     = input[2];         // Bottom left corner  (x-coordinate)
    let y0: f64     = input[3];         // Bottom left corner  (y-coordinate)
    let z0: f64     = input[4];         // Bottom left corner  (z-coordinate)
    let wd: f64     = input[5];         // Width (x)
    let lg: f64     = input[6];         // Length (y)
    let hg: f64     = input[7];         // Height (z)

    let rhol: f64   = input[8];         // Left density
    let rhor: f64   = input[9];         // Right density
    let pl: f64     = input[10];        // Left initial pressure
    let pr: f64     = input[11];        // Right initial pressure
    let xm: f64     = input[12];        // Discontinuity position
    
    let nxl: u32    = input[16] as u32; // Particle resolution in the x direction (left)
    let nxr: u32    = input[17] as u32; // Particle resolution in the x direction (right)

    let ul: f64     = pl/((gamma - 1.)*rhol);      // Left initial energy
    let ur: f64     = pr/((gamma - 1.)*rhor);      // Right initial energy
    
    let bxl: f64    = (xm-x0)/6.;       // Left boundary region
    let bxr: f64    = (x0+wd-xm)/6.;    // Right boundary region
    
    let mut particles : Vec<Particle> = Vec::new();

    // Left State
    partdistribution::init_dist_hcp(&mut particles, nxl, rhol, eta, xm-x0, lg, hg, x0, y0, z0);

    // Right State
    partdistribution::init_dist_hcp(&mut particles, nxr, rhor, eta, x0+wd-xm, lg, hg, xm, y0, z0);

    let n: usize    = particles.len();
    
    // Boundary particles
    for ii in 0..n {
        if particles[ii].x <= xm {
            particles[ii].u = ul;
            if particles[ii].x <= x0 + bxl {
                particles[ii].ptype = 1;
            }
        } else {
            particles[ii].u = ur;
            if particles[ii].x >= x0 + wd - bxr {
                particles[ii].ptype = 1;
            }
        }
    }
    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}