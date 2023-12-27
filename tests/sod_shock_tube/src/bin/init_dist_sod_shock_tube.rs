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
    let path = "./Sodtube/Ini_00.csv";
    let input_file = "./tests/sod_shock_tube/input";
    
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64 = input[0];         // Dimensionless constant specifying the smoothing length
    let gamma:f64= input[1];         // Heat capacity ratio
    let _d: i32  = input[2] as i32;  // Dimensions
    
    let x0:f64   = input[3];         // Bottom left corner  (x-coordinate)
    let y0:f64   = input[4];         // Bottom left corner  (y-coordinate)
    let z0:f64   = input[5];         // Bottom left corner  (z-coordinate)
    let wd:f64   = input[6];         // Width (x)
    let lg:f64   = input[7];         // Length (y)
    let hg:f64   = input[8];         // Heigth (z)

    let rhol:f64 = input[9];         // Left density
    let rhor:f64 = input[10];        // Right density
    let pl:f64   = input[11];        // Left initial pressure
    let pr:f64   = input[12];        // Right initial pressure
    let xm: f64  = input[13];        // Discontinuity position
    
    let nxl:u32  = input[18] as u32; // Particle resolution in the x direction (left)
    let nxr:u32  = input[19] as u32; // Particle resolution in the x direction (right)
    let nyl:u32  = input[20] as u32; // Particle resolution in the y direction (left)
    let nyr:u32  = input[21] as u32; // Particle resolution in the y direction (right)
    let nzl:u32  = input[22] as u32; // Particle resolution in the z direction (left)
    let nzr:u32  = input[23] as u32; // Particle resolution in the z direction (right)

    let ul:f64   = pl/((gamma - 1.)*rhol);      // Left initial energy
    let ur:f64   = pr/((gamma - 1.)*rhor);      // Right initial energy
    let n:u32  = nxl*nyl*nzl + nxr*nyr*nzr;   // Fluid particles
    
    let vol: f64 = wd*lg*hg;                    // Volumen
    let dm: f64  = 0.5*vol*(rhol+rhor)/n as f64;// Particles' mass
    let bxl: f64 = (xm-x0)/6.;
    let bxr: f64 = (x0+wd-xm)/6.;
    
    let mut particles :Vec<Particle> = Vec::new();

    // Left State
    partdistribution::init_dist_cubic(&mut particles, nxl, nyl, nzl, rhol, eta, xm-x0, lg, hg, x0, y0, z0, dm);

    // Right State
    partdistribution::init_dist_cubic(&mut particles, nxr, nyr, nzr, rhor, eta, x0+wd-xm, lg, hg, xm, y0, z0, dm);

    for ii in 0..n as usize {
        if particles[ii].x <= xm {
            particles[ii].u = ul;
            if particles[ii].x <= x0 + bxl {
                particles[ii].ptype = 1;
            }
        } else {
            particles[ii].u = ur;
            if particles[ii].x >= wd - bxr {
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