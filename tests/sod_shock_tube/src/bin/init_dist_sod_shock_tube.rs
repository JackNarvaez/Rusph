// Initial setting for the Sod Shock Tube Problem in 3D

use std::{
    error::Error,
    process,
};

use structures::Particle;
use sphfunctions::h_by_density;
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
    let n: u32   = nxl*nyl*nzl + nxr*nyr*nzr;   // Total number of particles
    
    let vol: f64 = wd*lg*hg;                    // Volumen
    let dm: f64  = 0.5*vol*(rhol+rhor)/n as f64;// Particles' mass
    
    let mut particles :Vec<Particle> = Vec::new();

    init_dist_sod_tube(&mut particles, nxl, nyl, nzl, nxr, nyr, nzr, rhol, rhor, ul, ur, xm, eta, wd, lg, hg, x0, y0, z0, dm);

    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn init_dist_sod_tube(particles: &mut Vec<Particle>, nxl: u32, nyl: u32, nzl: u32, nxr: u32, nyr: u32, nzr: u32, rhol: f64, rhor: f64,
                      ul: f64, ur: f64, xm: f64, eta: f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64, dm: f64){
    let dxl: f64 = (xm-x0)/nxl as f64;
    let dyl: f64 = lg/nyl as f64;
    let dzl: f64 = hg/nzl as f64;
    let dxr: f64 = (x0+wd-xm)/nxr as f64;
    let dyr: f64 = lg/nyr as f64;
    let dzr: f64 = hg/nzr as f64;
    let hpl: f64 = h_by_density(dm, rhol, eta);
    let hpr: f64 = h_by_density(dm, rhor, eta);

    // Left State
    for kk in 0..nzl {
        for jj in 0..nyl{
            for ii in 0..nxl{
                let xp: f64 = x0 + dxl*ii as f64;
                let yp: f64 = y0 + dyl*jj as f64;
                let zp: f64 = z0 + dzl*kk as f64;
                particles.push(Particle{x:xp, y:yp, z:zp,
                                        vx: 0.0, vy: 0.0, vz: 0.0,
                                        h:hpl, u: ul,
                                        ..Default::default()});
            }
        }
    }
    // Right State
    for kk in 0..nzr {
        for jj in 0..nyr{
            for ii in 0..nxr{
                let xp: f64 = xm + dxr*ii as f64;
                let yp: f64 = y0 + dyr*jj as f64;
                let zp: f64 = z0 + dzr*kk as f64;
                particles.push(Particle{x:xp, y:yp, z:zp,
                                        vx: 0.0, vy: 0.0, vz: 0.0,
                                        h:hpr, u: ur,
                                        ..Default::default()});
            }
        }
    }
}