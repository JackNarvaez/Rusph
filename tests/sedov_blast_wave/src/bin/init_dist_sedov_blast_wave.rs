// Initial setting for the Sedov Blast Wave  Problem in 3D

use std::{
    error::Error,
    process,
};

use structures::Particle;

use sphfunctions::f_cubic_kernel;
use datafunctions;

use partdistribution;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path = "./Sedov/Ini_00.csv";
    let input_file = "./tests/sedov_blast_wave/input";
    
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64 = input[0];         // Dimensionless constant specifying the smoothing length
    let _d: i32   = input[2] as i32; // Dimensions
    
    let x0:f64   = input[3];         // Bottom left corner  (x-coordinate)
    let y0:f64   = input[4];         // Bottom left corner  (y-coordinate)
    let z0:f64   = input[5];         // Bottom left corner  (z-coordinate)
    let wd:f64   = input[6];         // Width (x)
    let lg:f64   = input[7];         // Length (y)
    let hg:f64   = input[8];         // Heigth (z)

    let rho:f64  = input[9];         // Density
    let e0:f64   = input[10];        // Initial energy
    
    let nx:u32   = input[16] as u32; // Particle resolution in the x direction
    let ny:u32   = input[17] as u32; // Particle resolution in the y direction
    let nz:u32   = input[18] as u32; // Particle resolution in the z direction
    let mut n: u32   = nx*ny*nz;         // Total number of particles
    
    let rkern: f64 = 2.;             // Kernel radius
    
    let vol: f64 = wd*lg*hg;         // Volumen
    let mut dm: f64  = rho*vol/n as f64; // Particles' mass
    let h0: f64  = 2.*eta*(wd/nx as f64);
    
    let mut u_norm:f64 = 0.0;        // Normalization's constant for the initial energy
    let mut rad_part: Vec<usize> = Vec::new(); // Particles inside the initial sphere
    let mut particles :Vec<Particle> = Vec::new();

    partdistribution::init_dist_hcp(&mut particles, nx, ny, nz, rho, eta, wd, lg, hg, x0, y0, z0, dm);
    n = particles.len() as u32;
    dm = rho*vol/n as f64;
    init_dist_sedov(&mut particles, &mut rad_part, &mut u_norm, n as usize, rkern, h0, wd, lg, hg, f_cubic_kernel);
    norm_energy(&mut particles, &mut rad_part, e0, u_norm, dm);

    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn norm_energy(particles: &mut Vec<Particle>, rad_part: & Vec<usize>, e0: f64, u_norm: f64, dm: f64) {
    let u0: f64 = e0/(dm * u_norm);
    for ii in rad_part {
        particles[*ii].u *= u0;
    }
}

fn init_dist_sedov(particles: &mut Vec<Particle>, rad_part: &mut Vec<usize>, u_norm: &mut f64, n: usize, rkern: f64, h0: f64,
                   wd:f64, lg:f64, hg: f64, kernel: fn(f64) -> f64){
    let xc: f64 = wd/2.;
    let yc: f64 = lg/2.;
    let zc: f64 = hg/2.;

    let mut q: f64;
    let mut up: f64;
    for ii in 0..n {
        q = ((particles[ii].x-xc)*(particles[ii].x-xc) + (particles[ii].y-yc)*(particles[ii].y-yc) + (particles[ii].z-zc)*(particles[ii].z-zc)).sqrt()/ h0;
        up = 0.;
        if q<rkern {
            up = kernel(q);
            *u_norm += up;
            rad_part.push(ii);
        }
        particles[ii].u = up;
    }
}