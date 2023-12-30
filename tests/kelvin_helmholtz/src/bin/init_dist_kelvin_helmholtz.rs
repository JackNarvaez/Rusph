// Initial setting for the Kelvin Helmholtz problem in 3D

use std::{
    error::Error,
    process,
};

use datafunctions;
use partdistribution;

use structures::Particle;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path = "./Kelvinhelmholtz/Ini_00.csv";
    let input_file = "./tests/kelvin_helmholtz/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta:f64  = input[0];         // Dimensionless constant specifying the smoothing length
    let gamma:f64= input[1];         // Heat capacity ratio
    let _d: i32   = input[2] as i32; // Dimensions
    
    let x0:f64   = input[3];         // Bottom left corner  (x-coordinate)
    let y0:f64   = input[4];         // Bottom left corner  (y-coordinate)
    let z0:f64   = input[5];         // Bottom left corner  (z-coordinate)
    let wd:f64   = input[6];         // Width (x)
    let lg:f64   = input[7];         // Length (y)
    let hg:f64   = input[8];         // Heigth (z)
    let y1:f64   = input[9];         // Y-lower edge of fluid 2 
    let y2:f64   = input[10];        // Y-upper edge of fluid 2
    let rho1:f64 = input[11];        // Initial density fluid 1
    let rho2:f64 = input[12];        // Initial density fluid 2
    let vx1:f64  = input[13];        // Initial x velocity fluid 1
    let vx2:f64  = input[14];        // Initial x velocity fluid 2
    let p0:f64   = input[15];        // Initial pressure
    
    let nx:u32   = input[20] as u32;        // Particle resolution in the x direction
    let ny:u32   = input[21] as u32;        // Particle resolution in the y direction
    let nz:u32   = input[22] as u32;        // Particle resolution in the z direction
    
    let m:f64    = (rho2 *(y2-y1) + rho1*(lg-y2+y1))*wd*hg; // Total mass  #### Check THIS!!!!! ####
    let mut n:usize  = (3*nx*ny*nz) as usize;        // Total number of particles
    let dm: f64  = m/n as f64;        // Particles' mass
    
    let mut particles :Vec<Particle> = Vec::new();

    partdistribution::init_dist_hcp(&mut particles, nx, ny, nz, rho1, eta, wd, 0.25*lg, hg, x0, y0, z0, dm);
    partdistribution::init_dist_hcp(&mut particles, (1.25*nx as f64) as u32, ny, nz, rho2, eta, wd, 0.5*lg, hg, x0, y0+0.25*lg, z0, dm);
    partdistribution::init_dist_hcp(&mut particles, nx, ny, nz, rho1, eta, wd, 0.25*lg, hg, x0, y0+0.75*lg, z0, dm);

    n = particles.len();

    kh_init_setup(&mut particles, n,  0.25*lg, y1, y2, rho1, rho2, vx1, vx2, p0, gamma-1.);
    
    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn delta_vy(x: f64, y: f64, y1: f64, y2: f64) -> f64 {
    let w0: f64 = 0.1;
    let lamd: f64 = 0.5;
    let sigm: f64 = 0.05;
    let exp1: f64 = (y-y1)/sigm;
    let exp2: f64 = (y-y2)/sigm;
    let vy: f64 = w0*(2.*PI *x / lamd).sin()*((-exp1*exp1).exp()+(-exp2*exp2).exp());
    return vy;
}


fn kh_init_setup(particles: &mut Vec<Particle>, n: usize, lgmid: f64, y1: f64, y2: f64, rho1: f64, rho2: f64, vx1: f64, vx2: f64, p:f64, gamm1: f64) {
    let u1 = p/(gamm1*rho1);
    let u2 = p/(gamm1*rho2);
    let ym: f64 = 0.5*(y1+y2); 
    for ii in 0..n {
        // If Centre Zone
        if (particles[ii].y - ym).abs() <= lgmid {
            particles[ii].vx = vx2;
            particles[ii].u  = u2;
        } else {
            particles[ii].vx = vx1;
            particles[ii].u  = u1;
        }
        particles[ii].vy = delta_vy(particles[ii].x, particles[ii].y, y1, y2);
    }
}