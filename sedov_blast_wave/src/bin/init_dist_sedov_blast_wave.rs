// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

use structures::Particle;

use datafunctions;

fn main() -> Result<(), Box<dyn Error>> {

    let mut particles :Vec<Particle> = Vec::new();

    let path = "./Data/initial_distribution/sedov_blast_wave.csv";
    let input_file = "./sedov_blast_wave/input";

    let input: Vec<f64> = datafunctions::read_input(input_file);

    let nx:u32 = input[16] as u32; // Number of Particles in x direction
    let ny:u32 = input[17] as u32; // Number of Particles in y direction
    let nz:u32 = input[18] as u32; // Number of Particles in z direction
    let n: u32 = nx*ny*nz; // Total number of particles
    let wd:f64 = input[6]; // Width
    let lg:f64 = input[7]; // Longth
    let hg:f64 = input[8]; // Heigth
    let x0:f64 = input[3]; // Initial x
    let y0:f64 = input[4]; // Initial y
    let z0:f64 = input[5]; // Initial y
    let rho:f64 = input[10]; // Density
    let e0:f64 = input[11]; // Energy
    let mut u_norm:f64 = 0.0; // Normalization parameter
    let mut rad_part: Vec<usize> = Vec::new(); // Particles inside the initial sphere
    let eta: f64 = input [0];
    let rkern: f64 = 2.; // Kernel radius
    let d: i32 = input[1] as i32; // Dimensions
    let vol: f64 = wd*lg*hg; // Volumen
    let h0: f64 = eta*( vol / n as f64).powf(1./d as f64); // Smoothing length
    let dm: f64 = rho*vol/n as f64; // Particles' mass

    init_dist_sedov(&mut particles, &mut rad_part, &mut u_norm, nx, ny, nz, rho, rkern, h0, eta, d, wd, lg, hg, x0, y0, z0, dm, sphfunctions::f_cubic_kernel);
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

fn init_dist_sedov(particles: &mut Vec<Particle>, rad_part: &mut Vec<usize>, u_norm: &mut f64, nx: u32, ny: u32, nz: u32, rho: f64, radius: f64, h0: f64, eta: f64, d: i32, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64, dm: f64, kernel: fn(f64) -> f64){
    let dx: f64 = wd/nx as f64;
    let dy: f64 = lg/ny as f64;
    let dz: f64 = hg/nz as f64;
    let xc: f64 = wd/2.;
    let yc: f64 = lg/2.;
    let zc: f64 = hg/2.;
    let hp: f64 = sphfunctions::h_by_density(dm, rho, eta, d);
    for kk in 0..nz {
        for jj in 0..ny{
            for ii in 0..nx{
                let xp: f64 = x0 + dx*ii as f64;
                let yp: f64 = y0 + dy*jj as f64;
                let zp: f64 = z0 + dz*kk as f64;
                let mut up: f64 = 0.0;
                let r: f64 = ((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc)).sqrt()/ h0;
                if r <= radius {
                    up = kernel(r);
                    *u_norm += up;
                    rad_part.push(((kk*ny + jj)*nx + ii) as usize);
                }
                particles.push(Particle{x:xp, y:yp, z:zp,
                                        vx: 0.0, vy: 0.0, vz: 0.0,
                                        h:hp, u: up,
                                        ..Default::default()});
            }
        }
    }
}