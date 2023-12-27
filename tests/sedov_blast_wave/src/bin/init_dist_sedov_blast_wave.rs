// Initial setting for the Sedov Blast Wave  Problem in 3D

use std::{
    error::Error,
    process,
};

use structures::Particle;

use sphfunctions::{
    f_cubic_kernel,
    h_by_density
};
use datafunctions;

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
    let r0: f64  = input[11];        // Radius that contains e0.
    
    let nx:u32   = input[16] as u32; // Particle resolution in the x direction
    let ny:u32   = input[17] as u32; // Particle resolution in the y direction
    let nz:u32   = input[18] as u32; // Particle resolution in the z direction
    let n: u32   = nx*ny*nz;         // Total number of particles
    
    let rkern: f64 = 2.;             // Kernel radius
    
    let vol: f64 = wd*lg*hg;         // Volumen
    let dm: f64 = rho*vol/n as f64;  // Particles' mass
    
    let mut u_norm:f64 = 0.0;        // Normalization's constant for the initial energy
    let mut rad_part: Vec<usize> = Vec::new(); // Particles inside the initial sphere
    let mut particles :Vec<Particle> = Vec::new();


    init_dist_sedov(&mut particles, &mut rad_part, &mut u_norm, nx, ny, nz, rho, rkern, r0, eta, wd, lg, hg, x0, y0, z0, dm, f_cubic_kernel);
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

fn init_dist_sedov(particles: &mut Vec<Particle>, rad_part: &mut Vec<usize>, u_norm: &mut f64, nx: u32, ny: u32, nz: u32, rho: f64, rkern: f64, r0: f64,
                   eta: f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64, dm: f64, kernel: fn(f64) -> f64){
    let dx: f64 = wd/nx as f64;
    let dy: f64 = lg/ny as f64;
    let dz: f64 = hg/nz as f64;
    let xc: f64 = wd/2.;
    let yc: f64 = lg/2.;
    let zc: f64 = hg/2.;
    let hp: f64 = h_by_density(dm, rho, eta);

    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let xstart: f64 = x0 + 0.5*dx;
    let ystart: f64 = y0 + 0.5*dy;
    let zstart: f64 = z0 + 0.5*dz;

    for kk in 0..nz {
        zp = zstart + dz*kk as f64;
        for jj in 0..ny{
            yp = ystart + dy*jj as f64;
            for ii in 0..nx{
                xp = xstart + dx*ii as f64;
                let mut up: f64 = 0.0;
                let r: f64 = ((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc)).sqrt()/ r0;
                if r <= 1.{
                    up = kernel(r*rkern);
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