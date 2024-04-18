// Initial setting for the Toy Star Problem in 3D.

use std::{
    error::Error,
    process,
};

use structures::Particle;
use sphfunctions::h_by_density;
use datafunctions;
use partdistribution;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path: &str      = "./Accretiondisc/Ini_00.csv";
    let input_file: &str= "./tests/accretion_disc/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length
    let x0: f64     = input[2];         // Bottom left corner  (x-coordinate)
    let y0: f64     = input[3];         // Bottom left corner  (y-coordinate)
    let z0: f64     = input[4];         // Bottom left corner  (z-coordinate)
    let wd: f64     = input[5];         // Width (x)
    let lg: f64     = input[6];         // Length (y)
    let hg: f64     = input[7];         // Height (z)
    
    let rho: f64    = input[8];         // Density
    let m: f64      = input[9];         // Star's mass
    let r_disc: f64 = input[10];
    let nx: u32     = input[14] as u32; // Particle resolution
    
    let r_in: f64   = 0.5*r_disc;
    let r2_in: f64  = r_in*r_in;
    let r2_out: f64 = r_disc*r_disc;

    let x_c: f64 = x0 + 0.5*wd;
    let y_c: f64 = y0 + 0.5*lg;

    
    let mut particles :Vec<Particle> = Vec::new();    
    partdistribution::init_dist_hcp(&mut particles, nx, rho, eta, wd, lg, hg, x0, y0, z0);
    particles.retain(|particle| distance_center(&particle, x_c, y_c, r2_in, r2_out));

    let n: usize    = particles.len();
    let dm:f64      = PI*rho*hg*(r_disc-r_in)*(r_disc+r_in)/n as f64;
    let h: f64      = h_by_density(dm, rho, eta);
    let m_t: f64    = m + dm;
    for ii in 0..n {
        let omega: f64 = keplerian_velocity(&particles[ii], m_t, x_c, y_c);
        particles[ii].h  = h;
        particles[ii].vx = -particles[ii].y * omega;
        particles[ii].vy = particles[ii].x * omega;
        particles[ii].vz = 0.0;
        particles[ii].u  = 0.0;
    }

    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn distance_center(particle: &Particle, x0: f64, y0: f64, r2_in: f64, r2_out: f64) -> bool {
    let xtem: f64 = particle.x - x0;
    let ytem: f64 = particle.y - y0;
    let r2: f64  = xtem*xtem + ytem*ytem;
    return (r2_in < r2) && (r2 < r2_out);
}

fn keplerian_velocity(particle: &Particle, m_t: f64, x0: f64, y0: f64) -> f64 {
    let xtem: f64 = particle.x - x0;
    let ytem: f64 = particle.y - y0;
    let r2: f64   = xtem*xtem + ytem*ytem;
    return (m_t/(r2).powf(1.5)).sqrt();
}