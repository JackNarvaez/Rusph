// ------------------------------------------------------------------------- //
// Initial setup for a thin disc using uniform density profile               //
// ------------------------------------------------------------------------- //

use std::{
    error::Error,
    process,
};

use structures::Particle;
use sphfunctions::h_from_density;
use datafunctions;
use partdistribution;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path: &str      = "./Accretiondiscuniform/Ini_00.csv";
    let input_file: &str= "./tests/accretion_disc_uniform/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length   
    let m_star: f64 = input[3];         // Star's mass
    let m_disc: f64 = input[4];         // Disc's mass
    let r_in: f64   = input[5];        // Inner radius
    let r_out: f64  = input[6];        // Outer radius
    let h_disc: f64 = input[7];        // Hight
    let nx: u32     = input[11] as u32; // Particle resolution
    
    let r2_in: f64  = r_in*r_in;
    let r2_out: f64 = r_out*r_out;

    let x_c: f64    = 0.0;
    let y_c: f64    = 0.0;
    let z_c: f64    = 0.0;

    let wd: f64 = 2.0*r_out;            // Width (x)
    let lg: f64 = 2.0*r_out;            // Length (y)

    let x0: f64 = x_c - r_out;          // Bottom left corner  (x-coordinate)
    let y0: f64 = y_c - r_out;          // Bottom left corner  (y-coordinate)
    let z0: f64 = z_c - 0.5*h_disc;         // Bottom left corner  (z-coordinate)        

    let rho: f64    = m_disc/(PI*h_disc*(r_out-r_in)*(r_out+r_in));
    
    let mut particles :Vec<Particle> = Vec::new();    
    partdistribution::init_dist_hcp(&mut particles, nx, rho, eta, wd, lg, h_disc, x0, y0, z0);
    particles.retain(|particle| distance_center(&particle, x_c, y_c, r2_in, r2_out));

    let n: usize    = particles.len();
    let dm: f64     = m_disc/n as f64;
    let h: f64      = h_from_density(dm, rho, eta);
    let m_t: f64    = m_star + dm;
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