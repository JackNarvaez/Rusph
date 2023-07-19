// Initialize a set of particle to simulate a Toy Star system in 3D.

use std::{
    error::Error,
    process,
};

use csv::Writer;

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/toy_star.csv";
    let d: i32 = 2; // Dimension
    let eta: f64 = 1.2;
    let n: u32 = 1000; // Number of Particles
    let r: f64 = 0.75; // Star's radius
    let m: f64 = 2.0; // Star's mass
    let rho: f64 = 3. * m/(4.*PI*r*r*r); // Density
    let h: f64 = 0.1*eta*(m/(n as f64 * rho)).powf(1./d as f64); // Smoothing length
    let (x0, y0, z0) = (0.0, 0.0, 0.0); // Sphere's center
    let (vx0, vy0, vz0) = (0.0, 0.0, 0.0); // Sphere's velocity
    let u0: f64 = 0.0; // Sphere's energy
    if let Err(err) = init_random_circle(path, n, r, rho, h, x0, y0, z0, vx0, vy0, vz0, u0) {
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

fn init_random_circle(path: &str, n: u32, r:f64, rho:f64, h:f64, x0:f64, y0:f64, z0: f64, vx0:f64, vy0:f64, vz0: f64, u0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "u"])?;
    for _ii in 0..n{
        let r_i = r*(rng.gen::<f64>()).sqrt();
        let theta_i = 2.0*PI*rng.gen::<f64>();
        let phi_i = PI*rng.gen::<f64>();
        let x = r_i*theta_i.cos()*phi_i.sin() + x0;
        let y = r_i*theta_i.sin()*phi_i.sin() + y0;
        let z = r_i*phi_i.cos() + z0;
        wtr.write_record(&[x.to_string(), y.to_string(), z.to_string(),
                           vx0.to_string(), vy0.to_string(), vz0.to_string(),
                           h.to_string(), u0.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}
