// Initial setting for the Toy Star Problem in 3D.

use std::{
    error::Error,
    process,
};

use csv::Writer;

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use datafunctions;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path = "./Toystar/Ini_00.csv";
    let input_file = "./tests/toy_star/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta:f64  = input[0];         // Dimensionless constant specifying the smoothing length
    let d: i32   = input[2] as i32;  // Dimensions
    let m: f64   = input[4];         // Star's mass
    let r: f64   = input[5];         // Star's radius
    
    let x0:f64   = input[6];         // Star's center (x-coordinate)
    let y0:f64   = input[7];         // Star's center (y-coordinate)
    let z0:f64   = input[8];         // Star's center (z-coordinate)
    let vx0:f64  = input[9];         // Star's velocity (x-coordinate)
    let vy0:f64  = input[10];        // Star's velocity (y-coordinate)
    let vz0:f64  = input[11];        // Star's velocity (z-coordinate)
    let u0:f64   = input[12];        // Initial energy
    
    let n:u32    = input[17] as u32; // Resolution

    let rd:u32   = input[18] as u32; // Random initial distribution

    if rd == 1 {
        let rho:f64  = 3. * m/(4.*PI*r*r*r); // Density
        let nt: u32  = n*n*n;         // Total number of particles
        let h:f64    = eta*(m/(nt as f64 * rho)).powf(1./d as f64); // Smoothing length

        if let Err(err) = init_random_circle(path, nt, r, h, x0, y0, z0, vx0, vy0, vz0, u0) {
            println!("{}", err);
            process::exit(1);
        }
    } else {
        let h:f64    = eta*(2.*r/n as f64); // Smoothing length
        if let Err(err) = init_circle(path, n, r, h, x0, y0, z0, vx0, vy0, vz0, u0) {
            println!("{}", err);
            process::exit(1);
        }

    }
    Ok(())
}

fn init_random_circle(path: &str, n: u32, r:f64, h:f64, x0:f64, y0:f64, z0: f64, vx0:f64, vy0:f64, vz0: f64, u0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "u"])?;
    for _ii in 0..n{
        // let r_i = r*((rng.gen::<f64>()).powf(1./3.)).cbrt();
        let r_i = r*(rng.gen::<f64>()).cbrt();
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

pub fn init_circle(path: &str, n: u32, r:f64, h:f64, x0:f64, y0:f64, z0: f64, vx0:f64, vy0:f64, vz0: f64, u0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "rho"])?;
    let dl: f64 = 2.*r/n as f64; 
    let r2:f64 = r*r;
    for kk in 0..n+1{
        for jj in 0..n+1{
            for ii in 0..n+1{
                let x: f64 = (x0-r) + ii as f64*dl;
                let y: f64 = (y0-r) + jj as f64*dl;
                let z: f64 = (z0-r) + kk as f64*dl;
                if (x*x + y*y + z*z) < r2 {
                    wtr.write_record(&[x.to_string(), y.to_string(), z.to_string(),
                           vx0.to_string(), vy0.to_string(), vz0.to_string(),
                           h.to_string(), u0.to_string()])?;
                }
            }
        }
    }
    wtr.flush()?;
    Ok(())
}
