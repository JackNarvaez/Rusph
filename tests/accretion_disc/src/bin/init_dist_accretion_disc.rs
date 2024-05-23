// ------------------------------------------------------------------------- //
// Initial Setup for a 3D accretion Disc                                     //
// ------------------------------------------------------------------------- //

use std::{
    error::Error,
    process,
};

use structures::Particle;
use datafunctions;
use partdistribution;

const G: f64    = 1.0;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path: &str      = "./Accretiondisc/Ini_00.csv";
    let input_file: &str= "./tests/accretion_disc/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta: f64    = input[0];         // eta: dimensionless constant specifying the smoothing length
    let gamm: f64   = input[1];         // gamma: Heat capacity ratio
    let x_c: f64    = input[3];         // x_c: center (x-coordinate)
    let y_c: f64    = input[4];         // y_c: center (y-coordinate)
    let z_c: f64    = input[5];         // z_c: center (z-coordinate)
    let r_in: f64   = input[6];         // inner radius of the acc. disc
    let r_ref: f64  = input[7];         // reference radius of the acc. disc
    let r_out: f64  = input[8];         // outer radius of the acc. disc
    let m_dc: f64   = input[9];         // portion of the disc's mass w.r.t. the star mass
    let m_star: f64 = input[10];        // star's mass
    let p_index: f64= input[11];        // p index - density profile 
    let q_index: f64= input[12];        // q index - density profile 
    let h_r: f64    = input[13];        // H over r_ref
    
    let n: u32      = input[17] as u32; // Particle resolution
    
    let m_disc: f64 = m_dc*m_star;      // Disc's mass
    let dm: f64     = m_disc/n as f64;  // Particle's mass
    let vx0: f64    = 0.0;              // x velocity in CoM
    let vy0: f64    = 0.0;              // y velocity in CoM
    let vz0: f64    = 0.0;              // z velocity in CoM
    let nbins: usize= 10000;             // Number of bins for integration

    let m0_disc: f64= partdistribution::disc_mass(r_in, r_out, r_ref, p_index, 1.0, nbins);
    let sigma0: f64 = m_disc/m0_disc;

    let cs0: f64    = h_r*(G*m_star/r_ref).sqrt()*r_ref.powf(q_index);

    let mut particles :Vec<Particle> = Vec::new();

    partdistribution::init_dist_disc1(&mut particles, n, m_star, r_in, r_out, m_disc, p_index, q_index, r_ref, sigma0, cs0, eta, nbins);
    partdistribution::init_dist_disc_velocities(&mut particles, n, m_star, r_in, p_index, q_index, cs0, gamm);
    partdistribution::com_frame(&mut particles, n, dm, x_c, y_c, z_c, vx0, vy0, vz0);

    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}