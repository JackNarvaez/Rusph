// ------------------------------------------------------------------------- //
// 3D Accretion Disc                                                         //
// ------------------------------------------------------------------------- //

use std::{
    fs::File,
    io::Write,
    error::Error,
    time::Instant,
    process,
};

use structures::{
    Particle,
    Node,
    Pointer,
    Star,
};

use datafunctions;
use sphfunctions;

use tree_algorithm::BuildTree;
use std::f64::consts::PI;
//const G: f64 = 1.0;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source: &str   = "./Accretiondisc/Ini_00.csv";
    let input_file: &str    = "./tests/accretion_disc/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta: f64    = input[0];         // eta: dimensionless constant specifying the smoothing length
    let gamm: f64   = input[1];         // gamma: Heat capacity ratio
    let x_c: f64    = input[2];         // x_c: center (x-coordinate)
    let y_c: f64    = input[3];         // y_c: center (y-coordinate)
    let z_c: f64    = input[4];         // z_c: center (z-coordinate)
    //let r_in: f64   = input[5];         // inner radius of the acc. disc
    let r_out: f64  = input[6];         // outer radius of the acc. disc
    let m_dc: f64   = input[7];         // portion of the disc's mass w.r.t. the star mass
    let m_star: f64 = input[8];         // star's mass
    //let p_index: f64    = input[9];     // p index - density profile 
    //let q_index: f64    = input[10];    // q index - density profile 
    let h_r: f64    = input[11];        // H over r_ref
    
    let t0: f64     = input[12];        // Initial time
    let tf: f64     = input[13];        // Final time
    let dt_sav: f64 = input[14];        // Recording time step
    let n: usize    = input[15] as usize; // Particle resolution
    
    // Tree's parameters
    let s_: i32     = input[16] as i32; // Bucket size
    let alpha_: f64 = input[17];        // Fraction of the bucket size
    let beta_: f64  = input[18];        // Maximum ratio of cells with less than alpha*s particles
    
    let m_disc: f64 = m_dc*m_star;      // Disc's mass
    let dm: f64     = m_disc/n as f64;  // Particle's mass
    //let r_ref: f64  = r_in;             // Reference radius 

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.01*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    //---------------------------------------------------------------------------------------------
    let wd: f64 = 4.0*r_out;
    let lg: f64 = 4.0*r_out;
    let hg: f64 = 20.0*h_r*r_out;

    // Create Particles
    let x0: f64 = x_c - 0.5*wd;
    let y0: f64 = y_c - 0.5*lg;
    let z0: f64 = z_c - 0.5*hg;


    //let cs0: f64    = h_r*(G*m_star/r_ref).sqrt()*r_ref.powf(q_index);
    //let cs02: f64   = cs0*cs0;

    let mut particles :Vec<Particle> = Vec::new();
    let star: Star = Star{ m: m_star, x: x_c, y: y_c, z: z_c, ..Default::default()};
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    let mut t: f64  = t0;               // Time
    let mut it: u32 = 0;                // Time iterations

    // Save time evolution
    let mut time_file = File::create("./Accretiondisc/Time.txt").expect("creation failed"); // Save time steps
    
    //------------------------------------ kernel -------------------------------------------------
    let sigma: f64  = 1./(120.*PI);     // Normalization constant of kernel
    let rkern: f64  = 3.;               // Kernel radius
    //---------------------------------------------------------------------------------------------

    for ii in 0..n{
        particles[ii].rho = sphfunctions::density_from_h(dm, particles[ii].h, eta);
    }

    let mut tree: Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);

    //------------------------------------ Main Loop ----------------------------------------------
    let start = Instant::now();   // Runing time
    while t < tf {
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, eos_isothermal_disc, sound_speed_isothermal_disc, gamm,
                                       sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern, 
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_gravitation, &star, true,
                                       sphfunctions::periodic_boundary, xper, yper, zper, wd, lg, hg,  x0, y0, z0);
        dt = sphfunctions::time_step_bale(&particles, n, gamm, rkern, wd, lg, hg, &mut tree, s_, sound_speed_isothermal_disc);
        tree.restart(n);
        datafunctions::time_step(&mut t, &mut dt, dt_sav, &mut sav, &mut it_sav);
        println!("dt: {}\tt: {}", dt, t);
        if sav {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondisc/Ev_") + &(it_sav-2).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondisc/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// ------------------------------------------------------------------------- //
// Locally isothermal equation of state                                      //
// Returns the presure                                                       //
//      P = Cs^2 * rho                                                       //
// where Cs is the sound speed at R and rho is the density.                  //
// Lodato & Pringle (2007)                                                   //
// ------------------------------------------------------------------------- //
fn eos_isothermal_disc(
    rho:f64, _u:f64, _gamma:f64, x: f64, y: f64, z: f64
) -> f64 {
    0.0007905694150411133 * (x*x + y*y + z*z).powf(-0.25)*rho
}

// ------------------------------------------------------------------------- //
//Locally isothermal equation of state                                       //
// Returns the speed of sound                                                //
//      Cs = Cs_0 R^(-q)                                                     //
// where Cs_o is the sound speed at R_in, R = sqrt[x^2 + y^2 + z^2], and     //
// q is a constant index.                                                    //
// ------------------------------------------------------------------------- //
fn sound_speed_isothermal_disc(
    _rho: f64, _u:f64, _gamma: f64, x: f64, y: f64, z: f64
) -> f64 {
    (0.0007905694150411133 * (x*x + y*y + z*z).powf(-0.25)).sqrt()
}