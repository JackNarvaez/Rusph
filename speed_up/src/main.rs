// Speedup

use std::{
    error::Error,
    process,
    time::Instant,
};

use sphfunctions;
use datafunctions;

use tree_algorithm::BuildTree;

use structures::{
    Particle,
    Node,
    Pointer,
    Star,
};

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // File
    let path_source = "./Data/initial_distribution/hydro32_00020.csv";

    // Create Particles
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation parameters
    let it_tot = 2; // Total iterations
    let n : usize = particles.len(); // Number of particles

    //---------------------------------------------------------------------------------------------
    // Parameters
    let eta:f64     = 1.2;              // Dimensionless constant related to the ratio of smoothing length
    let wd:f64      = 1.;               // Domain's width
    let lg:f64      = 1.;               // Domain's large
    let hg:f64      = 1.;               // Domain's large
    let x0:f64      = 0.;               // x-coordinate of the bottom left corner
    let y0:f64      = 0.;               // y-coordinate of the bottom left corner
    let z0:f64      = 0.;               // y-coordinate of the bottom left corner
    let gamma:f64   = 1.;               // Gamma factor (heat capacity ratio)
    let dm:f64      = 3.0517578125e-05; // Particles' mass (32K Particles)
    // SBW let gamma:f64   = 5./3.;         // Gamma factor (heat capacity ratio)
    // SBW let rho0:f64    = 1.;             // Initial density
    // SBW let dm:f64      = rho0*wd*lg/n as f64;// Particles' mass
    let rkern:f64   = 2.;               // Kernel radius
    let sigma:f64   = 1./PI;            // Normalization's constant of kernel
    
    let s_:i32      = 10;
    let alpha_:f64  = 0.5;
    let beta_:f64   = 0.5;
    
    let mut dt:f64  = 0.0001;           // Time step
    let mut it:u32  = 0;                // Time iterations

    // Boundary conditions
    let xper:bool   = true;
    let yper:bool   = true;
    let zper:bool   = true;

    //---------------------------------------------------------------------------------------------
    
    for ii in 0..n {
        particles[ii].rho = sphfunctions::density_from_h(dm, particles[ii].h, eta);
    }
    
    let mut tree:Node= <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);

    //------------------------------------ Main Loop ----------------------------------------------
    let start       = Instant::now();   // Runing time
    let star : Star = Star{..Default::default()};
    while it < it_tot  {
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, eos_isothermal, sound_speed_isothermal, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, rkern,
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_null, &star, false,
                                       sphfunctions::periodic_boundary, xper, yper, zper, wd, lg, hg, x0, y0, z0);
        dt = sphfunctions::time_step_mon(&particles, n, gamma, rkern, wd, lg, hg, x0, y0, z0, &mut tree, s_, sound_speed_isothermal, xper, yper, zper);
        tree.restart(n);
        it += 1;
    }
    println!("{}", start.elapsed().as_secs());
    //---------------------------------------------------------------------------------------------

    Ok(())
}

// ------------------------------------------------------------------------- //
// The Polytropic equation of state.                                         //
// Returns the pressure                                                      //
//      P = rho^gamm                                                       //
// where rho is the density, K is the polytropic constant, and gamm is the   //
// adiabatic index.                                                          //
// ------------------------------------------------------------------------- //
fn eos_isothermal(
    rho:f64, _u:f64, gamma:f64, _x: f64, _y: f64, _z: f64
) -> f64 {
    rho.powf(gamma)
}

// ------------------------------------------------------------------------- //
// The Polytropic equation of state.                                         //
// Returns the speed of sound                                                //
//      Cs = sqrt[gamm * rho^(gamm-1)]                                     //
// where rho is the density, K is the polytropic constant, and gamm is the   //
// adiabatic index.                                                          //
// ------------------------------------------------------------------------- //
fn sound_speed_isothermal(
    rho: f64, _u:f64, gamma: f64, _x: f64, _y: f64, _z: f64
) -> f64 {
    (gamma*rho.powf(gamma-1.)).sqrt()
}