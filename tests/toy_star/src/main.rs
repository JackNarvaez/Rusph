// The Toy Star Problem in 3D

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
};

use datafunctions;
use sphfunctions;
use spfunc::gamma::*;

use tree_algorithm::BuildTree;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source: &str   = "./Toystar/Ini_00.csv";
    let input_file: &str    = "./tests/toy_star/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length
    let gamm: f64   = input[1];         // Heat capacity ratio
    let nu: f64     = input[2];         // Viscocity parameter
    let m: f64      = input[3];         // Star's mass
    let r: f64      = input[4];         // Star's radius
    
    let x0: f64     = input[5];         // Bottom left corner  (x-coordinate)
    let y0: f64     = input[6];         // Bottom left corner  (y-coordinate)
    let z0: f64     = input[7];         // Bottom left corner  (z-coordinate)
    
    let t0: f64     = input[12];        // Initial time
    let tf: f64     = input[13];        // Final time
    let dt_sav: f64 = input[14];        // Recording time step
    
    // Tree's parameters
    let s_: i32     = input[16] as i32; // Bucket size
    let alpha_: f64 = input[17];        // Fraction of the bucket size
    let beta_: f64  = input[18];        // Maximum ratio of cells with less than alpha*s particles

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.01*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    //---------------------------------------------------------------------------------------------

    // Create Particles
    let mut particles: Vec<Particle> = Vec::new();
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());


    let mut t: f64  = t0;               // Time
    let n: usize    = particles.len();
    let dm: f64     = m/n as f64;       // Particles' mass
    let lmbda: f64  = coeff_static_grav_potential(0.05, gamm, m, r);
    let mut it: u32 = 0;                // Time iterations
    // Save time evolution
    let mut time_file = File::create("./Toystar/Time.txt").expect("creation failed"); // Save time steps
    
    //------------------------------------ kernel -------------------------------------------------
    let sigma: f64  = 1./(120.*PI);     // Normalization constant of kernel
    let rkern: f64  = 3.;               // Kernel radius
    //---------------------------------------------------------------------------------------------

    for ii in 0..n {
        particles[ii].rho = sphfunctions::density_by_smoothing_length(dm, particles[ii].h, eta);
    }

    let mut tree: Node = <Node as BuildTree>::new(n as i32, x0-2.*r, y0-2.*r, z0-2.*r, 4.*r, 4.*r, 4.*r);

    //------------------------------------ Main Loop ----------------------------------------------
    let start = Instant::now();   // Runing time
    while t < tf {
        // In toy star, body forces depend on the particles' velocity.
        // Therefore, it is not straightforward to use the basic LF integrator.
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_polytropic, sphfunctions::sound_speed_ideal_gas, gamm,
                                                 sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern,
                                                 eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                                 sphfunctions::mon97_art_vis,
                                                 sphfunctions::body_forces_toy_star, nu, lmbda, true,
                                                 sphfunctions::periodic_boundary, xper, yper, zper, 4.*r, 4.*r, 4.*r, x0-2.*r, y0-2.*r, z0-2.*r);
        dt = sphfunctions::time_step_mon(&particles, n, gamm, rkern, 4.*r, 4.*r, 4.*r,  x0, y0, z0, &mut tree, s_, sphfunctions::sound_speed_polytropic, xper, yper, zper);
        tree.restart(n);
        datafunctions::time_step(&mut t, dt, dt_sav, &mut sav, &mut it_sav);
        if sav {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Toystar/Ev_") + &(it_sav-2).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Toystar/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// Coefficient of gravital force
fn coeff_static_grav_potential(k:f64, gamm:f64, m:f64, r:f64) -> f64 {
    let gamma_func: f64 = gamma(1.5 + gamm/(gamm-1.))/gamma(gamm/(gamm-1.));
    return 2.*k*(gamm/(gamm-1.))/PI.powf(1.5*(gamm-1.)) * (gamma_func* m/(r*r*r)).powf(gamm-1.)/(r*r);
}