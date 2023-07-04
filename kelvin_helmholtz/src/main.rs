// Sedov Blast Wave

use std::{
    error::Error,
    process,
    time::Instant,
};

use sphfunctions;

use tree_algorithm::{
    BuildTree,
};

use structures::{
    Particle,
    Node,
    Pointer,
};

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // File's information
    let path_source = "./Data/initial_distribution/kelvin_helmholtz.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // Initial time
    let tf:f64 = 1.0; // Final time
    let mut t:f64 = t0; // Time
    let n : usize = particles.len(); // Number of particles

    // System's parameters
    let eta :f64 = 1.2; // Dimensionless constant related to the ratio of smoothing length
    let d: i32 = 2; // Dimension of the system
    let gamma:f64 = 5./3.;  // Gamma factor (heat capacity ratio)
    let sigma :f64 = 10.0/(7.*PI); // Normalization's constant of kernel
    let w :f64 = 1.; // Domain's width
    let l :f64 = 1.; // Domain's large
    let x0: f64 = 0.; // x-coordinate of the bottom left corner
    let y0: f64 = 0.; // y-coordinate of the bottom left corner
    let rho0: f64 = 1.0; // Initial density
    let dm :f64 = rho0*w*l/n as f64; // Particles' mass
    let p0: f64 = 2.5; // Initial pressure

    kh_conf(&mut particles, n, ramp, p0, gamma);

    // Save initial information
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/initial.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.05;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as u32, x0-0.1*(x0.abs()+0.1), y0-0.1*(x0.abs()+0.1), l+0.2*(x0.abs()+0.1));
    
    let mut dt :f64 = 0.001; // Time step
    let mut it: u32 = 0; // Time iterations
    let it_save: u32 = 10; // Frequency of data saving

    // Main loop
    let start = Instant::now(); // Runing time
    while t < tf  {
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_gaussian_kernel, sphfunctions::dfdq_gaussian_kernel, sigma,
                                       d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon89_art_vis,
                                       sphfunctions::body_forces_null, 0.0, 0.0, false,
                                       sphfunctions::periodic_boundary, w, l, x0, y0);
        dt = sphfunctions::time_step_bale(&particles, n, gamma);
        t += dt;
        println!("{}", t);
        if (it%it_save) == 0 {
            if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/") + &(it/it_save).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);

    // Save final information
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// Setting initial configuration
fn ramp(y: f64, delta: f64) -> f64 {
    let f: f64 = 1./(1.+(2.*(y - 0.25)/delta).exp());
    let g: f64 = 1./(1.+(2.*(0.75 - y)/delta).exp());
    return (1.-f)*(1.-g);
}

fn kh_conf(particles: &mut Vec<Particle>, n: usize, rfn: fn(f64, f64) -> f64, p: f64, gamma: f64) {
    let rho1: f64 = 1.;
    let rho2: f64 = 2.;
    let v1: f64 = -0.5;
    let v2: f64 = 0.5;
    let delta: f64 = 0.25;
    for ii in 0..n{
        particles[ii].rho = rho1 + rfn(particles[ii].y, delta)*(rho2-rho1);
        particles[ii].u = p/((gamma - 1.)*particles[ii].rho);
        particles[ii].vx = v1 + rfn(particles[ii].y, delta)*(v2-v1);
        particles[ii].vy = 0.1*(2.*PI*particles[ii].x).sin();
    }
}