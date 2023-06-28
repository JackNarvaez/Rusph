// Be careful with the time of integration.
// This code is still unstable.

use std::{
    error::Error,
    process,
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
    let path_source: &str = "./Data/initial_distribution/ideal_gas.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 0.3; // final time
    let mut t:f64 = t0; // Time
    let n : usize = particles.len(); // Number of particles

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d: i32 = 2; // Dimension of the system
    let gamma:f64 = 5./3.;  // Polytropic index
    let sigma :f64 = 10.0/(7.*PI); // Normalization's constant of kernel
    let w :f64 = 1.; // width
    let l :f64 = 1.; // large
    let x0 = -0.5;
    let y0 = -0.5; 
    let dm :f64 = 1.; // Particles' mass

    for ii in 0..n{
        particles[ii].vx = -0.1;
    }

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as u32, x0-0.1, y0-0.1, 1.2*l);

    let mut dt :f64 = 0.004;
    let mut it: u32 = 0;
    while t < tf  {
        sphfunctions::euler_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                                 sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma,
                                                 d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                                 sphfunctions::body_forces_null, 0.0, 0.0, false,
                                                 sphfunctions::periodic_boundary, w, l, x0, y0);
        dt = sphfunctions::time_step_bale(&particles, n, gamma);
        t += dt;
        println!("{}", t);
        if (it%100)==0 {
            if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/ideal_gas/") + &(it/100).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully. /n Iterations: {}", it);
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/ideal_gas/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}