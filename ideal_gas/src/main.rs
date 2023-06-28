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

use rayon::prelude::*;

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
    let tf:f64 = 0.6; // initial time
    let mut t:f64 = t0; // Time
    let n : usize = particles.len(); // Number of particles

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d: u32 = 2; // Dimension of the system
    let gamma:f64 = 5./3.;  // Polytropic index
    let sigma :f64 = 10.0/(7.*PI); // Normalization's constant of kernel
    let w :f64 = 1.; // width
    let l :f64 = 1.; // large
    let dm :f64 = 1.; // Particles' mass

    for ii in 0..particles.len(){
        particles[ii].vx = -0.1;
    }

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(particles.len() as u32, -0.6, -0.6, 1.2);

    let mut dt :f64 = 0.004;
    let mut it: u32 = 0;
    while t < tf  {
        tree.build_tree(d, s_, alpha_, beta_, &particles, 1.0e-02);
        sphfunctions::smoothing_length(&mut particles, dm, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_, n, particles_ptr);
        sphfunctions::accelerations(&mut particles, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_,n, particles_ptr);
        dt = sphfunctions::time_step_bale(&particles, n, gamma);
        particles.par_iter_mut().for_each(|particle|{
            sphfunctions::euler_integrator(particle, dt);
            sphfunctions::periodic_boundary(particle, w, l, -0.5, -0.5);
            // Initialize variables to zero
            particle.ax = 0.;
            particle.ay = 0.;
            particle.divv = 0.;
            particle.du = 0.;
        });
        tree.restart(n);
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