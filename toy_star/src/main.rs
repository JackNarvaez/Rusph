// Toy Star system in 2D

use std::{
    error::Error,
    time::Instant,
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
    let path_source: &str = "./Data/initial_distribution/toy_star_2D.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 15.; // initial time
    let mut t:f64 = t0; // Time
    let n: usize = particles.len(); // Number of particles 

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d: u32 = 2; // Dimension of the system
    let gamma:f64 = 2.0;  // Polytropic index
    let m:f64 = 2.0; // Star's mass
    let r:f64 = 0.75; // Star's radius
    let x0:f64 = 0.; // Star's center
    let y0:f64 = 0.; // Star's center
    let nu:f64 = 1.0; // Viscocity parameter
    let lmbda: f64 = sphfunctions::coeff_static_grav_potential(0.25, gamma, m, r);
    let sigma :f64 = 10.0/(7.*PI); //  Normalization's constant of kernel
    let dm:f64 = m/n as f64; // mass of each particle

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as u32, x0-2.*r, y0-2.*r, 4.*r);

    // Main loop
    let start: Instant = Instant::now();
    let mut dt :f64 = 0.004;
    let mut it: u32 = 0;
    while t < tf {
        tree.build_tree(d, s_, alpha_, beta_, &particles, 1.0e-02);
        sphfunctions::smoothing_length(&mut particles, dm, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_, n, particles_ptr);
        sphfunctions::accelerations(&mut particles, dm, sphfunctions::eos_polytropic, sphfunctions::sound_speed_ideal_gas, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_, n, particles_ptr);
        dt = sphfunctions::time_step_bale(&particles, n, gamma);
        particles.par_iter_mut().for_each(|particle|{
            sphfunctions::body_forces_toy_star(particle, nu, lmbda);
            sphfunctions::euler_integrator(particle, dt);
            // Initialize variables to zero
            particle.ax = 0.;
            particle.ay = 0.;
            particle.divv = 0.;
            particle.du = 0.;
        });
        tree.restart(n);
        t += dt;
        println!("{}", t);
        if (it%100) == 0 {
            if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/toy_star/") + &(it/100).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully. /n Iterations: {} /n Time: {} s", it, start.elapsed().as_secs());
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/toy_star/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}