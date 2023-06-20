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
    let path_source = "./Data/initial_distribution/toy_star_2D.csv";
    let path_result = "./Data/results/toy_star_2D.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 10.; // initial time
    let dt :f64 = 0.004; // time step
    let t_iter :u32 = ((tf-t0)/dt) as u32; // time iterations
    let n : usize = particles.len(); // Number of particles 
    println!("{}", t_iter);

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d = 2; // Dimension of the system
    let k:f64 = 0.05; // Pressure constant
    let gamma:f64 = 1.0;  // Polytropic index
    let m:f64 = 2.0; // Star's mass
    let r:f64 = 0.75; // Star's radius
    let x0:f64 = 0.; // Star's center
    let y0:f64 = 0.; // Star's center
    let nu:f64 = 1.0; // Viscocity parameter
    let lmbda = sphfunctions::coeff_static_grav_potential(k, gamma, m, r);
    let sigma :f64 = 10.0/(7.*PI); //  Normalization's constant of kernel
    let dm:f64 = m/n as f64; // mass of each particle

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as u32, x0-2.*r, y0-2.*r, 4.*r);

    // Main loop
    let start = Instant::now();
    for tt in 0..t_iter {
        tree.build_tree(d, s_, alpha_, beta_, &particles, 1.0e-02);

        sphfunctions::smoothing_length(&mut particles, dm, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_, n, particles_ptr);

        sphfunctions::accelerations(&mut particles, dm, sphfunctions::eos_polytropic, k, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_, n, particles_ptr);

        particles.par_iter_mut().for_each(|particle|{
            sphfunctions::body_forces_toy_star(particle, nu, lmbda);
            sphfunctions::euler_integrator(particle, dt);
            // Initialize variables to zero
            particle.ax = 0.;
            particle.ay = 0.;
            particle.dh = 0.;
            particle.du = 0.;
        });
        tree.restart();
        println!{"{}", tt};
    }
    println!("{} s", start.elapsed().as_secs());
    if let Err(err) = sphfunctions::save_data(path_result, &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}