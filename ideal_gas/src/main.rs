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
};

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // File's information
    let path_source = "./Data/initial_distribution/ideal_gas.csv";
    let path_result = "./Data/results/ideal_gas.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 0.6; // initial time
    let dt :f64 = 0.004; // time step
    let t_iter :u32 = ((tf-t0)/dt) as u32; // time steps
    let n : usize = particles.len(); // Number of particles 
    println!("{}", t_iter);

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d = 2; // Dimension of the system
    let k:f64 = 0.1; // Pressure constant
    let gamma:f64 = 5./3.;  // Polytropic index
    let sigma :f64 = 10.0/(7.*PI); // 
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
    let mut tree : Node = <Node as BuildTree>::new(particles.len() as u32, -0.1, -0.1, 1.2);

    for tt in 0..t_iter {
        tree.build_tree(d, s_, alpha_, beta_, &particles, 1.0e-02);
        let iter_values: Vec<(f64, f64)> = sphfunctions::smoothing_length(&mut particles, dm, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_, n);
        for ii in 0..n{
            particles[ii].h = iter_values[ii].0;
            particles[ii].rho = iter_values[ii].1;
        }
        sphfunctions::accelerations(&mut particles, dm, sphfunctions::eos_ideal_gas, k, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_,n);
        for ii in 0..particles.len(){
            sphfunctions::euler_integrator(&mut particles[ii], dt);
            sphfunctions::periodic_boundary(&mut particles[ii], w, l);
            // Initialize variables to zero
            particles[ii].ax = 0.;
            particles[ii].ay = 0.;
            particles[ii].dh = 0.;
            particles[ii].du = 0.;
        }
        tree.restart();
        println!{"{}", tt};
    }
    if let Err(err) = sphfunctions::save_data(path_result, &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}