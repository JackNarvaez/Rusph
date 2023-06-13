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
    let path_source = "./Data/initial_distribution/toy_star_2D.csv";
    let path_result = "./Data/results/toy_star_2D.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 1.; // initial time
    let dt :f64 = 0.004; // time step
    let t_iter :u32 = ((tf-t0)/dt) as u32; // time steps
    println!("{}", t_iter);

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d = 2; // Dimension of the system
    let k:f64 = 0.1; // Pressure constant
    let gamma:f64 = 1.0;  // Polytropic index
    let m:f64 = 2.0; // Star's mass
    let r:f64 = 0.75; // Star's radius
    let x0:f64 = 0.; // Star's center
    let y0:f64 = 0.; // Star's center
    let nu:f64 = 1.0; // 1.0Viscocity parameter
    let lmbda = sphfunctions::coeff_static_grav_potential(k, gamma, m, r);
    let sigma :f64 = 10.0/(7.*PI); // 

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(particles.len() as u32, x0-2.*r, y0-2.*r, 4.*r);

    for tt in 0..t_iter {
        tree.build_tree(d, s_, alpha_, beta_, &particles);
        sphfunctions::smoothing_length(&mut particles, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_);
        sphfunctions::accelerations(&mut particles, sphfunctions::eos_polytropic, k, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_);
        for ii in 0..particles.len(){
            sphfunctions::body_forces_toy_star(&mut particles[ii], nu, lmbda);
            sphfunctions::euler_integrator(&mut particles[ii], dt);

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