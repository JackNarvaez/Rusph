// Sedov Blast Wave

use std::{
    error::Error,
    process,
    time::Instant,
};

use sphfunctions;

use tree_algorithm::BuildTree;

use structures::{
    Particle,
    Node,
    Pointer,
};

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // File's information
    let path_source = "./Data/initial_distribution/sedov_blast_wave.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let it_tot = 10; // Total iterations
    let n : usize = particles.len(); // Number of particles

    // System's parameters
    let eta :f64 = 1.2; // Dimensionless constant related to the ratio of smoothing length
    let d: i32 = 2; // Dimension of the system
    let gamma:f64 = 5./3.;  // Gamma factor (heat capacity ratio)
    let sigma :f64 = 10.0/(7.*PI); // Normalization's constant of kernel
    let wd :f64 = 1.; // Domain's width
    let lg :f64 = 1.; // Domain's large
    let hg :f64 = 1.; // Domain's large
    let x0: f64 = 0.; // x-coordinate of the bottom left corner
    let y0: f64 = 0.; // y-coordinate of the bottom left corner
    let z0: f64 = 0.; // y-coordinate of the bottom left corner
    let rho0: f64 = 1.0; // Initial density
    let dm :f64 = rho0*wd*lg/n as f64; // Particles' mass
    let rkern: f64 = 2.;

    for ii in 0..n {
        particles[ii].rho = sphfunctions::density_by_smoothing_length(dm, particles[ii].h, eta, d);
    }

    // Tree's parameters
    let s_ : i32 = 10;
    let alpha_ : f64 = 0.05;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, lg);
    
    let mut dt :f64 = 0.001; // Time step
    let mut it: u32 = 0; // Time iterations

    // Main loop
    let start = Instant::now(); // Runing time
    while it < it_tot  {
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, rkern,
                                       d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_null, 0.0, 0.0, false,
                                       sphfunctions::periodic_boundary, wd, lg, hg, x0, y0, z0);
        dt = sphfunctions::time_step_mon(&particles, n, gamma, rkern, d, wd, lg, hg, &mut tree, s_);
        tree.restart(n);
        it += 1;
    }
    println!("{}", start.elapsed().as_secs());
    Ok(())
}