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
    let path_source = "./Data/initial_distribution/sedov_blast_wave.csv";
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // Initial time
    let tf:f64 = 0.1; // Final time
    let mut t:f64 = t0; // Time
    let n : usize = particles.len(); // Number of particles

    // System's parameters
    let eta :f64 = 1.2; // Dimensionless constant related to the ratio of smoothing length
    let d: i32 = 2; // Dimension of the system
    let gamma:f64 = 5./3.;  // Gamma factor (heat capacity ratio)
    let sigma :f64 = 10.0/(7.*PI); // Normalization's constant of kernel
    let w :f64 = 1.; // Domain's width
    let l :f64 = 1.; // Domain's large
    let x0: f64 = -0.5; // x-coordinate of the bottom left corner
    let y0: f64 = -0.5; // y-coordinate of the bottom left corner
    let rho0: f64 = 1.0; // Initial density
    let dm :f64 = rho0*w*l/n as f64; // Particles' mass
    let h0: f64 = 2.*eta*(w*l / n as f64).sqrt(); // Initial radius of Sedov's wave

    sedov_conf(&mut particles, n, h0, 1.0, sphfunctions::f_cubic_kernel, sigma);

    // Save initial information
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/sedov_blast_wave/initial.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.05;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as u32, x0-0.1*x0.abs(), y0-0.1*x0.abs(), l+0.2*x0.abs());
    
    let mut dt :f64 = 0.001; // Time step
    let mut it: u32 = 0; // Time iterations
    let it_save: u32 = 100; // Frequency of data saving

    // Main loop
    let start = Instant::now(); // Runing time
    while t < tf  {
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma,
                                       d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon89_art_vis,
                                       sphfunctions::body_forces_null, 0.0, 0.0, false,
                                       sphfunctions::periodic_boundary, w, l, x0, y0);
        dt = sphfunctions::time_step_bale(&particles, n, gamma);
        t += dt;
        println!("{}", t);
        if (it%it_save) == 0 {
            if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/sedov_blast_wave/") + &(it/it_save).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);

    // Save final information
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/sedov_blast_wave/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// Setting initial configuration
fn sedov_conf(particles: &mut Vec<Particle>, n: usize, h0: f64, radius: f64, kernel: fn(f64) -> f64, sigma: f64) {
    let mut init: Vec<usize> = Vec::new();
    let mut sum: f64 = 0.0;
    for ii in 0..n{
        let r = ((particles[ii].x * particles[ii].x + particles[ii].y * particles[ii].y).sqrt())/h0;
        if r <= radius {
            init.push(ii);
            let energy = sigma*kernel(r);
            particles[ii].u = energy;
            sum += energy;
        } else{
            particles[ii].u = 0.0;
        }
    }
    // Normalize thermal energy
    for ii in init {
        particles[ii].u /= sum;
    }
}