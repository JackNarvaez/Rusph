// Toy Star system in 2D

use std::{
    fs::File,
    io::Write,
    error::Error,
    time::Instant,
    process,
};

use sphfunctions;

use tree_algorithm::BuildTree;

use structures::{
    Particle,
    Node,
    Pointer,
};

use std::f64::consts::PI;

use datafunctions;

fn main() -> Result<(), Box<dyn Error>> {

    // File's information
    let path_source: &str = "./Data/initial_distribution/toy_star.csv";
    let input_file = "./toy_star/input";
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = input[13]; // Initial time
    let tf:f64 = input[14]; // Final time
    let mut t:f64 = t0; // Time
    let n: usize = particles.len(); // Number of particles 
    let mut time_file = File::create("./Data/results/toy_star/Time.txt").expect("creation failed"); // Save time steps

    // System's parameters.
    let eta :f64 = input[0]; // Dimensionless constant related to the ratio of smoothing length
    let d: i32 = input[1] as i32; // Dimension of the system
    let gamma:f64 = input[2];  // Gamma factor (heat capacity ratio)
    let nu:f64 = input[3]; // Viscocity parameter
    let m:f64 = input[4]; // Star's mass
    let r:f64 = input[5]; // Star's radius
    let (x0, y0, z0) = (input[6], input[7], input[8]); // Sphere's center
    let lmbda: f64 = sphfunctions::coeff_static_grav_potential(0.05, gamma, m, r, d);
    let sigma :f64 = 21./(16.*PI);//1./PI; //  Normalization's constant of kernel
    let rkern: f64 = 2.; // Kernel radius
    let dm:f64 = m/n as f64; // Particles' mass

    // Tree's parameters
    let s_ : i32 = input[18] as i32;
    let alpha_ : f64 = input[19];
    let beta_ : f64 = input[20];
    let mut tree : Node = <Node as BuildTree>::new(n as i32, x0-2.*r, y0-2.*r, z0-2.*r, 4.*r, 4.*r, 4.*r);

    let mut dt :f64 = input[15]; // Time step
    let mut it: u32 = 0; // Time iterations
    let it_save: u32 = input[16] as u32; // Frequency of data saving

    time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
    if let Err(err) = datafunctions::save_data(&(String::from("./Data/results/toy_star/initial.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }

    // Main loop
    let start: Instant = Instant::now();
    while t < tf {
        // In toy star, body forces depend on the particles' velocity.
        // Therefore, it is not straightforward to use the basic LF integrator.
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_polytropic, sphfunctions::sound_speed_ideal_gas, gamma,
                                                 sphfunctions::dwdh, sphfunctions::f_wendland_kernel, sphfunctions::dfdq_wendland_kernel, sigma, rkern,
                                                 d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                                 sphfunctions::mon97_art_vis,
                                                 sphfunctions::body_forces_toy_star, nu, lmbda, true,
                                                 sphfunctions::periodic_boundary, 4.*r, 4.*r, 4.*r, x0-2.*r, y0-2.*r, z0-2.*r);
        dt = sphfunctions::time_step_mon(&particles, n, gamma, rkern, d, 4.*r, 4.*r, 4.*r,  x0, y0, z0, &mut tree, s_, sphfunctions::sound_speed_polytropic);
        tree.restart(n);
        t += dt;
        if (it%it_save) == 0 {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data(&(String::from("./Data/results/toy_star/") + &(it/it_save).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    
    // Save final information
    if let Err(err) = datafunctions::save_data(&(String::from("./Data/results/toy_star/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}