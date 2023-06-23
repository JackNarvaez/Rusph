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
    Star,
};

use std::f64::consts::PI;

use rayon::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {

    // Binary stars system

    let star1: Star = Star{m: 4.*PI*PI, x: 0., y: 0., vx: 0., vy: 0., ax: 0., ay:0.};
    let star2: Star = Star{m: 4.*PI*PI*0.08, x: 1., y: 0., vx: 0., vy: 0., ax: 0., ay:0.};

    let p_orb: f64 = 1000.; // Orbital Period [seg]
    let a: f64 = (p_orb*p_orb*(star1.m + star2.m)/(4.*PI*PI)).cbrt(); // Orbital separation;

    let cm_x: f64 = (star1.m*star1.x + star2.m*star2.x)/(star1.m+star2.m); // Center of mass x
    let cm_y: f64 = (star1.m*star1.y + star2.m*star2.y)/(star1.m+star2.m); // Center of mass y

    let r_12: f64 = ((star2.x-star1.x)*(star2.x-star1.x) + (star2.y-star1.y)*(star2.y-star1.y)).sqrt(); // Distance between the two stars
    let ang_12: f64 = (star2.y-star1.y).atan2(star2.x-star1.x); // Angle

    let l1: f64 = r_12*(1.-(star2.m/(3.*(star1.m + star2.m))).cbrt()); // Lagrange point L1

    let roche_x: f64 = l1 * ang_12.cos() + cm_x;
    let roche_y: f64 = l1 * ang_12.sin() + cm_y;

    let w_orb: f64 = ((star1.m + star2.m)/(a*a*a)).sqrt(); // Angular velocity

    // Disk's particles information

    let n : usize = 100; // Number of particles

    // Injection rate

    let n_e: usize = 10;
    let dl_e: f64 = 0.01*a;
    let h_ini: f64 = 0.02*a;

    
    // File's information
    let path_result = "./Data/results/accretion_disk_2D.csv";
    let mut particles :Vec<Particle> = Vec::new();
    let particles_ptr = Pointer(particles.as_mut_ptr());

    // Simulation's parameters
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 1.; // initial time
    let mut t:f64 = t0; // Time

    // System's parameters
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d = 2; // Dimension of the system
    let k:f64 = 0.05; // Pressure constant
    let gamma:f64 = 1.01;  // Polytropic index

    let sigma :f64 = 10.0/(7.*PI); //  Normalization's constant of kernel
    let dm:f64 = star2.m/(4.*PI*PI * n as f64); // mass of each particle

    // Tree's parameters
    let s_ : u32 = 10;
    let alpha_ : f64 = 0.5;
    let beta_ : f64 = 0.5;
    
    // Main loop
    let start = Instant::now();
    let mut dt :f64 = 0.004;
    let mut it: u32 = 0;
    let mut nt: usize = 0;
    println!("{} {}", h_ini, dl_e);
    for ii in 0..5 {
        nt += n_e;
        sphfunctions::inject_particles(&mut particles, roche_x-(ii as f64)*dl_e, roche_y, dl_e, n_e, h_ini);
    }

    let mut tree : Node = <Node as BuildTree>::new(5*n_e as u32, -a, -a, 4.*a);
    while t < tf {
        sphfunctions::inject_particles(&mut particles, roche_x, roche_y, dl_e, n_e, h_ini);
        nt += n_e;
        if let Err(err) = sphfunctions::save_data(path_result, &particles){
            println!("{}", err);
            process::exit(1);
        }
        tree.build_tree(d, s_, alpha_, beta_, &particles, 1.0e-02);
        sphfunctions::smoothing_length(&mut particles, dm, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, 1e-03, 100, dt, &tree, s_, nt, particles_ptr);
        sphfunctions::accelerations(&mut particles, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, k, gamma, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d as i32, &tree, s_, nt, particles_ptr);
        dt = sphfunctions::time_step_bale(&particles, nt, gamma, k);
        println!("{} {}", particles.len(), dt);
        for ii in 0..nt {
            sphfunctions::body_forces_grav_2obj(&mut particles[ii], &star1, &star2, w_orb);
        }
        particles.par_iter_mut().for_each(|particle|{
            //sphfunctions::body_forces_grav_2obj(particle, unsafe {*{star1_ptr}.0.add(0)}, unsafe { *{star2_ptr}.0.add(0)}, w_orb);
            sphfunctions::euler_integrator(particle, dt);
            // Initialize variables to zero
            particle.ax = 0.;
            particle.ay = 0.;
            particle.divv = 0.;
            particle.du = 0.;
        });
        tree.restart(nt);
        t += dt;
        it += 1;
    }
    println!("Simulation run successfully. /n Iterations: {} /n Time: {} s", it, start.elapsed().as_secs());
    if let Err(err) = sphfunctions::save_data(path_result, &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}