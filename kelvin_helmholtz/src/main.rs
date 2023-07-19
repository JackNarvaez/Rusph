// Sedov Blast Wave

use std::{
    fs::File,
    io::Write,
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

    let mut particles :Vec<Particle> = Vec::new();
    
    // Simulation's parameters
    let t0:f64 = 0.0; // Initial time
    let tf:f64 = 1.5; // Final time
    let mut t:f64 = t0; // Time
    let mut time_file = File::create("./Data/results/kelvin_helmholtz/Time.txt").expect("creation failed"); // Save time steps
    
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
    let y1: f64 = 0.25; // Y-lower edge of fluid 2 
    let y2: f64 = 0.75; // Y-upper edge of fluid 2
    let nx1: usize = 60; 
    let nx2: usize = 100; 
    let ny1: usize = 40; 
    let ny2: usize = 96; 
    let rho1: f64 = 1.0; // Initial density fluid 1
    let rho2: f64 = 2.0; // Initial density fluid 2
    let vx1: f64 = -0.5; // Initial x velocity fluid 1
    let vx2: f64 = 0.5; // Initial x velocity fluid 2
    let p0: f64 = 2.5; // Initial pressure

    let m: f64 = rho1 *(y2-y1)*w + rho2*(l-y2+y1)*w;
    let dm: f64 = m / (2*nx1*ny1 + nx2*ny2) as f64; // Particles' mass

    kh_init_setup(&mut particles, nx1, nx2, ny1, ny2, w, l, x0, y0, y1, y2, rho1, rho2, vx1, vx2, p0, gamma, eta, dm, d);
    
    let particles_ptr = Pointer(particles.as_mut_ptr());

    let n: usize = particles.len(); // Number of particles

    
    // Save initial information
    time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/initial.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    
    // Tree's parameters
    let s_ : i32 = 10;
    let alpha_ : f64 = 0.05;
    let beta_ : f64 = 0.5;
    let mut tree : Node = <Node as BuildTree>::new(n as i32, x0, y0, l);
    
    let mut dt :f64 = 0.0001; // Time step
    let mut it: u32 = 0; // Time iterations
    let it_save: u32 = 10; // Frequency of data saving
    
    // Main loop
    let start = Instant::now(); // Runing time
    while t < tf  {
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma,
                                       d, eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_null, 0.0, 0.0, false,
                                       sphfunctions::periodic_boundary, w, l, x0, y0);
        dt = sphfunctions::time_step_bale(&particles, n, gamma, d, w, l, &mut tree, s_);
        tree.restart(n);
        t += dt;
        if (it%it_save) == 0 {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/") + &(it/it_save).to_string() + &".csv"), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);

    // Save final information
    time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/results/kelvin_helmholtz/final.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

fn delta_vy(x: f64, y: f64, y1: f64, y2: f64) -> f64 {
    let a: f64 = 0.025;
    let ld: f64 = 1./6.;
    if (y-y1).abs() < a {
        return a*(-2.*PI*(x+0.5)/ld).sin();
    } else if (y-y2).abs() < a {
        return a*(2.*PI*(x+0.5)/ld).sin();
    } else {
        return 0.0;
    }
} 

fn kh_init_setup(particles: &mut Vec<Particle>, nx1: usize, nx2: usize, ny1: usize, ny2: usize, w: f64, l: f64, x0: f64, y0: f64, y1: f64, y2: f64, rho1: f64, rho2: f64, vx1: f64, vx2: f64, p:f64, gamma: f64, eta: f64, dm: f64, d: i32) {
    let dx1 = w/nx1 as f64;
    let dx2 = w/nx2 as f64;
    let dy1 = 0.5*(l-y2+y1)/ny1 as f64;
    let dy2 = (y2-y1)/ny2 as f64;
    let h01: f64 = eta *(dm/rho1).powf(1./d as f64);
    let h02: f64 = eta *(dm/rho2).powf(1./d as f64);
    for ii in 0..nx1 {
        for jj in 0..ny1 {
            particles.push(Particle{x:x0+(dx1*ii as f64), y:y0+(dy1*jj as f64),
                                    h:h01, rho:rho1,
                                    vx: vx1, vy: delta_vy(x0+(dx1*ii as f64), y0+(dy1*jj as f64), y1, y2),
                                    u: p/((gamma - 1.)*rho1),
                                    ..Default::default()});
            particles.push(Particle{x:x0+(dx1*ii as f64), y:y0+y2+(dy1*jj as f64),
                                    h:h01, rho:rho1,
                                    vx: vx1, vy: delta_vy(x0+(dx1*ii as f64), y0+y2+(dy1*jj as f64), y1, y2),
                                    u: p/((gamma - 1.)*rho1),
                                    ..Default::default()});
        }
    }
    for ii in 0..nx2 {
        for jj in 0..ny2 {
            particles.push(Particle{x:x0+(dx2*ii as f64), y:y0+y1+(dy2*jj as f64),
                                    h:h02, rho:rho2,
                                    vx: vx2, vy: delta_vy(x0+(dx2*ii as f64), y0+y1+(dy2*jj as f64), y1, y2),
                                    u: p/((gamma - 1.)*rho2),
                                    ..Default::default()});
        }
    }
}
