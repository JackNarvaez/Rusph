// The Sod Shock Tube Problem in 3D

use std::{
    fs::File,
    io::Write,
    error::Error,
    process,
    time::Instant,
};

use sphfunctions;
use datafunctions;

use tree_algorithm::BuildTree;

use structures::{
    Particle,
    Node,
    Pointer,
};

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source = "./Sodtube/Ini_00.csv";
    let input_file = "./tests/sod_shock_tube/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta:f64     = input[0];         // Dimensionless constant specifying the smoothing length
    let gamma:f64   = input[1];         // Heat capacity ratio
    let _d:i32      = input[2] as i32;  // Dimensions
    
    let x0:f64      = input[3];         // Bottom left corner  (x-coordinate)
    let y0:f64      = input[4];         // Bottom left corner  (y-coordinate)
    let z0:f64      = input[5];         // Bottom left corner  (z-coordinate)
    let wd:f64      = input[6];         // Width (x)
    let lg:f64      = input[7];         // Length (y)
    let hg:f64      = input[8];         // Heigth (z)
    
    let rhol:f64    = input[9];         // Left density
    let rhor:f64    = input[10];        // Right density
    
    let t0:f64      = input[14];        // Initial time
    let tf:f64      = input[15];        // Final time
    let mut dt:f64  = input[16];        // Initial time step
    let it_save:u32 = input[17] as u32; // Frequency of data saving
    
    // Tree's parameters
    let s_:i32      = input[24] as i32; // Bucket size
    let alpha_:f64  = input[25];        // Fraction of the bucket size
    let beta_:f64   = input[26];        // Maximum ratio of cells with less than alpha*s particles
    
    //---------------------------------------------------------------------------------------------
    
    // Create particles
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr = Pointer(particles.as_mut_ptr());
    

    let mut t:f64   = t0;               // Time
    let n : usize   = particles.len();  // Number of particles
    let dm: f64     = 0.5*(wd*lg*hg)*(rhol+rhor)/n as f64;// Particles' mass
    let mut it: u32 = 0;                // Time iterations

    // Save time evolution
    let mut time_file = File::create("./Sodtube/Time.txt").expect("creation failed");
    
    //------------------------------------ kernel -------------------------------------------------
    let sigma :f64  = 1./(120.*PI);     // Normalization constant of kernel
    let rkern: f64  = 3.;               // Kernel radius
    //---------------------------------------------------------------------------------------------
    
    for ii in 0..n {
        particles[ii].rho = sphfunctions::density_by_smoothing_length(dm, particles[ii].h, eta);
    }
    
    let mut tree : Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);

    //------------------------------------ Main Loop ----------------------------------------------
    let start       = Instant::now();   // Runing time
    while t < tf  {
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamma,
                                       sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern,
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_null, 0.0, 0.0, false,
                                       sphfunctions::periodic_boundary, wd, lg, hg, x0, y0, z0);
        dt = sphfunctions::time_step_mon(&particles, n, gamma, rkern, wd, lg, hg, x0, y0, z0, &mut tree, s_, sphfunctions::sound_speed_ideal_gas_u);
        tree.restart(n);
        t += dt;
        if (it%it_save) == 0 {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Sodtube/Ev_") + &(it/it_save).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Sodtube/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}