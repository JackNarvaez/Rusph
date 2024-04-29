// ------------------------------------------------------------------------- //
// A thin accretion disc with uniform density profiles under the             //
// gravitational force of a central massive star                             //
// ------------------------------------------------------------------------- //

use std::{
    fs::File,
    io::Write,
    error::Error,
    time::Instant,
    process,
};

use structures::{
    Particle,
    Node,
    Pointer,
    Star,
};

use datafunctions;
use sphfunctions;

use tree_algorithm::BuildTree;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source: &str   = "./Accretiondiscuniform/Ini_00.csv";
    let input_file: &str    = "./tests/accretion_disc_uniform/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length
    let gamm: f64   = input[1];         // Heat capacity ratio

    let x0: f64     = input[2];         // Bottom left corner  (x-coordinate)
    let y0: f64     = input[3];         // Bottom left corner  (y-coordinate)
    let z0: f64     = input[4];         // Bottom left corner  (z-coordinate)
    let wd: f64     = input[5];         // Width (x)
    let lg: f64     = input[6];         // Length (y)
    let hg: f64     = input[7];         // Height (z)
    
    let rho: f64    = input[8];         // Density
    let m: f64      = input[9];         // Star's mass
    let r_disc: f64 = input[10];
    let r_in: f64   = 0.5*r_disc;
    
    let t0: f64     = input[11];        // Initial time
    let tf: f64     = input[12];        // Final time
    let dt_sav: f64 = input[13];        // Recording time step
    
    // Tree's parameters
    let s_: i32     = input[15] as i32; // Bucket size
    let alpha_: f64 = input[16];        // Fraction of the bucket size
    let beta_: f64  = input[17];        // Maximum ratio of cells with less than alpha*s particles

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.1*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    //---------------------------------------------------------------------------------------------

    // Create Particles
    let x_c: f64 = x0 + 0.5*wd;
    let y_c: f64 = y0 + 0.5*lg;
    let z_c: f64 = z0 + 0.5*hg;

    
    let mut particles :Vec<Particle> = Vec::new();
    let star: Star = Star{ m: m, x: x_c, y: y_c, z: z_c, ..Default::default()};
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());


    let mut t: f64  = t0;               // Time
    let n: usize    = particles.len();
    let dm:f64      = PI*rho*hg*(r_disc-r_in)*(r_disc+r_in)/n as f64;
    //let lmbda: f64  = coeff_static_grav_potential(0.05, gamm, m, r);
    let mut it: u32 = 0;                // Time iterations
    // Save time evolution
    let mut time_file = File::create("./Accretiondiscuniform/Time.txt").expect("creation failed"); // Save time steps
    
    //------------------------------------ kernel -------------------------------------------------
    let sigma: f64  = 1./(120.*PI);     // Normalization constant of kernel
    let rkern: f64  = 3.;               // Kernel radius
    //---------------------------------------------------------------------------------------------

    for ii in 0..n {
        particles[ii].rho = sphfunctions::density_from_h(dm, particles[ii].h, eta);
    }

    let mut tree: Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);

    //------------------------------------ Main Loop ----------------------------------------------
    let start = Instant::now();   // Runing time
    while t < tf {
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamm,
                                       sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern, 
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_gravitation, &star, true,
                                       sphfunctions::periodic_boundary, xper, yper, zper, wd, lg, hg,  x0, y0, z0);
        dt = sphfunctions::time_step_bale(&particles, n, gamm, rkern, wd, lg, hg, &mut tree, s_, sphfunctions::sound_speed_ideal_gas);
        tree.restart(n);
        datafunctions::time_step(&mut t, dt, dt_sav, &mut sav, &mut it_sav);
        if sav {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondiscuniform/Ev_") + &(it_sav-2).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondiscuniform/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}
