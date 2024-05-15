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
    let eos_t: bool = input[2] != 0.0;  // EoS (0=isoth[No u]; 1=adiab[u])
    let m_star: f64 = input[3];         // Star's mass
    let m_disc: f64 = input[4];         // Disc's mass
    let r_in: f64   = input[5];         // Inner radius
    let r_out: f64  = input[6];         // Outer radius

    let t0: f64     = input[8];         // Initial time
    let tf: f64     = input[9];         // Final time
    let dt_sav: f64 = input[10];        // Recording time step
    
    // Tree's parameters
    let s_: i32     = input[12] as i32; // Bucket size
    let alpha_: f64 = input[13];        // Fraction of the bucket size
    let beta_: f64  = input[14];        // Maximum ratio of cells with less than alpha*s particles
    
    let x_c: f64    = 0.0;
    let y_c: f64    = 0.0;
    let z_c: f64    = 0.0;

    let wd: f64 = 3.0*r_out;            // Width (x)
    let lg: f64 = 3.0*r_out;            // Length (y)
    let hg: f64 = 1.0*r_out;           // Height (z)

    let x0: f64 = x_c - 0.5*wd;         // Bottom left corner  (x-coordinate)
    let y0: f64 = y_c - 0.5*lg;         // Bottom left corner  (y-coordinate)
    let z0: f64 = z_c - 0.5*hg;         // Bottom left corner  (z-coordinate)        

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.1*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    //---------------------------------------------------------------------------------------------

    let hacc: f64 = 0.5*r_in;
    let facc: f64 = 0.8;

    let mut particles :Vec<Particle> = Vec::new();
    let mut star: Star = Star{ m: m_star, x: x_c, y: y_c, z: z_c, hacc:hacc, facc: facc, ..Default::default()};
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());


    let mut t: f64  = t0;               // Time
    let mut n: usize= particles.len();
    let dm:f64      = m_disc/n as f64;
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
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, eos_t, sphfunctions::eos_ideal_gas, sphfunctions::sound_speed_ideal_gas, gamm,
                                       sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern, 
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::mon97_art_vis,
                                       sphfunctions::body_forces_gravitation, &star, true,
                                       sphfunctions::periodic_boundary, xper, yper, zper, wd, lg, hg,  x0, y0, z0);
        sphfunctions::accretion_boundary(&mut star, &mut particles, dm, &mut n, & tree, s_, wd, lg, hg, x0, y0, z0, xper, yper, zper);
        sphfunctions::star_integrator(&mut star, dt);
        dt = sphfunctions::time_step_bale(&particles, n, gamm, rkern, wd, lg, hg, &mut tree, s_, sphfunctions::sound_speed_ideal_gas);
        tree.restart(n);
        datafunctions::time_step(&mut t, &mut dt, dt_sav, &mut sav, &mut it_sav);
        println!("dt: {:.4}\tt: {:.4}\tn:{}", dt, t, n);
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