// ------------------------------------------------------------------------- //
// The Toy Star Problem in 3D                                                //
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
use spfunc::gamma::*;

use tree_algorithm::BuildTree;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source: &str   = "./Toystar/Ini_00.csv";
    let input_file: &str    = "./tests/toy_star/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);
    
    let eta: f64    = input[0];         // Dimensionless constant specifying the smoothing length
    let gamm: f64   = input[1];         // Heat capacity ratio
    let k: f64      = input[2];         // Constant coefficient [EoS]
    let eos_t: bool = input[3] != 0.0;  // EoS (0=isoth[No u]; 1=adiab[u])
    let nu: f64     = input[4];         // Viscocity parameter
    let m_star: f64 = input[5];         // Star's mass
    let r: f64      = input[6];         // Star's radius
    
    let x_c: f64    = input[7];         // Bottom left corner  (x-coordinate)
    let y_c: f64    = input[8];         // Bottom left corner  (y-coordinate)
    let z_c: f64    = input[9];         // Bottom left corner  (z-coordinate)
    
    let t0: f64     = input[14];        // Initial time
    let tf: f64     = input[15];        // Final time
    let dt_sav: f64 = input[16];        // Recording time step
    
    // Tree's parameters
    let s_: i32     = input[18] as i32; // Bucket size
    let alpha_: f64 = input[19];        // Fraction of the bucket size
    let beta_: f64  = input[20];        // Maximum ratio of cells with less than alpha*s particles

    let mut wd: f64 = 3.0*r;            // Bottom left corner  (x-coordinate)
    let mut lg: f64 = 3.0*r;            // Bottom left corner  (y-coordinate)
    let mut hg: f64 = 3.0*r;            // Bottom left corner  (z-coordinate)

    let mut x0: f64 = x_c - 0.5*wd;
    let mut y0: f64 = y_c - 0.5*lg;
    let mut z0: f64 = z_c - 0.5*hg;

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.01*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    let lmbda: f64  = coeff_static_grav_potential(k, gamm, m_star, r);
    //---------------------------------------------------------------------------------------------
    
    // Create Particles
    let mut particles: Vec<Particle> = Vec::new();
    
    let star : Star = Star{m: m_star, x: x_c, y: y_c, z: z_c, hacc: nu, facc: lmbda, ..Default::default()};
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    let mut t: f64  = t0;               // Time
    let n: usize    = particles.len();
    let dm: f64     = star.m/n as f64;  // Particles' mass
    let mut it: u32 = 0;                // Time iterations
    // Save time evolution
    let mut time_file = File::create("./Toystar/Time.txt").expect("creation failed"); // Save time steps
    
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
        sphfunctions::velocity_verlet_integrator(&mut particles, dt, dm, eos_t, sphfunctions::eos_polytropic, sphfunctions::sound_speed_polytropic, gamm, k,
                                                 sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern,
                                                 eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                                 sphfunctions::mon97_art_vis,
                                                 sphfunctions::body_forces_toy_star, &star, true,
                                                 sphfunctions::none_boundary, xper, yper, zper, wd, lg, hg, x0, y0, z0);
        dt = sphfunctions::time_step_mon(&particles, n, gamm, k, rkern, wd, lg, hg,  x0, y0, z0, &mut tree, s_, sphfunctions::sound_speed_polytropic, xper, yper, zper);
        tree.restart(n);
        sphfunctions::open_boundary(&particles, &mut wd, &mut lg, &mut hg, &mut x0, &mut y0, &mut z0);
        datafunctions::time_step(&mut t, &mut dt, dt_sav, &mut sav, &mut it_sav);
        println!("dt: {:.4}\tt: {:.4}", dt, t);
        if sav {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Toystar/Ev_") + &(it_sav-2).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        }
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Toystar/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// Coefficient of gravital force
fn coeff_static_grav_potential(k:f64, gamm:f64, m:f64, r:f64) -> f64 {
    let gamma_func: f64 = gamma(1.5 + gamm/(gamm-1.))/gamma(gamm/(gamm-1.));
    return 2.*k*(gamm/(gamm-1.))/PI.powf(1.5*(gamm-1.)) * (gamma_func* m/(r*r*r)).powf(gamm-1.)/(r*r);
}