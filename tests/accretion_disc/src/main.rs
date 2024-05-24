// ------------------------------------------------------------------------- //
// 3D Accretion Disc                                                         //
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

use tree_algorithm::BuildTree;
use std::f64::consts::PI;

const G: f64 = 1.0;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path_source: &str   = "./Accretiondisc/Ini_00.csv";
    let input_file: &str    = "./tests/accretion_disc/input";

    //---------------------------------------------------------------------------------------------
    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta: f64    = input[0];         // eta: dimensionless constant specifying the smoothing length
    let eos_t: bool = input[2] != 0.0;  // EoS (0=isoth[No u]; 1=adiab[u])
    let x_c: f64    = input[3];         // x_c: center (x-coordinate)
    let y_c: f64    = input[4];         // y_c: center (y-coordinate)
    let z_c: f64    = input[5];         // z_c: center (z-coordinate)
    let r_ref: f64  = input[7];         // reference radius of the acc. disc
    let r_out: f64  = input[8];         // outer radius of the acc. disc
    let m_dc: f64   = input[9];         // portion of the disc's mass w.r.t. the star mass
    let m_star: f64 = input[10];         // star's mass
    let q_index: f64= input[12];        // q index - density profile 
    let h_r: f64    = input[13];        // H over r_ref
    
    let t0: f64     = input[14];        // Initial time
    let tf: f64     = input[15];        // Final time
    let dt_sav: f64 = input[16];        // Recording time step
    let mut n: usize= input[17] as usize; // Particle resolution
    
    // Tree's parameters
    let s_: i32     = input[18] as i32; // Bucket size
    let alpha_: f64 = input[19];        // Fraction of the bucket size
    let beta_: f64  = input[20];        // Maximum ratio of cells with less than alpha*s particles
    
    let m_disc: f64 = m_dc*m_star;      // Disc's mass
    let dm: f64     = m_disc/n as f64;  // Particle's mass

    // Boundary conditions
    let xper: bool  = false;
    let yper: bool  = false;
    let zper: bool  = false;

    let mut dt: f64     = 0.01*dt_sav;  // Initial time step
    let mut sav: bool   = false;        // Save data
    let mut it_sav: u32 = 1;            // Save data iteration
    
    //---------------------------------------------------------------------------------------------
    let wd: f64 = 4.0*r_out;
    let lg: f64 = 4.0*r_out;
    let hg: f64 = 20.0*h_r*r_out;

    // Create Particles
    let x0: f64 = x_c - 0.5*wd;
    let y0: f64 = y_c - 0.5*lg;
    let z0: f64 = z_c - 0.5*hg;

    let hacc: f64 = 1.0;
    let facc: f64 = 0.8;

    let cs0: f64    = h_r*(G*m_star/r_ref).sqrt()*r_ref.powf(q_index);
    let cs02: f64   = cs0*cs0;

    let mut particles :Vec<Particle> = Vec::new();
    let mut star: Star = Star{ m: m_star, x: x_c, y: y_c, z: z_c, hacc:hacc, facc: facc, ..Default::default()};
    if let Err(err) = datafunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let particles_ptr: Pointer = Pointer(particles.as_mut_ptr());

    let mut t: f64  = t0;               // Time
    let mut it: u32 = 0;                // Time iterations

    // Save time evolution
    let mut time_file = File::create("./Accretiondisc/Time.txt").expect("creation failed"); // Save time steps
    
    //------------------------------------ kernel -------------------------------------------------
    let sigma: f64  = 1./(120.*PI);     // Normalization constant of kernel
    let rkern: f64  = 3.;               // Kernel radius
    //---------------------------------------------------------------------------------------------

    for ii in 0..n{
        particles[ii].rho = sphfunctions::density_from_h(dm, particles[ii].h, eta);
    }

    let mut tree: Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);

    //------------------------------------ Main Loop ----------------------------------------------
    let start = Instant::now();   // Runing time
    while t < tf {
        sphfunctions::predictor_kdk_integrator(&mut particles, dt, dm, eos_t, sphfunctions::eos_isothermal_disc, sphfunctions::sound_speed_isothermal_disc, q_index, cs02,
                                       sphfunctions::dwdh, sphfunctions::f_quintic_kernel, sphfunctions::dfdq_quintic_kernel, sigma, rkern, 
                                       eta, &mut tree, s_, alpha_, beta_, n, particles_ptr,
                                       sphfunctions::lodatoprice10_art_vis,
                                       sphfunctions::body_forces_gravitation, &star, true,
                                       sphfunctions::none_boundary, xper, yper, zper, wd, lg, hg,  x0, y0, z0);
        sphfunctions::accretion_boundary(&mut star, &mut particles, dm, &mut n, & tree, s_, wd, lg, hg, x0, y0, z0, xper, yper, zper);
        sphfunctions::star_integrator(&mut star, dt);
        tree.restart(n);
        tree.build_tree(s_, alpha_, beta_, &particles, 1.0e-02);
        dt = sphfunctions::time_step_mon(&particles, n, q_index, cs02, rkern, wd, lg, hg, x0, y0, z0, &mut tree, s_, sphfunctions::sound_speed_isothermal_disc, xper, yper, zper);
        tree.restart(n);
        datafunctions::time_step(&mut t, &mut dt, dt_sav, &mut sav, &mut it_sav);
        println!("dt: {:.4}\tt: {:.4}\tn:{}", dt, t, n);
        if sav {
            time_file.write((t.to_string() + &"\n").as_bytes()).expect("write failed");
            if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondisc/Ev_") + &(it_sav-2).to_string()), &particles){
                println!("{}", err);
                process::exit(1);
            }
        } 
        it += 1;
    }
    println!("Simulation run successfully.\n Time {} s.\n Iterations: {}.", start.elapsed().as_secs(), it);
    //---------------------------------------------------------------------------------------------

    // Save final information
    if let Err(err) = datafunctions::save_data_bin(&(String::from("./Accretiondisc/Fin_00")), &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}

// ------------------------------------------------------------------------- //
//Locally isothermal equation of state                                       //
// Returns the speed of sound                                                //
//      Cs = Cs_0 R^(-q)     0.0007905694150411133                                                  //
