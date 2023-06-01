use std::{
    error::Error,
    process,
};

use sphfunctions;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {
    let path_source = "./Data/initial_distribution/ToyStar2D.csv";
    let path_result = "./Data/results/ToyStar2D.csv";
    let mut particles :Vec<sphfunctions::Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let sigma :f64 = 1.0/PI; // 
    let d = 3; // Dimension of the system
    let k:f64 = 0.1; // Pressure constant
    let gamma:f64 = 1.0;  // Polytropic index
    let m:f64 = 2.0; // Star's mass
    let r:f64 = 0.75; // Star's radius
    let nu:f64 = 1.0; // Viscocity parameter
    let lmbda = sphfunctions::coeff_static_grav_potential(k, gamma, m, r);
    let dt :f64 = 0.001;
    let t_iter :u32 = 10;

    for tt in 0..t_iter {
        sphfunctions::smoothing_length(&mut particles, eta, sphfunctions::cubic_kernel, sphfunctions::dwdq_cubic_kernel, sigma, d, 1e-04, 100);
        sphfunctions::accelerations(&mut particles, sphfunctions::eos_polytropic, k, gamma, sphfunctions::dwdh, sphfunctions::cubic_kernel, sphfunctions::dwdq_cubic_kernel, nu, lmbda, sigma, d);
        for ii in 0..particles.len(){
            sphfunctions::euler_integrator(&mut particles[ii], dt);
        }
        println!{"{}", tt};
    }
    if let Err(err) = sphfunctions::save_data(path_result, &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}