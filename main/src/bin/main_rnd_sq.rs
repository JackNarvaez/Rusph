use std::{
    error::Error,
    process,
};

use sphfunctions;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {
    let path_source = "./Data/initial_distribution/Rand_Square.csv";
    let path_result = "./Data/results/Rand_Square.csv";
    let mut particles :Vec<sphfunctions::Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let t0:f64 = 0.0; // initial time
    let tf:f64 = 1.; // initial time
    let dt :f64 = 0.004; // time step
    let t_iter :u32 = ((tf-t0)/dt) as u32; // time steps
    println!("{}", t_iter);
    let eta :f64 = 1.2; // dimensionless constant related to the ratio of smoothing length
    let d = 2; // Dimension of the system
    let k:f64 = 0.1; // Pressure constant
    let gamma:f64 = 5./3.;  // Polytropic index
    let sigma :f64 = 10.0/(7.*PI); // 

    let w :f64 = 1.; // width
    let l :f64 = 1.; // large

    for ii in 0..particles.len(){
        particles[ii].vx = -0.1;
    }

    for tt in 0..t_iter {
        sphfunctions::smoothing_length(&mut particles, eta, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d, 1e-03, 100, dt);
        sphfunctions::accelerations(&mut particles, sphfunctions::eos_ideal_gas, k, gamma, sphfunctions::dwdh, sphfunctions::f_cubic_kernel, sphfunctions::dfdq_cubic_kernel, sigma, d);
        for ii in 0..particles.len(){
            sphfunctions::euler_integrator(&mut particles[ii], dt);
            sphfunctions::periodic_boundary(&mut particles[ii], w, l);
            // Initialize variables to zero
            particles[ii].ax = 0.;
            particles[ii].ay = 0.;
            particles[ii].dh = 0.;
            particles[ii].du = 0.;
        }
        println!{"{}", tt};
    }
    if let Err(err) = sphfunctions::save_data(path_result, &particles){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}