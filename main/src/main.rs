use std::{
    error::Error,
    process,
};

use sphfunctions;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {
    let path_source = "./Data/initial_distribution/UniformDist.csv";
    let path_result = "./Data/results/UniformDist.csv";
    let mut particles :Vec<sphfunctions::Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path_source, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let it :u32 = 100;
    let eta :f64 = 1.2;
    let sigma :f64 = 1.0/PI;
    let nu = 3;
    let gamma :f64 = 7.0;
    let k :f64 = 1.2;

    let dt :f64 = 0.001;
    let n :u32 = 2;

    for tt in 0..n {
        sphfunctions::smoothing_length(&mut particles, eta, sphfunctions::cubic_kernel, sphfunctions::dwdq_cubic_kernel, sigma, nu, 1e-04, it);
        sphfunctions::accelerations(&mut particles, sphfunctions::eos_polytropic, k, gamma, sphfunctions::dwdh, sphfunctions::cubic_kernel, sphfunctions::dwdq_cubic_kernel, nu, sigma);
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