use std::{
    error::Error,
    process,
};

use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/initial_distribution/UniformDist.csv";
    let n:u32 = 20;
    let m:f64 = 1.0;
    let rho:f64 = 1.0;
    if let Err(err) = sphfunctions::init_square(path, n, m, rho){
        println!("{}", err);
        process::exit(1);
    }
    Ok(())
}