use std::{
    error::Error,
    f64,
};

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use csv::Writer;

use structures::Particle;

use sphfunctions::h_by_density;

// -------- Write data --------

pub fn init_random_square(path: &str, n: u32, h:f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "rho"])?;
    for _ii in 0..n{
        let x: f64 = wd*rng.gen::<f64>();
        let y: f64 = lg*rng.gen::<f64>();
        let z: f64 = hg*rng.gen::<f64>();
            wtr.write_record(&[(x0 + x).to_string(), (y0 + y).to_string(), (z0+z).to_string(),
                               String::from("0.0"), String::from("0.0"), String::from("0.0"),
                               h.to_string(), String::from("0.0")])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_square(path: &str, n: u32, h:f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "rho"])?;
    let dl: f64 = ((wd*lg*hg)/n as f64).powf(1./3.);
    let nx: u32 = (wd/dl).ceil() as u32;  
    let ny: u32 = (lg/dl).ceil() as u32;  
    let nz: u32 = (hg/dl).ceil() as u32;  
    for kk in 0..nz{
        for jj in 0..ny{
            for ii in 0..nx{
                let x: f64 = x0 + ii as f64*dl;
                let y: f64 = y0 + jj as f64*dl;
                let z: f64 = z0 + kk as f64*dl;
                wtr.write_record(&[(x0 + x).to_string(), (y0 + y).to_string(), (z0+z).to_string(),
                               String::from("0.0"), String::from("0.0"), String::from("0.0"),
                               h.to_string(), String::from("0.0")])?;
            }
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_dist_cubic(particles: &mut Vec<Particle>, nx: u32, ny: u32, nz: u32, rho: f64,
                   eta: f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64, dm: f64){
    let dx: f64 = wd/nx as f64;
    let dy: f64 = lg/ny as f64;
    let dz: f64 = hg/nz as f64;
    let hp: f64 = h_by_density(dm, rho, eta);

    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let xstart: f64 = x0 + 0.5*dx;
    let ystart: f64 = y0 + 0.5*dy;
    let zstart: f64 = z0 + 0.5*dz;

    for kk in 0..nz {
        zp = zstart + dz*kk as f64;
        for jj in 0..ny{
            yp = ystart + dy*jj as f64;
            for ii in 0..nx{
                xp = xstart + dx*ii as f64;
                particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
                                        ..Default::default()});
            }
        }
    }
}