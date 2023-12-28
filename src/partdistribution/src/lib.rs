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

pub fn init_dist_hcp(particles: &mut Vec<Particle>, nx: u32, _ny: u32, _nz: u32, rho: f64,
                   eta: f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64, _dm: f64){
    // Based on the Phantom implementation
    let mut dx: f64 = wd/nx as f64;
    let mut dy: f64 = (0.75_f64).sqrt()*dx;
    let mut dz: f64 = (2./3.0_f64).sqrt()*dx;

    let dxp: f64= 0.5*dx;
    let dyp: f64= dy/3.;

    let nxnew: u32 = 2_u32*(((wd/dx)as u32 + 1_u32)/2);
    let nynew: u32 = 2_u32*(((lg/dy)as u32 + 1_u32)/2);
    let nznew: u32 = 2_u32*(((hg/dz)as u32 + 1_u32)/2);

    let nnew: u32 = nxnew*nynew*nznew;

    dx = wd/nxnew as f64;
    dy = lg/nynew as f64;
    dz = hg/nznew as f64;

    let dmp: f64  = rho*wd*lg*hg/nnew as f64;

    let hp: f64 = h_by_density(dmp, rho, eta);

    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let mut xstart: f64;
    let mut ystart: f64;
    let zstart: f64 = z0 + 0.5*dz;

    for kk in 0..nznew {
        zp = zstart + dz*kk as f64;
        for jj in 0..nynew{
            ystart = y0 + 0.5*dyp;
            xstart = x0 + 0.5*dxp;
            if (kk%2) == 0 {
                ystart += dyp;
                if (jj%2) == 1 {
                    xstart += dxp;
                }
            } else {
                if (jj%2) == 0 {
                    xstart += dxp;
                }
            }
            yp = ystart + dy*jj as f64;
            for ii in 0..nxnew{
                xp = xstart + dx*ii as f64;
                particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
                                        ..Default::default()});
            }
        }
    }
}