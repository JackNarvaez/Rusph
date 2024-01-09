use std::f64;

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use structures::Particle;

use sphfunctions::h_by_density;

// -------- Uniform distributions --------

pub fn init_dist_random(
    particles: &mut Vec<Particle>, nx: u32, rho: f64, eta: f64,
    wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    let nnew: u32 = nx*nx*nx;
    let dm  : f64 = rho*wd*lg*hg/nnew as f64;
    let hp  : f64 = h_by_density(dm, rho, eta);
    
    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let mut rng = Pcg64::seed_from_u64(SEED);

    for _ii in 0..nnew {
        xp = x0 + wd*rng.gen::<f64>();
        yp = y0 + lg*rng.gen::<f64>();
        zp = z0 + hg*rng.gen::<f64>();
        particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
            ..Default::default()});
    }
}

pub fn init_dist_cubic(
    particles: &mut Vec<Particle>, nx: u32, rho: f64, eta: f64,
    wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    let dx: f64 = wd/nx as f64;
    let ny: u32 = (lg/dx) as u32;
    let nz: u32 = (hg/dx) as u32;

    let nnew: u32 = nx*ny*nz;
    let dm  : f64 = rho*wd*lg*hg/nnew as f64;
    let hp  : f64 = h_by_density(dm, rho, eta);

    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let xstart: f64 = x0 + 0.5*dx;
    let ystart: f64 = y0 + 0.5*dx;
    let zstart: f64 = z0 + 0.5*dx;

    for kk in 0..nz {
        zp = zstart + dx*kk as f64;
        for jj in 0..ny {
            yp = ystart + dx*jj as f64;
            for ii in 0..nx {
                xp = xstart + dx*ii as f64;
                particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
                                        ..Default::default()});
            }
        }
    }
}

pub fn init_dist_hcp(
    particles: &mut Vec<Particle>, nx: u32, rho: f64, eta: f64,
    wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    let mut dx: f64 = wd/nx as f64;
    let mut dy: f64 = (0.75_f64).sqrt()*dx;
    let mut dz: f64 = (2./3.0_f64).sqrt()*dx;

    let dxp: f64= 0.5*dx;
    let dyp: f64= dy/3.;

    let nxnew: u32 = 2_u32*(((wd/dx)as u32 + 1_u32)/2);
    let nynew: u32 = 2_u32*(((lg/dy)as u32 + 1_u32)/2);
    let nznew: u32 = 2_u32*(((hg/dz)as u32 + 1_u32)/2);

    let nnew : u32 = nxnew*nynew*nznew;

    dx = wd/nxnew as f64;
    dy = lg/nynew as f64;
    dz = hg/nznew as f64;

    let dmp : f64  = rho*wd*lg*hg/nnew as f64;

    let hp: f64 = h_by_density(dmp, rho, eta);

    let mut xp: f64;
    let mut yp: f64;
    let mut zp: f64;

    let mut xstart: f64;
    let mut ystart: f64;
    let zstart: f64 = z0 + 0.5*dz;

    for kk in 0..nznew {
        zp = zstart + dz*kk as f64;
        for jj in 0..nynew {
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
            for ii in 0..nxnew {
                xp = xstart + dx*ii as f64;
                particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
                                        ..Default::default()});
            }
        }
    }
}