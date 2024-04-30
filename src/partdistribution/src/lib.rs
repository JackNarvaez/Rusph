use std::f64;

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;
const G: f64    = 1.0;

use structures::Particle;

use sphfunctions::h_from_density;

use std::f64::consts::PI;

// -------- Uniform distributions --------

pub fn init_dist_random(
    particles: &mut Vec<Particle>, nx: u32, rho: f64, eta: f64,
    wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    let nnew: u32 = nx*nx*nx;
    let dm  : f64 = rho*wd*lg*hg/nnew as f64;
    let hp  : f64 = h_from_density(dm, rho, eta);
    
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
    let hp  : f64 = h_from_density(dm, rho, eta);

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

    let hp: f64 = h_from_density(dmp, rho, eta);

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

// -------- Accretion disc distributions --------

pub fn sigma_profile(
    r: f64, p: f64, r_ref:f64, r_in: f64
) -> f64 {
    (r/r_ref).powf(-p)*(1.0 -(r_in/r).sqrt())
}

pub fn cs_disc(
    r:f64, r_in: f64, q: f64, cs0: f64
) -> f64 {
    cs0*(r/r_in).powf(-q)
}

pub fn disc_mass(
    r_in: f64, r_out: f64, r_ref: f64,
    p_index:f64, sigm0: f64, nbins: usize
) -> f64 {
    let mut r   : f64;
    let mut sigm: f64;
    let mut dm  : f64;
    let r_wd    : f64 = r_out - r_in;
    let dr      : f64 = r_wd/nbins as f64;
    let mut m_disc: f64 = 0.0;
    for ii in 0..nbins {
        r = r_in + dr*ii as f64;
        sigm = sigm0*sigma_profile(r, p_index, r_ref, r_in);
        dm = 2.0*PI*r*sigm*dr;
        m_disc += dm;
    }
    return m_disc;
}

pub fn init_dist_disc1(
    particles: &mut Vec<Particle>, n: u32, m_star: f64, r_in: f64, r_out: f64, m_disc: f64,
    p_index: f64, q_index: f64, r_ref: f64, sigma0: f64, cs0: f64,
    eta: f64, nbins: usize, x_c: f64, y_c: f64, z_c: f64
) {
    let mut r   : f64 = 0.0;
    let mut z   : f64 = 0.0;
    let mut phi : f64;
    let mut xp  : f64;
    let mut yp  : f64;
    let mut zp  : f64;
    let mut zmin: f64;
    let mut zmax: f64;
    let mut z_wd: f64;
    let mut hp  : f64;
    let mut f   : f64;
    let mut f_ii: f64;
    let mut omeg: f64;
    let mut cs  : f64;
    let mut h2  : f64;
    let mut s2h2: f64;
    let mut fz_max: f64;
    
    let mut f_max   : f64 = 0.0;
    let mut sigm: f64 = 0.0;
    let mut rho : f64 = 0.0;
    let r_wd    : f64 = r_out - r_in;
    let dr      : f64 = r_wd/nbins as f64;
    let dm      : f64 = m_disc/n as f64;

    // Calculate f_max = r*SIGMA(r)
    for ii in 0..nbins {
        r = r_in + dr*ii as f64;
        f_ii = r * sigma0 * sigma_profile(r, p_index, r_ref, r_in);
        if f_ii > f_max {
            f_max = f_ii;
        }
    }

    let mut rng = Pcg64::seed_from_u64(SEED);

    for _ii in 0..n {
        // find r
        // the rejection method
        f = 0.0;
        f_ii = 1.0;
        while f_ii > f {
            r   = r_in + r_wd*rng.gen::<f64>();             // r
            f_ii= f_max*rng.gen::<f64>();                   // u2
            f   = r*sigma0*sigma_profile(r, p_index, r_ref, r_in); // f = r*SIGMA(r)
            sigm= f/r;
        }
        // find z
        // the rejection method
        cs      = cs_disc(r, r_in, q_index,cs0);
        omeg    = (G*m_star/r.powi(3)).sqrt();
        h2      = cs/omeg;
        s2h2    = 2.0_f64.sqrt()*h2;

        zmin    = -3.0_f64.sqrt()*s2h2;
        zmax    = 3.0_f64.sqrt()*s2h2;
        z_wd    = zmax - zmin;

        fz_max  = sigm/(s2h2 * PI.sqrt());
        f = 0.0;
        f_ii = 1.0;
        while f_ii > f {
            z   = zmin + z_wd*rng.gen::<f64>();         // z
            f_ii= fz_max*rng.gen::<f64>();               // u3
            f = sigm*(-(z/(s2h2)).powi(2)).exp()/(s2h2 * PI.sqrt());  // rho(x, y, z)
            rho= f;
        }

        // Calculate cartesian coordinates
        phi = 2.0*PI*rng.gen::<f64>();
        xp = x_c + r*phi.cos();
        yp = y_c + r*phi.sin();
        zp = z_c + z;
        hp  = h_from_density(dm, rho, eta);
        particles.push(Particle{x:xp, y:yp, z:zp, h:hp,
            ..Default::default()});
    }   
}

pub fn init_dist_disc_velocities(
    particles: &mut Vec<Particle>, n: u32, m_star: f64, r_in: f64,
    p_index: f64, q_index: f64, cs0: f64, gamm: f64
) {
    let mut kepl: f64;
    let mut r   : f64;
    let mut phi : f64;
    let mut cs  : f64;
    let mut cs2 : f64;
    let mut f_p : f64;
    let mut vphi: f64;
    
    for ii in 0..n as usize {
        r = (particles[ii].x * particles[ii].x + particles[ii].y*particles[ii].y).sqrt();
        kepl = G*m_star/r;
        phi = (particles[ii].y).atan2(particles[ii].x);
        cs = cs_disc(r, r_in, q_index, cs0);
        cs2 = cs*cs;
        f_p = -cs2*(1.5+p_index+q_index);
        vphi = (kepl + f_p).sqrt();
        particles[ii].vx = -vphi* phi.sin();
        particles[ii].vy = vphi* phi.cos();
        particles[ii].vz = 0.0_f64;
        particles[ii].u  = cs2/((gamm-1.)*gamm);
    }
}

pub fn com_frame(
    particles: &mut Vec<Particle>, n: u32, dm: f64,
    x_c: f64, y_c:f64, z_c:f64, vx0: f64, vy0: f64, vz0: f64
) {
    let mut xcm     : f64 = 0.0;
    let mut ycm     : f64 = 0.0;
    let mut zcm     : f64 = 0.0;
    let mut vxcm    : f64 = 0.0;
    let mut vycm    : f64 = 0.0;
    let mut vzcm    : f64 = 0.0;
    let mut totmass : f64 = 0.0;
    let invtmass: f64;

    for ii in 0..n as usize {
        totmass += dm;
        xcm += particles[ii].x * dm;
        ycm += particles[ii].y * dm;
        zcm += particles[ii].z * dm;
        vxcm += particles[ii].vx * dm;
        vycm += particles[ii].vy * dm;
        vzcm += particles[ii].vz * dm;
    }

    invtmass= 1./totmass;
    xcm = xcm*invtmass;
    ycm = ycm*invtmass;
    zcm = zcm*invtmass;
    vxcm = vxcm*invtmass;
    vycm = vycm*invtmass;
    vzcm = vzcm*invtmass;

    for ii in 0..n as usize {
        particles[ii].x += x_c - xcm;  
        particles[ii].y += y_c - ycm;  
        particles[ii].z += z_c - zcm;  
        particles[ii].vx += vx0 - vxcm;  
        particles[ii].vy += vy0 - vycm;  
        particles[ii].vz += vz0 - vzcm;  
    }
}