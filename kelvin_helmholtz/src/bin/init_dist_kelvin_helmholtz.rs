// Initialize a set of particles inside a square

use std::{
    error::Error,
    process,
};

use sphfunctions;

use structures::Particle;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    let path = "./Data/initial_distribution/kelvin_helmholtz.csv";

    let mut particles :Vec<Particle> = Vec::new();

    // System's parameters
    let gamma:f64 = 5./3.;  // Gamma factor (heat capacity ratio)
    let eta :f64 = 1.2; // Dimensionless constant related to the ratio of smoothing length
    let d: i32 = 2; // Dimension of the system
    let wd :f64 = 1.; // Domain's width
    let lg :f64 = 1.; // Domain's large
    let hg :f64 = 1.; // Domain's large
    let x0: f64 = 0.; // x-coordinate of the bottom left corner
    let y0: f64 = 0.; // y-coordinate of the bottom left corner
    let z0: f64 = 0.; // y-coordinate of the bottom left corner
    let y1: f64 = 0.25; // Y-lower edge of fluid 2 
    let y2: f64 = 0.75; // Y-upper edge of fluid 2
    let nx: usize = 60; // Number of particles in x direction
    let ny: usize = 100; // Number of particles in region 1 in y direction // TOTAL NY = ny + (rho2/rho1)*ny
    let nz: usize = 40; // Number of particles in z direction

    let rho1: f64 = 1.0; // Initial density fluid 1
    let rho2: f64 = 2.0; // Initial density fluid 2
    let vx1: f64 = -0.5; // Initial x velocity fluid 1
    let vx2: f64 = 0.5; // Initial x velocity fluid 2
    let p0: f64 = 2.5; // Initial pressure

    let m: f64 = (rho1 *(y2-y1) + rho2*(lg-y2+y1))*wd*hg;
    let n: usize = 3*nx*ny*nz;
    let dm: f64 = m / n as f64; // Particles' mass

    kh_init_setup(&mut particles, nx, ny, nz, wd, lg, hg, x0, y0, z0, y1, y2, rho1, rho2, vx1, vx2, p0, gamma, eta, dm, d);

    if let Err(err) = sphfunctions::save_data(&(String::from("./Data/initial_distribution/kelvin_helmholtz.csv")), &particles){
        println!("{}", err);
        process::exit(1);
    }

    Ok(())
}

fn delta_vy(x: f64, y: f64, y1: f64, y2: f64) -> f64 {
    let a: f64 = 0.025;
    let ld: f64 = 1./6.;
    if (y-y1).abs() < a {
        return a*(-2.*PI*(x+0.5)/ld).sin();
    } else if (y-y2).abs() < a {
        return a*(2.*PI*(x+0.5)/ld).sin();
    } else {
        return 0.0;
    }
} 


fn kh_init_setup(particles: &mut Vec<Particle>, nx: usize, ny: usize, nz: usize, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, y1: f64, y2: f64, rho1: f64, rho2: f64, vx1: f64, vx2: f64, p:f64, gamma: f64, eta: f64, dm: f64, d: i32) {
    let dx = wd/nx as f64;
    let dy = 0.5*lg/ny as f64;
    let dy2 = 0.5*dy;
    let dz = hg/nz as f64;
    let h01: f64 = sphfunctions::h_by_density(dm, rho1, eta, d);
    let h02: f64 = sphfunctions::h_by_density(dm, rho2, eta, d);
    for kk in 0..nz{
        for ii in 0..nx {
            for jj in 0..(ny/2) {
                particles.push(Particle{x:x0+(dx*ii as f64), y:y0+(dy*jj as f64), z:z0+(dz*kk as f64),
                                        vx: vx1, vy: delta_vy(x0+(dx*ii as f64), y0+(dy*jj as f64), y1, y2),
                                        h:h01, u: p/((gamma - 1.)*rho1),
                                        ..Default::default()});
                particles.push(Particle{x:x0+(dx*ii as f64), y:y0+y2+(dy*jj as f64), z:z0+(dz*kk as f64),
                                        vx: vx1, vy: delta_vy(x0+(dx*ii as f64), y0+y2+(dy*jj as f64), y1, y2),
                                        h:h01, u: p/((gamma - 1.)*rho1),
                                        ..Default::default()});
            }
            for jj in 0..(2*ny) {
                particles.push(Particle{x:x0+(dx*ii as f64), y:y0+y1+(dy2*jj as f64), z:z0+(dz*kk as f64),
                                        vx: vx2, vy: delta_vy(x0+(dx*ii as f64), y0+y1+(dy2*jj as f64), y1, y2),
                                        h:h02, u: p/((gamma - 1.)*rho2),
                                        ..Default::default()});
            }
        }
    }
}