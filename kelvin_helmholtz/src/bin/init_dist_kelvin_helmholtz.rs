// Initial setting for the Kelvin Helmholtz problem in 3D

use std::{
    error::Error,
    process,
};

use datafunctions;

use structures::Particle;

use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn Error>> {

    // Files
    let path = "./Data/initial_distribution/kelvin_helmholtz.csv";
    let input_file = "./sedov_blast_wave/input";

    // Parameters
    let input: Vec<f64> = datafunctions::read_input(input_file);

    let eta:f64  = input[0];         // Dimensionless constant specifying the smoothing length
    let gamma:f64= input[1];         // Heat capacity ratio
    let d: i32   = input[2] as i32;  // Dimensions
    
    let x0:f64   = input[3];         // Bottom left corner  (x-coordinate)
    let y0:f64   = input[4];         // Bottom left corner  (y-coordinate)
    let z0:f64   = input[5];         // Bottom left corner  (z-coordinate)
    let wd:f64   = input[6];         // Width (x)
    let lg:f64   = input[7];         // Length (y)
    let hg:f64   = input[8];         // Heigth (z)
    let y1:f64   = input[9];         // Y-lower edge of fluid 2 
    let y2:f64   = input[10];        // Y-upper edge of fluid 2
    let rho1:f64 = input[11];        // Initial density fluid 1
    let rho2:f64 = input[12];        // Initial density fluid 2
    let vx1:f64  = input[13];        // Initial x velocity fluid 1
    let vx2:f64  = input[14];        // Initial x velocity fluid 2
    let p0:f64   = input[15];        // Initial pressure
    
    let nx:usize = input[20] as usize;// Particle resolution in the x direction
    let ny:usize = input[21] as usize;// Particle resolution in the y direction
    let nz:usize = input[22] as usize;// Particle resolution in the z direction
    

    let m:f64    = (rho1 *(y2-y1) + rho2*(lg-y2+y1))*wd*hg; // Total mass
    let n:usize  = 3*nx*ny*nz;        // Total number of particles
    let dm: f64  = m/n as f64;        // Particles' mass
    
    let mut particles :Vec<Particle> = Vec::new();


    kh_init_setup(&mut particles, nx, ny, nz, wd, lg, hg, x0, y0, z0, y1, y2, rho1, rho2, vx1, vx2, p0, gamma, eta, dm, d);
    
    if let Err(err) = datafunctions::save_data(path, &particles){
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