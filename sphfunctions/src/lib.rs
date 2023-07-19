use std::{
    error::Error,
    f64,
};

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use csv::Writer;

use rayon::prelude::*;

use std::f64::consts::PI;

use tree_algorithm::{
    FindNeighbors,
    BuildTree,
};

use structures::{
    Particle,
    Node,
    Pointer,
    Star,
};

// -------- Write data --------

pub fn init_random_square(path: &str, n: u32, h:f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "z", "h", "rho"])?;
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

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "u"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].x.to_string(), particles[ii].y.to_string(), particles[ii].z.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(), particles[ii].vz.to_string(),
                           particles[ii].h.to_string(), particles[ii].u.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

// -------- Read data --------

pub fn read_data(path: &str, particles: &mut Vec<Particle>) -> Result<(), Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;
    for result in rdr.records() {
        let record = result?;
        particles.push(Particle{x:(&record[0]).parse::<f64>().unwrap(), y:(&record[1]).parse::<f64>().unwrap(), z:(&record[2]).parse::<f64>().unwrap(),
                                vx:(&record[3]).parse::<f64>().unwrap(), vy:(&record[4]).parse::<f64>().unwrap(), vz:(&record[5]).parse::<f64>().unwrap(),
                                h:(&record[6]).parse::<f64>().unwrap(), u:(&record[7]).parse::<f64>().unwrap(),
                                ..Default::default()});
    }
    Ok(())
}

// -------- Basic vector functions --------

// Periodic relative distance
pub fn periodic_rel_vector(p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64) -> (f64, f64, f64) {
    
    let mut x_temp: f64 = p1.x - p2.x;
    let mut y_temp: f64 = p1.y - p2.y;
    let mut z_temp: f64 = p1.z - p2.z;

    if x_temp.abs() > wd-2.*eps {
        if x_temp > 0. {
            x_temp -= wd;
        } else {
            x_temp += wd;
        }
    }
    if y_temp.abs() > lg-2.*eps {
        if y_temp > 0. {
            y_temp -= lg;
        } else {
            y_temp += lg;
        }
    }
    if z_temp.abs() > hg-2.*eps {
        if z_temp > 0. {
            z_temp -= hg;
        } else {
            z_temp += hg;
        }
    }
    return (x_temp, y_temp, z_temp);
}

// Periodic Distance
pub fn periodic_norm(p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64) -> f64 { 
    let (x_temp, y_temp, z_temp) = periodic_rel_vector(p1, p2, wd, lg, hg, eps);
    return (x_temp*x_temp + y_temp*y_temp + z_temp*z_temp).sqrt();
}


// -------- Kernel function --------

// Cubic Kernel
pub fn f_cubic_kernel(q:f64) -> f64 {
    let mut f:f64 = 0.;
    if q < 1. {
        f = 1. - 0.75*q*q*(2.0 - q); 
    } else if q < 2.{
        f = 0.25*(2.-q).powi(3);
    }
    f
}

// Derivative of cubic kernel
pub fn dfdq_cubic_kernel(q:f64) -> f64 {
    let mut df:f64 = 0.;
    if q < 1. {
        df = (2.25*q-3.)*q;
    } else if q < 2.{
        df = -3.*(0.25*q*q-q+1.);
    }
    df
}

// Gaussian Kernel
pub fn f_gaussian_kernel(q:f64) -> f64 {
    let mut f:f64 = 0.;
    if q < 2. {
        f = (-q*q).exp();
    }
    f
}

// Derivative of Gaussian Kernel
pub fn dfdq_gaussian_kernel(q:f64) -> f64 {
    let mut f:f64 = 0.;
    if q < 2. {
        f = -2.0*q*(-q*q).exp();
    }
    f
}

// Derivative of kernel w.r.t the smoothing length
pub fn dwdh(q: f64, f: fn(f64) -> f64, df: fn(f64) -> f64, d:i32) -> f64 {
    (d as f64) *f(q) + q*df(q)
}


// -------- Kernel approximations --------

// Kernel approximation of density
pub fn density_kernel(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, sigma:f64, rkern: f64, d:i32, f: fn(f64)->f64, wd: f64, lg: f64, hg: f64) -> f64 {
    let mut rho :f64 = 0.0;
    for jj in neigh_particles{
        let r = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h);
        rho += f(r/h);
    }
    rho * dm * sigma / h.powi(d)
}

// Density calculated by smoothing function
pub fn density_by_smoothing_length(m:f64, h:f64, eta:f64, d:i32) -> f64{
    m*(eta/h).powi(d)
}

// Density calculated by smoothing function
pub fn h_by_density(m:f64, rho:f64, eta:f64, d:i32) -> f64{
    eta*(m/rho).powf(1./d as f64)
}

// Omega operator
pub fn omega(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64, d:i32, wd: f64, lg: f64, hg: f64) -> f64{
    let mut omeg :f64 = 0.0;
    for jj in neigh_particles {
        let q = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*h)/h;
        omeg -= dwdh_(q, f, dfdq, d);
    }
    omeg *= dm*sigma/(h.powi(d)*rho*(d as f64));
    return omeg + 1.;
}


// -------- Root solver --------

// -- Newton-Raphson iterator --

// function and derivative of function
pub fn f_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, d:i32, wd: f64, lg: f64, hg: f64) -> (f64 , f64) {
    let rho_kernel = density_kernel(particles, ii, neigh_particles, dm, h, sigma, rkern, d, f, wd, lg, hg);
    let rho_h = density_by_smoothing_length(dm, h, eta, d);
    let f_h = rho_h - rho_kernel;
    let omeg = omega(particles, ii, neigh_particles, dm, h, rho_h, dwdh, f, dfdq, sigma, rkern, d, wd, lg, hg);
    let df = -(d as f64)*rho_h*omeg/ h;
    (f_h, df)
}

// Calculate a new value of 'h'
fn nr_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h_old: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, d:i32, wd: f64, lg: f64, hg: f64) -> f64 {
    let (f, df) = f_iter(particles, ii, neigh_particles, dm, h_old, eta, f, dfdq, sigma, rkern, d, wd, lg, hg);
    h_old - f / df
}

// Newton raphson solver to find the value of 'h' for particle 'ii'
pub fn newton_raphson(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, d:i32, tol: f64, it: u32, tree: &Node, s_: i32, wd: f64, lg: f64, hg: f64) -> (f64, Vec<usize>) {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    let mut neighbors: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighboring particles
        neighbors.clear();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, wd, lg, hg, h_old);
        // Obtain h_new
        h_new = nr_iter(particles, ii, &neighbors, dm, h_old, eta, f, dfdq, sigma, rkern, d, wd, lg, hg);
        if ((h_new - h_old)/h_old).abs() <=  tol {
            i = it + 2;
        } else{
            i += 1;
            h_old = h_new;
        }
    }
    if i == it+1 {
        println!("PROBLEM NEWTHON RAPHSON METHOD: Solution not found for particle {}. Initial guess was h = {}.", ii, h_guess);
        (0.0, neighbors)
    } else{
        (h_new, neighbors)
    }
}

// -- Bisection iterator --

// Bisection solver to find the value of 'h' for particle 'ii'
pub fn bisection(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, sigma:f64, rkern: f64, d:i32, tol: f64, it: u32, tree: &Node, s_: i32, wd: f64, lg: f64, hg: f64) -> (f64, Vec<usize>) {
    let mut h_left :f64 = 0.1*h_guess;
    let mut h_right :f64 = 0.2;
    let mut h_mid = 0.5*(h_left + h_right);
    let mut i : u32 = 1;
    let mut neighbors_left: Vec<usize> = Vec::new();
    let mut neighbors_right: Vec<usize> = Vec::new();
    let mut neighbors_mid: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighboring particles
        neighbors_left.clear();
        neighbors_right.clear();
        neighbors_mid.clear();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_left, wd, lg, hg, h_left);
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_right, wd, lg, hg, h_right);
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_mid, wd, lg, hg, h_mid);

        let f_left = density_by_smoothing_length(dm, h_left, eta, d) - density_kernel(particles, ii, &neighbors_left, dm, h_left, sigma, rkern, d, f, wd, lg, hg);
        let f_right = density_by_smoothing_length(dm, h_right, eta, d) - density_kernel(particles, ii, &neighbors_right, dm, h_right, sigma, rkern, d, f, wd, lg, hg);
        let f_mid = density_by_smoothing_length(dm, h_mid, eta, d) - density_kernel(particles, ii, &neighbors_mid, dm, h_mid, sigma, rkern, d, f, wd, lg, hg);

        if f_right.signum() == f_left.signum() {
            i = it + 2;
        }

        if ((h_right - h_left)/h_mid).abs() <=  tol || f_mid.abs() < 0.001*tol  {
            i = it + 2;
        } else{
            i += 1;
            if f_mid.signum() == f_left.signum() {
                h_left = h_mid;
            } else {
                h_right = h_mid;
            }
            h_mid = 0.5*(h_left + h_right);
        }
    }
    if i == it+1 {
        println!("PROBLEM BISECTION METHOD: Solution not found for particle {}. Final value was h = {}.", ii, h_mid);
        (0.0, neighbors_mid)
    } else{
        (h_mid, neighbors_mid)
    }
}


// -------- Smoothing length --------

// Calculate the smoothing function for each particle in a given time.
pub fn smoothing_length(particles: &mut Vec<Particle>, dm:f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, d:i32, tol: f64, it: u32, dt:f64, tree: &Node, s_: i32, n: usize, ptr : Pointer, wd: f64, lg: f64, hg: f64){
    (0..n).into_par_iter().for_each(|ii| {
        let (mut h_new, mut neighbors) = newton_raphson(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/(d as f64)), eta, f, dfdq, sigma, rkern, d, tol, it, tree, s_, wd, lg, hg);
        let particle = unsafe { &mut *{ptr}.0.add(ii)};
        if h_new != 0.0 {
            // If h is not found, then keep it constant in time.
            particle.h = h_new;
        } else {
            (h_new, neighbors) = bisection(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/(d as f64)), eta, f, sigma, rkern, d, tol, it, tree, s_, wd, lg, hg);
            if h_new != 0.0 {
                // If h is not found, then keep it constant in time.
                particle.h = h_new;
            } else {
                neighbors.clear();
                tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, wd, lg, hg, particle.h);
            }
        }
        particle.rho = density_kernel(particles, ii, &neighbors, dm, particle.h, sigma, rkern, d, f, wd, lg, hg);
    });
}


// -------- Equations of state --------

// -- Toy Star 2D --

// Polytropic equation
pub fn eos_polytropic(rho:f64, _k:f64, gamma:f64) -> f64 {
    let k: f64 = 0.05;
    k * rho.powf(gamma)
}

// Coefficient of gravital force
pub fn coeff_static_grav_potential(k:f64, gamma:f64, m:f64, r:f64) -> f64 {
    2.0*k*PI.powf(1.-gamma) * (m*(gamma/(gamma - 1.))/(r*r)).powf(gamma)/m
}

// -- Ideal Gas --

pub fn eos_ideal_gas(rho:f64, u:f64, gamma:f64) -> f64 {
    (gamma-1.)*rho*u
}

pub fn thermal_energy(rho:f64, p:f64, gamma:f64) -> f64 {
    p/((gamma-1.)*rho)
}

// Sound speed for the ideal gas eos
pub fn sound_speed_ideal_gas(rho:f64, p:f64, gamma:f64) -> f64 {
    (gamma*p/rho).sqrt()
}


// -------- Artificial Viscosity --------

// Monaghan (1989): "Standard" SPH viscous term
pub fn mon89_art_vis(r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, h_i: f64, h_j: f64, rho_i: f64, rho_j: f64) -> (f64, f64) {
    if dot_r_v <= 0. {
        // Mean values
        let cs_mean :f64 = 0.5*(cs_i+cs_j);
        let h_mean :f64 = 0.5*(h_i+h_j);
        let rho_mean :f64 = 0.5*(rho_i+rho_j);
    
        // Parameters
        let alpha :f64 = 1.0;
        let beta :f64 = 2.0;
        let eps :f64 = 0.01;
        let nu_visc :f64 = h_mean*dot_r_v/(r_ij*r_ij+eps*h_mean*h_mean);
        let dvdt :f64 = (-alpha*cs_mean+beta*nu_visc)*nu_visc/rho_mean;
    
        return (dvdt, 0.5*dvdt*dot_r_v);
    } else {
        return (0.0, 0.0);
    }
}

// Monaghan (1997): AV by Rieman solvers
pub fn mon97_art_vis(r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, _h_i: f64, _h_j: f64, rho_i: f64, rho_j: f64) -> (f64, f64) {
    if dot_r_v <= 0. {
        // Parameters
        let alpha: f64 = 1.0;
        let beta: f64 = 2.0;

        let v_sig:f64 = 0.5*alpha*(cs_i + cs_j - beta*dot_r_v/r_ij);
        let rho_mean :f64 = 0.5*(rho_i+rho_j);
        let dvdt :f64 = -v_sig*dot_r_v/(r_ij*rho_mean);
        let dudt :f64 = 0.5*dvdt*dot_r_v;
        return (dvdt, dudt);
    } else {
        return (0.0, 0.0);
    }
}

// Price (2008): Thermal conductivity switches
pub fn price08_therm_cond(p_i: f64, p_j: f64, rho_i: f64, rho_j: f64, u_i: f64, u_j: f64) -> f64 {
    // Parameters
    let alpha_u: f64 = 1.0;

    let rho_mean :f64 = 0.5*(rho_i+rho_j);
    let v_sig_u:f64 = ((p_i - p_j).abs()/rho_mean).sqrt();
    let dudt :f64 = alpha_u*v_sig_u*(u_i-u_j)/rho_mean;

    return dudt;
}


// -------- Dynamic Equations --------

// Internal forces: Due to pressure gradient and AV
pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, x_rel: f64, y_rel: f64, z_rel: f64, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64, art_visc: f64) -> (f64, f64, f64) {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb + 0.5*art_visc*(grad_ha+grad_hb);
    (-acc*x_rel, -acc*y_rel, -acc*z_rel)
}

// No body forces
pub fn body_forces_null(_particles: &mut Particle, _nu: f64, _lmbda: f64) {
}

// Body forces for a toy star in 2D
pub fn body_forces_toy_star(particle: &mut Particle, nu: f64, lmbda: f64) {
    particle.ax -= nu * particle.vx + lmbda*particle.x;
    particle.ay -= nu * particle.vy + lmbda*particle.y; 
    particle.az -= nu * particle.vz + lmbda*particle.z; 
}

// Gravitational Force due to two massive objects
pub fn body_forces_grav_2obj(particle: &mut Particle, m1: & Star, m2: & Star, omega: f64) {
    let x_1: f64 = particle.x - m1.x;
    let y_1: f64 = particle.y - m1.y;
    let z_1: f64 = particle.z - m1.z;
    let x_2: f64 = particle.x - m2.x;
    let y_2: f64 = particle.y - m2.y;
    let z_2: f64 = particle.z - m2.z;
    let mr_1: f64 = m1.m * (x_1*x_1 + y_1*y_1 + z_1*z_1).powf(-1.5);
    let mr_2: f64 = m2.m * (x_2*x_2 + y_2*y_2 + z_2*z_2).powf(-1.5);

    // Gravitational Force
    particle.ax -= mr_1 * x_1  + mr_2 * x_2;
    particle.ay -= mr_1 * y_1  + mr_2 * y_2;
    particle.az -= mr_1 * z_1  + mr_2 * z_2;

    // Coriolis and Centripetal Forces
    particle.ax += omega * (omega * particle.x + 2.*particle.vy);
    particle.ay += omega * (omega * particle.y - 2.*particle.vx);
}

// Calculate acceleration for each particle in the system
pub fn accelerations(particles: &mut Vec<Particle>, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                     dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                     d:i32, tree: &Node, s_: i32, n: usize, ptr : Pointer, wd: f64, lg: f64, hg: f64,
                     artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                     body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool){
    
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, wd, lg, hg, particles[ii].h);
        return neighbors;
    }).collect();

    (0..n).into_par_iter().for_each(move |ii| {

        // Pointer to iith-particle
        let particle_i = unsafe { &mut *{ptr}.0.add(ii)};

        // Initialize variables to zero
        particle_i.ax = 0.;
        particle_i.ay = 0.;
        particle_i.az = 0.;
        particle_i.divv = 0.;
        particle_i.du = 0.;

        let p_i = eos(particles[ii].rho, particles[ii].u, gamma);
        let cs_i = cs(particles[ii].rho, p_i, gamma);
        let omeg_i = omega(particles, ii, &neighbors[ii], dm, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, rkern, d, wd, lg, hg);
        
        for jj in 0..n {
            if ii != jj {
                let p_j = eos(particles[jj].rho, particles[jj].u, gamma);
                let cs_j = cs(particles[jj].rho, p_j, gamma);
                let omeg_j = omega(particles, jj, &neighbors[jj], dm, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, rkern, d, wd, lg, hg);
                
                let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[jj], wd, lg, hg, rkern*particles[ii].h);
                let r_ij = (x_rel*x_rel + y_rel*y_rel).sqrt();
                let grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*(particles[ii].h).powi(d+1));
                let grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*(particles[jj].h).powi(d+1));
    
                // Velocity dot position
                let dot_r_v = (particles[ii].vx-particles[jj].vx)*x_rel
                             +(particles[ii].vy-particles[jj].vy)*y_rel
                             +(particles[ii].vz-particles[jj].vz)*z_rel;
    
                // Artificial viscosity
                let (art_visc_mom, art_visc_ene) = artificial_viscosity(r_ij, dot_r_v, cs_i, cs_j, particles[ii].h, particles[jj].h, particles[ii].rho, particles[jj].rho);
    
                // Artificial thermal conductivity
                let art_therm_cond: f64 = price08_therm_cond(p_i, p_j, particles[ii].rho, particles[jj].rho, particles[ii].u, particles[jj].u);
    
                // Acceleration
                let (f_ij_x, f_ij_y, f_ij_z) = acceleration_ab(&particles[ii], &particles[jj], x_rel, y_rel, z_rel, p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc_mom);
                particle_i.ax += dm * f_ij_x;
                particle_i.ay += dm * f_ij_y;
                particle_i.az += dm * f_ij_z;
                
                // Divergence of v per unit of mass
                let div_vel :f64 = grad_hi*dot_r_v / (omeg_i*particles[ii].rho);
                particle_i.divv -= dm*div_vel;
                
                // Thermal change
                particle_i.du += dm * ((p_i/particles[ii].rho)*div_vel + 0.5*(art_visc_ene + art_therm_cond*r_ij)*(grad_hi/omeg_i+grad_hj/omeg_j));
            }
        }
            
        // Body forces
        if bf {
            body_forces(particle_i, nu, lmbda);
        }
    });
}


// -------- Time integrator --------
pub fn euler_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                        dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                        d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                        artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                        body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                        boundary: fn(&mut Vec<Particle>, f64, f64,f64, f64, f64, f64), wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, d, 1e-03, 50, dt, tree, s_, n, ptr, wd, lg, hg);
    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, d, tree, s_, n, ptr, wd, lg, hg, artificial_viscosity, body_forces, nu, lmbda, bf);
    particles.par_iter_mut().for_each(|particle|{
        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
        particle.z += dt * particle.vz;
        particle.vx += dt * particle.ax;
        particle.vy += dt * particle.ay;
        particle.vz += dt * particle.az;
        particle.u += dt * particle.du;
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
}

// Velocity Verlet integrator
pub fn velocity_verlet_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                                  d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        particle.vz += 0.5 * dt * particle.az;

        particle.u += 0.5 * dt * particle.du;

        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
        particle.z += dt * particle.vz;
    });

    boundary(particles, wd, lg, hg, x0, y0, z0);

    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, d, 1e-03, 50, dt, tree, s_, n, ptr, wd, lg, hg);

    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, d, tree, s_, n, ptr, wd, lg, hg, artificial_viscosity, body_forces, nu, lmbda, bf);
    particles.par_iter_mut().for_each(|particle|{
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        particle.vz += 0.5 * dt * particle.az;
        particle.u += 0.5 * dt * particle.du;
    });
}

// Velocity Verlet integrator
pub fn predictor_kdk_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                                  d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        // println!("1 : {:?}", particle);
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        particle.vz += 0.5 * dt * particle.az;
        
        particle.u += 0.5 * dt * particle.du;
        
        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
        particle.z += dt * particle.vz;
        
        // Predictor
        
        particle.vx_star = particle.vx;
        particle.vy_star = particle.vy;
        particle.vz_star = particle.vz;
        particle.u_star = particle.u;
        
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        particle.vz += 0.5 * dt * particle.az;
        
        particle.u += 0.5 * dt * particle.du;
    });

    boundary(particles, wd, lg, hg, x0, y0, z0);

    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, d, 1e-03, 50, dt, tree, s_, n, ptr, wd, lg, hg);

    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, d, tree, s_, n, ptr, wd, lg, hg, artificial_viscosity, body_forces, nu, lmbda, bf);
    
    particles.par_iter_mut().for_each(|particle|{
        particle.vx = particle.vx_star + 0.5 * dt * particle.ax;
        particle.vy = particle.vy_star + 0.5 * dt * particle.ay;
        particle.vz = particle.vz_star + 0.5 * dt * particle.az;
        particle.u = particle.u_star + 0.5 * dt * particle.du;
    });
}


// -------- Boundary conditions --------

// Periodic Boundary Conditions
pub fn periodic_boundary(particles: &mut Vec<Particle>, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0: f64){
    // We assume that the domain's system is a rectangular box.
    particles.par_iter_mut().for_each(|particle|{
        if particle.x > (wd+x0) {
            particle.x -= wd;
        } else if particle.x < x0 {
            particle.x += wd;
        }
        if particle.y > (lg + y0) {
            particle.y -= lg;
        } else if particle.y < y0 {
            particle.y += lg;
        }
        if particle.z > (hg + z0) {
            particle.z -= hg;
        } else if particle.z < z0 {
            particle.z += hg;
        }
    });
}

// -------- Timestepping Criteria --------

// Bate at al. (1995). CFL criterion
pub fn cfl_dt(h: f64, cs: f64, div_v:f64, alpha:f64, beta: f64) -> f64{
    if div_v < 0. {
        return 0.3*h / (cs + h*div_v.abs() + 1.2*(alpha*cs + beta*h*div_v.abs()));
    } else {
        return 0.3*h / (cs + h*div_v.abs());
    }
}

// MOnaghan (1989) Force conditon
pub fn force_dt(h: f64, a: f64, f: f64) -> f64 {
    f*(h/a).sqrt()
}

// Timestepping Criteria Cossins P. J. (2010)
pub fn time_step_bale_toy_star(particles: & Vec<Particle>, n: usize, gamma: f64, _rkern: f64, _d: i32, _wd: f64, _lg: f64, _hg: f64, _tree: &mut Node, _s_: i32) -> f64{
    let k: f64 = 0.05;
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let cs: f64 = (gamma*k*(particles[ii].rho).powf(gamma-1.)).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Cossins P. J. (2010)
pub fn time_step_bale(particles: & Vec<Particle>, n: usize, gamma: f64, _rkern: f64, _d: i32, _wd: f64, _lg: f64, _hg: f64, _tree: &mut Node, _s_: i32) -> f64{
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let cs: f64 = (gamma*(gamma-1.)*particles[ii].u).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Monaghan (1997)
pub fn time_step_mon_toy_star(particles: & Vec<Particle>, n: usize, gamma: f64, rkern: f64, d: i32, wd: f64, lg: f64, hg: f64, tree: &mut Node, s_: i32) -> f64{
    let k: f64 = 0.05;
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, wd, lg, hg, particles[ii].h);
        return neighbors;
    }).collect();
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let alpha: f64 = 1.;
        let beta: f64 = 2.;
        let mut v_sig:f64 = 0.0;
        let cs_i: f64 = (gamma*k*(particles[ii].rho).powf(gamma-1.)).sqrt();
        for jj in &neighbors[ii] {
            let cs_j: f64 = (gamma*k*(particles[*jj].rho).powf(gamma-1.)).sqrt();

            // Velocity dot position
            let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h);
            let r_ij = (x_rel*x_rel + y_rel*y_rel + z_rel*z_rel).sqrt();
            let dot_r_v = (particles[ii].vx-particles[*jj].vx)*x_rel
                         +(particles[ii].vy-particles[*jj].vy)*y_rel
                         +(particles[ii].vz-particles[*jj].vz)*z_rel;
            
            if dot_r_v < 0. {
                let v_sig_ij = alpha*(cs_i+cs_j - beta*(dot_r_v/r_ij));
                if v_sig_ij > v_sig {
                    v_sig = v_sig_ij;
                }
            }
        }
        let a_norm: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a_norm, 0.25);
        let dt_cfl: f64 = 0.3*particles[ii].h / v_sig;
    return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Monaghan (1997)
pub fn time_step_mon(particles: & Vec<Particle>, n: usize, gamma: f64, rkern: f64, d: i32, wd: f64, lg: f64, hg: f64, tree: &mut Node, s_: i32) -> f64{
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, wd, lg, hg, particles[ii].h);
        return neighbors;
    }).collect();
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let alpha: f64 = 1.;
        let beta: f64 = 2.;
        let mut v_sig:f64 = 0.0;
        let cs_i: f64 = (gamma*(gamma-1.)*particles[ii].u).sqrt();
        for jj in &neighbors[ii] {
            let cs_j: f64 = (gamma*(gamma-1.)*particles[*jj].u).sqrt();

            // Velocity dot position
            let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h);
            let r_ij = (x_rel*x_rel + y_rel*y_rel + z_rel*z_rel).sqrt();
            let dot_r_v = (particles[ii].vx-particles[*jj].vx)*x_rel
                         +(particles[ii].vy-particles[*jj].vy)*y_rel
                         +(particles[ii].vz-particles[*jj].vz)*z_rel;
            
            if dot_r_v < 0. {
                let v_sig_ij = alpha*(cs_i+cs_j - beta*(dot_r_v/r_ij));
                if v_sig_ij > v_sig {
                    v_sig = v_sig_ij;
                }
            }
        }
        let a_norm: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a_norm, 0.25);
        let dt_cfl: f64 = 0.3*particles[ii].h / v_sig;
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}