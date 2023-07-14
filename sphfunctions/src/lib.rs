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

pub fn init_square(path: &str, n: u32, rho:f64, h:f64, w:f64, l:f64, x0: f64, y0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let dx = (w*l / n as f64).sqrt();
    let nx :i64 = (w/dx) as i64;
    let ny :i64 = (l/dx) as i64;
    wtr.write_record(&["x", "y", "h", "rho"])?;
    for jj in 0..ny{
        for ii in 0..nx{
            wtr.write_record(&[(x0 + dx*ii as f64).to_string(), (y0 + dx*jj as f64).to_string(), h.to_string(), rho.to_string()])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_square(path: &str, n: u32, rho:f64, h:f64, w:f64, l:f64, x0: f64, y0: f64)-> Result<(), Box<dyn Error>>{
   let mut wtr = Writer::from_path(path)?;
   let mut rng = Pcg64::seed_from_u64(SEED);
   wtr.write_record(&["x", "y", "h", "rho"])?;
   for _ii in 0..n{
       let x = w*rng.gen::<f64>();
       let y = l*rng.gen::<f64>();
       wtr.write_record(&[(x+x0).to_string(), (y+y0).to_string(), h.to_string(), rho.to_string()])?;
   }
   wtr.flush()?;
   Ok(())
}

pub fn init_random_circle(path: &str, n: u32, r:f64, rho:f64, h:f64, x0:f64, y0:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "h", "rho"])?;
    for _ii in 0..n{
        let r_i = r*(rng.gen::<f64>()).sqrt();
        let theta_i = 2.0*PI*rng.gen::<f64>();
        let x = r_i*theta_i.cos() + x0;
        let y = r_i*theta_i.sin() + y0;
        wtr.write_record(&[x.to_string(), y.to_string(), h.to_string(), rho.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "vx", "vy", "u", "h", "rho"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].x.to_string(), particles[ii].y.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(),
                           particles[ii].u.to_string(), particles[ii].h.to_string(), particles[ii].rho.to_string()])?;
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
        particles.push(Particle{x:(&record[0]).parse::<f64>().unwrap(), y:(&record[1]).parse::<f64>().unwrap(),
                                h:(&record[2]).parse::<f64>().unwrap(), rho:(&record[3]).parse::<f64>().unwrap(),
                                ..Default::default()});
    }
    Ok(())
}


// -------- Injecting Particles --------

// Create Particle --------
pub fn create_particle(particles: &mut Vec<Particle>, x_p: f64, y_p: f64, vx_p: f64, vy_p: f64, h_p:f64) {
    particles.push(Particle{
    x: x_p,
    y: y_p,
    vx: vx_p,
    vy: vy_p,
    vx_star: 0.0,
    vy_star: 0.0,
    ax: 0.,
    ay: 0.,
    h: h_p,
    rho: 1.0,
    divv: 0.0,
    u: 0.0,
    u_star: 0.0, 
    du: 0.0});
}

// Replace Particle --------
pub fn replace_particle(particle: &mut Particle, x_p: f64, y_p: f64, vx_p: f64, vy_p: f64, h: f64) {
    particle.x = x_p;
    particle.y = y_p;
    particle.vx = vx_p;
    particle.vy = vy_p;
    particle.ax = 0.;
    particle.ay = 0.;
    particle.vx_star = 0.;
    particle.vy_star = 0.;
    particle.h = h;
    particle.rho = 1.0;
    particle.divv = 0.0;
    particle.u = 0.0;
    particle.u_star = 0.0;
    particle.du = 0.0;
}

// Inject set of particles
pub fn inject_particles(particles: &mut Vec<Particle>, x: f64, y: f64, dx: f64, n_e: usize, h:f64) {
    let x0: f64 = x - 0.5 * dx * n_e as f64;
    for ii in 0..n_e{
        create_particle(particles, x0, y+ (ii as f64) * dx, 0., 0., h);
    }
}


// -------- Basic vector functions --------

// Euclidean distance
pub fn euclidean_norm(p1: &Particle, p2: &Particle) -> f64 {
    let sum :f64 = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    sum.sqrt()
}

// Periodic relative distance
pub fn periodic_rel_vector(p1: &Particle, p2: &Particle, w: f64, l: f64, eps: f64) -> (f64, f64) {
    
    let mut x_temp: f64 = p1.x - p2.x;
    let mut y_temp: f64 = p1.y - p2.y;

    if x_temp.abs() > w-2.*eps {
        if x_temp > 0. {
            x_temp -= w;
        } else {
            x_temp += w;
        }
    }
    if y_temp.abs() > l-2.*eps {
        if y_temp > 0. {
            y_temp -= l;
        } else {
            y_temp += l;
        }
    }
    return (x_temp, y_temp);
}

// Periodic Distance
pub fn periodic_norm(p1: &Particle, p2: &Particle, w: f64, l: f64, eps: f64) -> f64 { 
    let (x_temp, y_temp) = periodic_rel_vector(p1, p2, w, l, eps);
    return (x_temp*x_temp + y_temp*y_temp).sqrt();
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
pub fn density_kernel(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, sigma:f64, d:i32, f: fn(f64)->f64, w: f64, l: f64) -> f64 {
    let mut rho :f64 = 0.0;
    for jj in neigh_particles{
        let r = periodic_norm(&particles[ii], &particles[*jj], w, l, 2.*particles[ii].h);
        rho += f(r/h);
    }
    rho * dm * sigma / h.powi(d)
}

// Density calculated by smoothing function
pub fn density_by_smoothing_length(m:f64, h:f64, eta:f64, d:i32) -> f64{
    m*(eta/h).powi(d)
}

// Omega operator
pub fn omega(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, d:i32, w: f64, l: f64) -> f64{
    let mut omeg :f64 = 0.0;
    for jj in neigh_particles {
        let q = periodic_norm(&particles[ii], &particles[*jj], w, l, 2.*h)/h;
        omeg -= dwdh_(q, f, dfdq, d);
    }
    omeg *= dm*sigma/(h.powi(d)*rho*(d as f64));
    return omeg + 1.;
}


// -------- Root solver --------

// -- Newton-Raphson iterator --

// function and derivative of function
pub fn f_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, w: f64, l: f64) -> (f64 , f64) {
    let rho_kernel = density_kernel(particles, ii, neigh_particles, dm, h, sigma, d, f, w, l);
    let rho_h = density_by_smoothing_length(dm, h, eta, d);
    let f_h = rho_h - rho_kernel;
    let omeg = omega(particles, ii, neigh_particles, dm, h, rho_h, dwdh, f, dfdq, sigma, d, w, l);
    let df = -(d as f64)*rho_h*omeg/ h;
    (f_h, df)
}

// Calculate a new value of 'h'
fn nr_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h_old: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, w: f64, l: f64) -> f64 {
    let (f, df) = f_iter(particles, ii, neigh_particles, dm, h_old, eta, f, dfdq, sigma, d, w, l);
    h_old - f / df
}

// Newton raphson solver to find the value of 'h' for particle 'ii'
pub fn newton_raphson(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, tree: &Node, s_: i32, w: f64, l: f64) -> (f64, Vec<usize>) {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    let mut neighbors: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighboring particles
        neighbors.clear();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, w, l, h_old);
        // Obtain h_new
        h_new = nr_iter(particles, ii, &neighbors, dm, h_old, eta, f, dfdq, sigma, d, w, l);
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
pub fn bisection(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, tree: &Node, s_: i32, w: f64, l: f64) -> (f64, Vec<usize>) {
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
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_left, w, l, h_left);
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_right, w, l, h_right);
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors_mid, w, l, h_mid);

        let f_left = density_by_smoothing_length(dm, h_left, eta, d) - density_kernel(particles, ii, &neighbors_left, dm, h_left, sigma, d, f, w, l);
        let f_right = density_by_smoothing_length(dm, h_right, eta, d) - density_kernel(particles, ii, &neighbors_right, dm, h_right, sigma, d, f, w, l);
        let f_mid = density_by_smoothing_length(dm, h_mid, eta, d) - density_kernel(particles, ii, &neighbors_mid, dm, h_mid, sigma, d, f, w, l);

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
pub fn smoothing_length(particles: &mut Vec<Particle>, dm:f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, dt:f64, tree: &Node, s_: i32, n: usize, ptr : Pointer, w: f64, l: f64){
    (0..n).into_par_iter().for_each(|ii| {
        let (mut h_new, mut neighbors) = newton_raphson(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/(d as f64)), eta, f, dfdq, sigma, d, tol, it, tree, s_, w, l);
        let particle = unsafe { &mut *{ptr}.0.add(ii)};
        if h_new != 0.0 {
            // If h is not found, then keep it constant in time.
            particle.h = h_new;
        } else {
            (h_new, neighbors) = bisection(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/(d as f64)), eta, f, sigma, d, tol, it, tree, s_, w, l);
            if h_new != 0.0 {
                // If h is not found, then keep it constant in time.
                particle.h = h_new;
            } else {
                neighbors.clear();
                tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, w, l, particle.h);
            }
        }
        particle.rho = density_kernel(particles, ii, &neighbors, dm, particle.h, sigma, d, f, w, l);
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
    let dudt :f64 = -alpha_u*v_sig_u*(u_i-u_j)/rho_mean;

    return dudt;
}


// -------- Dynamic Equations --------

// Internal forces: Due to pressure gradient and AV
pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64, art_visc: f64, w: f64, l: f64) -> Vec<f64> {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb + 0.5*art_visc*(grad_ha+grad_hb);
    let (x_rel, y_rel) = periodic_rel_vector(particle_a, particle_b, w, l, 2.*particle_a.h);
    vec![-acc*x_rel, -acc*y_rel]
}

// No body forces
pub fn body_forces_null(_particles: &mut Particle, _nu: f64, _lmbda: f64) {
}

// Body forces for a toy star in 2D
pub fn body_forces_toy_star(particle: &mut Particle, nu: f64, lmbda: f64) {
    particle.ax -= nu * particle.vx + lmbda*particle.x;
    particle.ay -= nu * particle.vy + lmbda*particle.y; 
}

// Gravitational Force due to two massive objects
pub fn body_forces_grav_2obj(particle: &mut Particle, m1: & Star, m2: & Star, omega: f64) {
    let x_1: f64 = particle.x - m1.x;
    let y_1: f64 = particle.y - m1.y;
    let x_2: f64 = particle.x - m2.x;
    let y_2: f64 = particle.y - m2.y;
    let mr_1: f64 = m1.m * (x_1*x_1 + y_1*y_1).powf(-1.5);
    let mr_2: f64 = m2.m * (x_2*x_2 + y_2*y_2).powf(-1.5);

    // Gravitational Force
    particle.ax -= mr_1 * x_1  + mr_2 * x_2;
    particle.ay -= mr_1 * y_1  + mr_2 * y_2;

    // Coriolis and Centripetal Forces
    particle.ax += omega * (omega * particle.x + 2.*particle.vy);
    particle.ay += omega * (omega * particle.y - 2.*particle.vx);
}

// Calculate acceleration for each particle in the system
pub fn accelerations(particles: &mut Vec<Particle>, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                     dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64,
                     d:i32, tree: &Node, s_: i32, n: usize, ptr : Pointer, w: f64, l: f64,
                     artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                     body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool){
    
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, w, l, particles[ii].h);
        return neighbors;
    }).collect();

    (0..n).into_par_iter().for_each(move |ii| {

        // Pointer to iith-particle
        let particle_i = unsafe { &mut *{ptr}.0.add(ii)};

        // Initialize variables to zero
        particle_i.ax = 0.;
        particle_i.ay = 0.;
        particle_i.divv = 0.;
        particle_i.du = 0.;

        let p_i = eos(particles[ii].rho, particles[ii].u, gamma);
        let cs_i = cs(particles[ii].rho, p_i, gamma);
        let omeg_i = omega(particles, ii, &neighbors[ii], dm, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, d, w, l);
        
        for jj in 0..n {
            if ii != jj {
                let p_j = eos(particles[jj].rho, particles[jj].u, gamma);
                let cs_j = cs(particles[jj].rho, p_j, gamma);
                let omeg_j = omega(particles, jj, &neighbors[jj], dm, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, d, w, l);
                
                let (x_rel, y_rel) = periodic_rel_vector(&particles[ii], &particles[jj], w, l, 2.*particles[ii].h);
                let r_ij = (x_rel*x_rel + y_rel*y_rel).sqrt();
                let grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*(particles[ii].h).powi(d+1));
                let grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*(particles[jj].h).powi(d+1));
    
                // Velocity dot position
                let dot_r_v = (particles[ii].vx-particles[jj].vx)*x_rel
                             +(particles[ii].vy-particles[jj].vy)*y_rel;
    
                // Artificial viscosity
                let (art_visc_mom, art_visc_ene) = artificial_viscosity(r_ij, dot_r_v, cs_i, cs_j, particles[ii].h, particles[jj].h, particles[ii].rho, particles[jj].rho);
    
                // Artificial thermal conductivity
                let art_therm_cond: f64 = price08_therm_cond(p_i, p_j, particles[ii].rho, particles[jj].rho, particles[ii].u, particles[jj].u);
    
                // Acceleration
                let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc_mom, w, l);
                particle_i.ax += dm *f_ij[0];
                particle_i.ay += dm *f_ij[1];
                
                // Divergence of v per unit of mass
                let div_vel :f64 = grad_hi*dot_r_v / (omeg_i*particles[ii].rho);
                particle_i.divv -= dm*div_vel;
                
                // Thermal change
                particle_i.du += dm * ((p_i/particles[ii].rho)*div_vel + 0.5*(art_visc_ene + 0.*art_therm_cond*r_ij)*(grad_hi+grad_hj));
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
                        dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64,
                        d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                        artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                        body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                        boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64), w: f64, l: f64, x0: f64, y0: f64) {
    
    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, d, 1e-03, 50, dt, tree, s_, n, ptr, w, l);
    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, d, tree, s_, n, ptr, w, l, artificial_viscosity, body_forces, nu, lmbda, bf);
    particles.par_iter_mut().for_each(|particle|{
        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
        particle.vx += dt * particle.ax;
        particle.vy += dt * particle.ay;
        particle.u += dt * particle.du;
    });
    boundary(particles, w, l, x0, y0);
}

// Velocity Verlet integrator
pub fn velocity_verlet_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64,
                                  d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64), w: f64, l: f64, x0: f64, y0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;

        particle.u += 0.5 * dt * particle.du;

        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
    });

    boundary(particles, w, l, x0, y0);

    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, d, 1e-03, 50, dt, tree, s_, n, ptr, w, l);

    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, d, tree, s_, n, ptr, w, l, artificial_viscosity, body_forces, nu, lmbda, bf);
    particles.par_iter_mut().for_each(|particle|{
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        particle.u += 0.5 * dt * particle.du;
    });
}

// Velocity Verlet integrator
pub fn predictor_kdk_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64,
                                  d:i32, eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, f64, f64), nu:f64, lmbda: f64, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64), w: f64, l: f64, x0: f64, y0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;

        particle.u += 0.5 * dt * particle.du;

        particle.x += dt * particle.vx;
        particle.y += dt * particle.vy;
        
        // Predictor
        
        particle.vx_star = particle.vx;
        particle.vy_star = particle.vy;
        particle.u_star = particle.u;
        
        particle.vx += 0.5 * dt * particle.ax;
        particle.vy += 0.5 * dt * particle.ay;
        
        particle.u += 0.5 * dt * particle.du;
    });

    boundary(particles, w, l, x0, y0);

    tree.build_tree(d as u32, s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, d, 1e-03, 50, dt, tree, s_, n, ptr, w, l);

    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, d, tree, s_, n, ptr, w, l, artificial_viscosity, body_forces, nu, lmbda, bf);
    
    particles.par_iter_mut().for_each(|particle|{
        particle.vx = particle.vx_star + 0.5 * dt * particle.ax;
        particle.vy = particle.vy_star + 0.5 * dt * particle.ay;
        particle.u = particle.u_star + 0.5 * dt * particle.du;
    });
}


// -------- Boundary conditions --------

// Periodic Boundary Conditions
pub fn periodic_boundary(particles: &mut Vec<Particle>, w: f64, h: f64, x0:f64, y0: f64){
    // We assume that the domain's system is a rectangular box.
    particles.par_iter_mut().for_each(|particle|{
        if particle.x > (w+x0) {
            particle.x -= w;
        } else if particle.x < x0 {
            particle.x += w;
        }
        if particle.y > (h + y0) {
            particle.y -= h;
        } else if particle.y < y0 {
            particle.y += h;
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
pub fn time_step_bale_toy_star(particles: & Vec<Particle>, n: usize, gamma: f64, _d: i32, _w: f64, _l: f64, _tree: &mut Node, _s_: i32) -> f64{
    let k: f64 = 0.05;
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay).sqrt();
        let cs: f64 = (gamma*k*(particles[ii].rho).powf(gamma-1.)).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Cossins P. J. (2010)
pub fn time_step_bale(particles: & Vec<Particle>, n: usize, gamma: f64, _d: i32, _w: f64, _l: f64, _tree: &mut Node, _s_: i32) -> f64{
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay).sqrt();
        let cs: f64 = (gamma*(gamma-1.)*particles[ii].u).sqrt();
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Monaghan (1997)
pub fn time_step_mon_toy_star(particles: & Vec<Particle>, n: usize, gamma: f64, d: i32, w: f64, l: f64, tree: &mut Node, s_: i32) -> f64{
    let k: f64 = 0.05;
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, w, l, particles[ii].h);
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
            let (x_rel, y_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], w, l, 2.*particles[ii].h);
            let r_ij = (x_rel*x_rel + y_rel*y_rel).sqrt();
            let dot_r_v = (particles[ii].vx-particles[*jj].vx)*x_rel
            +(particles[ii].vy-particles[*jj].vy)*y_rel;
            
            if dot_r_v < 0. {
                let v_sig_ij = alpha*(cs_i+cs_j - beta*(dot_r_v/r_ij));
                if v_sig_ij > v_sig {
                    v_sig = v_sig_ij;
                }
            }
        }
        let dt_a: f64 = force_dt(particles[ii].h, (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay).sqrt(), 0.25);
        let dt_cfl: f64 = 0.3*particles[ii].h / v_sig;
    return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Monaghan (1997)
pub fn time_step_mon(particles: & Vec<Particle>, n: usize, gamma: f64, d: i32, w: f64, l: f64, tree: &mut Node, s_: i32) -> f64{
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors, w, l, particles[ii].h);
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
            let (x_rel, y_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], w, l, 2.*particles[ii].h);
            let r_ij = (x_rel*x_rel + y_rel*y_rel).sqrt();
            let dot_r_v = (particles[ii].vx-particles[*jj].vx)*x_rel
            +(particles[ii].vy-particles[*jj].vy)*y_rel;
            
            if dot_r_v < 0. {
                let v_sig_ij = alpha*(cs_i+cs_j - beta*(dot_r_v/r_ij));
                if v_sig_ij > v_sig {
                    v_sig = v_sig_ij;
                }
            }
        }
        let dt_a: f64 = force_dt(particles[ii].h, (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay).sqrt(), 0.25);
        let dt_cfl: f64 = 0.3*particles[ii].h / v_sig;
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}