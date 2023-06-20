use std::{
    error::Error,
    f64,
};

use rand::{thread_rng, Rng};

use csv::Writer;

use rayon::prelude::*;

use std::f64::consts::PI;

use tree_algorithm::{
    FindNeighbors,
};

use structures::{
    Particle,
    Node,
    Pointer,
};

// -------- Write data --------

pub fn init_square(path: &str, n: u32, rho:f64, h:f64, w:f64, l:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let dx = (w*l / n as f64).sqrt();
    let nx :i64 = (w/dx) as i64;
    let ny :i64 = (l/dx) as i64;
    wtr.write_record(&["x", "y", "h", "rho"])?;
    for jj in 0..ny{
        for ii in 0..nx{
            wtr.write_record(&[(dx*ii as f64).to_string(), (dx*jj as f64).to_string(), h.to_string(), rho.to_string()])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_square(path: &str, n: u32, rho:f64, h:f64, w:f64, l:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = thread_rng();
    wtr.write_record(&["x", "y", "h", "rho"])?;
    for _ii in 0..n{
        let x = rng.gen_range(0.0f64, w);
        let y = rng.gen_range(0.0f64, l);
        wtr.write_record(&[x.to_string(), y.to_string(), h.to_string(), rho.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_circle(path: &str, n: u32, r:f64, rho:f64, h:f64, x0:f64, y0:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = thread_rng();
    wtr.write_record(&["x", "y", "h", "rho"])?;
    for _ii in 0..n{
        let r_i = r*(rng.gen_range(0.0f64, 1.0f64)).sqrt();
        let theta_i = 2.0*PI*rng.gen_range(0.0f64, 1.0f64);
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


// -------- Basic vector functions --------

// Euclidean distance
pub fn euclidean_norm(p1: &Particle, p2: &Particle) -> f64 {
    let sum :f64 = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    sum.sqrt()
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
pub fn density_kernel(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, sigma:f64, d:i32, f: fn(f64)->f64) -> f64 {
    let mut rho :f64 = 0.0;
    for jj in neigh_particles{
        let r = euclidean_norm(&particles[ii], &particles[*jj]);
        rho += f(r/h);
    }
    rho * dm * sigma / h.powi(d)
    //let rho :f64 = neigh_particles.par_iter().map(|jj| f(euclidean_norm(&particles[ii], &particles[*jj])/h)).sum(); // Parallel calculation
}

// Density calculated by smoothing function
pub fn density_by_smoothing_length(m:f64, h:f64, eta:f64, d:i32) -> f64{
    let vol = eta/h;
    m*vol.powi(d)
}

// Omega operator
pub fn omega(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, d:i32) -> f64{
    let mut omeg :f64 = 0.0;
    for jj in neigh_particles {
        let q = euclidean_norm(&particles[ii], &particles[*jj])/h;
        omeg -= dwdh_(q, f, dfdq, d);
    }
    omeg *= dm*sigma/(h.powi(d)*rho*(d as f64));
    omeg + 1.
}


// -------- Root solver --------

// -- Newton-Raphson iterator --

// function and derivative of function
pub fn f_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32) -> (f64 , f64) {
    let rho_kernel = density_kernel(particles, ii, neigh_particles, dm, h, sigma, d, f);
    let rho_h = density_by_smoothing_length(dm, h, eta, d);
    let f_h = rho_h - rho_kernel;
    let omeg = omega(particles, ii, neigh_particles, dm, h, rho_kernel, dwdh, f, dfdq, sigma, d);
    let df = -(d as f64)*rho_h*omeg/ h;
    (f_h, df)
}

// Calculate a new value of 'h'
fn nr_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h_old: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32) -> f64 {
    let (f, df) = f_iter(particles, ii, neigh_particles, dm, h_old, eta, f, dfdq, sigma, d);
    (h_old - f / df).abs()
}

// Newton raphson solver to find the value of 'h' for particle 'ii'
pub fn newton_raphson(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, tree: &Node, s_: u32) -> (f64, Vec<usize>) {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    let mut neighbors: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighboring particles
        neighbors.clear();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors);
        // Obtain h_new
        h_new = nr_iter(particles, ii, &neighbors, dm, h_old, eta, f, dfdq, sigma, d);
        if (h_new - h_old).abs() <=  tol {
            i = it + 2;
        } else{
            i += 1;
            h_old = h_new;
        }
    }
    if i == it+1 {
        (0.0, neighbors)
    } else{
        (h_new, neighbors)
    }
}


// -------- Smoothing length --------

// Calculate the smoothing function for each particle in a given time.
pub fn smoothing_length(particles: &mut Vec<Particle>, dm:f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, dt:f64, tree: &Node, s_: u32, n: usize, ptr : Pointer){
    (0..n).into_par_iter().for_each(move |ii| {
        let (h_new, neighbors) = newton_raphson(ii, particles, dm, particles[ii].h+dt*particles[ii].dh, eta, f, dfdq, sigma, d, tol, it, tree, s_);
        let particle = unsafe { &mut *{ptr}.0.add(ii)};
        if h_new != 0.0 {
            // If h is not found, then keep it constant in time.
            particle.h = h_new;
        }
        particle.rho = density_kernel(particles, ii, &neighbors, dm, particle.h, sigma, d, f);
    });
}


// -------- Equations of state --------

// -- Toy Star 2D --

// Polytropic equation
pub fn eos_polytropic(rho:f64, k:f64, gamma:f64) -> f64 {
    k * rho.powf(1.+1./gamma)
}

// Coefficient of gravital force
pub fn coeff_static_grav_potential(k:f64, gamma:f64, m:f64, r:f64) -> f64 {
    2.0*k/(PI.powf(1./gamma)) * (m*(1.+gamma)/(r*r)).powf(1.+1./gamma)/m
}

// -- Ideal Gas --

pub fn eos_ideal_gas(rho:f64, k:f64, gamma:f64) -> f64 {
    k*rho.powf(gamma)
}

pub fn thermal_energy(rho:f64, p:f64, gamma:f64) -> f64 {
    p/((gamma-1.)*rho)
}


// -------- Dynamic Equations --------

// Force due to the pressure's gradient
pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64, art_visc: f64) -> Vec<f64> {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb + art_visc;
    vec![-acc*(particle_a.x - particle_b.x), -acc*(particle_a.y - particle_b.y)]
}

// Body forces for a toy star in 2D
pub fn body_forces_toy_star(particle: &mut Particle, nu: f64, lmbda: f64) {
    particle.ax -= nu * particle.vx + lmbda*particle.x;
    particle.ay -= nu * particle.vy + lmbda*particle.y; 
}

// Calculate acceleration for each particle in the system
pub fn accelerations(particles: &mut Vec<Particle>, dm:f64, eos: fn(f64, f64, f64)->f64, k:f64, gamma:f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, d:i32, tree: &Node, s_: u32, n: usize, ptr : Pointer){
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, d as f64, s_, particles, &mut neighbors);
        return neighbors;
    }).collect();

    (0..n).into_par_iter().for_each(move |ii| {
        let p_i = eos(particles[ii].rho, k, gamma);
        let omeg_i = omega(particles, ii, &neighbors[ii], dm, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, d);
        // Pointer to iith-particle
        let particle_i = unsafe { &mut *{ptr}.0.add(ii)};
        for jj in (ii+1)..n {
            let p_j = eos(particles[jj].rho, k, gamma);
            let omeg_j = omega(particles, jj, &neighbors[jj], dm, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, d);
            let r_ij = euclidean_norm(&particles[ii], &particles[jj]);
            let grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*(particles[ii].h).powi(d+1));
            let grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*(particles[jj].h).powi(d+1));

            // Divergence of velocity
            let dot_r_v = (particles[ii].vx-particles[jj].vx)*(particles[ii].x-particles[jj].x)
                         +(particles[ii].vy-particles[jj].vy)*(particles[ii].y-particles[jj].y);
            
            // Artificial viscosity
            let alpha :f64 = 1.0;
            let eps :f64 = 0.01;
            let h_mean = 0.5*(particles[ii].h+particles[jj].h);
            let nu_visc = alpha*2.0*h_mean/(particles[ii].rho+particles[jj].rho);
            let mut art_visc = 0.0;
            if dot_r_v < 0.0 {
                art_visc = -nu_visc*dot_r_v/(r_ij*r_ij+eps*h_mean*h_mean);
            }

            // Pointer to jjth-particle
            let particle_j = unsafe { &mut *{ptr}.0.add(jj)};

            // Acceleration
            let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc);
            particle_i.ax += dm *f_ij[0];
            particle_i.ay += dm *f_ij[1];
            particle_j.ax -= dm *f_ij[0];
            particle_j.ay -= dm *f_ij[1];
            
            // Smoothing length change
            let div_vel :f64 = grad_hi*dot_r_v;
            particle_i.dh += particles[ii].h*div_vel/(d as f64);
            particle_j.dh += particles[jj].h*div_vel/(d as f64);
            
            // Thermal change
            particle_i.du = dm*p_i / (omeg_i*particles[ii].rho*particles[ii].rho) * div_vel;
            particle_j.du = dm*p_j / (omeg_j*particles[jj].rho*particles[jj].rho) * div_vel;
        }
    });
}

// -------- Time integrator --------

// Euler-Raphson method
pub fn euler_integrator(particle: &mut Particle, dt: f64) {
    particle.x += dt * particle.vx;
    particle.y += dt * particle.vy;
    particle.vx += dt * particle.ax;
    particle.vy += dt * particle.ay;
    particle.u += dt * particle.du;
}

// -------- Boundary conditions --------

// Periodic Boundary Conditions
pub fn periodic_boundary(particle: &mut Particle, w: f64, h: f64){
    // We assume that the domain's system is a rectangular box.
    if particle.x > w {
        particle.x -= w;
    } else if particle.x < 0.0 {
        particle.x += w;
    }
    if particle.y > h {
        particle.y -= h;
    } else if particle.y < 0.0 {
        particle.y += h;
    }
}