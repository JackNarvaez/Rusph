use std::f64;

use rayon::prelude::*;

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

// -------- Basic vector functions --------

// Periodic relative distance
pub fn periodic_rel_vector(p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> (f64, f64, f64) {
    
    let mut x_temp: f64 = p1.x - p2.x;
    let mut y_temp: f64 = p1.y - p2.y;
    let mut z_temp: f64 = p1.z - p2.z;
    let twoeps: f64 = 2.*eps;
    if xperiodic && x_temp.abs() > wd-twoeps {
        if x_temp > 0. {
            x_temp -= wd;
        } else {
            x_temp += wd;
        }
    }
    if yperiodic && y_temp.abs() > lg-twoeps {
        if y_temp > 0. {
            y_temp -= lg;
        } else {
            y_temp += lg;
        }
    }
    if zperiodic && z_temp.abs() > hg-twoeps {
        if z_temp > 0. {
            z_temp -= hg;
        } else {
            z_temp += hg;
        }
    }
    return (x_temp, y_temp, z_temp);
}

// Periodic Distance
pub fn periodic_norm(p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> f64 { 
    let (x_temp, y_temp, z_temp) = periodic_rel_vector(p1, p2, wd, lg, hg, eps, xperiodic, yperiodic, zperiodic);
    return (x_temp*x_temp + y_temp*y_temp + z_temp*z_temp).sqrt();
}


// -------- Kernel functions --------

// *** B-Spline Kernels *** //

// Cubic Kernel
pub fn f_cubic_kernel(q:f64) -> f64 {
    if q < 1. {
        return 1. + 0.75*q*q*(q-2.); 
    } else if q < 2.{
        let f1: f64 = 2.-q;
        return 0.25*f1*f1*f1;
    } else {
        return 0.;
    }
}

// Derivative of cubic kernel
pub fn dfdq_cubic_kernel(q:f64) -> f64 {
    if q < 1. {
        return (2.25*q-3.)*q;
    } else if q < 2.{
        return -3.*(0.25*q*q-q+1.);
    } else {
        return 0.;
    }
}

// Quintic Kernel
pub fn f_quintic_kernel(q:f64) -> f64 {
    let f1: f64 = 3.-q;
    let f1sq: f64 = f1*f1;
    if q < 1. {
        let f2: f64 = 2.-q;
        let f3: f64 = 1.-q;
        let f2sq: f64 = f2*f2;
        let f3sq: f64 = f3*f3;
        return f1sq*f1sq*f1-6.*f2sq*f2sq*f2+15.*f3sq*f3sq*f3; 
    } else if q < 2.{
        let f2: f64 = 2.-q;
        let f2sq: f64 = f2*f2;
        return f1sq*f1sq*f1-6.*f2sq*f2sq*f2;
    } else if q < 3. {
        return f1sq*f1sq*f1;
    } else {
        return 0.;
    }
}

// Derivative of quintic kernel
pub fn dfdq_quintic_kernel(q:f64) -> f64 {
    let f1: f64 = 3.-q;
    let f1sq: f64 = f1*f1;
    if q < 1. {
        let f2: f64 = 2.-q;
        let f3: f64 = 1.-q;
        let f2sq: f64 = f2*f2;
        let f3sq: f64 = f3*f3;
        return 5.*(-f1sq*f1sq+6.*f2sq*f2sq-15.*f3sq*f3sq); 
    } else if q < 2.{
        let f2: f64 = 2.-q;
        let f2sq: f64 = f2*f2;
        return 5.*(-f1sq*f1sq+6.*f2sq*f2sq);
    } else if q < 3. {
        return -5.*f1sq*f1sq;
    } else {
        return 0.;
    }
}

// *** Gaussian Kernels *** //

// Gaussian Kernel
pub fn f_gaussian_kernel(q:f64) -> f64 {
    if q < 3. {
        return (-q*q).exp();
    } else {
        return 0.;
    }
}

// Derivative of Gaussian Kernel
pub fn dfdq_gaussian_kernel(q:f64) -> f64 {
    if q < 3. {
        return -2.0*q*(-q*q).exp();
    } else {
        return 0.;
    }
}

// *** Wendland Kernels *** //

// C2 Wendland kernel
pub fn f_c2wendland_kernel(q:f64) -> f64 {
    if q < 2. {
        let f1: f64 = 1.-0.5*q;
        return f1*f1*f1*f1*(2.*q+1.);
    } else {
        return 0.;
    }
}

// Derivative of C2 Wendland kernel
pub fn dfdq_c2wendland_kernel(q:f64) -> f64 {
    if q < 2.{
        let f1: f64 = 1.-0.5*q;
        return -5.*q*f1*f1*f1;
    } else {
        return 0.;
    }
}

// ----------------------------------

// Derivative of kernel w.r.t the smoothing length
pub fn dwdh(q: f64, f: fn(f64) -> f64, df: fn(f64) -> f64) -> f64 {
    3. *f(q) + q*df(q)
}


// -------- Kernel approximations --------

// Kernel approximation of density
pub fn density_kernel(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, sigma:f64, rkern: f64, f: fn(f64)->f64, wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> f64 {
    let mut rho :f64 = 0.0;
    for jj in neigh_particles{
        let r = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
        rho += f(r/h);
    }
    rho * dm * sigma / (h*h*h)
}

// Density calculated from smoothing length
pub fn density_by_smoothing_length(dm:f64, h:f64, eta:f64) -> f64{
    let vol: f64 = eta/h;
    dm*vol*vol*vol
}

// Smoothing length calculated from density number
pub fn h_by_density(dm:f64, rho:f64, eta:f64) -> f64{
    eta*(dm/rho).cbrt()
}

// Omega operator
pub fn omega(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64, wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> f64{
    let mut omeg :f64 = 0.0;
    for jj in neigh_particles {
        let q = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*h, xperiodic, yperiodic, zperiodic)/h;
        omeg -= dwdh_(q, f, dfdq);
    }
    omeg *= dm*sigma/(3.*h*h*h*rho);
    return omeg + 1.;
}


// -------- Root solver --------

// -- Newton-Raphson iterator --

// function and derivative of function
pub fn f_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> (f64 , f64) {
    let rho_kernel = density_kernel(particles, ii, neigh_particles, dm, h, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    let rho_h = density_by_smoothing_length(dm, h, eta);
    let f_h = rho_h - rho_kernel;
    let omeg = omega(particles, ii, neigh_particles, dm, h, rho_h, dwdh, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    let df = -3.*rho_h*omeg/ h;
    (f_h, df)
}

// Calculate a new value of 'h'
fn nr_iter(particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h_old: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> f64 {
    let (f, df) = f_iter(particles, ii, neigh_particles, dm, h_old, eta, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    h_old - f / df
}

// Newton raphson solver to find the value of 'h' for particle 'ii'
pub fn newton_raphson(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, tol: f64, it: u32, tree: &Node, s_: i32, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> (f64, Vec<usize>) {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    let mut neighbors: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighboring particles
        neighbors.clear();
        tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, h_old*rkern, xperiodic, yperiodic, zperiodic);
        // Obtain h_new
        h_new = nr_iter(particles, ii, &neighbors, dm, h_old, eta, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        if h_new > 1.2*particles[ii].h {
            h_new = 1.2*particles[ii].h;
        } else if h_new < 0.8*particles[ii].h {
            h_new = 0.8*particles[ii].h;
        }
        if ((h_new - h_old)/h_old).abs() <=  tol {
            i = it + 2;
        } else if h_new < 0. || h_new > wd {
            i = it + 1;
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
pub fn bisection(ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64, f: fn(f64) -> f64, sigma:f64, rkern: f64, tol: f64, it: u32, tree: &Node, s_: i32, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> (f64, Vec<usize>) {
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
        tree.find_neighbors(ii, s_, particles, &mut neighbors_left, wd, lg, hg, x0, y0, z0, h_left*rkern, xperiodic, yperiodic, zperiodic);
        tree.find_neighbors(ii, s_, particles, &mut neighbors_right, wd, lg, hg, x0, y0, z0, h_right*rkern, xperiodic, yperiodic, zperiodic);
        tree.find_neighbors(ii, s_, particles, &mut neighbors_mid, wd, lg, hg, x0, y0, z0, h_mid*rkern, xperiodic, yperiodic, zperiodic);

        let f_left = density_by_smoothing_length(dm, h_left, eta) - density_kernel(particles, ii, &neighbors_left, dm, h_left, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        let f_right = density_by_smoothing_length(dm, h_right, eta) - density_kernel(particles, ii, &neighbors_right, dm, h_right, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        let f_mid = density_by_smoothing_length(dm, h_mid, eta) - density_kernel(particles, ii, &neighbors_mid, dm, h_mid, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);

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
pub fn smoothing_length(particles: &mut Vec<Particle>, dm:f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, tol: f64, it: u32, dt:f64, tree: &Node, s_: i32, n: usize, ptr : Pointer, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool){
    (0..n).into_par_iter().for_each(|ii| {
        if particles[ii].ptype==0 {
            let (mut h_new, mut neighbors) = newton_raphson(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/3.), eta, f, dfdq, sigma, rkern, tol, it, tree, s_, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
            let particle = unsafe { &mut *{ptr}.0.add(ii)};
            if h_new != 0.0 {
                // If h is not found, then keep it constant in time.
                particle.h = h_new;
            } else {
                (h_new, neighbors) = bisection(ii, particles, dm, particles[ii].h*(1.+dt*particles[ii].divv/3.), eta, f, sigma, rkern, tol, it, tree, s_, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
                if h_new != 0.0 {
                    // If h is not found, then keep it constant in time.
                    particle.h = h_new;
                } else {
                    neighbors.clear();
                    tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particle.h*rkern, xperiodic, yperiodic, zperiodic);
                }
            }
         particle.rho = density_kernel(particles, ii, &neighbors, dm, particle.h, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        }
    });
}


// -------- Equations of state --------

// -- Toy Star 2D --

// Polytropic equation
pub fn eos_polytropic(rho:f64, _k:f64, gamma:f64) -> f64 {
    let k: f64 = 0.05;
    k * rho.powf(gamma)
}

pub fn sound_speed_polytropic(rho: f64, gamma: f64) -> f64 {
    let k: f64 = 0.05;
    (gamma*k*rho.powf(gamma-1.)).sqrt()
}

// -- Ideal Gas --

pub fn eos_ideal_gas(rho:f64, u:f64, gamma:f64) -> f64 {
    (gamma-1.)*rho*u
}

// Sound speed for the ideal gas eos
pub fn sound_speed_ideal_gas_u(u:f64, gamma:f64) -> f64 {
    ((gamma-1.)*gamma*u).sqrt()
}

pub fn sound_speed_ideal_gas(rho:f64, p:f64, gamma:f64) -> f64 {
    (gamma*p/rho).sqrt()
}

// -- Isotheraml EoS

pub fn eos_isothermal(rho:f64, _u:f64, gamma:f64) -> f64 {
    rho.powi(gamma as i32)
}

pub fn sound_speed_isothermal(rho:f64, _u:f64, gamma:f64) -> f64 {
    gamma*rho.powi((gamma-1.) as i32)
}

pub fn sound_speed_isothermal_dt(rho:f64, gamma:f64) -> f64 {
    gamma*rho.powi((gamma-1.) as i32)
}


// -------- Artificial Viscosity --------

// Monaghan (1989): "Standard" SPH viscous term
pub fn mon89_art_vis(r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, h_i: f64, h_j: f64, rho_mean: f64) -> (f64, f64) {
    if dot_r_v <= 0. {
        // Mean values
        let cs_mean :f64 = 0.5*(cs_i+cs_j);
        let h_mean :f64 = 0.5*(h_i+h_j);
    
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
pub fn mon97_art_vis(r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, _h_i: f64, _h_j: f64, rho_mean: f64) -> (f64, f64) {
    if dot_r_v <= 0. {
        // Parameters
        let alpha: f64 = 1.0;
        let beta: f64 = 2.0;

        let v_sig:f64 = 0.5*alpha*(cs_i + cs_j - beta*dot_r_v/r_ij);
        let dvdt :f64 = -v_sig*dot_r_v/(r_ij*rho_mean);
        let dudt :f64 = 0.5*dvdt*dot_r_v;
        return (dvdt, dudt);
    } else {
        return (0.0, 0.0);
    }
}

// Price (2008): Thermal conductivity switches
pub fn price08_therm_cond(p_i: f64, p_j: f64, rho_mean: f64, u_i: f64, u_j: f64) -> f64 {
    // Parameters
    let alpha_u: f64 = 1.;

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
pub fn body_forces_null(_particles: &mut Particle, _m_star: & Star) {
}

// Body forces for a toy star in 2D
pub fn body_forces_toy_star(particle: &mut Particle, _m_star: &Star) {
    let nu: f64 = 1.0;
    let lmbda: f64 = 0.5030082152040148;
    particle.ax -= nu * particle.vx + lmbda*particle.x;
    particle.ay -= nu * particle.vy + lmbda*particle.y; 
    particle.az -= nu * particle.vz + lmbda*particle.z; 
}

// Newtonian Gravitational Force
pub fn body_forces_gravitation(particle: &mut Particle, m_star: & Star) {
    let x_r: f64 = particle.x - m_star.x;
    let y_r: f64 = particle.y - m_star.y;
    let z_r: f64 = particle.z - m_star.z;
    let f_grav: f64 = m_star.m * (x_r*x_r + y_r*y_r + z_r*z_r+0.000625).powf(-1.5);
    particle.ax += -f_grav*x_r;
    particle.ay += -f_grav*y_r;
    particle.az += -f_grav*z_r;
}

// Calculate acceleration for each particle in the system
pub fn accelerations(particles: &mut Vec<Particle>, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                     dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                     tree: &Node, s_: i32, n: usize, ptr : Pointer, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64,
                     artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                     body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool, xperiodic: bool, yperiodic:bool, zperiodic:bool){
    
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[ii].h*rkern, xperiodic, yperiodic, zperiodic);
        return neighbors;
    }).collect();
    (0..n).into_par_iter().for_each(move |ii| {
        if particles[ii].ptype==0 {

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
            let omeg_i = omega(particles, ii, &neighbors[ii], dm, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
            for jj in 0..n {
                if ii != jj {
                    let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
                    let r_ij = (x_rel*x_rel + y_rel*y_rel+ z_rel*z_rel).sqrt();

                    let mut grad_hi = 0.0;
                    let mut grad_hj = 0.0;
                    if r_ij <= rkern*particles[ii].h {
                        let hisq = particles[ii].h*particles[ii].h;
                        grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*hisq*hisq);
                    }
                    if r_ij <= rkern*particles[jj].h {
                        let hjsq: f64 = particles[jj].h*particles[jj].h;
                        grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*hjsq*hjsq);
                    }
                    if grad_hi != 0. || grad_hj != 0.0 {
                        let p_j = eos(particles[jj].rho, particles[jj].u, gamma);
                        let cs_j = cs(particles[jj].rho, p_j, gamma);
                        let omeg_j = omega(particles, jj, &neighbors[jj], dm, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);


                        // Velocity dot position
                        let dot_r_v = (particles[ii].vx-particles[jj].vx)*x_rel
                        +(particles[ii].vy-particles[jj].vy)*y_rel
                        +(particles[ii].vz-particles[jj].vz)*z_rel;

                        // Mean density
                        let rho_mean: f64 = 0.5 * (particles[ii].rho + particles[jj].rho);
                        // Artificial viscosity
                        let (art_visc_mom, art_visc_ene) = artificial_viscosity(r_ij, dot_r_v, cs_i, cs_j, rho_mean, particles[ii].rho, particles[jj].rho);

                        // Artificial thermal conductivity
                        let art_therm_cond: f64 = price08_therm_cond(p_i, p_j, rho_mean, particles[ii].u, particles[jj].u);

                        // Acceleration
                        let (f_ij_x, f_ij_y, f_ij_z) = acceleration_ab(&particles[ii], &particles[jj], x_rel, y_rel, z_rel, p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc_mom);
                        particle_i.ax += dm * f_ij_x;
                        particle_i.ay += dm * f_ij_y;
                        particle_i.az += dm * f_ij_z;

                        // Divergence of v per unit of mass
                        let div_vel :f64 = grad_hi*dot_r_v / (omeg_i*particles[ii].rho);
                        particle_i.divv -= dm*div_vel;

                        // Thermal change
                        particle_i.du += dm * (0.5*(art_visc_ene + art_therm_cond*r_ij)*(grad_hi/omeg_i+grad_hj/omeg_j));
                    }
                }
            }
            particle_i.du -= (p_i/particles[ii].rho)*particle_i.divv;
            // Body forces
            if bf {
                body_forces(particle_i, star);
            }
        }
    });
}


// -------- Time integrator --------
pub fn euler_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                        dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                        eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                        artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                        body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
                        boundary: fn(&mut Vec<Particle>, f64, f64,f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 10, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.x += dt * particle.vx;
            particle.y += dt * particle.vy;
            particle.z += dt * particle.vz;
            particle.vx += dt * particle.ax;
            particle.vy += dt * particle.ay;
            particle.vz += dt * particle.az;
            particle.u += dt * particle.du;
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
}

// Velocity Verlet integrator
pub fn velocity_verlet_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                                  eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;

            particle.u += 0.5 * dt * particle.du;

            particle.x += dt * particle.vx;
            particle.y += dt * particle.vy;
            particle.z += dt * particle.vz;
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
    
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 30, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);

    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;
            particle.u += 0.5 * dt * particle.du;
        }
    });
}

// Velocity Verlet integrator
pub fn predictor_kdk_integrator(particles: &mut Vec<Particle>, dt:f64, dm:f64, eos: fn(f64, f64, f64)->f64, cs: fn(f64, f64, f64)->f64, gamma:f64,
                                  dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
                                  eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
                                  artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
                                  body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
                                  boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64) {
    
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
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
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 30, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
    accelerations(particles, dm, eos, cs, gamma, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx = particle.vx_star + 0.5 * dt * particle.ax;
            particle.vy = particle.vy_star + 0.5 * dt * particle.ay;
            particle.vz = particle.vz_star + 0.5 * dt * particle.az;
            particle.u = particle.u_star + 0.5 * dt * particle.du;
        }
    });
}


// -------- Boundary conditions --------

// Periodic Boundary Conditions
pub fn periodic_boundary(particles: &mut Vec<Particle>, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0: f64){
    // We assume that the domain's system is a rectangular box.
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            if particle.x >= (wd+x0) {
                particle.x -= wd;
            } else if particle.x < x0 {
                particle.x += wd;
            }
            if particle.y >= (lg + y0) {
                particle.y -= lg;
            } else if particle.y < y0 {
                particle.y += lg;
            }
            if particle.z >= (hg + z0) {
                particle.z -= hg;
            } else if particle.z < z0 {
                particle.z += hg;
            }
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

// Monaghan (1989) Force conditon
pub fn force_dt(h: f64, a: f64, f: f64) -> f64 {
    f*(h/a).sqrt()
}


// Timestepping Criteria Cossins P. J. (2010)
pub fn time_step_bale(particles: & Vec<Particle>, n: usize, gamma: f64, _rkern: f64, _wd: f64, _lg: f64, _hg: f64,
                      _tree: &mut Node, _s_: i32, cs: fn(f64, f64) -> f64) -> f64 {

    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let cs: f64 = cs(particles[ii].u, gamma);
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// Timestepping Criteria Monaghan (1997)
pub fn time_step_mon(particles: & Vec<Particle>, n: usize, gamma: f64, rkern: f64, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64,
                     tree: &mut Node, s_: i32, cs: fn(f64, f64) -> f64, xperiodic: bool, yperiodic:bool, zperiodic:bool) -> f64 {
                        
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[ii].h*rkern, xperiodic, yperiodic, zperiodic);
        return neighbors;
    }).collect();
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        let alpha: f64 = 1.;
        let beta: f64 = 2.;
        let mut v_sig:f64 = 0.0;
        let cs_i: f64 = cs(particles[ii].u, gamma);
        for jj in &neighbors[ii] {
            let cs_j: f64 = cs(particles[*jj].u, gamma);

            // Velocity dot position
            let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
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