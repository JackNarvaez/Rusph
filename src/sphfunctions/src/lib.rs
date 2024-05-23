// ------------------------------------------------------------------------- //
// Main functions used by Rusph.                                             //
// It includes:                                                              //
//      Vector functions                                                     //
//      SPH functions                                                        //
//      External forces                                                      //
//      Force calculation                                                    //
//      Time stepping                                                        //
// ------------------------------------------------------------------------- //

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


// !!!---------------------- Basic Vector Functions ---------------------!!! //

// ------------------------------------------------------------------------- //
// Periodic relative distance:                                               //
// Returns                                                                   //
//      r(p1)-r(p2) = r_12 = (x_12, y_12, z_12),                             //
// the relative distance between particles p1 and p2.                        //
// Periodic boundary conditions are included using booleans rperiodic.       //
// ------------------------------------------------------------------------- //
pub fn periodic_rel_vector(
    p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64,
    xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> (f64, f64, f64) {
    
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

// ------------------------------------------------------------------------- //
// Norm of periodic distance:                                                //
// Returns                                                                   //
//      |r(p1)-r(p2)| = |r_12| = sqrt(x_12^2 + y_12^2 + z_12^2),             //
// the norm of relative distance between particles p1 and p2.                //
// Periodic boundary conditions are included using booleans rperiodic.       //
// ------------------------------------------------------------------------- //
pub fn periodic_norm(
    p1: &Particle, p2: &Particle, wd: f64, lg: f64, hg: f64, eps: f64,
    xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> f64 { 
    let (x_temp, y_temp, z_temp) = periodic_rel_vector(p1, p2, wd, lg, hg, eps, xperiodic, yperiodic, zperiodic);
    return (x_temp*x_temp + y_temp*y_temp + z_temp*z_temp).sqrt();
}

// ------------------------------------------------------------------------- //
// Norm of non-periodic distance:                                            //
// Returns                                                                   //
//      |r(star)-r(p)| = |r_12| = sqrt(x_12^2 + y_12^2 + z_12^2),            //
// the norm of relative distance between star and particle p.                //
// ------------------------------------------------------------------------- //
pub fn distance_star(
    star: & Star, p: &Particle
) -> f64 { 
    let dx: f64 = star.x - p.x;
    let dy: f64 = star.y - p.y;
    let dz: f64 = star.z - p.z;
    return (dx*dx + dy*dy + dz*dz).sqrt();
}


// !!!------------------------- Kernel Functions ------------------------!!! //

// ***------------------------- B-Spline Kernels ------------------------*** //

// ------------------------------------------------------------------------- //
// Cubic kernel:                                                             //
// Returns the M4 cubic spline.                                              //
// Monaghan & Lattanzio (1985)                                               //
// ------------------------------------------------------------------------- //
pub fn f_cubic_kernel(
    q:f64
) -> f64 {
    if q < 1. {
        return 1. + 0.75*q*q*(q-2.); 
    } else if q < 2.{
        let f1: f64 = 2.-q;
        return 0.25*f1*f1*f1;
    } else {
        return 0.;
    }
}

// ------------------------------------------------------------------------- //
// Derivative of the cubic kernel:                                           //
// Returns df(q)/dq for the M4 cubic spline.                                 //
// Monaghan & Lattanzio (1985)                                               //
// ------------------------------------------------------------------------- //
pub fn dfdq_cubic_kernel(
    q:f64
) -> f64 {
    if q < 1. {
        return (2.25*q-3.)*q;
    } else if q < 2.{
        return -3.*(0.25*q*q-q+1.);
    } else {
        return 0.;
    }
}

// ------------------------------------------------------------------------- //
// Quintic kernel:                                                           //
// Returns the M6 quintic spline.                                            //
// ------------------------------------------------------------------------- //
pub fn f_quintic_kernel(
    q:f64
) -> f64 {
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

// ------------------------------------------------------------------------- //
// Derivative of the quintic kernel:                                         //
// Returns df(q)/dq for the M6 quintic spline.                               //
// ------------------------------------------------------------------------- //
pub fn dfdq_quintic_kernel(
    q:f64
) -> f64 {
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

// ***------------------------- Gaussian Kernels ------------------------*** //

// ------------------------------------------------------------------------- //
// Gaussian kernel:                                                          //
// Returns the Gaussian kernel with support radius q < 3.                    //
// Liu (2010).                                                               //
// ------------------------------------------------------------------------- //
pub fn f_gaussian_kernel(
    q:f64
) -> f64 {
    if q < 3. {
        return (-q*q).exp();
    } else {
        return 0.;
    }
}

// ------------------------------------------------------------------------- //
// Derivative of the Gaussian kernel:                                        //
// Returns df(q)/dq for the Gaussian kernel with support radius q < 3.       //
// ------------------------------------------------------------------------- //
pub fn dfdq_gaussian_kernel(
    q:f64
) -> f64 {
    if q < 3. {
        return -2.0*q*(-q*q).exp();
    } else {
        return 0.;
    }
}

// ***------------------------- Wendland Kernels ------------------------*** //

// ------------------------------------------------------------------------- //
// C2 Wendland kernel:                                                       //
// Returns the C2 Wendland kernel scaled to a radius of 2h.                  //
// Wendland (1995)                                                           //
// ------------------------------------------------------------------------- //
pub fn f_c2wendland_kernel(
    q:f64
) -> f64 {
    if q < 2. {
        let f1: f64 = 1.-0.5*q;
        return f1*f1*f1*f1*(2.*q+1.);
    } else {
        return 0.;
    }
}

// ------------------------------------------------------------------------- //
// Derivative of te C2 Wendland kernel:                                      //
// Returns Returns df(q)/dq for the C2 Wendland kernel.                      //
// Wendland (1995)                                                           //
// ------------------------------------------------------------------------- //
pub fn dfdq_c2wendland_kernel(
    q:f64
) -> f64 {
    if q < 2.{
        let f1: f64 = 1.-0.5*q;
        return -5.*q*f1*f1*f1;
    } else {
        return 0.;
    }
}


// !!!-------------------- Kernel Derivative w.r.t. h -------------------!!! //

// ------------------------------------------------------------------------- //
// Returns                                                                   //
//      dW/dh ~ 3f(q)+qf'(q),                                                //
// the derivative of kernel function w.r.t the smoothing length.             //
// Here, f(q) is the kernel and f'(q) = df/dq.                               //
// ------------------------------------------------------------------------- //
pub fn dwdh(
    q: f64, f: fn(f64) -> f64, df: fn(f64) -> f64
) -> f64 {
    3. *f(q) + q*df(q)
}


// !!!----------------------------- ρ and h -----------------------------!!! //

// ------------------------------------------------------------------------- //
// Returns                                                                   //
//      rho(h) = m * (eta/h)^3,                                              //
// the density calculated from the smoothing length.                         //
// ------------------------------------------------------------------------- //
pub fn density_from_h(
    dm:f64, h:f64, eta:f64
) -> f64{
    let vol: f64 = eta/h;
    dm*vol*vol*vol
}

// ------------------------------------------------------------------------- //
// Returns                                                                   //
//      h(rho) = eta * (m/rho)^(1/3),                                        //
// the smoothing length calculated from the density number.                  //
// ------------------------------------------------------------------------- //
pub fn h_from_density(
    dm:f64, rho:f64, eta:f64
) -> f64{
    eta*(dm/rho).cbrt()
}


// !!!---------------------- Kernel Approximations ----------------------!!! //

// ------------------------------------------------------------------------- //
// Returns                                                                   //
//      rho_a = SUM_b m_b W(r_ab, h_a),                                      //
// the usual SPH density sum.                                                //
// ------------------------------------------------------------------------- //
pub fn density_kernel(
    particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64,
    sigma:f64, rkern: f64, f: fn(f64)->f64, wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> f64 {
    let mut rho :f64 = 0.0;
    for jj in neigh_particles{
        let r: f64 = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
        rho += f(r/h);
    }
    rho * dm * sigma / (h*h*h)
}

// ------------------------------------------------------------------------- //
// Returns the Omega operator, which is related to the gradient of the       //
// smoothing length.                                                         //
// Monaghan (2002)                                                           //
// ------------------------------------------------------------------------- //
pub fn omega(
    particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, rho: f64,
    dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
    wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> f64{
    let mut omeg :f64 = 0.0;
    for jj in neigh_particles {
        let q: f64 = periodic_norm(&particles[ii], &particles[*jj], wd, lg, hg, rkern*h, xperiodic, yperiodic, zperiodic)/h;
        omeg -= dwdh_(q, f, dfdq);
    }
    omeg *= dm*sigma/(3.*h*h*h*rho);
    if omeg <= -1.0 {
        omeg = 0.0;
    }
    return omeg + 1.;
}


// !!!--------------------- Root-Finding Algorithms ---------------------!!! //

// ***--------------------- Newton-Raphson iterator ---------------------*** //

// ------------------------------------------------------------------------- //
// Returns                                                                   //
//      f(h_a) = rho(h_a) - rho_sph(h_a)                                     //
//      df(h_a) = df(h_a)/dh_a = -3*rho*OMEGA/h_a                            //
// the function to find the root and its derivative.                         //
// ------------------------------------------------------------------------- //
pub fn f_iter(
    particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h: f64, eta:f64,
    f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64,
    wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> (f64 , f64) {
    let rho_kernel: f64 = density_kernel(particles, ii, neigh_particles, dm, h, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    let rho_h: f64  = density_from_h(dm, h, eta);
    let f_h: f64    = rho_h - rho_kernel;
    let omeg: f64   = omega(particles, ii, neigh_particles, dm, h, rho_h, dwdh, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    let df: f64     = -3.*rho_h*omeg/ h;
    (f_h, df)
}

// ------------------------------------------------------------------------- //
// Iterates the smoothing function.                                          //
// Returns                                                                   //
//      h_new = h_old - f´(h_old)/f(h_old)                                   //
// the new value of h after one iteration.                                   //
// ------------------------------------------------------------------------- //
fn nr_iter(
    particles: & Vec<Particle>, ii:usize, neigh_particles: & Vec<usize>, dm:f64, h_old: f64, eta:f64,
    f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64,
    wd: f64, lg: f64, hg: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> f64 {
    let (f_i, df_i) = f_iter(particles, ii, neigh_particles, dm, h_old, eta, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    h_old - f_i / df_i
}

// ------------------------------------------------------------------------- //
// Newton Raphson solver: finds the value of 'h' for the iith- particle.     //
// Returns                                                                   //
//      h_ii                                                                 //
// the smoothing length of th iith-particle. It iterates the function:       //
//      f(h) = rho(h) - rho_sph(h)                                           //
// ------------------------------------------------------------------------- //
pub fn newton_raphson(
    ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64,
    f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64, tol: f64, it: u32, tree: &Node, s_: i32,
    wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> (f64, Vec<usize>) {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    let mut neighbors: Vec<usize> = Vec::new();
    while i <= it {
        // Searching neighbouring particles
        neighbors.clear();
        tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, h_old*rkern, xperiodic, yperiodic, zperiodic);
        // Obtain h_new
        h_new = nr_iter(particles, ii, &neighbors, dm, h_old, eta, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        
        // Restrict result to [0.8h_old, 1.2h_old]
        if h_new > 1.2*particles[ii].h {
            h_new = 1.2*particles[ii].h;
        } else if h_new < 0.8*particles[ii].h {
            h_new = 0.8*particles[ii].h;
        }

        if ((h_new - h_old)/h_old).abs() <=  tol {
            // Success
            i = it + 2;
        } else if h_new < 0. || h_new > wd {
            // No root found
            i = it + 1;
        } else{
            // Continue
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

// ***----------------------- Bisection iterator ------------------------*** //

// ------------------------------------------------------------------------- //
// Bisection solver: finds the value of 'h' for the iith- particle.          //
// Returns                                                                   //
//      h_ii                                                                 //
// the smoothing length of th iith-particle. It iterates the function:       //
//      f(h) = rho(h) - rho_sph(h)                                           //
// ------------------------------------------------------------------------- //
pub fn bisection(
    ii: usize, particles: & Vec<Particle>, dm:f64, h_guess: f64, eta:f64,
    f: fn(f64) -> f64, sigma:f64, rkern: f64, tol: f64, it: u32, tree: &Node, s_: i32,
    wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> (f64, Vec<usize>) {
    let mut h_left :f64 = 0.6*h_guess;
    let mut h_right :f64= 1.667*h_guess;
    let mut h_mid: f64  = 0.5*(h_left+h_right);
    let mut i : u32     = 1;

    // f_left
    let mut neighbors_left: Vec<usize> = Vec::new();
    neighbors_left.clear();
    tree.find_neighbors(ii, s_, particles, &mut neighbors_left, wd, lg, hg, x0, y0, z0, h_left*rkern, xperiodic, yperiodic, zperiodic);
    let mut f_left: f64 = density_from_h(dm, h_left, eta) - density_kernel(particles, ii, &neighbors_left, dm, h_left, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);

    // f_right
    let mut neighbors_mid: Vec<usize> = Vec::new();
    neighbors_mid.clear();
    tree.find_neighbors(ii, s_, particles, &mut neighbors_mid, wd, lg, hg, x0, y0, z0, h_right*rkern, xperiodic, yperiodic, zperiodic);
    let mut f_mid: f64  = density_from_h(dm, h_right, eta) - density_kernel(particles, ii, &neighbors_mid, dm, h_right, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
    
    if f_mid*f_left > 0.0 {
        i = it + 1;
    }

    while i <= it {
        h_mid = 0.5*(h_left+h_right);

        // f_middle
        neighbors_mid.clear();
        tree.find_neighbors(ii, s_, particles, &mut neighbors_mid, wd, lg, hg, x0, y0, z0, h_mid*rkern, xperiodic, yperiodic, zperiodic);
        f_mid  = density_from_h(dm, h_mid, eta) - density_kernel(particles, ii, &neighbors_mid, dm, h_mid, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);

        if ((h_right - h_left)/h_mid).abs() <=  tol  {
            i = it + 2;
        } else{
            i += 1;
            if f_mid*f_left < 0.0 {
                h_right = h_mid;
            } else {
                h_left = h_mid;
                f_left = f_mid;
            }
        }
    }
    if i == it+1 {
        println!("PROBLEM BISECTION METHOD: Solution not found for particle {}. Final value was h = {}.", ii, h_mid);
        (0.0, neighbors_mid)
    } else{
        (h_mid, neighbors_mid)
    }
}


// !!!------------------ Calculating Smoothing Length -------------------!!! //

// ------------------------------------------------------------------------- //
// Calculate the smoothing length (h_ii) for every active gas particle       //
// It first uses the Newton Raphson solver. Then, if not root found, it      //
// uses bisection solver. Finally, if not root found, it keeps h constant    //
// ------------------------------------------------------------------------- //
pub fn smoothing_length(
    particles: &mut Vec<Particle>, dm:f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, rkern: f64,
    tol: f64, it: u32, dt:f64, tree: &Node, s_: i32, n: usize, ptr : Pointer,
    wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
){
    (0..n).into_par_iter().for_each(|ii| {
        if particles[ii].ptype==0 {
            let h_guess: f64 = particles[ii].h*(1.+dt*particles[ii].divv/3.);
            let (mut h_new, mut neighbors) = newton_raphson(ii, particles, dm, h_guess, eta, f, dfdq, sigma, rkern, tol, it, tree, s_, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
            let particle = unsafe { &mut *{ptr}.0.add(ii)};
            if h_new != 0.0 {
                particle.h = h_new;
            } else {
                (h_new, neighbors) = bisection(ii, particles, dm, h_guess, eta, f, sigma, rkern, tol, it, tree, s_, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
                if h_new != 0.0 {
                    particle.h = h_new;
                } else {
                    // h is not found, then keep it constant in time.
                    neighbors.clear();
                    tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particle.h*rkern, xperiodic, yperiodic, zperiodic);
                }
            }
         particle.rho = density_kernel(particles, ii, &neighbors, dm, particle.h, sigma, rkern, f, wd, lg, hg, xperiodic, yperiodic, zperiodic);
        }
    });
}

// !!!----------------------- Equations of State ------------------------!!! //

// ***------------------------- Polytropic EoS --------------------------*** //

// ------------------------------------------------------------------------- //
// The Polytropic equation of state.                                         //
// Returns the pressure                                                      //
//      P = K rho^gamm                                                       //
// where rho is the density, K is the polytropic constant, and gamm is the   //
// adiabatic index.                                                          //
// ------------------------------------------------------------------------- //
pub fn eos_polytropic(
    rho:f64, _u:f64, gamma:f64, _x: f64, _y: f64, _z: f64, coeff: f64
) -> f64 {
    coeff * rho.powf(gamma)
}

// ------------------------------------------------------------------------- //
// The Polytropic equation of state.                                         //
// Returns the speed of sound                                                //
//      Cs = sqrt[gamm * K rho^(gamm-1)]                                     //
// where rho is the density, K is the polytropic constant, and gamm is the   //
// adiabatic index.                                                          //
// ------------------------------------------------------------------------- //
pub fn sound_speed_polytropic(
    rho: f64, _u:f64, gamma: f64, _x: f64, _y: f64, _z: f64, coeff: f64
) -> f64 {
    (gamma*coeff*rho.powf(gamma-1.)).sqrt()
}


// ***---------------------------- Ideal Gas ----------------------------*** //

// ------------------------------------------------------------------------- //
// The equation of state for an ideal gas.                                   //
// Returns the pressure                                                      //
//      P = (gamm-1)*rho*u                                                   //
// where rho is the density, u is the internal energy, and gamm is the       //
// adiabatic index.                                                          //
// ------------------------------------------------------------------------- //
pub fn eos_ideal_gas(
    rho:f64, u:f64, gamma:f64, _x: f64, _y: f64, _z: f64, _coeff: f64
) -> f64 {
    (gamma-1.)*rho*u
}

// ------------------------------------------------------------------------- //
// The equation of state for an ideal gas.                                   //
// Returns the speed of sound                                                //
//      Cs = sqrt[gamm*(gamm-1)*u]                                           //
// where u is the internal energy and gamm is the adiabatic index.           //
// ------------------------------------------------------------------------- //
pub fn sound_speed_ideal_gas(
    _rho: f64, u:f64, gamma:f64, _x: f64, _y: f64, _z: f64, _coeff: f64
) -> f64 {
    ((gamma-1.)*gamma*u).sqrt()
}


// ***----------------------- Locally Isothermal ------------------------*** //

// ------------------------------------------------------------------------- //
// Locally isothermal equation of state                                      //
// Returns the presure                                                       //
//      P = Cs^2 * rho                                                       //
// where Cs is the sound speed at R and rho is the density.                  //
// Lodato & Pringle (2007)                                                   //
// ------------------------------------------------------------------------- //
pub fn eos_isothermal_disc(
    rho:f64, _u:f64, gamma:f64, x: f64, y: f64, z: f64, cs02: f64
) -> f64 {
    cs02 * (x*x + y*y + z*z).powf(-gamma)*rho
}

// ------------------------------------------------------------------------- //
//Locally isothermal equation of state                                       //
// Returns the speed of sound                                                //
//      Cs = Cs_0 R^(-q)                                                     //
// where Cs_o is the sound speed at R_in, R = sqrt[x^2 + y^2 + z^2], and     //
// q is a constant index.                                                    //
// ------------------------------------------------------------------------- //
pub fn sound_speed_isothermal_disc(
    _rho: f64, _u:f64, gamma: f64, x: f64, y: f64, z: f64, cs02: f64
) -> f64 {
    (cs02*(x*x + y*y + z*z).powf(-gamma)).sqrt()
}

// !!!---------------------- Artificial Viscosity -----------------------!!! //

// ------------------------------------------------------------------------- //
// Artificial Viscosity proposed by Monaghan (1992)                          //
// Returns:                                                                  //
//      [PI_ab, 0.5*PI_ab*(V_ab \cdot r_ab)] if V_ab \cdot r_ab < 0          //
//      [0.0, 0.0]                           if V_ab \cdot r_ab > 0          //
// i.e., the contribution from artificial viscosity to the momentum and      //
// energy equations. Note that this term have to be scaled by the kernel     //
// gradient and the particle mass.                                           //
// Note: dot_r_v = \vec_{v}_ab \cdot \vec{r}_ab. r_ab is not the unit vector //
// ------------------------------------------------------------------------- //
pub fn mon92_art_vis(
    r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, h_i: f64, h_j: f64, rho_mean: f64
) -> (f64, f64) {
    if dot_r_v <= 0. {
        //particles approaching
        // Mean values
        let cs_mean: f64= 0.5*(cs_i+cs_j);
        let h_mean: f64 = 0.5*(h_i+h_j);
    
        // Parameters
        let alpha: f64  = 1.0;
        let beta: f64   = 2.0;
        let eps: f64    = 0.01;
        let nu_visc: f64= h_mean*dot_r_v/(r_ij*r_ij+eps*h_mean*h_mean);
        let dvdt: f64   = (-alpha*cs_mean+beta*nu_visc)*nu_visc/rho_mean;
        let dudt :f64   = 0.5*dvdt*dot_r_v;
    
        return (dvdt, dudt);
    } else {
        // particles receding
        return (0.0, 0.0);
    }
}

// ------------------------------------------------------------------------- //
// Artificial Viscosity proposed by Monaghan (1997) - Riemann Solvers        //
// Returns:                                                                  //
//      [PI_ab, 0.5*PI_ab*(V_ab \cdot r_ab)] if V_ab \cdot r_ab < 0          //
//      [0.0, 0.0]                           if V_ab \cdot r_ab > 0          //
// i.e., the contribution from artificial viscosity to the momentum and      //
// energy equations. Note that this term have to be scaled by the kernel     //
// gradient and the particle mass.                                           //
// Note: dot_r_v = \vec_{v}_ab \cdot \vec{r}_ab. r_ab is not the unit vector //
// ------------------------------------------------------------------------- //
pub fn mon97_art_vis(
    r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, _h_i: f64, _h_j: f64, rho_mean: f64
) -> (f64, f64) {
    if dot_r_v <= 0. {
        //particles approaching
        // Parameters
        let alpha: f64 = 1.0;
        let beta: f64 = 2.0;

        let v_sig:f64 = 0.5*alpha*(cs_i + cs_j - beta*dot_r_v/r_ij);
        let dvdt :f64 = -v_sig*dot_r_v/(r_ij*rho_mean);
        let dudt :f64 = 0.5*dvdt*dot_r_v;

        return (dvdt, dudt);
    } else {
        // particles receding
        return (0.0, 0.0);
    }
}

// ------------------------------------------------------------------------- //
// Artificial Viscosity proposed by Lodato & Price (2010) - Accretion Discs  //
// Returns:                                                                  //
//      [PI_ab, 0.5*PI_ab*(V_ab \cdot r_ab)] if V_ab \cdot r_ab < 0          //
//      [0.0, 0.0]                           if V_ab \cdot r_ab > 0          //
// i.e., the contribution from artificial viscosity to the momentum and      //
// energy equations. Note that this term have to be scaled by the kernel     //
// gradient and the particle mass.                                           //
// Note: dot_r_v = \vec_{v}_ab \cdot \vec{r}_ab. r_ab is not the unit vector //
// ------------------------------------------------------------------------- //
pub fn lodatoprice10_art_vis(
    r_ij: f64, dot_r_v: f64, cs_i: f64, cs_j: f64, h_i: f64, h_j: f64, rho_mean: f64
) -> (f64, f64) {
    let alpha: f64 = 1.0;
    let h_mean: f64 = 0.5*(h_i+h_j);
    let v_sig: f64;
    if dot_r_v <= 0. {
        // particles approaching
        // Parameters
        let beta: f64 = 2.0;

        v_sig = 0.5*alpha*(cs_i + cs_j - beta*dot_r_v/r_ij);
    } else {
        // particles receding
        // Parameters
        v_sig = 0.5*alpha*(cs_i + cs_j);
    }
    let dvdt: f64 = -h_mean*v_sig*dot_r_v/(r_ij*rho_mean*r_ij);
    let dudt: f64 = 0.5*dvdt*dot_r_v;
    return (dvdt, dudt);
}


// !!!-------------------- Artificial conductivity ----------------------!!! //

// ------------------------------------------------------------------------- //
// Artificial conductivity proposed by Price (2008)                          //
// Returns:                                                                  //
//      du/dt                                                                //
// i.e., the contribution from thermal conductivity to the energy equation.  //
// Note that this term have to be scaled by the kernel gradient and the      //
// particle mass.                                                            //
// ------------------------------------------------------------------------- //
pub fn price08_therm_cond(
    p_i: f64, p_j: f64, rho_mean: f64, u_i: f64, u_j: f64
) -> f64 {
    // Parameters
    let alpha_u: f64 = 1.;

    let v_sig_u: f64 = ((p_i - p_j).abs()/rho_mean).sqrt();
    let dudt: f64    = alpha_u*v_sig_u*(u_i-u_j)/rho_mean;

    return dudt;
}


// !!!-------------------- Fluid Dynamics Equations ---------------------!!! //

// ------------------------------------------------------------------------- //
// The momentum equation for particle_a due to particle_b                    //
// Returns:                                                                  //
//      (acc_x, acc_y, acc_z) per unit of mass                               //
// i.e., the particle_a's acceleration per unit of mass.                     //
// It includes the pressure gradient and AV terms                            //
// ------------------------------------------------------------------------- //
pub fn acceleration_ab(
    particle_a: &Particle, particle_b: &Particle, x_rel: f64, y_rel: f64, z_rel: f64,
    p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64, art_visc: f64
) -> (f64, f64, f64) {
    let acc: f64 = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb + 0.5*art_visc*(grad_ha+grad_hb);
    (-acc*x_rel, -acc*y_rel, -acc*z_rel)
}

// ***------------------------- External Forces -------------------------*** //

// ------------------------------------------------------------------------- //
// No external forces are applied.                                           //
// Nothing is returned                                                       //
// ------------------------------------------------------------------------- //
pub fn body_forces_null(
    _particles: &mut Particle, _m_star: & Star
) {
}

// ------------------------------------------------------------------------- //
// Toy Star model, Monaghan & Price (2004):                                  //
// Computes external forces:                                                 //
//      a_i = -nu * v_i - lambda * x_i                                       //
// where nu is a damping coefficient for the equilibrium state, and the      //
// lambda term is the simpliﬁed gravitational term.                          //
// For this case, we use some features of star to save the nu and lmbda      //
// parameters. star.hacc = nu; star.facc = lmbda.                            //
// This is a temporal implementation. To be fixed in future versions.        //
// ------------------------------------------------------------------------- //
pub fn body_forces_toy_star(
    particle: &mut Particle, star: &Star
) { 
    particle.ax -= star.hacc * particle.vx + star.facc*particle.x;
    particle.ay -= star.hacc * particle.vy + star.facc*particle.y; 
    particle.az -= star.hacc * particle.vz + star.facc*particle.z; 
}

// ------------------------------------------------------------------------- //
// Newtonian gravitational force due to a sink particle (star):              //
// Computes external forces:                                                 //
//      a_i = -G*M/(r^2 + e^2) \hat{r}                                       //
// where G is the gravitational constant, M is the star mass, r the distance //
// to the star and e is the softening term.                                  //
// ------------------------------------------------------------------------- //
pub fn body_forces_gravitation(
    particle: &mut Particle, m_star: & Star
) {
    let x_r: f64 = particle.x - m_star.x;
    let y_r: f64 = particle.y - m_star.y;
    let z_r: f64 = particle.z - m_star.z;
    let f_grav: f64 = m_star.m * (x_r*x_r + y_r*y_r + z_r*z_r+0.000625).powf(-1.5);
    particle.ax += -f_grav*x_r;
    particle.ay += -f_grav*y_r;
    particle.az += -f_grav*z_r;
}


// ***-------------------------- HD Equations ---------------------------*** //

// ------------------------------------------------------------------------- //
// Solve the SPH hydrodynamical equations for a given time step.             //
// Updates the state of particles:
//      acceleration: ax, ay, az                                             //
//      div(v)                                                               //
//      Delta u (du)                                                         //
// ------------------------------------------------------------------------- //
pub fn accelerations(
    particles: &mut Vec<Particle>, dm:f64, eos_type: bool, eos: fn(f64, f64, f64, f64, f64, f64, f64)->f64, cs: fn(f64, f64, f64, f64, f64, f64, f64)->f64, gamma: f64, coeff: f64,
    dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
    tree: &Node, s_: i32, n: usize, ptr : Pointer, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64,
    artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
    body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool, xperiodic: bool, yperiodic:bool, zperiodic:bool
) {
    // Find every neighbour of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[ii].h*rkern, xperiodic, yperiodic, zperiodic);
        return neighbors;
    }).collect();
    // Update particles state
    (0..n).into_par_iter().for_each(move |ii| {
        if particles[ii].ptype==0 {

            // Pointer to iith-particle
            let particle_i = unsafe { &mut *{ptr}.0.add(ii)};

            // Initialize variables to zero
            particle_i.ax   = 0.;
            particle_i.ay   = 0.;
            particle_i.az   = 0.;
            particle_i.divv = 0.;
            particle_i.du   = 0.;

            let p_i: f64    = eos(particles[ii].rho, particles[ii].u, gamma, particles[ii].x, particles[ii].y, particles[ii].z, coeff);
            let cs_i: f64   = cs(particles[ii].rho, particles[ii].u, gamma, particles[ii].x, particles[ii].y, particles[ii].z, coeff);
            let omeg_i: f64 = omega(particles, ii, &neighbors[ii], dm, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);
            
            for jj in 0..n {
                if ii != jj && particles[jj].ptype != 2 {
                    let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
                    let r_ij: f64 = (x_rel*x_rel + y_rel*y_rel+ z_rel*z_rel).sqrt();

                    let mut grad_hi: f64 = 0.0;
                    let mut grad_hj: f64 = 0.0;
                    if r_ij <= rkern*particles[ii].h {
                        let hisq: f64 = particles[ii].h*particles[ii].h;
                        grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*hisq*hisq);
                    }
                    if r_ij <= rkern*particles[jj].h {
                        let hjsq: f64 = particles[jj].h*particles[jj].h;
                        grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*hjsq*hjsq);
                    }
                    if grad_hi != 0. || grad_hj != 0.0 {
                        let p_j: f64    = eos(particles[jj].rho, particles[jj].u, gamma, particles[jj].x, particles[jj].y, particles[jj].z, coeff);
                        let cs_j: f64   = cs(particles[jj].rho, particles[jj].u, gamma, particles[jj].x, particles[jj].y, particles[jj].z, coeff);
                        let omeg_j: f64 = omega(particles, jj, &neighbors[jj], dm, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, rkern, wd, lg, hg, xperiodic, yperiodic, zperiodic);

                        // Velocity dot position
                        let dot_r_v = (particles[ii].vx-particles[jj].vx)*x_rel
                                          +(particles[ii].vy-particles[jj].vy)*y_rel
                                          +(particles[ii].vz-particles[jj].vz)*z_rel;

                        // Mean density
                        let rho_mean: f64 = 0.5 * (particles[ii].rho + particles[jj].rho);
                        
                        // Artificial viscosity
                        let (art_visc_mom, art_visc_ene) = artificial_viscosity(r_ij, dot_r_v, cs_i, cs_j, rho_mean, particles[ii].rho, particles[jj].rho);

                        // Acceleration
                        let (f_ij_x, f_ij_y, f_ij_z) = acceleration_ab(&particles[ii], &particles[jj], x_rel, y_rel, z_rel, p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc_mom);
                        particle_i.ax += dm * f_ij_x;
                        particle_i.ay += dm * f_ij_y;
                        particle_i.az += dm * f_ij_z;
                        
                        // Divergence of v per unit of mass
                        let div_vel :f64 = grad_hi*dot_r_v / (omeg_i*particles[ii].rho);
                        particle_i.divv -= dm*div_vel;
                        
                        // Internal energy change
                        if eos_type {
                            // Artificial thermal conductivity
                            let art_therm_cond: f64 = price08_therm_cond(p_i, p_j, rho_mean, particles[ii].u, particles[jj].u);
                            particle_i.du += dm * (0.5*(art_visc_ene + art_therm_cond*r_ij)*(grad_hi/omeg_i+grad_hj/omeg_j));
                        }
                    }
                }
            }
            if eos_type {
                particle_i.du -= (p_i/particles[ii].rho)*particle_i.divv;
            }
            // Body forces
            if bf {
                body_forces(particle_i, star);
            }
        }
    });
}


// ***------------------------- Sink Evolution --------------------------*** //

// ------------------------------------------------------------------------- //
// The Kick-Drift-Kick (KDK) integrator (second order)                       //
// Updates the star's state by                                               //
//      velocities (half a step),                                            //
//      positions  (full step),                                              //
//      velocities (half a step)                                             //
// ------------------------------------------------------------------------- //
pub fn star_integrator(
    star: &mut Star, dt: f64
) {
    star.vx += 0.5 * dt * star.ax;
    star.vy += 0.5 * dt * star.ay;
    star.vz += 0.5 * dt * star.az;
    
    star.x  += dt * star.vx;
    star.y  += dt * star.vy;
    star.z  += dt * star.vz;
    
    star.vx += 0.5 * dt * star.ax;
    star.vy += 0.5 * dt * star.ay;
    star.vz += 0.5 * dt * star.az;
}


// !!!------------------------- Time Integrator -------------------------!!! //

// ------------------------------------------------------------------------- //
// Euler Integrator (first order)                                            //
// Updates system's state for one time step:                                 //
//      f(t + dt) = f(t) + dt * f'(t)                                        //
// ------------------------------------------------------------------------- //
pub fn euler_integrator(
    particles: &mut Vec<Particle>, dt:f64, dm:f64, eos_type: bool, eos: fn(f64, f64, f64, f64, f64, f64, f64)->f64, cs: fn(f64, f64, f64, f64, f64, f64, f64)->f64, gamma:f64, coeff: f64,
    dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
    eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
    artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
    body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
    boundary: fn(&mut Vec<Particle>, f64, f64,f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 30, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
    accelerations(particles, dm, eos_type, eos, cs, gamma, coeff, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.x  += dt * particle.vx;
            particle.y  += dt * particle.vy;
            particle.z  += dt * particle.vz;
            particle.vx += dt * particle.ax;
            particle.vy += dt * particle.ay;
            particle.vz += dt * particle.az;
            particle.u  += dt * particle.du;
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
}

// ------------------------------------------------------------------------- //
// Velocity Verlet integrator (second order)                                 //
// Updates system's state for one time step.                                 //
// Verlet (1967)                                                             //
// ------------------------------------------------------------------------- //
pub fn velocity_verlet_integrator(
    particles: &mut Vec<Particle>, dt:f64, dm:f64, eos_type: bool, eos: fn(f64, f64, f64, f64, f64, f64, f64)->f64, cs: fn(f64, f64, f64, f64, f64, f64, f64)->f64, gamma:f64, coeff: f64,
    dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
    eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
    artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
    body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
    boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64
) {    
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;

            particle.u  += 0.5 * dt * particle.du;

            particle.x  += dt * particle.vx;
            particle.y  += dt * particle.vy;
            particle.z  += dt * particle.vz;
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
    
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 30, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);

    accelerations(particles, dm, eos_type, eos, cs, gamma, coeff, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;
            particle.u  += 0.5 * dt * particle.du;
        }
    });
}

// ------------------------------------------------------------------------- //
// The Kick-Drift-Kick (KDK) integrator (second order)                       //
// Updates the system's state by                                             //
//      velocities (half a step),                                            //
//      positions  (full step),                                              //
//      velocities (half a step)                                             //
// ------------------------------------------------------------------------- //
pub fn predictor_kdk_integrator(
    particles: &mut Vec<Particle>, dt:f64, dm:f64, eos_type: bool, eos: fn(f64, f64, f64, f64, f64, f64, f64)->f64, cs: fn(f64, f64, f64, f64, f64, f64, f64)->f64, gamma:f64, coeff: f64,
    dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, rkern: f64,
    eta: f64, tree: &mut Node, s_: i32, alpha_: f64, beta_:f64, n: usize, ptr : Pointer,
    artificial_viscosity: fn(f64, f64, f64, f64, f64, f64, f64) -> (f64, f64),
    body_forces: fn(&mut Particle, &Star), star: &Star, bf: bool,
    boundary: fn(&mut Vec<Particle>, f64, f64, f64, f64, f64, f64), xperiodic: bool, yperiodic:bool, zperiodic:bool, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64
) {
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;
            
            particle.u  += 0.5 * dt * particle.du;
            
            particle.x  += dt * particle.vx;
            particle.y  += dt * particle.vy;
            particle.z  += dt * particle.vz;
            
            // Predictor 
            particle.vx_star = particle.vx;
            particle.vy_star = particle.vy;
            particle.vz_star = particle.vz;
            particle.u_star = particle.u;
            
            particle.vx += 0.5 * dt * particle.ax;
            particle.vy += 0.5 * dt * particle.ay;
            particle.vz += 0.5 * dt * particle.az;
            
            particle.u  += 0.5 * dt * particle.du;
        }
    });
    boundary(particles, wd, lg, hg, x0, y0, z0);
    tree.build_tree(s_, alpha_, beta_, particles, 1.0e-02);
    smoothing_length(particles, dm, eta, f, dfdq, sigma, rkern, 1e-03, 30, dt, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, xperiodic, yperiodic, zperiodic);
    accelerations(particles, dm, eos_type, eos, cs, gamma, coeff, dwdh_, f, dfdq, sigma, rkern, tree, s_, n, ptr, wd, lg, hg, x0, y0, z0, artificial_viscosity, body_forces, star, bf, xperiodic, yperiodic, zperiodic);
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            particle.vx = particle.vx_star + 0.5 * dt * particle.ax;
            particle.vy = particle.vy_star + 0.5 * dt * particle.ay;
            particle.vz = particle.vz_star + 0.5 * dt * particle.az;
            particle.u = particle.u_star + 0.5 * dt * particle.du;
        }
    });
}


// ***----------------------- Boundary conditions -----------------------*** //

// ------------------------------------------------------------------------- //
// None Boundary Condition                                                   //
// ------------------------------------------------------------------------- //
pub fn none_boundary(
    _particles: &mut Vec<Particle>, _wd: f64, _lg: f64, _hg: f64, _x0:f64, _y0: f64, _z0: f64
){
}

// ------------------------------------------------------------------------- //
// Periodic Boundary Conditions                                              //
// ------------------------------------------------------------------------- //
pub fn periodic_boundary(
    particles: &mut Vec<Particle>, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0: f64
){
    // We assume that the system's domain is a rectangular box.
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

// ------------------------------------------------------------------------- //
// Open Boundary Conditions                                              //
// ------------------------------------------------------------------------- //
pub fn open_boundary(
    particles: & Vec<Particle>, wd: &mut f64, lg: &mut f64, hg: &mut f64, x0: &mut f64, y0: &mut f64, z0: &mut f64
){
    let xmin: f64 = particles.iter().fold(*x0, |a, part| a.min(part.x));
    let ymin: f64 = particles.iter().fold(*y0, |a, part| a.min(part.y));
    let zmin: f64 = particles.iter().fold(*z0, |a, part| a.min(part.z));
    let xmax: f64 = particles.iter().fold(*x0+*wd, |a, part| a.max(part.x));
    let ymax: f64 = particles.iter().fold(*y0+*lg, |a, part| a.max(part.y));
    let zmax: f64 = particles.iter().fold(*z0+*hg, |a, part| a.max(part.z));

    if xmin < *x0 { *x0 = xmin - 0.1*xmin.abs();}
    if ymin < *y0 { *y0 = ymin - 0.1*ymin.abs();}
    if zmin < *z0 { *z0 = zmin - 0.1*zmin.abs();}
    if xmax > *x0+*wd { *wd = (xmax-*x0) + 0.1*(xmax-*x0).abs();}
    if ymax > *y0+*lg { *lg = (ymax-*y0) + 0.1*(ymax-*y0).abs();}
    if zmax > *z0+*hg { *hg = (zmax-*z0) + 0.1*(zmax-*z0).abs();}
}

// ------------------------------------------------------------------------- //
// Non-Periodic Boundary Conditions                                              //
// ------------------------------------------------------------------------- //
pub fn box_boundary(
    particles: &mut Vec<Particle>, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0: f64
){
    // We assume that the system's domain is a rectangular box.
    particles.par_iter_mut().for_each(|particle|{
        if particle.ptype==0 {
            if (particle.x >= (wd+x0)) || (particle.x < x0) || (particle.y >= (lg + y0)) || (particle.y < y0)
                || (particle.z >= (hg + z0)) || (particle.y < z0){
                particle.ptype = 2;
            }
        }
    });
}

// ------------------------------------------------------------------------- //
// Accretion onto sink particles                                             //
// ------------------------------------------------------------------------- //
pub fn accretion_boundary(
    star: &mut Star, particles: &mut Vec<Particle>, dm:f64, n: &mut usize, tree: & Node, s_: i32, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0: f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
){
    for ii in 0..*n {
        if particles[ii].ptype == 2 {
            particles.remove(ii);
            *n -= 1;
        }
    }
    star.ax = 0.0;
    star.ay = 0.0;
    star.az = 0.0;

    let mut neighbors: Vec<usize> = Vec::new();
    tree.find_neighbors_star(&star, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, star.hacc, xperiodic, yperiodic, zperiodic);
    
    for ii in neighbors {
        let dx: f64 = star.x - particles[ii].x;
        let dy: f64 = star.y - particles[ii].y;
        let dz: f64 = star.z - particles[ii].z;
        let r2: f64 = dx*dx + dy*dy + dz*dz;
        let dvx: f64= star.vx - particles[ii].vx;
        let dvy: f64= star.vy - particles[ii].vy;
        let dvz: f64= star.vz - particles[ii].vz;
        let v2: f64 = dvx*dvx + dvy*dvy + dvz*dvz;
        let rds: f64= star.facc*star.hacc;
        let mut acc: bool = false;
        if r2 < rds*rds {
            acc = true;
        } else if r2 <  star.hacc*star.hacc {
            // 1. |Lai| < |Lacc|
            let rdotv: f64  = dx*dvx + dy*dvy + dz*dvz;
            let dl2: f64    = r2*v2 - rdotv*rdotv;
            let dlacc2: f64 = star.hacc*star.m;
            if dl2 < dlacc2 {
                // 2. e < 0
                let e: f64  = 0.5 * v2 - star.m/r2.sqrt();
                if e < 0.0 {
                    acc = true;
                }
            }
        }
        if acc {
            star.m += dm;
            let inv_mtot: f64 = 1.0/star.m;
            star.x = (particles[ii].x*dm + star.x*star.m)*inv_mtot;
            star.y = (particles[ii].y*dm + star.y*star.m)*inv_mtot;
            star.z = (particles[ii].z*dm + star.z*star.m)*inv_mtot;
            star.vx = (particles[ii].vx*dm + star.vx*star.m)*inv_mtot;
            star.vy = (particles[ii].vy*dm + star.vy*star.m)*inv_mtot;
            star.vz = (particles[ii].vz*dm + star.vz*star.m)*inv_mtot;
            star.ax = (particles[ii].ax*dm + star.ax*star.m)*inv_mtot;
            star.ay = (particles[ii].ay*dm + star.ay*star.m)*inv_mtot;
            star.az = (particles[ii].az*dm + star.az*star.m)*inv_mtot;
            particles[ii].ptype = 2;
            particles.remove(ii);
            *n -= 1;
        }
    }
}


// !!!---------------------- Timestepping Criteria ----------------------!!! //

// ------------------------------------------------------------------------- //
// The Courant-Friedrichs-Lewy (CFL) condition.                              //
// Returns the maximum stable timestep,                                      //
//      dt_a = 0.3 * h_a / v_sig_a.                                          //
//  Bate at al. (1995).                                                      //
// ------------------------------------------------------------------------- //
pub fn cfl_dt(
    h: f64, cs: f64, div_v:f64, alpha:f64, beta: f64
) -> f64{
    if div_v < 0. {
        return 0.3*h / (cs + h*div_v.abs() + 1.2*(alpha*cs + beta*h*div_v.abs()));
    } else {
        return 0.3*h / (cs + h*div_v.abs());
    }
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// The force conditon.                                                       //
// Returns the maximum stable timestep from accelerations,                   //
//      dt_a = f * sqrt[h_a / |acc_a|].                                      //
//  Monaghan (1989).                                                         //
// ------------------------------------------------------------------------- //
pub fn force_dt(
    h: f64, a: f64, f: f64
) -> f64 {
    f*(h/a).sqrt()
}

// ------------------------------------------------------------------------- //
// Timestepping Criteria by Cossins P. J. (2010)                             //
// Returns the minimum time step between the CFL and the force conditions.   //
// ------------------------------------------------------------------------- //
pub fn time_step_bale(
    particles: & Vec<Particle>, n: usize, gamma: f64, coeff: f64, _rkern: f64, _wd: f64, _lg: f64, _hg: f64,
    _tree: &mut Node, _s_: i32, cs: fn(f64, f64, f64, f64, f64, f64, f64) -> f64
) -> f64 {
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        if particles[ii].ptype == 0 {
        let a: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let cs: f64 = cs(particles[ii].rho, particles[ii].u, gamma, particles[ii].x, particles[ii].y, particles[ii].z, coeff);
        let dt_a: f64 = force_dt(particles[ii].h, a, 0.3);
        let dt_cfl: f64 = cfl_dt(particles[ii].h, cs, particles[ii].divv, 1., 2.);
        return (dt_a).min(dt_cfl);}
        else {
            return f64::INFINITY;
        }
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}

// ------------------------------------------------------------------------- //
// Timestepping Criteria by Monaghan (1997)                                  //
// Returns the minimum time step between the CFL and the force conditions.   //
// ------------------------------------------------------------------------- //
pub fn time_step_mon(
    particles: & Vec<Particle>, n: usize, gamma: f64, coeff: f64, rkern: f64, wd: f64, lg: f64, hg: f64, x0: f64, y0: f64, z0: f64,
    tree: &mut Node, s_: i32, cs: fn(f64, f64, f64, f64, f64, f64, f64) -> f64, xperiodic: bool, yperiodic:bool, zperiodic:bool
) -> f64 {
    // Find every neighbor of every particle.
    let neighbors: Vec<Vec<usize>> = (0..n).into_par_iter().map(|ii: usize| {
        let mut neighbors: Vec<usize> = Vec::new();
        if particles[ii].ptype == 0 {
            tree.find_neighbors(ii, s_, particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[ii].h*rkern, xperiodic, yperiodic, zperiodic);
        }
        return neighbors;
    }).collect();
    let dts :Vec<f64> = (0..n).into_par_iter().map(|ii| -> f64 {
        if particles[ii].ptype == 0 {
        let alpha: f64  = 1.;
        let beta: f64   = 2.;
        let mut v_sig:f64 = 0.0;
        let cs_i: f64 = cs(particles[ii].rho, particles[ii].u, gamma, particles[ii].x, particles[ii].y, particles[ii].z, coeff);
        for jj in &neighbors[ii] {
            let cs_j: f64 = cs(particles[*jj].rho, particles[*jj].u, gamma, particles[*jj].x, particles[*jj].y, particles[*jj].z, coeff);

            // Velocity dot position
            let (x_rel, y_rel, z_rel) = periodic_rel_vector(&particles[ii], &particles[*jj], wd, lg, hg, rkern*particles[ii].h, xperiodic, yperiodic, zperiodic);
            let r_ij: f64 = (x_rel*x_rel + y_rel*y_rel + z_rel*z_rel).sqrt();
            let dot_r_v: f64 = (particles[ii].vx-particles[*jj].vx)*x_rel
                              +(particles[ii].vy-particles[*jj].vy)*y_rel
                              +(particles[ii].vz-particles[*jj].vz)*z_rel;
            
            if dot_r_v < 0. {
                let v_sig_ij: f64 = alpha*(cs_i+cs_j - beta*(dot_r_v/r_ij));
                if v_sig_ij > v_sig {
                    v_sig = v_sig_ij;
                }
            }
        }
        let a_norm: f64 = (particles[ii].ax*particles[ii].ax + particles[ii].ay*particles[ii].ay + particles[ii].az*particles[ii].az).sqrt();
        let dt_a: f64   = force_dt(particles[ii].h, a_norm, 0.25);
        let dt_cfl: f64 = 0.3*particles[ii].h / v_sig;
        return (dt_a).min(dt_cfl);
        } else {
            return f64::INFINITY;
        }
    }).collect();
    dts.iter().fold(f64::INFINITY, |a, &b| a.min(b))
}