use std::{
    error::Error,
    f64,
};

use rand::{thread_rng, Rng};

use csv::Writer;

use std::f64::consts::PI;

pub struct Particle {
	pub m: f64,
	pub x: f64,
    pub y: f64,
    pub h: f64,
    pub dh: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    pub ax: f64,
    pub ay: f64,
    pub u: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            m: 1.,
            x: 0.,
            y: 0.,
            h: 0.1,
            dh: 0.0,
            rho: 1.0,
            vx: 0.0,
            vy: 0.0,
            ax: 0.0,
            ay: 0.0,
            u: 0.0,
            du: 0.0,
        }
    }
}


// Write data
pub fn init_square(path: &str, n: u32, m:f64, rho:f64, h:f64, w:f64, l:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let dx = (w*l / n as f64).sqrt();
    let nx :i64 = (w/dx) as i64;
    let ny :i64 = (l/dx) as i64;
    wtr.write_record(&["m", "x", "y", "h", "rho"])?;
    for jj in 0..ny{
        for ii in 0..nx{
            wtr.write_record(&[m.to_string(), (dx*ii as f64).to_string(), (dx*jj as f64).to_string(), h.to_string(), rho.to_string()])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_square(path: &str, n: u32, m:f64, rho:f64, h:f64, w:f64, l:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = thread_rng();
    wtr.write_record(&["m", "x", "y", "h", "rho"])?;
    for _ii in 0..n{
        let x = rng.gen_range(0.0f64, w);
        let y = rng.gen_range(0.0f64, l);
        wtr.write_record(&[m.to_string(), x.to_string(), y.to_string(), h.to_string(), rho.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_circle(path: &str, n: u32, r:f64, m:f64, rho:f64, h:f64, x0:f64, y0:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = thread_rng();
    wtr.write_record(&["m", "x", "y", "h", "rho"])?;
    for _ii in 0..n{
        let r_i = r*(rng.gen_range(0.0f64, 1.0f64)).sqrt();
        let theta_i = 2.0*PI*rng.gen_range(0.0f64, 1.0f64);
        let x = r_i*theta_i.cos() + x0;
        let y = r_i*theta_i.sin() + y0;
        wtr.write_record(&[m.to_string(), x.to_string(), y.to_string(), h.to_string(), rho.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["m", "x", "y", "vx", "vy", "u", "h", "rho"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].m.to_string(), particles[ii].x.to_string(), particles[ii].y.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(),
                           particles[ii].u.to_string(), particles[ii].h.to_string(), particles[ii].rho.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

// Read data
pub fn read_data(path: &str, particles: &mut Vec<Particle>) -> Result<(), Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;
    for result in rdr.records() {
        let record = result?;
        particles.push(Particle{m:(&record[0]).parse::<f64>().unwrap(), x:(&record[1]).parse::<f64>().unwrap(), y:(&record[2]).parse::<f64>().unwrap(),
                                h:(&record[3]).parse::<f64>().unwrap(), rho:(&record[4]).parse::<f64>().unwrap(), ..Default::default()});
    }
    Ok(())
}


// Basic vector functions
pub fn euclidean_norm(p1: &Particle, p2: &Particle) -> f64 {
    let sum :f64 = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    sum.sqrt()
}

// Kernel

// Cubic Kernel
pub fn f_cubic_kernel(q:f64) -> f64 {
    let mut f:f64 = 0.;
    if q < 1. {
        f = 1. - 1.5*q*q + 0.75*q*q*q;
    } else if q < 2.{
        f = 0.25*(2.-q).powi(3);
    }
    f
}

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

pub fn dfdq_gaussian_kernel(q:f64) -> f64 {
    let mut f:f64 = 0.;
    if q < 2. {
        f = -2.0*q*(-q*q).exp();
    }
    f
}

pub fn dwdh(q: f64, f: fn(f64) -> f64, df: fn(f64) -> f64, d:i32) -> f64 {
    (d as f64) *f(q) + q*df(q)
}

pub fn density_kernel(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64, sigma:f64, d:i32, f: fn(f64)->f64) -> f64 {
    let mut rho :f64 = 0.0;
    for ii in 0..neigh_particles.len(){
        let r = euclidean_norm(&particle_a, &neigh_particles[ii]);
        rho += f(r/h);
    }
    rho * &particle_a.m * sigma / h.powi(d)
}

pub fn density_by_smoothing_length(m:f64, h:f64, eta:f64, d:i32) -> f64{
    let vol = eta/h;
    m*vol.powi(d)
}

// Iterations
pub fn omega(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, d:i32) -> f64{
    let n = neigh_particles.len();
    let mut omeg :f64 = 0.0;
    for ii in 0..n {
        let q = euclidean_norm(&particle_a, &neigh_particles[ii])/h;
        omeg -= dwdh_(q, f, dfdq, d);
    }
    omeg *= &particle_a.m*sigma/(h.powi(d)*rho*(d as f64));
    omeg + 1.
}

pub fn f_iter(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32) -> (f64 , f64) {
    let rho_kernel = density_kernel(particle_a, neigh_particles, h, sigma, d, f);
    let rho_h = density_by_smoothing_length(particle_a.m, h, eta, d);
    let f_h = rho_h - rho_kernel;
    let omeg = omega(particle_a, neigh_particles, h, rho_kernel, dwdh, f, dfdq, sigma, d);
    //println!("{} = {} - {}    {} {}", f_h, rho_h, rho_kernel, h, omeg);
    let df = -(d as f64)*rho_h*omeg/ h;
    (f_h, df)
}

fn nr_iter(particle_a: &Particle, neigh_particles: & Vec<Particle>, h_old: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32) -> f64 {
    let (f, df) = f_iter(particle_a, neigh_particles, h_old, eta, f, dfdq, sigma, d);
    (h_old - f / df).abs()
}

pub fn newton_raphson(particle_a: &Particle, particles: & Vec<Particle>, h_guess: f64, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32) -> f64 {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    while i <= it {
        // Look for neighboring particles
        h_new = nr_iter(particle_a, particles, h_old, eta, f, dfdq, sigma, d);
        if (h_new - h_old).abs() <=  tol {
            i = it + 2;
        } else{
            i += 1;
            h_old = h_new;
        }
    }
    if i == it+1 {
        0.0
    } else{
        h_new
    }
}

pub fn smoothing_length(particles: &mut Vec<Particle>, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32, dt:f64){
    for ii in 0..particles.len(){
        // Look for neighboring particles
        let h_new = newton_raphson(&particles[ii], particles, particles[ii].h+dt*particles[ii].dh, eta, f, dfdq, sigma, d, tol, it);
        if h_new != 0.0 {
            particles[ii].h = h_new;
        }
        particles[ii].rho = density_kernel(&particles[ii], &particles, particles[ii].h, sigma, d, f);
    }
}

// Equation of state

//  --- Toy Star 2D ---
pub fn eos_polytropic(rho:f64, k:f64, gamma:f64) -> f64 {
    k * rho.powf(1.+1./gamma)
}

pub fn coeff_static_grav_potential(k:f64, gamma:f64, m:f64, r:f64) -> f64 {
    2.0*k/(PI.powf(1./gamma)) * (m*(1.+gamma)/(r*r)).powf(1.+1./gamma)/m
}

// --- Ideal Gases ---
pub fn eos_ideal_gas(rho:f64, k:f64, gamma:f64) -> f64 {
    k*rho.powf(gamma)
}

pub fn thermal_energy(rho:f64, p:f64, gamma:f64) -> f64 {
    p/((gamma-1.)*rho)
}


// Dynamic Equations

// --- Force due to the gradient of pressure ---
pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64, art_visc: f64) -> Vec<f64> {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb + art_visc;
    vec![-acc*(particle_a.x - particle_b.x), -acc*(particle_a.y - particle_b.y)]
}

// --- Body Forces ---
pub fn body_forces_toy_star(particle: &mut Particle, nu: f64, lmbda: f64) {
    particle.ax -= nu * particle.vx + lmbda*particle.x;
    particle.ay -= nu * particle.vy + lmbda*particle.y; 
}

pub fn accelerations(particles: &mut Vec<Particle>, eos: fn(f64, f64, f64)->f64, k:f64, gamma:f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma: f64, d:i32){
    let n = particles.len();
    for ii in 0..n {
        let p_i = eos(particles[ii].rho, k, gamma);
        let omeg_i = omega(&particles[ii], particles, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, d);
        for jj in (ii+1)..n {
            let p_j = eos(particles[jj].rho, k, gamma);
            let omeg_j = omega(&particles[jj], particles, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, d);
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
            // Acceleration
            let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj, art_visc);
            particles[ii].dh += grad_hi*dot_r_v;

            particles[ii].ax += particles[jj].m *f_ij[0];
            particles[ii].ay += particles[jj].m *f_ij[1];
            particles[jj].ax -= particles[ii].m *f_ij[0];
            particles[jj].ay -= particles[ii].m *f_ij[1];
        }

        // Thermal change
        particles[ii].du = particles[ii].m*p_i / (omeg_i*particles[ii].rho*particles[ii].rho) * particles[ii].dh;
        // Smoothing length change
        particles[ii].dh *= -particles[ii].m * particles[ii].h/ (omeg_i*particles[ii].rho*d as f64);
    }
}


// Time integrator
pub fn euler_integrator(particle: &mut Particle, dt: f64) {
    particle.x += dt * particle.vx;
    particle.y += dt * particle.vy;
    particle.vx += dt * particle.ax;
    particle.vy += dt * particle.ay;
    particle.u += dt * particle.du;
}

// Boundary conditions

// --- Periodic Boundary Conditions ---
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