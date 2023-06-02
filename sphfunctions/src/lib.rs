use std::{
    error::Error,
};

use rand::{thread_rng, Rng};

use csv::Writer;

use std::f64::consts::PI;

pub struct Particle {
	pub m: f64,
	pub x: f64,
    pub y: f64,
    //pub z: f64,
    pub h: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    //pub vz: f64,
    pub ax: f64,
    pub ay: f64,
    //pub az: f64,
    pub u: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            m: 1.,
            x: 0.,
            y: 0.,
            //z: 0.,
            h: 0.1,
            rho: 1.0,
            vx: 0.0,
            vy: 0.0,
            //vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            //az: 0.0,
            u: 0.0,
            du: 0.0,
        }
    }
}


// Write data
pub fn init_square(path: &str, n: u32, m:f64, rho:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let h:f64 = 1./n as f64;
    let offset = 0.5 / n as f64;
    wtr.write_record(&["m", "x", "y", "h", "rho"])?;
    for ii in 0..n{
        for jj in 0..n{
            wtr.write_record(&[m.to_string(), (ii as f64 / n as f64 + offset).to_string(), (jj as f64 / n as f64 + offset).to_string(), h.to_string(), rho.to_string()])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_random_circle(path: &str, n:u32, r:f64, m:f64, rho:f64, h:f64, x0:f64, y0:f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    //let mut g = rand::seed_from_f64(seed);

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
pub fn sust_vec(vec1: & Vec<f64>, vec2: & Vec<f64>) -> Vec<f64> {
    vec1.into_iter().zip(vec2).map(|(a, b)| a - b).collect()
}

pub fn rel_distance(p1: &Particle, p2: &Particle) -> Vec<f64> {
    vec![p1.x - p2.x, p1.y - p2.y]
}

pub fn euclidean_norm(p1: &Particle, p2: &Particle) -> f64 {
    let vec = rel_distance(p1, p2);
    let mut sum :f64 = 0.0;
    for x in vec.iter() {
        sum += x*x;
    }
    sum.sqrt()
}


// Kernel
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
    //println!("{} {}", f_h, h);
    let omeg = omega(particle_a, neigh_particles, h, rho_kernel, dwdh, f, dfdq, sigma, d);
    let df = -(d
         as f64)*rho_h*omeg/ h;
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

pub fn smoothing_length(particles: &mut Vec<Particle>, eta:f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, sigma:f64, d:i32, tol: f64, it: u32){
    for ii in 0..particles.len(){
        // Look for neighboring particles
        particles[ii].h = newton_raphson(&particles[ii], particles, particles[ii].h, eta, f, dfdq, sigma, d, tol, it);
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


pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64) -> Vec<f64> {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha + p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb;
    vec![-acc*(particle_a.x - particle_b.x), -acc*(particle_a.y - particle_b.y)]
}

pub fn body_forces_toy_star(x: f64, y: f64, vx: f64, vy: f64, nu: f64, lmbda: f64) -> Vec<f64> {
    vec![-nu * vx - lmbda*x, -nu * vy - lmbda*y]  
}

pub fn accelerations(particles: &mut Vec<Particle>, eos: fn(f64, f64, f64)->f64, k:f64, gamma:f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, f: fn(f64) -> f64, dfdq: fn(f64) -> f64, nu:f64, lmbda:f64, sigma: f64, d:i32){
    let n = particles.len();
    for ii in 0..(n) {
        let p_i = eos(particles[ii].rho, k, gamma);
        let omeg_i = omega(&particles[ii], particles, particles[ii].h, particles[ii].rho, dwdh_, f, dfdq, sigma, d);
        let mut dudt = 0.0;
        //println!("{} {} {}", ii, p_i, particles[ii].h);
        for jj in (ii+1)..n {
            let p_j = eos(particles[jj].rho, k, gamma);
            let omeg_j = omega(&particles[jj], particles, particles[jj].h, particles[jj].rho, dwdh_, f, dfdq, sigma, d);
            let r_ij = euclidean_norm(&particles[ii], &particles[jj]);
            let grad_hi = dfdq(r_ij/particles[ii].h)*sigma/(r_ij*(particles[ii].h).powi(d+1));
            let grad_hj = dfdq(r_ij/particles[jj].h)*sigma/(r_ij*(particles[jj].h).powi(d+1));
            //println!("{}   {}   {}   {}   {}   {}   {}", ii, jj, p_j, omeg_j, r_ij, grad_hi, grad_hj);
            // Acceleration
            let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj);
            //let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, 1., 1., grad_hi, grad_hi);
            // Thermal change
            let dot_r_v = (particles[ii].vx-particles[jj].vx)*(particles[ii].x-particles[jj].x)
                         +(particles[ii].vy-particles[jj].vy)*(particles[ii].y-particles[jj].y);
            dudt += particles[jj].m*grad_hi*dot_r_v;

            particles[ii].ax += particles[jj].m *f_ij[0];
            particles[ii].ay += particles[jj].m *f_ij[1];
            particles[jj].ax -= particles[ii].m *f_ij[0];
            particles[jj].ay -= particles[ii].m *f_ij[1];
        }
        let body_forces = body_forces_toy_star(particles[ii].x, particles[ii].y, particles[ii].vx, particles[ii].vy, nu, lmbda);
        particles[ii].ax += body_forces[0];
        particles[ii].ay += body_forces[1];
        particles[ii].du = p_i / (omeg_i*particles[ii].h*particles[ii].h) *dudt;
    }
}

pub fn euler_integrator(particle: &mut Particle, dt: f64) {
    particle.x += dt * particle.vx;
    particle.y += dt * particle.vy;
    particle.vx += dt * particle.ax;
    particle.vy += dt * particle.ay;
    particle.u += dt * particle.du;
}