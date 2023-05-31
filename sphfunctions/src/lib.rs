use std::{
    error::Error,
};

//use rand::Rng;

use csv::Writer;

//use std::f64::consts::PI;

pub struct Particle {
	pub m: f64,
	pub x: f64,
    pub y: f64,
    pub z: f64,
    pub h: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub ax: f64,
    pub ay: f64,
    pub az: f64,
    pub u: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            m: 1.,
            x: 0.,
            y: 0.,
            z: 0.,
            h: 0.1,
            rho: 1.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
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
    wtr.write_record(&["m", "x", "y", "z", "h", "rho"])?;
    for ii in 0..n{
        for jj in 0..n{
            wtr.write_record(&[m.to_string(), (ii as f64 / n as f64 + offset).to_string(), (jj as f64 / n as f64 + offset).to_string(), (jj as f64 / n as f64 + offset).to_string(), h.to_string(), rho.to_string()])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

//pub fn init_random_circle(path: &str, n:u32, r:f64, x0:f64, y0:f64, z0:f64, seed:u32){
//    let mut wtr = Writer::from_path(path)?;
//    let mut g = rand::seed_from_f64(1234.0);
//    let rad = rng.gen_range(0.0 .. 1.0);
    //theta = 2*pi*rand(Float64, N)
    //x = r.*cos.(theta) .+ x0
    //y = r.*sin.(theta) .+ y0
//    wtr.write_record(&["m", "x", "y", "z", "h", "rho"])?;
//    for ii in 0..n{
//        let r_i = r*(rng.gen_range(0.0 .. 1.0)).sqrt();
//        let theta_i = 2.0*PI*rng.gen_range(0.0 .. 1.0);
        //let x = r_i*theta_i.cos() + x0;
        //let y = r_i*theta_i.sin() + y0;
        //let z = r_i*theta_i.cos() + z0;
//        wtr.write_record(&["1.0", x.to_string(), y.to_string(), z.to_string()])?;
//    }
//    wtr.flush()?;
//    Ok(())
//}

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["m", "x", "y", "z", "vx", "vy", "vz", "u", "h", "rho"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].m.to_string(), particles[ii].x.to_string(), particles[ii].y.to_string(), particles[ii].z.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(), particles[ii].vz.to_string(),
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
        particles.push(Particle{m:(&record[0]).parse::<f64>().unwrap(), x:(&record[1]).parse::<f64>().unwrap(), y:(&record[2]).parse::<f64>().unwrap(), z:(&record[3]).parse::<f64>().unwrap(),
                                h:(&record[4]).parse::<f64>().unwrap(), rho:(&record[5]).parse::<f64>().unwrap(), ..Default::default()});
    }
    Ok(())
}

// Basic vector functions
pub fn sust_vec(vec1: & Vec<f64>, vec2: & Vec<f64>) -> Vec<f64> {
    vec1.into_iter().zip(vec2).map(|(a, b)| a - b).collect()
}

pub fn rel_distance(p1: &Particle, p2: &Particle) -> Vec<f64> {
    vec![p1.x - p2.x, p1.y - p2.y, p1.z - p2.z]
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
pub fn cubic_kernel(q:f64) -> f64 {
    let mut w:f64 = 0.;
    if q < 1. {
        w = 1. - 1.5*q*q + 0.75*q*q*q;
    } else if q < 2.{
        w = 0.25*(2.-q).powi(3);
    }
    w
}

pub fn dwdq_cubic_kernel(q:f64) -> f64 {
    let mut dw:f64 = 0.;
    if q < 1. {
        dw = (2.25*q-3.)*q;
    } else if q < 2.{
        dw = -3.*(0.25*q*q-q+1.);
    }
    dw
}

pub fn dwdh(q: f64, w: fn(f64) -> f64, dw: fn(f64) -> f64, nu:i32) -> f64 {
    nu as f64 *w(q)+q*dw(q)
}

// Density
pub fn density_local(m:f64, pos:f64, h:f64, kernel: fn(f64) -> f64) -> f64{
    let q = pos/h;
    let rho = m * kernel(q);
    rho
}

pub fn density_kernel(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64) -> f64 {
    let mut rho :f64 = 0.0;
    for ii in 0..neigh_particles.len(){
        let r = euclidean_norm(&particle_a, &neigh_particles[ii]);
        rho += density_local((&neigh_particles[ii]).m, r, h, cubic_kernel);
    }
    rho
}

pub fn density_by_smoothing_length(m:f64, h:f64, eta:f64, nu:i32) -> f64{
    let vol = eta/h;
    m*vol.powi(nu)
}

// Iterations
pub fn omega(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64, rho: f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, nu:i32, sigma: f64) -> f64{
    let n = neigh_particles.len();
    let mut omeg :f64 = 0.0;
    for ii in 0..n{
        let q = euclidean_norm(&particle_a, &neigh_particles[ii])/h;
        let k = dwdh_(q, w, dwdq, nu);
        omeg -= (&neigh_particles[ii]).m * k;
    }
    omeg *= sigma/(h.powi(nu)*rho*nu as f64);
    omeg + 1.
}

pub fn f_iter(particle_a: &Particle, neigh_particles: & Vec<Particle>, h: f64, eta:f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, sigma:f64, nu:i32) -> (f64 , f64) {
    let rho_kernel = density_kernel(particle_a, neigh_particles, h);
    let rho_h = density_by_smoothing_length(particle_a.m, h, eta, nu);
    let f = rho_h - rho_kernel;
    let omeg = omega(particle_a, neigh_particles, h, rho_kernel, dwdh, w, dwdq, nu, sigma);
    let df = -(nu as f64)*rho_h*omeg/ h;
    (f, df)
}

fn nr_iter(particle_a: &Particle, neigh_particles: & Vec<Particle>, h_old: f64, eta:f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, sigma:f64, nu:i32) -> f64 {
    let (f, df) = f_iter(particle_a, neigh_particles, h_old, eta, w, dwdq, sigma, nu);
    (h_old - f / df).abs()
}

pub fn newton_raphson(particle_a: &Particle, particles: & Vec<Particle>, h_guess: f64, eta:f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, sigma:f64, nu:i32, tol: f64, it: u32) -> f64 {
    let mut h_new :f64 = 0.0;
    let mut h_old :f64 = h_guess;
    let mut i : u32 = 1;
    while i <= it {
        // Look for neighboring particles
        h_new = nr_iter(particle_a, particles, h_old, eta, w, dwdq, sigma, nu);
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

pub fn smoothing_length(particles: &mut Vec<Particle>, eta:f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, sigma:f64, nu:i32, tol: f64, it: u32){
    for ii in 0..particles.len(){
        // Look for neighboring particles
        particles[ii].h = newton_raphson(&particles[ii], particles, particles[ii].h, eta, w, dwdq, sigma, nu, tol, it);
    }
}

// Equation of state
pub fn eos_ideal_gas(rho:f64, k:f64, gamma:f64) -> f64 {
    k*rho.powf(gamma)
}

pub fn eos_polytropic(rho:f64, k:f64, gamma:f64) -> f64 {
    k * rho.powf(1.+1./gamma)
}

pub fn thermal_energy(rho:f64, p:f64, gamma:f64) -> f64 {
    p/((gamma-1.)*rho)
}

pub fn acceleration_ab(particle_a: &Particle, particle_b: &Particle, p_a: f64, p_b: f64, omeg_a: f64, omeg_b: f64, grad_ha: f64, grad_hb: f64) -> Vec<f64> {
    let acc = p_a/(omeg_a*particle_a.rho*particle_a.rho)*grad_ha - p_b/(omeg_b*particle_b.rho*particle_b.rho) * grad_hb;
    vec![-acc*(particle_a.x - particle_b.x), -acc*(particle_a.y - particle_b.y), -acc*(particle_a.z - particle_b.z)]
}

pub fn accelerations(particles: &mut Vec<Particle>, eos: fn(f64, f64, f64)->f64, k:f64, gamma:f64, dwdh_: fn(f64, fn(f64) -> f64, fn(f64) -> f64, i32) -> f64, w: fn(f64) -> f64, dwdq: fn(f64) -> f64, nu:i32, sigma: f64){
    let n = particles.len();
    for ii in 0..n {
        let p_i = eos(particles[ii].rho, k, gamma);
        let omeg_i = omega(&particles[ii], particles, particles[ii].h, particles[ii].rho, dwdh_, w, dwdq, nu, sigma);
        let mut dudt = 0.0;
        for jj in ii+1..n {
            let p_j = eos(particles[jj].rho, k, gamma);
            let omeg_j = omega(&particles[jj], particles, particles[jj].h, particles[jj].rho, dwdh_, w, dwdq, nu, sigma);
            let r_ij = euclidean_norm(&particles[ii], &particles[jj]);
            let grad_hi = dwdq(r_ij/particles[ii].h)/(r_ij*particles[ii].h);
            let grad_hj = dwdq(r_ij/particles[jj].h)/(r_ij*particles[jj].h);

            // Acceleration
            let f_ij = acceleration_ab(&particles[ii], &particles[jj], p_i, p_j, omeg_i, omeg_j, grad_hi, grad_hj);

            // Thermal change
            let dot_r_v = (particles[ii].vx-particles[jj].vx)*(particles[ii].x-particles[jj].x)
                         +(particles[ii].vy-particles[jj].vy)*(particles[ii].y-particles[jj].y)
                         +(particles[ii].vz-particles[jj].vz)*(particles[ii].z-particles[jj].z);
            dudt += particles[jj].m*grad_hi*dot_r_v;

            particles[ii].ax += particles[jj].m *f_ij[0];
            particles[ii].ay += particles[jj].m *f_ij[1];
            particles[ii].az += particles[jj].m *f_ij[2];
            particles[jj].ax -= particles[ii].m *f_ij[0];
            particles[jj].ay -= particles[ii].m *f_ij[1];
            particles[jj].az -= particles[ii].m *f_ij[2];
        }
        particles[ii].du = p_i / (omeg_i*particles[ii].h*particles[ii].h) *dudt;
    }
}

pub fn euler_integrator(particle: &mut Particle, dt: f64) {
    particle.x += dt * particle.vx;
    particle.y += dt * particle.vy;
    particle.z += dt * particle.vz;
    particle.vx += dt * particle.ax;
    particle.vy += dt * particle.ay;
    particle.vz += dt * particle.az;
    particle.u += dt * particle.du;
}