// Tree algorithm: Tree builder and neighbors finder 

use std::{
    error::Error,
    process,
    time::Instant,
};

use tree_algorithm::{
    BuildTree,
    FindNeighbors,
    save_tree,
    save_neighbors
};

use structures::{
    Particle,
    Node,
};

use datafunctions;
use partdistribution;
use sphfunctions;

fn main() -> Result<(), Box<dyn Error>> {
    let path: &str          = "./FindNeigh/set_particles.csv";
    let path_tree: &str     = "./FindNeigh/set_tree.csv";
    let path_neighbors: &str= "./FindNeigh/set_neighbors.csv";
    let nx: u32     = 16;   // Particle resolution
    let x0: f64     = -0.5; // Bottom left corner  (x-coordinate)
    let y0: f64     = -0.5; // Bottom left corner  (y-coordinate)
    let z0: f64     = -0.5; // Bottom left corner  (z-coordinate)
    let wd: f64     = 1.;   // Width (x)
    let lg: f64     = 1.;   // Length (y)
    let hg: f64     = 1.;   // Height (z)
    let rho:f64     = 1.;   // Density
    let rkern: f64  = 2.;   // Kernel radius
    let eta: f64    = 1.2;  // Dimensionless constant specifying the smoothing length

    let xper: bool   = true;
    let yper: bool   = true;
    let zper: bool   = true;

    let mut particles :Vec<Particle> = Vec::new();
    
    // Initialize system
    partdistribution::init_dist_hcp(&mut particles, nx, rho, eta, wd, lg, hg, x0, y0, z0);
    if let Err(err) = datafunctions::save_data(path, &particles){
        println!("{}", err);
        process::exit(1);
    }

    // Read data
    particles.clear();
    if let Err(err) = datafunctions::read_data(path, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }

    let n: u32  = particles.len() as u32;
    let dm: f64 = rho*wd*lg*hg/n as f64;
    let h: f64  = sphfunctions::h_by_density(dm, rho, eta);

    // Tree parameters
    let s: i32      = 4;
    let alpha: f64  = 0.5;
    let beta: f64   = 0.5;
    
    // Tree builder
    let start1 = Instant::now();
    let mut root : Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd, lg, hg);
    root.build_tree(s, alpha, beta, &particles, 0.1*h);
    println!("Tree Builder: {} s", start1.elapsed().as_secs());
    save_tree(path_tree, &root);

    // Neighbors finder
    let start2 = Instant::now();
    let mut neighbors: Vec<usize> = Vec::new();
    for p in 0..n as usize {
        root.find_neighbors(p, s, &particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[p].h*rkern, xper, yper, zper);
        neighbors = Vec::new();
    }
    let p: usize = 124;
    root.find_neighbors(p, s, &particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[p].h*rkern, xper, yper, zper);
    
    println!("Neighbors Finder: {} s", start2.elapsed().as_secs());
    save_neighbors(path_neighbors, p, & neighbors);
    Ok(())
}