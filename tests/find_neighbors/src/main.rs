// Tree algorithm: Tree builder and neighbors finder 
use std::{
    error::Error,
    process,
    time::Instant,
};

use datafunctions;
use partdistribution;
use sphfunctions;

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

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./FindNeigh/set_particles.csv";
    let path_tree = "./FindNeigh/set_tree.csv";
    let path_neighbors = "./FindNeigh/set_neighbors.csv";
    let mut n:u32 = 16; // Number of Particles
    let x0:f64 = -0.5; // circle's center
    let y0:f64 = -0.5; // circle's center
    let z0:f64 = -0.5; // circle's center
    let wd:f64 = 1.; // circle's center
    let lg:f64 = 1.; // circle's center
    let hg:f64 = 1.; // circle's center
    let rho:f64= 1.; // Smoothing length
    let rkern :f64 = 2.; // Smoothing length
    let eta:f64= 1.2;

    let xper:bool = true;
    let yper:bool = true;
    let zper:bool = true;

    let mut particles :Vec<Particle> = Vec::new();
    
    // Initialize system
    partdistribution::init_dist_hcp(&mut particles, n, n, n, rho, eta, wd, lg, hg, x0, y0, z0, 0.1);
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
    n = particles.len() as u32;
    let dm: f64 = 1.0/n as f64;
    let h: f64 = sphfunctions::h_by_density(dm, rho, eta);

    // Tree parameters
    let s : i32 = 4;
    let alpha : f64 = 0.5;
    let beta : f64 = 0.5;
    
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