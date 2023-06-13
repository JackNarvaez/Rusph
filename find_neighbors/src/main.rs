use std::{
    error::Error,
    process,
};

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
    let path = "./Data/tree_algorithm/set_particles.csv";
    let path_tree = "./Data/tree_algorithm/set_tree.csv";
    let path_neighbors = "./Data/tree_algorithm/set_neighbors.csv";
    let n:u32 = 1000; // Number of Particles
    let dm:f64 = 0.1; // particle's mass
    let x0:f64 = 0.; // circle's center
    let y0:f64 = 0.; // circle's center
    let r:f64 = 0.75; // radius
    let rho:f64 = 1.0; // density
    let h = 0.1; // Smoothing length
    if let Err(err) = sphfunctions::init_random_circle(path, n, r, dm, rho, h, x0, y0){
        println!("{}", err);
        process::exit(1);
    }
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let k : u32 = 2;
    let s : u32 = 10;
    let alpha : f64 = 0.5;
    let beta : f64 = 0.5;
    let mut root : Node = <Node as BuildTree>::new(n, x0-2.*r, y0-2.*r, 4.*r);
    root.build_tree(k, s, alpha, beta, &particles, 0.1*h);
    save_tree(path_tree, &root);
    let mut neighbors: Vec<usize> = Vec::new();
    let p: usize = 1;
    root.find_neighbors(p, k as f64, s, &particles, &mut neighbors);
    save_neighbors(path_neighbors, p, & neighbors);
    Ok(())
}