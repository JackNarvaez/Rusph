use std::{
    error::Error,
    process,
};

use sphfunctions;
use tree_algorithm::{
    Node,
    BuildTree,
    FindNeighbors,
    save_tree,
    save_neighbors,
};

fn main() -> Result<(), Box<dyn Error>> {
    let path = "./Data/tree_algorithm/set_particles.csv";
    let path_tree = "./Data/tree_algorithm/set_tree.csv";
    let path_neighbors = "./Data/tree_algorithm/set_neighbors.csv";
    let n:u32 = 1000; // Number of Particles
    let dm:f64 = 0.1; // particle's mass
    let w:f64 = 1.; // width
    let l:f64 = 1.; // large
    let rho:f64 = 1.0; // density
    let h = 0.1; // Smoothing length
    if let Err(err) = sphfunctions::init_random_square(path, n, dm, rho, h, w, l){
        println!("{}", err);
        process::exit(1);
    }
    let mut particles :Vec<sphfunctions::Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    let k : u32 = 2;
    let s : u32 = 10;
    let alpha : f64 = 0.5;
    let beta : f64 = 0.5;
    let mut root : Node = <Node as BuildTree>::new(n);
    root.build_tree(k, s, alpha, beta, &particles);
    save_tree(path_tree, &root);
    let mut neighbors: Vec<usize> = Vec::new();
    let p: usize = 1;
    root.find_neighbors(p, k as f64, s, &particles, &mut neighbors);
    save_neighbors(path_neighbors, p, & neighbors);
    Ok(())
}