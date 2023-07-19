// Tree algorithm: Tree builder and neighbors finder 

use std::{
    error::Error,
    process,
    time::Instant,
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
    let n:u32 = 10000; // Number of Particles
    let x0:f64 = 0.; // circle's center
    let y0:f64 = 0.; // circle's center
    let z0:f64 = 0.; // circle's center
    let wd:f64 = 1.; // circle's center
    let lg:f64 = 1.; // circle's center
    let hg:f64 = 1.; // circle's center
    let h :f64 = 0.1; // Smoothing length
    // Initialize system
    if let Err(err) = sphfunctions::init_random_square(path, n, h, wd, lg, hg, x0, y0, z0){//init_random_circle(path, n, r, rho, h, x0, y0){
        println!("{}", err);
        process::exit(1);
    }
    // Read data
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = sphfunctions::read_data(path, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }
    // Tree parameters
    let k : u32 = 2; // Dimension
    let s : i32 = 10;
    let alpha : f64 = 0.5;
    let beta : f64 = 0.5;

    // Tree builder
    let start1 = Instant::now();
    let mut root : Node = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd);
    for _ii in 0..10000 {
        root = <Node as BuildTree>::new(n as i32, x0, y0, z0, wd);
        root.build_tree(k, s, alpha, beta, &particles, 0.1*h);
    }
    println!("Tree Builder: {} s", start1.elapsed().as_secs());
    save_tree(path_tree, &root);

    // Neighbors finder
    let start2 = Instant::now();
    let mut neighbors: Vec<usize> = Vec::new();
    let p: usize = 0;
    for _ii in 0..100000 {
        neighbors = Vec::new();

        root.find_neighbors(p, k as f64, s, &particles, &mut neighbors, wd, lg, hg, particles[p].h);
    }
    println!("Neighbors Finder: {} s", start2.elapsed().as_secs());
    save_neighbors(path_neighbors, p, & neighbors);
    Ok(())
}