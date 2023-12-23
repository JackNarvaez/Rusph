// Tree algorithm: Tree builder and neighbors finder 

use std::{
    error::Error,
    process,
    time::Instant,
};

use datafunctions;

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
    let n:u32 = 32*32*32; // Number of Particles
    let x0:f64 = 0.; // circle's center
    let y0:f64 = 0.; // circle's center
    let z0:f64 = 0.; // circle's center
    let wd:f64 = 1.; // circle's center
    let lg:f64 = 1.; // circle's center
    let hg:f64 = 1.; // circle's center
    let h :f64 = 0.1; // Smoothing length
    let rkern :f64 = 2.; // Smoothing length

    let xper:bool = true;
    let yper:bool = true;
    let zper:bool = true;


    // Initialize system
    if let Err(err) = datafunctions::init_square(path, n, h, wd, lg, hg, x0, y0, z0){
        println!("{}", err);
        process::exit(1);
    }
    // Read data
    let mut particles :Vec<Particle> = Vec::new();
    if let Err(err) = datafunctions::read_data(path, &mut particles) {
        println!("{}", err);
        process::exit(1);
    }

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
    let p: usize = 10000;
    root.find_neighbors(p, s, &particles, &mut neighbors, wd, lg, hg, x0, y0, z0, particles[p].h*rkern, xper, yper, zper);
    
    println!("Neighbors Finder: {} s", start2.elapsed().as_secs());
    save_neighbors(path_neighbors, p, & neighbors);
    Ok(())
}