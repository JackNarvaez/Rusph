use std::{
    error::Error,
    f64,
    u32,
};

pub struct Nodes {
	pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub depth: u32,
    pub n: u32,
    //pub children: Vec<address>,
}

pub struct Leaves {
	pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub depth: u32,
    pub n: u32,
    pub particles: Vec<u32>,
}

pub struct Tree {
    Nodes root;
}

impl build_tree for Tree {
    fn build_tree {
        
    }
}

fn branching_factor(k: f64, n: f64, s:f64) -> u32 {
    ((n/s).powf(1./k)).ceil() as u32
}

fn 2d_to_1d(i: u32; j :u32; b: u32) -> u32 {
    i + j*b
}

fn add_particle(i:f64, cell: Cell) {
    cell.particles.append(i);
}

fn distribute_particle(x: f64, y: f64){

}

fn build_tree(cell: Cell, i: u32, k: u32, s: u32) {
    let redistribution :bool = True;
    while redistribution {
        let b = branching_factor(k as f64, cell.n as f64, s as f64);
        let sub_cells = b.pow(k);
        let sub_cells = 
        for p in cell.particles {
            let x_p = ((particles[p].x - cell.xmin) / cell.side * b as f64).floor() as u32;
            let y_p = ((particles[p].y - cell.ymin) / cell.side * b as f64).floor() as u32;
            let j = 2d_to_1d(x_p;y_p; b);

        }
    }
}