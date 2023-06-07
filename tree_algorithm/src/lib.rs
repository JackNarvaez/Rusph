use std::{
    error::Error,
    f64,
};

pub struct Cell {
	pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub particles: Vec<u32>,
    pub n: u32,
}

fn branching_factor(k: f64, n: f64, s:f64) -> f64 {
    ((n/s).powf(1./k)).ceil()
}

fn 2d_to_1d(i: f64; j :f64; b: f64) -> f64 {
    i + j*b
}

fn add_particle(i:f64, cell: Cell) {
    cell.particles.append(i);
}

fn distribute_particle(x: f64, y: f64){

}

fn build_tree(cell: Cell, i: u32) {
    let redistribution :bool = True;
    while redistribution {
        let b = branching_factor(k, cell.n, s);
        let sub_cells = b.powf(k) as u32;
        for p in cell.particles {
            let x_p = ((particles[p].x - cell.xmin) / side * b ).floor();
            let y_p = ((particles[p].y - cell.ymin) / side * b ).floor();
            let j = 2d_to_1d(x_p;y_p; b);
        }
    }
}