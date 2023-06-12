use std::{
    io,
    f64,
    u32,
};

use csv::Writer;

use sphfunctions;

#[derive(Debug)]
pub struct Node {
	pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub depth: u32,
    pub n: u32,
    pub branches: u32,
    pub leave: bool,
    pub children: Vec<Node>,
    pub particles: Vec<usize>,
}

impl Default for Node {
    fn default() -> Node {
        Node {
            xmin: 0.0,
            ymin: 0.0,
            side: 1.,
            id: 0,
            depth: 0,
            n: 0,
            branches: 0,
            leave: false,
            children: Vec::new(),
            particles: Vec::new(),
        }
    }
}

pub trait BuildTree {

    fn new(n_p: u32) -> Node;

    fn branching_factor(& self, k: f64, s:f64) -> u32;

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Node;

    fn create_sub_cells(&mut self, b: u32);

    fn delete_sub_cells(&mut self);

    fn delete_particles(&mut self);

    fn distribution_ratio(&self, limit: u32, b: u32) -> f64;

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64, particles: & Vec<sphfunctions::Particle>);
}

impl BuildTree for Node {

    fn new(n_p: u32) -> Node {
        Node{n: n_p, particles: (0..n_p as usize).collect(), ..Default::default()}
    }
    fn branching_factor(& self, k: f64, s:f64) -> u32 {
        ((self.n as f64 /s).powf(1./k)).ceil() as u32
    }

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Node {
        Node{xmin: self.xmin + dx*(j%b) as f64,
            ymin: self.ymin + dx*(j/b) as f64,
            side: dx,
            id: j,
            depth: self.depth +1,
            n: 0,
            branches: 0,
            leave: true,
            children: Vec::new(),
            particles: Vec::new(),
        }
    }

    fn create_sub_cells(&mut self, b: u32) {
        let dx = self.side / b as f64;
        for ii in 0..self.branches {
            self.children.push(self.create_child(ii, b, dx));
        }
    }

    fn delete_sub_cells(&mut self) {
        self.children.clear();
    }

    fn delete_particles(&mut self) {
        self.particles.clear();
    }

    fn distribution_ratio(&self, limit: u32, b:u32) -> f64 {
        let mut r :f64 = 0.0;
        for child in &self.children{
            if child.n <= limit {
                r += 1.0;
            }
        }
        r / (b.pow(self.depth) as f64)
    }

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64, particles: & Vec<sphfunctions::Particle>) {
        let mut redistribution :bool = true;
        let smallest_cell = 1.0e-2;
        let mut b = self.branching_factor(k as f64, s as f64);
        while redistribution {
            self.branches = b.pow(k);
            self.create_sub_cells(b);
            for p in &self.particles {
                let x_p = ((particles[*p].x - self.xmin) / self.side * b as f64).floor();
                let y_p = ((particles[*p].y - self.ymin) / self.side * b as f64).floor();
                let j :usize = (x_p + y_p*b as f64) as usize;
                add_particle(&mut self.children[j], *p);
            }
            let r = self.distribution_ratio((alpha * s as f64) as u32, b);
            if r >= beta {
                b = b/2;
                self.delete_sub_cells();
            } else {
                self.delete_particles();
                redistribution = false;
            }
        }
        for child in &mut self.children {
            if (child.n > s) && (child.side > smallest_cell) {
                child.build_tree(k, s, alpha, beta, particles);
            } else {
                child.leave = true;
            }
        }
    }
}

pub trait FindNeighbors {
    fn range_neigh(&self, x_p: f64, y_p: f64, h: f64, b: f64) -> (u32, u32, u32, u32);

    fn children_in_range(&self, xmin: u32, xmax: u32, ymin: u32, ymax:u32, b:u32) -> Vec<usize>;

    fn find_neighbors(&mut self, p: usize, k: f64, s: u32, particles: & Vec<sphfunctions::Particle>, neighbors_of_p: &mut Vec<usize>);
}

impl FindNeighbors for Node {

    fn range_neigh(&self, x_p: f64, y_p: f64, h: f64, b: f64) -> (u32, u32, u32, u32){
        let factor : f64 =  b/self.side;
        let x_min = (((x_p - 2.0*h) - self.xmin) * factor).floor() as u32;
        let x_max = (((x_p + 2.0*h) - self.xmin) * factor).floor() as u32;
        let y_min = (((y_p - 2.0*h) - self.ymin) * factor).floor() as u32;
        let y_max = (((y_p + 2.0*h) - self.ymin) * factor).floor() as u32;
        (x_min, x_max, y_min, y_max)
    }

    fn children_in_range(&self, xmin: u32, xmax: u32, ymin: u32, ymax:u32, b:u32) -> Vec<usize>{
        let mut neighbors : Vec<usize> = Vec::new();
        for child in &self.children{
            if (child.id%b >= xmin) && (child.id/b >= ymin) {
                if (child.id%b <= xmax) && (child.id/b <= ymax) {
                    neighbors.push(child.id as usize);
                }
            }
        }
        neighbors
    }

    fn find_neighbors(&mut self, p: usize, k: f64, s: u32, particles: & Vec<sphfunctions::Particle>, neighbors_of_p: &mut Vec<usize>) {
        let b = (self.branches as f64).powf(1./k) as u32;
        let (x_min, x_max, y_min, y_max) = self.range_neigh(particles[p].x, particles[p].y, particles[p].h, b as f64);
        let neighbors = self.children_in_range(x_min, x_max, y_min, y_max, b);
        for ii in neighbors {
            if self.children[ii].n <= s {
                for q in &self.children[ii].particles {
                    if sphfunctions::euclidean_norm(&particles[p], &particles[*q]) <= 2.0*particles[p].h {
                        neighbors_of_p.push(*q);
                    }
                }
            } else {
                self.children[ii].find_neighbors(p, k, s, particles, neighbors_of_p);
            }
        }
    }
}

fn add_particle(cell: &mut Node, i: usize) {
    cell.particles.push(i);
    cell.n += 1;
}

pub fn save_tree(path: &str, tree: & Node){
    let mut wtr = (Writer::from_path(path)).expect("REASON");
    wtr.write_record(&["x_min", "y_min", "side", "depth", "n"]).expect("Couldn't write data");
    save_child(&mut wtr, tree);
}

fn save_child<W: io::Write>(wtr: &mut Writer<W>, tree: & Node){
    wtr.write_record(&[tree.xmin.to_string(), tree.ymin.to_string(), tree.side.to_string(), tree.depth.to_string(), tree.n.to_string()]).expect("Couldn't write data");
    for child in &tree.children{
        save_child(wtr, child);
    }
}

pub fn save_neighbors(path: &str, p: usize, neighbors: & Vec<usize>){
    let mut wtr = (Writer::from_path(path)).expect("REASON");
    wtr.write_record(&["p"]).expect("Couldn't write data");
    wtr.write_record(&[p.to_string()]).expect("Couldn't write data");
    for ii in neighbors {
        wtr.write_record(&[ii.to_string()]).expect("Couldn't write data");
    }
}