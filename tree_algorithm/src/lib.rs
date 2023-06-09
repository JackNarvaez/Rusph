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

    fn distribution_ratio(&self, limit: u32) -> f64;

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
        println!("dx  {} {}", dx, self.id);
        for ii in 0..self.branches {
            self.children.push(self.create_child(ii, b, dx));
        }
    }

    fn delete_sub_cells(&mut self) {
        self.children.clear();// = Vec::new();
    }

    fn delete_particles(&mut self) {
        self.particles.clear();//= Vec::new();
    }

    fn distribution_ratio(&self, limit: u32) -> f64 {
        let mut r :f64 = 0.0;
        for child in &self.children{
            if child.n <= limit {
                r += 1.0;
            }
        }
        r / self.branches as f64
    }

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64, particles: & Vec<sphfunctions::Particle>) {
        let mut redistribution :bool = true;
        let smallest_cell = 1.0e-2;
        let mut b = self.branching_factor(k as f64, s as f64);
        while redistribution {
            self.branches = b.pow(k);
            println!("{} {} {}",self.depth, b, self.side);
            self.create_sub_cells(b);
            for p in &self.particles {
                let x_p = ((particles[*p].x - self.xmin) / self.side * b as f64).floor();
                let y_p = ((particles[*p].y - self.ymin) / self.side * b as f64).floor();
                let j :usize = (x_p + y_p*b as f64) as usize;
                add_particle(&mut self.children[j], *p);
            }
            let r = self.distribution_ratio((alpha * s as f64) as u32);
            println!("{}", r);
            if r < beta {
                b = b/2;
                //println!("{} {}", self.id, self.n);
                self.delete_sub_cells();
                //println!("Children  {}  {:?}", self.depth, self.children)
            } else {
                //println!("yep {} {}", self.side, self.branches);
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

fn add_particle(cell: &mut Node, i: usize) {
    cell.particles.push(i);
    cell.n += 1;
}

pub fn save_tree(path: &str, tree: & Node){
    let mut wtr = (Writer::from_path(path)).expect("REASON");
    wtr.write_record(&["x_min", "y_min", "side", "depth", "n"]);
    save_child(&mut wtr, tree);
}

fn save_child<W: io::Write>(wtr: &mut Writer<W>, tree: & Node){
    wtr.write_record(&[tree.xmin.to_string(), tree.ymin.to_string(), tree.side.to_string(), tree.depth.to_string(), tree.n.to_string()]);
    for child in &tree.children{
        save_child(wtr, child);
    }
}