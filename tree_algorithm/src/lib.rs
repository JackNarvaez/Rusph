use std::{
    error::Error,
    f64,
    u32,
};

pub struct Node {
	pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub depth: u32,
    pub n: u32,
    pub branches: u32,
    pub leave: bool,
    pub children: Vec<Nodes>,
    pub particles: Vec<u32>,
}
trait build_tree {
    fn branching_factor(& self, k: f64, s:f64) -> u32;

    fn 2d_to_1d(i: u32; j :u32; b: u32) -> u32;
    
    fn add_particle(&self, i:f64);

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Box<dyn Node>;

    fn create_sub_cells(&mut self, b: u32);

    fn delete_sub_cells(&mut self);

    fn distribution_ratio(&self, limit: f64) -> f64;

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64);
}

impl build_tree for Node {

    fn branching_factor(& self, k: f64, s:f64) -> u32 {
        ((self.n/s).powf(1./k)).ceil() as u32
    }

    fn 2d_to_1d(i: u32; j :u32; b: u32) -> u32 {
        i + j*b
    }
    
    fn add_particle(&self, i:f64) {
        self.particles.push(i);
        self.n += 1;
    }

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Box<dyn Node> {
        return Box::new(Node{xmin: self.xmin + dx*j.div_floor(b),
                             ymin: self.ymin + dx*(j%b),
                             side: dx,
                             id: j,
                             depth: self.depth +1,
                             n: 0,
                             branches: 0,
                             leave: true,
                             children: Vec::new(),
                             particles: Vec::new(),
        });
    }

    fn create_sub_cells(&mut self, b: u32) {
        let dx = self.side / b as f64; 
        for ii in 0..self.branches {
            self.children.push(create_child(self, ii, b, dx));
        }
    }

    fn delete_sub_cells(&mut self) {
        self.children = Vec::new();
    }

    fn delete_particles(&mut self) {
        self.particles = Vec::new();
    }

    fn distribution_ratio(&self, limit: f64) -> f64 {
        let r :f64 = 0.0;
        for child in self.children{
            if child.n <= limit {
                r += child.n;
            }
        }
        r / self.branches as f64
    }

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64) {
        let redistribution :bool = True;
        let smallest_cell = 1.0e-5;
        let mut b = branching_factor(self, k as f64, s as f64);
        while redistribution {
            self.branches = b.pow(k);
            create_sub_cells(self, b);
            for p in self.particles {
                let x_p = ((particles[p].x - self.xmin) / self.side * b as f64).floor() as u32;
                let y_p = ((particles[p].y - self.ymin) / self.side * b as f64).floor() as u32;
                let j = 2d_to_1d(x_p;y_p; b);
                add_particle(self.children[j], p);
            }
            let r = distribution_ratio(self, alpha * s as f64);
            if r > beta {
                b /= 2;
                delete_sub_cells(self);
            } else {
                delete_particles(self);
                redistribution = false;
            }
        }
        for child in self.children {
            if (child.n > s) && (child.side > smallest_cell) {
                build_tree(self, k, s, alpha, beta);
            } else {
                child.leave = true;
            }
        }

    }
}