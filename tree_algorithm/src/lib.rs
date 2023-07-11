use std::{
    io,
    f64,
    u32,
};

use csv::Writer;

use structures::{
    Particle,
    Node,
};

pub trait BuildTree {

    fn new(n_p: u32, x0: f64, y0: f64, l: f64) -> Node;

    fn branching_factor(& self, k: f64, s:f64) -> u32;

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Node;

    fn create_sub_cells(&mut self, b: u32);

    fn delete_sub_cells(&mut self);

    fn delete_particles(&mut self);

    fn distribution_ratio(&self, limit: u32, b: u32) -> f64;

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64);

    fn restart(&mut self, n: usize);
}

use rayon::prelude::*;

impl BuildTree for Node {

    fn new(n_p: u32, x0: f64, y0: f64, l: f64) -> Node {
        Node{n: n_p, particles: (0..n_p as usize).collect(), 
             xmin: x0,
             ymin: y0,
             side: l,
             ..Default::default()}
    }
    
    fn branching_factor(& self, k: f64, s:f64) -> u32 {
        ((self.n as f64 /s).powf(1./k)).ceil() as u32
    }

    fn create_child(&self, j: u32, b: u32, dx: f64) -> Node {
        Node{xmin: self.xmin + dx*(j%b) as f64,
            ymin: self.ymin + dx*(j/b) as f64,
            side: dx,
            id: j,
            depth: self.depth + 1,
            n: 0,
            branches: 0,
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

    fn build_tree(&mut self, k: u32, s: u32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64) {
        let mut redistribution :bool = true;
        let mut b = self.branching_factor(k as f64, s as f64);
        while redistribution {
            self.branches = b.pow(k);
            self.create_sub_cells(b);
            for p in &self.particles {
                let mut x_p = ((particles[*p].x - self.xmin) / self.side * b as f64).floor() as u32;
                if x_p == b {
                    x_p -= 1;
                }
                let mut y_p = ((particles[*p].y - self.ymin) / self.side * b as f64).floor() as u32;
                if y_p == b {
                    y_p -= 1;
                }
                let j :usize = (x_p + y_p * b) as usize;
                add_particle(&mut self.children[j], *p);
            }
            let r = self.distribution_ratio((alpha * s as f64) as u32, b);
            if r >= beta {
                b = b/2;
                self.delete_sub_cells();
            } else {
                redistribution = false;
                self.delete_particles();
            }
        }
        (self.children).par_iter_mut().for_each(|child| {
            if (child.n > s) && (child.side > smallest_cell) {
                child.build_tree(k, s, alpha, beta, particles, smallest_cell);
            }
        });
    }

    fn restart(&mut self, n: usize) {
        self.branches = 0;
        self.particles = (0..n).collect();
        self.delete_sub_cells();
    }
}

pub trait FindNeighbors {
    fn range_neigh(&self, x_p: f64, y_p: f64, h: f64, b: u32) -> (u32, u32, u32, u32);

    fn children_in_range(&self, xmin: u32, xmax: u32, ymin: u32, ymax:u32, b:u32) -> Vec<usize>;

    fn find_neighbors(& self, p: usize, k: f64, s: u32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, x_side: f64, y_side: f64, h: f64);
}

impl FindNeighbors for Node {

    fn range_neigh(&self, x_p: f64, y_p: f64, h: f64, b: u32) -> (u32, u32, u32, u32){
        let factor : f64 =  b as f64 /self.side;
        let x_min = (((x_p - 2.0*h) - self.xmin) * factor).floor() as i32;
        let x_max = (((x_p + 2.0*h) - self.xmin) * factor).floor() as i32;
        let y_min = (((y_p - 2.0*h) - self.ymin) * factor).floor() as i32;
        let y_max = (((y_p + 2.0*h) - self.ymin) * factor).floor() as i32;
        ((x_min).rem_euclid(b.try_into().unwrap()) as u32, (x_max).rem_euclid(b.try_into().unwrap()) as u32,
         (y_min).rem_euclid(b.try_into().unwrap()) as u32, (y_max).rem_euclid(b.try_into().unwrap()) as u32)
    }

    fn children_in_range(&self, xmin: u32, xmax: u32, ymin: u32, ymax:u32, b:u32) -> Vec<usize>{
        let mut neighbors : Vec<usize> = Vec::new();
        for child in &self.children{
            let x_id: u32 = child.id%b;
            let y_id: u32 = child.id/b;
            let index: u32 = index_range(xmin, ymin, xmax, ymax);

            if index == 0 { 
                if (x_id >= xmin) && (y_id >= ymin) {
                    if (x_id <= xmax) && (y_id <= ymax) {
                        neighbors.push(child.id as usize);
                    }
                }
            } else if index == 1 {
                if (x_id >= xmin) && (x_id <= xmax) {
                    if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                        neighbors.push(child.id as usize);
                    }
                }
            } else if index == 2 {
                if (y_id >= ymin) && (y_id <= ymax) {
                    if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                        neighbors.push(child.id as usize);
                    }
                }
            } else {
                if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                    if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                        neighbors.push(child.id as usize);
                    }
                }
            }
        }
        neighbors
    }

    fn find_neighbors(& self, p: usize, k: f64, s: u32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, x_side: f64, y_side:f64, h: f64) {
        let b = (self.branches as f64).powf(1./k) as u32;
        let (x_min, x_max, y_min, y_max) = self.range_neigh(particles[p].x, particles[p].y, h, b);
        // println!("{} - {} {} {} {}", b, x_min, x_max, y_min, y_max);
        let neighbors = self.children_in_range(x_min, x_max, y_min, y_max, b);
        for ii in neighbors {
            if self.children[ii].n <= s {
                for q in &self.children[ii].particles {
                    if periodic_norm(particles[p].x, particles[*q].x, particles[p].y, particles[*q].y, x_side, y_side, 2.*h) <= 4.0*h*h {
                        neighbors_of_p.push(*q);
                    }
                }
            } else {
                self.children[ii].find_neighbors(p, k, s, particles, neighbors_of_p, x_side, y_side, h);
            }
        }
    }
}

fn index_range(xmin: u32, ymin: u32, xmax: u32, ymax: u32) -> u32 {
    if xmin < xmax {
        if ymin < ymax {
            return 0;
        } else {
            return 1;
        }
    } else {
        if ymin < ymax {
            return 2;
        } else {
            return 3;
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


// Periodic Distance
pub fn periodic_norm(x1: f64, x2: f64, y1: f64, y2: f64, w: f64, h: f64, eps: f64) -> f64 {
    
    let mut x_temp: f64 = x1 - x2;
    let mut y_temp: f64 = y1 - y2;

    if x_temp.abs() > w-2.*eps {
        if x_temp > 0. {
            x_temp -= w;
        } else {
            x_temp += w;
        }
    }
    if y_temp.abs() > h-2.*eps {
        if y_temp > 0. {
            y_temp -= h;
        } else {
            y_temp += h;
        }
    }
    return x_temp*x_temp + y_temp*y_temp;
}