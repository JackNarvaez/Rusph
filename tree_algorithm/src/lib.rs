use std::{
    io,
    f64,
    i32,
};

use csv::Writer;

use structures::{
    Particle,
    Node,
};

use rayon::prelude::*;

pub trait BuildTree {

    fn new(n_p: i32, x0: f64, y0: f64, z0: f64, wd: f64, lg: f64, hg: f64) -> Node;

    fn branching_factor(& self, k: f64, s:f64) -> i32;

    fn create_child(&self, j: i32, b: i32, dx: f64, dy: f64, dz: f64) -> Node;

    fn create_sub_cells(&mut self, b: i32);

    fn delete_sub_cells(&mut self);

    fn delete_particles(&mut self);

    fn distribution_ratio(&self, limit: i32) -> f64;

    fn build_tree(&mut self, k: u32, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64);

    fn build_octtree(&mut self, k: u32, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64);

    fn restart(&mut self, n: usize);
}

impl BuildTree for Node {

    fn new(n_p: i32, x0: f64, y0: f64, z0: f64, wd: f64, lg: f64, hg: f64) -> Node {
        Node{n: n_p, particles: (0..n_p as usize).collect(), 
             xmin: x0,
             ymin: y0,
             zmin: z0,
             sidex: wd,
             sidey: lg,
             sidez: hg,
             ..Default::default()}
    }
    
    fn branching_factor(& self, k: f64, s:f64) -> i32 {
        ((self.n as f64 /s).powf(1./k)).ceil() as i32
    }

    fn create_child(&self, j: i32, b: i32, dx: f64, dy: f64, dz: f64) -> Node {
        let tem: i32 = j/b;
        Node{xmin: self.xmin + dx * (j % b) as f64,
            ymin: self.ymin + dy * (tem % b) as f64,
            zmin: self.zmin + dz * (tem / b) as f64,
            sidex: dx,
            sidey: dy,
            sidez: dz,
            id: j,
            depth: self.depth + 1,
            n: 0,
            branches: 0,
            children: Vec::new(),
            particles: Vec::new(),
        }
    }

    fn create_sub_cells(&mut self, b: i32) {
        let dx: f64 = self.sidex / (b as f64);
        let dy: f64 = self.sidey / (b as f64);
        let dz: f64 = self.sidez / (b as f64);
        for ii in 0..self.branches {
            self.children.push(self.create_child(ii, b, dx, dy, dz));
        }
    }

    fn delete_sub_cells(&mut self) {
        self.children.clear();
    }

    fn delete_particles(&mut self) {
        self.particles.clear();
    }

    fn distribution_ratio(&self, limit: i32) -> f64 {
        let mut r :f64 = 0.0;
        for child in &self.children{
            if child.n <= limit {
                r += 1.0;
            }
        }
        r / (self.children).len() as f64
    }

    fn build_tree(&mut self, k: u32, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64) {
        let mut redistribution :bool = true;
        let mut b: i32 = self.branching_factor(k as f64, s as f64);
        while redistribution {
            self.branches = b.pow(k);
            self.create_sub_cells(b);
            for p in &self.particles {
                let mut x_p:i32 = ((particles[*p].x - self.xmin) / self.sidex * b as f64).floor() as i32;
                if x_p == b {
                    println!("ALERT!!");
                    // x_p -= 1;
                }
                let mut y_p: i32 = ((particles[*p].y - self.ymin) / self.sidey * b as f64).floor() as i32;
                if y_p == b {
                    println!("ALERT!!");
                    // y_p -= 1;
                }
                let mut z_p: i32 = ((particles[*p].z - self.zmin) / self.sidez * b as f64).floor() as i32;
                if z_p == b {
                    println!("ALERT!!");
                    // z_p -= 1;
                }
                let j :usize = (x_p + y_p * b + z_p * b * b) as usize;
                add_particle(&mut self.children[j], *p);
            }
            let mut r: f64 = 0.0;
            if b != 2 {
                r = self.distribution_ratio((alpha * s as f64) as i32);
            }
            if r > beta {
                if b > 4 {
                    b = b/2;
                } else {
                    b = 2;
                }
                self.delete_sub_cells();
            } else {
                redistribution = false;
                self.delete_particles();
            }
        }
        (self.children).par_iter_mut().for_each(|child| {
            if (child.n > s) && (child.sidex > smallest_cell) && (child.sidey > smallest_cell) && (child.sidez > smallest_cell) {
                child.build_tree(k, s, alpha, beta, particles, smallest_cell);
            }
        });
    }

    fn build_octtree(&mut self, k: u32, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64) {
        let b:i32 = 2;
        self.branches = b.pow(k);
        self.create_sub_cells(b);
        for p in &self.particles {
            let mut x_p: i32 = ((particles[*p].x - self.xmin) / self.sidex * b as f64).floor() as i32;
            if x_p == b {
                x_p -= 1;
            }
            let mut y_p: i32 = ((particles[*p].y - self.ymin) / self.sidey * b as f64).floor() as i32;
            if y_p == b {
                y_p -= 1;
            }
            let mut z_p: i32 = ((particles[*p].z - self.zmin) / self.sidez * b as f64).floor() as i32;
            if z_p == b {
                z_p -= 1;
            }
            let j :usize = (x_p + y_p * b + z_p * b * b) as usize;
            add_particle(&mut self.children[j], *p);
        }
        self.delete_particles();
        (self.children).par_iter_mut().for_each(|child| {
            if (child.n > s) && (child.sidex > smallest_cell) && (child.sidey > smallest_cell) && (child.sidez > smallest_cell) {
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
    fn range_neigh(&self, x_p: f64, y_p: f64, z_p: f64, h: f64, b: i32, rkern: f64) -> (i32, i32, i32, i32, i32, i32);

    fn children_in_range(&self, xmin: i32, xmax: i32, ymin: i32, ymax:i32, zmin: i32, zmax:i32, b:i32) -> Vec<usize>;

    fn find_neighbors(& self, p: usize, k: f64, s: i32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, wd: f64, lg: f64, hg: f64, h: f64, rkern: f64);
}

impl FindNeighbors for Node {

    fn range_neigh(&self, x_p: f64, y_p: f64, z_p: f64, h: f64, b: i32, rkern: f64) -> (i32, i32, i32, i32, i32, i32){
        let factorx : f64 =  b as f64 /self.sidex;
        let factory : f64 =  b as f64 /self.sidey;
        let factorz : f64 =  b as f64 /self.sidez;
        let x_min: i32 = (((x_p - rkern*h) - self.xmin) * factorx).floor() as i32;
        let x_max: i32 = (((x_p + rkern*h) - self.xmin) * factorx).floor() as i32;
        let y_min: i32 = (((y_p - rkern*h) - self.ymin) * factory).floor() as i32;
        let y_max: i32 = (((y_p + rkern*h) - self.ymin) * factory).floor() as i32;
        let z_min: i32 = (((z_p - rkern*h) - self.zmin) * factorz).floor() as i32;
        let z_max: i32 = (((z_p + rkern*h) - self.zmin) * factorz).floor() as i32;
        let n: i32 = 4;
        (x_min-n, x_max+n, y_min-n, y_max+n, z_min-n, z_max+n)     
        // ((x_min).rem_euclid(b), (x_max).rem_euclid(b),
        //  (y_min).rem_euclid(b), (y_max).rem_euclid(b),
        //  (z_min).rem_euclid(b), (z_max).rem_euclid(b))
    }

    fn children_in_range(&self, xmin: i32, xmax: i32, ymin: i32, ymax:i32, zmin: i32, zmax:i32, b:i32) -> Vec<usize>{
        let mut neighbors : Vec<usize> = Vec::new();
        for child in &self.children{
            let tem: i32 = child.id/b;
            let x_id: i32 = child.id%b;
            let y_id: i32 = tem%b;
            let z_id: i32 = tem/b;
            let index: i32 = index_range(xmin, xmax, ymin, ymax, zmin, zmax);
            println!("{} {} {} {} {} {} {}", index, xmin, xmax, ymin, ymax, zmin, zmax);
            if index == 0 { 
                if (x_id >= xmin) && (x_id <= xmax) {
                    if (y_id >= ymin) && (y_id <= ymax) {
                        if (z_id >= zmin) && (z_id <= zmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 1 { 
                if (x_id >= xmin) && (x_id <= xmax) {
                    if (y_id >= ymin) && (y_id <= ymax) {
                        if ((zmin <= z_id)&&(z_id <= b)) || (z_id <= zmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 2 {
                if (x_id >= xmin) && (x_id <= xmax) {
                    if (z_id >= zmin) && (z_id <= zmax) {
                        if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 3 {
                if (x_id >= xmin) && (x_id <= xmax) {
                    if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                        if ((zmin <= z_id)&&(z_id <= b)) || (z_id <= zmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 4 { 
                if (y_id >= ymin) && (y_id <= ymax) {
                    if (z_id >= zmin) && (z_id <= zmax) {
                        if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 5 { 
                if (y_id >= ymin) && (y_id <= ymax) {
                    if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                        if ((zmin <= z_id)&&(z_id <= b)) || (z_id <= zmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 6 {
                if (z_id >= zmin) && (z_id <= zmax) {
                    if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                        if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else if index == 7 {
                if ((xmin <= x_id)&&(x_id <= b)) || (x_id <= xmax) {
                    if ((ymin <= y_id)&&(y_id <= b)) || (y_id <= ymax) {
                        if ((zmin <= z_id)&&(z_id <= b)) || (z_id <= zmax) {
                            neighbors.push(child.id as usize);
                        }
                    }
                }
            } else {
                println!("ERROR: Finding Index in Neighbors' Finder ")
            }
        }
        neighbors
    }

    fn find_neighbors(& self, p: usize, k: f64, s: i32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, wd: f64, lg:f64, hg:f64, h: f64, rkern: f64) {
        let b: i32 = (self.branches as f64).powf(1./k) as i32;
        let (x_min, x_max, y_min, y_max, z_min, z_max) = self.range_neigh(particles[p].x, particles[p].y, particles[p].z, h, b as i32, rkern);
        let neighbors = self.children_in_range(x_min, x_max, y_min, y_max, z_min, z_max, b);
        for ii in neighbors {
            if self.children[ii].n <= s {
                for q in &self.children[ii].particles {
                    if periodic_norm(particles[p].x, particles[*q].x, particles[p].y, particles[*q].y, particles[p].z, particles[*q].z, wd, lg, hg, rkern*h) <= rkern*rkern*h*h {
                        neighbors_of_p.push(*q);
                    }
                }
            } else {
                self.children[ii].find_neighbors(p, k, s, particles, neighbors_of_p, wd, lg, hg, h, rkern);
            }
        }
    }
}

fn index_range(xmin: i32, xmax: i32, ymin: i32, ymax: i32, zmin: i32, zmax: i32) -> i32 {
    if xmin <= xmax {
        if ymin <= ymax {
            if zmin <= zmax {
                return 0;
            } else {
                return 1;
            }
        } else {
            if zmin < zmax {
                return 2;
            } else {
                return 3;
            }
        }
    } else {
        if ymin < ymax {
            if zmin < zmax {
                return 4;
            } else {
                return 5;
            }
        } else {
            if zmin < zmax {
                return 6;
            } else {
                return 7;
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
    wtr.write_record(&["x_min", "y_min", "z_min", "sidex", "sidey", "sidez", "depth", "n"]).expect("Couldn't write data");
    save_child(&mut wtr, tree);
}

fn save_child<W: io::Write>(wtr: &mut Writer<W>, tree: & Node){
    wtr.write_record(&[tree.xmin.to_string(), tree.ymin.to_string(), tree.zmin.to_string(), tree.sidex.to_string(), tree.sidey.to_string(), tree.sidez.to_string(), tree.depth.to_string(), tree.n.to_string()]).expect("Couldn't write data");
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
pub fn periodic_norm(x1: f64, x2: f64, y1: f64, y2: f64, z1: f64, z2: f64, wd: f64, lg: f64, hg: f64, eps: f64) -> f64 {
    
    let mut x_temp: f64 = x1 - x2;
    let mut y_temp: f64 = y1 - y2;
    let mut z_temp: f64 = z1 - z2;

    if x_temp.abs() > wd-2.*eps {
        if x_temp > 0. {
            x_temp -= wd;
        } else {
            x_temp += wd;
        }
    }
    if y_temp.abs() > lg-2.*eps {
        if y_temp > 0. {
            y_temp -= lg;
        } else {
            y_temp += lg;
        }
    }
    if z_temp.abs() > hg-2.*eps {
        if z_temp > 0. {
            z_temp -= hg;
        } else {
            z_temp += hg;
        }
    }
    return x_temp*x_temp + y_temp*y_temp + z_temp*z_temp;
}