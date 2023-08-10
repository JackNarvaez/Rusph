use std::{
    io,
    f64,
    i32,
};

use num_integer::Roots;

use csv::Writer;

use structures::{
    Particle,
    Node,
};

use rayon::prelude::*;

pub trait BuildTree {

    fn new(n_p: i32, x0: f64, y0: f64, z0: f64, wd: f64, lg: f64, hg: f64) -> Node;

    fn branching_factor(& self, s:f64) -> i32;

    fn create_child(&self, j: i32, b: i32, dx: f64, dy: f64, dz: f64) -> Node;

    fn create_sub_cells(&mut self, b: i32);

    fn add_particle(&mut self, i: usize);

    fn delete_sub_cells(&mut self);

    fn delete_particles(&mut self);

    fn distribution_ratio(&self, limit: i32) -> f64;

    fn build_tree(&mut self, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64);

    fn build_octtree(&mut self, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64);

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
             id: 0,
             depth: 0,
             branches: 0,
             children: Vec::new()}
    }
    
    fn branching_factor(& self, s:f64) -> i32 {
        let b:i32 = ((self.n as f64 /s).cbrt()).ceil() as i32;
        if b==0 {
            println!("ERROOOOOOORRR, {}", (self.n as f64 /s).cbrt());
        }
        b
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

    fn add_particle(&mut self, i: usize) {
        self.particles.push(i);
        self.n += 1;
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

    fn build_tree(&mut self, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64) {
        let mut redistribution :bool = true;
        let mut b: i32 = self.branching_factor(s as f64);
        while redistribution {
            if b==0 {
                println!("ERROOOOOOORRR, {}", b);
            }
            self.branches = b*b*b;
            self.create_sub_cells(b);
            for p in &self.particles {
                let mut x_p:i32 = ((particles[*p].x - self.xmin) / self.sidex * b as f64).floor() as i32;
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
                let j :usize = (x_p + (y_p  + z_p * b) * b) as usize;
                self.children[j].add_particle(*p);
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
            if b==0 {
                println!("ERROOOOOOORRR, {}", b);
            }
        }
        (self.children).par_iter_mut().for_each(|child| {
            if (child.n > s) && (child.sidex > smallest_cell) && (child.sidey > smallest_cell) && (child.sidez > smallest_cell) {
                child.build_tree(s, alpha, beta, particles, smallest_cell);
            }
        });
    }

    fn build_octtree(&mut self, s: i32, alpha: f64, beta: f64, particles: & Vec<Particle>, smallest_cell: f64) {
        let b:i32 = 2;
        self.branches = b*b*b;
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
            let j :usize = (x_p + (y_p + z_p * b) * b) as usize;
            self.children[j].add_particle(*p);
        }
        self.delete_particles();
        (self.children).par_iter_mut().for_each(|child| {
            if (child.n > s) && (child.sidex > smallest_cell) && (child.sidey > smallest_cell) && (child.sidez > smallest_cell) {
                child.build_tree(s, alpha, beta, particles, smallest_cell);
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
    fn range_neigh(&self, x_p: f64, y_p: f64, z_p: f64, b: i32, hrkern: f64, x0: f64, y0: f64, z0: f64, wd: f64, lg: f64, hg: f64,) -> Vec<usize>;

    fn print_particles(&self);

    fn find_neighbors(& self, p: usize, s: i32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, wd: f64, lg: f64, hg: f64, x0:f64, y0: f64, z0:f64, hrkern: f64);
}

impl FindNeighbors for Node {

    fn range_neigh(&self, x_p: f64, y_p: f64, z_p: f64, b: i32, hrkern: f64, x0: f64, y0: f64, z0: f64, wd: f64, lg: f64, hg: f64,) -> Vec<usize>{
        let factorx : f64 =  b as f64 /self.sidex;
        let factory : f64 =  b as f64 /self.sidey;
        let factorz : f64 =  b as f64 /self.sidez;
        
        let (xlow, xup) =  limits(x_p-hrkern, x_p+hrkern, x0, wd, self.xmin, self.depth);
        let (ylow, yup) =  limits(y_p-hrkern, y_p+hrkern, y0, lg, self.ymin, self.depth);
        let (zlow, zup) =  limits(z_p-hrkern, z_p+hrkern,  z0, hg, self.zmin, self.depth);

        let mut x_min: i32 = ((xlow - self.xmin) * factorx).floor() as i32;
        let mut x_max: i32 = ((xup - self.xmin) * factorx).floor() as i32;
        let mut y_min: i32 = ((ylow - self.ymin) * factory).floor() as i32;
        let mut y_max: i32 = ((yup - self.ymin) * factory).floor() as i32;
        let mut z_min: i32 = ((zlow - self.zmin) * factorz).floor() as i32;
        let mut z_max: i32 = ((zup - self.zmin) * factorz).floor() as i32;

        let mut neighbors : Vec<usize> = Vec::new();

        if self.depth != 0 {
            let up: i32 = b-1;
            set_limits(&mut x_min, 0, up);
            set_limits(&mut x_max, 0, up);
            set_limits(&mut y_min, 0, up);
            set_limits(&mut y_max, 0, up);
            set_limits(&mut z_min, 0, up);
            set_limits(&mut z_max, 0, up);
            
            for kk in z_min..z_max+1{
                for jj in y_min..y_max+1{
                    for ii in x_min..x_max+1{
                        neighbors.push((ii + (jj + kk * b) * b) as usize);
                    }
                }
            }
            return neighbors;
        } else {
            for kk in z_min..z_max+1{
                for jj in y_min..y_max+1{
                    for ii in x_min..x_max+1{
                        neighbors.push((ii.rem_euclid(b)+ (jj.rem_euclid(b)+kk.rem_euclid(b)*b)*b) as usize);
                    }
                }
            }
            return neighbors;
        }
    }

    fn print_particles(&self) {
        if self.children.len() == 0 {
            println!("{} {}", self.depth, self.id);
            for ii in &self.particles {
                println!("p:{}.", ii);
                }
        } else {
            println!("{} {}", self.depth, self.id);
            for child in &self.children{
                child.print_particles();
            }
        }
    }

    fn find_neighbors(& self, p: usize, s: i32, particles: & Vec<Particle>, neighbors_of_p: &mut Vec<usize>, wd: f64, lg:f64, hg:f64, x0:f64, y0:f64, z0:f64, hrkern: f64) {
        let b: i32 = (self.branches).cbrt();
        let cell_neighbors = self.range_neigh(particles[p].x, particles[p].y, particles[p].z, b as i32, hrkern, x0, y0, z0, wd, lg, hg);
        for ii in cell_neighbors {
            if self.children[ii].branches == 0 {
                for q in &self.children[ii].particles {
                    let norm: f64 = sq_periodic_norm(particles[p].x, particles[*q].x, particles[p].y, particles[*q].y, particles[p].z, particles[*q].z, wd, lg, hg, hrkern);
                    if norm <= hrkern*hrkern {
                        neighbors_of_p.push(*q);
                    }
                }
            } else {
                self.children[ii].find_neighbors(p, s, particles, neighbors_of_p, wd, lg, hg, x0, y0, z0, hrkern);
            }
        }
    }

}


// Save information about tree and neighbours
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
pub fn sq_periodic_norm(x1: f64, x2: f64, y1: f64, y2: f64, z1: f64, z2: f64, wd: f64, lg: f64, hg: f64, eps: f64) -> f64 {
    
    let mut x_temp: f64 = x1 - x2;
    let mut y_temp: f64 = y1 - y2;
    let mut z_temp: f64 = z1 - z2;

    let diam: f64 = 2.*eps;

    if x_temp.abs() > wd-diam {
        if x_temp > 0. {
            x_temp -= wd;
        } else {
            x_temp += wd;
        }
    }
    if y_temp.abs() > lg-diam {
        if y_temp > 0. {
            y_temp -= lg;
        } else {
            y_temp += lg;
        }
    }
    if z_temp.abs() > hg-diam {
        if z_temp > 0. {
            z_temp -= hg;
        } else {
            z_temp += hg;
        }
    }
    return x_temp*x_temp + y_temp*y_temp + z_temp*z_temp;
}

fn limits(low_lim: f64, up_lim: f64, l0: f64, l: f64, lmin: f64, depth: i32) -> (f64, f64) {
    if depth == 0 {
        return  (low_lim, up_lim);
    } else {
        if low_lim < l0 {
            if lmin < (l0 + 0.5 * l) {
                return (0., up_lim);
            } else {
                return (low_lim + l, l0 + l);
            }
        } else if up_lim > (l0 + l) {
            if lmin < (l0 + 0.5 * l) {
                return (0., up_lim-l);
            } else {
                return (low_lim, l0 + l);
            }
        } else {
            return (low_lim, up_lim);
        }
    }
}

fn set_limits(x: &mut i32, low: i32, up: i32) {
    if *x < low {
        *x = low;
    } else if *x > up {
        *x = up;
    }
}