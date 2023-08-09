use std::{
    io::{BufRead, BufReader},
    fs::File,
    error::Error,
    f64,
};

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;

const SEED: u64 = 123;

use csv::Writer;

use structures::Particle;

// -------- Write data --------

pub fn init_random_square(path: &str, n: u32, h:f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    let mut rng = Pcg64::seed_from_u64(SEED);
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "rho"])?;
    for _ii in 0..n{
        let x: f64 = wd*rng.gen::<f64>();
        let y: f64 = lg*rng.gen::<f64>();
        let z: f64 = hg*rng.gen::<f64>();
            wtr.write_record(&[(x0 + x).to_string(), (y0 + y).to_string(), (z0+z).to_string(),
                               String::from("0.0"), String::from("0.0"), String::from("0.0"),
                               h.to_string(), String::from("0.0")])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn init_square(path: &str, n: u32, h:f64, wd:f64, lg:f64, hg: f64, x0: f64, y0: f64, z0: f64)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "rho"])?;
    let dl: f64 = ((wd*lg*hg)/n as f64).powf(1./3.);
    let nx: u32 = (wd/dl).ceil() as u32;  
    let ny: u32 = (lg/dl).ceil() as u32;  
    let nz: u32 = (hg/dl).ceil() as u32;  
    for kk in 0..nz{
        for jj in 0..ny{
            for ii in 0..nx{
                let x: f64 = x0 + ii as f64*dl;
                let y: f64 = y0 + jj as f64*dl;
                let z: f64 = z0 + kk as f64*dl;
                wtr.write_record(&[(x0 + x).to_string(), (y0 + y).to_string(), (z0+z).to_string(),
                               String::from("0.0"), String::from("0.0"), String::from("0.0"),
                               h.to_string(), String::from("0.0")])?;
            }
        }
    }
    wtr.flush()?;
    Ok(())
}

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h", "u"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].x.to_string(), particles[ii].y.to_string(), particles[ii].z.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(), particles[ii].vz.to_string(),
                           particles[ii].h.to_string(), particles[ii].u.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn save_data_iso(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["x", "y", "z", "vx", "vy", "vz", "h"])?;
    for ii in 0..particles.len() {
        wtr.write_record(&[particles[ii].x.to_string(), particles[ii].y.to_string(), particles[ii].z.to_string(),
                           particles[ii].vx.to_string(), particles[ii].vy.to_string(), particles[ii].vz.to_string(),
                           particles[ii].h.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

// -------- Read data --------

pub fn read_data(path: &str, particles: &mut Vec<Particle>) -> Result<(), Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;
    for result in rdr.records() {
        let record = result?;
        particles.push(Particle{x:(&record[0]).parse::<f64>().unwrap(), y:(&record[1]).parse::<f64>().unwrap(), z:(&record[2]).parse::<f64>().unwrap(),
                                vx:(&record[3]).parse::<f64>().unwrap(), vy:(&record[4]).parse::<f64>().unwrap(), vz:(&record[5]).parse::<f64>().unwrap(),
                                h:(&record[6]).parse::<f64>().unwrap(), u:(&record[7]).parse::<f64>().unwrap(),
                                ..Default::default()});
    }
    Ok(())
}

pub fn read_data_iso(path: &str, particles: &mut Vec<Particle>) -> Result<(), Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;
    for result in rdr.records() {
        let record = result?;
        particles.push(Particle{x:(&record[0]).parse::<f64>().unwrap(), y:(&record[1]).parse::<f64>().unwrap(), z:(&record[2]).parse::<f64>().unwrap(),
                                vx:(&record[3]).parse::<f64>().unwrap(), vy:(&record[4]).parse::<f64>().unwrap(), vz:(&record[5]).parse::<f64>().unwrap(),
                                h:(&record[6]).parse::<f64>().unwrap(),
                                ..Default::default()});
    }
    Ok(())
}

pub fn read_input(path: &str) -> Vec<f64> {
    let reader = BufReader::new(File::open(path).expect("Cannot open input file"));
    let mut input = Vec::new();
    for line in reader.lines() {
        for word in line.unwrap().split_whitespace() {
            if word == "#" {
                break;
            } else {
                input.push(word.parse().unwrap());
            }
        }
    }
    return input;
}