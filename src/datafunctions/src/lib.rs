use std::{
    io::{BufRead, BufReader, prelude::*},
    fs::File,
    error::Error,
    f64,
};

use csv::Writer;

use structures::Particle;

pub fn save_data(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = Writer::from_path(path)?;
    wtr.write_record(&["ptype", "x", "y", "z", "vx", "vy", "vz", "h", "u"])?;
    for particle in particles {
        wtr.write_record(&[particle.ptype.to_string(), particle.x.to_string(), particle.y.to_string(), particle.z.to_string(),
                           particle.vx.to_string(), particle.vy.to_string(), particle.vz.to_string(),
                           particle.h.to_string(), particle.u.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn save_data_bin(path: &str, particles: & Vec<Particle>)-> Result<(), Box<dyn Error>>{
    let mut wtr = File::create(path)?;
    for particle in particles {
        wtr.write_all(&particle.x.to_le_bytes())?;
        wtr.write_all(&particle.y.to_le_bytes())?;
        wtr.write_all(&particle.z.to_le_bytes())?;
        wtr.write_all(&particle.vx.to_le_bytes())?;
        wtr.write_all(&particle.vy.to_le_bytes())?;
        wtr.write_all(&particle.vz.to_le_bytes())?;
        wtr.write_all(&particle.h.to_le_bytes())?;
        wtr.write_all(&particle.u.to_le_bytes())?;
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
        particles.push(Particle{ptype:(&record[0]).parse::<u8>().unwrap(), x:(&record[1]).parse::<f64>().unwrap(), y:(&record[2]).parse::<f64>().unwrap(), z:(&record[3]).parse::<f64>().unwrap(),
                                vx:(&record[4]).parse::<f64>().unwrap(), vy:(&record[5]).parse::<f64>().unwrap(), vz:(&record[6]).parse::<f64>().unwrap(),
                                h:(&record[7]).parse::<f64>().unwrap(), u:(&record[8]).parse::<f64>().unwrap(),
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