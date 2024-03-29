#[derive(Debug)]
pub struct Particle {
    pub ptype: u8,
	pub x: f64,
    pub y: f64,
    pub z: f64,
    pub h: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub vx_star: f64,
    pub vy_star: f64,
    pub vz_star: f64,
    pub divv: f64,
    pub ax: f64,
    pub ay: f64,
    pub az: f64,
    pub u: f64,
    pub u_star: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            ptype: 0,
            x: 0.,
            y: 0.,
            z: 0.,
            h: 0.1,
            rho: 1.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            vx_star: 0.0,
            vy_star: 0.0,
            vz_star: 0.0,
            divv: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
            u: 1.0,
            u_star: 0.0,
            du: 0.0,
        }
    }
}

pub struct Node {
    pub xmin: f64,
    pub ymin: f64,
    pub zmin: f64,
    pub sidex: f64,
    pub sidey: f64,
    pub sidez: f64,
    pub id: i32,
    pub depth: i32,
    pub n: i32,
    pub branches: i32,
    pub children: Vec<Node>,
    pub particles: Vec<usize>,
}

#[derive(Copy, Clone)]
pub struct Pointer(pub *mut Particle);
unsafe impl Send for Pointer {}
unsafe impl Sync for Pointer {}

pub struct Star {
    pub m: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub ax: f64,
    pub ay: f64,
    pub az: f64,
}

impl Default for Star {
    fn default() -> Star {
        Star {
            m: 0.,
            x: 0.,
            y: 0.,
            z: 0.,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
        }
    }
}