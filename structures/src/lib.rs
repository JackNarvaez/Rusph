#[derive(Debug)]
pub struct Particle {
	pub x: f64,
    pub y: f64,
    pub h: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    pub vx_star: f64,
    pub vy_star: f64,
    pub divv: f64,
    pub ax: f64,
    pub ay: f64,
    pub u: f64,
    pub u_star: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            x: 0.,
            y: 0.,
            h: 0.1,
            rho: 1.0,
            vx: 0.0,
            vy_star: 0.0,
            vx_star: 0.0,
            vy: 0.0,
            divv: 0.0,
            ax: 0.0,
            ay: 0.0,
            u: 1.0,
            u_star: 0.0,
            du: 0.0,
        }
    }
}

pub struct Node {
    pub xmin: f64,
    pub ymin: f64,
    pub side: f64,
    pub id: u32,
    pub depth: u32,
    pub n: u32,
    pub branches: u32,
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
            children: Vec::new(),
            particles: Vec::new(),
        }
    }
}

#[derive(Copy, Clone)]
pub struct Pointer(pub *mut Particle);
unsafe impl Send for Pointer {}
unsafe impl Sync for Pointer {}

pub struct Star {
    pub m: f64,
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub ax: f64,
    pub ay: f64,
}

impl Default for Star {
    fn default() -> Star {
        Star {
            m: 0.,
            x: 0.,
            y: 0.,
            vx: 0.0,
            vy: 0.0,
            ax: 0.0,
            ay: 0.0,
        }
    }
}