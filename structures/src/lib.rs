#[derive(Debug)]
pub struct Particle {
	pub x: f64,
    pub y: f64,
    pub h: f64,
    pub dh: f64,
    pub rho: f64,
    pub vx: f64,
    pub vy: f64,
    pub ax: f64,
    pub ay: f64,
    pub u: f64,
    pub du: f64,
}

impl Default for Particle {
    fn default() -> Particle {
        Particle {
            x: 0.,
            y: 0.,
            h: 0.1,
            dh: 0.0,
            rho: 1.0,
            vx: 0.0,
            vy: 0.0,
            ax: 0.0,
            ay: 0.0,
            u: 0.0,
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