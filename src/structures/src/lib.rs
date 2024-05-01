// ------------------------------------------------------------------------- //
// SPH STRUCTURES                                                            //
// ------------------------------------------------------------------------- //


// ------------------------------------------------------------------------- //
// Particle: SPH particle                                                    //
//  ptype:      u8  - 0: gas; 1: boundary; 2: dead                           //
//  rho:        f64 - density                                                //
//  x, y, z:    f64 - position                                               //
//  vx, vy, vz: f64 - velocity                                               //
//  v_xyz_star: f64 - auxiliar velocity for time integrator                  //
//  divv:       f64 - divergence of velocity div(v)                          //
//  ax, ay, az: f64 - acceleration                                           //
//  u:          f64 - internal energy                                        //
//  u_star:     f64 - auxiliar internal energy for time integrator           //
//  du:         f64 - internal energy change                                 //
// ------------------------------------------------------------------------- //
#[derive(Debug)]
pub struct Particle {
    pub ptype: u8,
    pub rho: f64,
	pub x: f64,
    pub y: f64,
    pub z: f64,
    pub h: f64,
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
            rho: 1.0,
            x: 0.,
            y: 0.,
            z: 0.,
            h: 0.1,
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

#[derive(Copy, Clone)]
pub struct Pointer(pub *mut Particle);
unsafe impl Send for Pointer {}
unsafe impl Sync for Pointer {}

// ------------------------------------------------------------------------- //
// Node: Tree node                                                           //
//  xmin, ymin, zmin:       f64 - bottom left corner                         //
//  sidex, sidey, sidez:    f64 - sides                                      //
//  id:                     i32 - id                                         //
//  depth:                  i32 - generation depth                           //
//  n:                      i32 - number of particles                        //
//  branches:               i32 - number of children                         //
//  chidren:                Vec - array of children                          //
// particles:               Vec - array of particles' id's                   //
// ------------------------------------------------------------------------- //
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

// ------------------------------------------------------------------------- //
// Star                                                                      //
//  m:          f64 - mass                                                   //
//  x, y, z:    f64 - position                                               //
//  vx, vy, vz: f64 - velocity                                               //
//  ax, ay, az: f64 - acceleration                                           //
// ------------------------------------------------------------------------- //
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