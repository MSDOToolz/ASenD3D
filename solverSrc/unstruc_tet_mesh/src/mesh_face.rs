use crate::constants::*;


#[derive(Clone)]
pub struct MeshFace {
    pub nodes : [usize; 3],
    pub elements : [usize; 2],
    pub norm_dir : [f64; 3],
    pub proj_dist : f64,
}

impl MeshFace {
    pub fn new() -> MeshFace {
        MeshFace {
            nodes : [0; 3],
            elements : [MAX_INT; 2],
            norm_dir : [0.0; 3],
            proj_dist : 0.0,
        } 
    }
}


