use crate::constants::MAX_INT;
use crate::diff_doub::*;

#[derive(Clone)]
pub struct FluidFace {
    pub glob_nodes : [usize; 3],
    pub area : DiffDoub1,
    pub normal : [DiffDoub1; 3],
    pub on_surf : bool,
    pub twin_id : usize,
    pub host_el : usize,
}

impl FluidFace {
    pub fn new() -> FluidFace {
        FluidFace {
            glob_nodes : [MAX_INT; 3],
            area : DiffDoub1::new(),
            normal : [DiffDoub1::new(); 3],
            on_surf : true,
            twin_id : MAX_INT,
            host_el : MAX_INT,
        }
    }
}

#[derive(Clone)]
pub struct CellData {
    pub v_grad : [DiffDoub1; 9],
    pub t_grad : [DiffDoub1; 3],
}

impl CellData {

    pub fn new() -> CellData {
        CellData {
            v_grad : [DiffDoub1::new(); 9],
            t_grad : [DiffDoub1::new(); 3],
        }
    }
}

pub mod fluid_face_meth;