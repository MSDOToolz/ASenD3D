use crate::constants::MAX_INT;
use crate::diff_doub::*;

#[derive(Clone)]
pub struct Face {
    pub glob_nodes : [usize; 3],
    pub loc_nodes : [usize; 3],
    pub on_surf : bool,
    pub twin_id : usize,
    pub host_cell : usize,
}

impl Face {
    pub fn new() -> Face {
        Face {
            glob_nodes : [MAX_INT; 3],
            loc_nodes : [MAX_INT; 3],
            on_surf : true,
            twin_id : MAX_INT,
            host_cell : MAX_INT,
        }
    }
}

#[derive(Clone)]
pub struct CellData {
    pub den : DiffDoub1,
    pub vel : [DiffDoub1; 3],
    pub v_rel : [DiffDoub1; 3],
    pub temp : DiffDoub1,
    pub turb : DiffDoub1,
    pub v_grad : [DiffDoub1; 9],
    pub t_grad : [DiffDoub1; 3],
}

impl CellData {

    pub fn new() -> CellData {
        CellData {
            den : DiffDoub1::new(),
            vel : [DiffDoub1::new(); 3],
            v_rel : [DiffDoub1::new(); 3],
            temp : DiffDoub1::new(),
            turb : DiffDoub1::new(),
            v_grad : [DiffDoub1::new(); 9],
            t_grad : [DiffDoub1::new(); 3],
        }
    }
}

#[derive(Clone)]
pub struct FaceData {
    pub area : DiffDoub1,
    pub normal : [DiffDoub1; 3],
    pub den : DiffDoub1,
    pub vel : [DiffDoub1; 3],
    pub v_rel : [DiffDoub1; 3],
    pub temp : DiffDoub1,
    pub turb : DiffDoub1,
    pub v_grad : [DiffDoub1; 9],
    pub t_grad : [DiffDoub1; 3],
}

impl FaceData {
    pub fn new() -> FaceData {
        FaceData {
            area : DiffDoub1::new(),
            normal : [DiffDoub1::new(); 3],
            den : DiffDoub1::new(),
            vel : [DiffDoub1::new(); 3],
            v_rel : [DiffDoub1::new(); 3],
            temp : DiffDoub1::new(),
            turb : DiffDoub1::new(),
            v_grad : [DiffDoub1::new(); 9],
            t_grad : [DiffDoub1::new(); 3],
        }
    }
}

pub mod face_meth;