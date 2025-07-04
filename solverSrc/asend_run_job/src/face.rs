use crate::constants::MAX_INT;
use crate::diff_doub::*;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Face {
    pub num_nds : usize,
    pub loc_nodes : [usize; 6],
    pub glob_nodes : [usize; 6],
    pub int_pts : Vec<f64>,
    pub num_ip : usize,
    pub on_surf : bool,
    pub twin_id : usize,
    pub host_el : usize,
}

impl Face {
    pub fn new() -> Face {
        Face {
            num_nds : 0usize,
            loc_nodes : [MAX_INT; 6],
            glob_nodes : [MAX_INT; 6],
            int_pts : Vec::new(),
            num_ip : 0usize,
            on_surf : true,
            twin_id : MAX_INT,
            host_el : MAX_INT,
        }
    }
}

#[derive(Clone)]
pub struct FacePtList {
    pub fc_list : LinkedList<usize>,
}

impl FacePtList {
    pub fn new() -> FacePtList {
        FacePtList {
            fc_list : LinkedList::new(),
        }
    }
}

#[derive(Clone)]
pub struct CentData {
    pub v_grad_dfd0 : [DiffDoub0; 9],
    pub v_grad_dfd1 : [DiffDoub1; 9],
    pub t_grad_dfd0 : [DiffDoub0; 3],
    pub t_grad_dfd1 : [DiffDoub1; 3],
}

impl CentData {
    pub fn new() -> CentData {
        CentData {
            v_grad_dfd0 : [DiffDoub0::new(); 9],
            v_grad_dfd1 : [DiffDoub1::new(); 9],
            t_grad_dfd0 : [DiffDoub0::new(); 3],
            t_grad_dfd1 : [DiffDoub1::new(); 3],
        }
    }
}