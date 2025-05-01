use crate::constants::MAX_INT;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Face {
    pub num_nds : usize,
    pub loc_nodes : [usize; 6],
    pub glob_nodes : [usize; 6],
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