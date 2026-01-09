use crate::cpp_str::*;
use crate::list_ent::*;
use crate::constants::*;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Load {
    pub this_type : CppStr,
    pub cell_set : CppStr,
    pub set_pt : usize,
    pub load : LinkedList<QuadFloat>,
    pub center : [f64; 3],
    pub axis : [f64; 3],
    pub angular_vel : f64,
    pub active_time : [f64; 2],
}

impl Load {
    pub fn new() -> Load {
        Load {
            this_type : CppStr::new(),
            cell_set : CppStr::new(),
            set_pt : MAX_INT,
            load : LinkedList::new(),
            center : [0f64; 3],
            axis : [0f64; 3],
            angular_vel : 0f64,
            active_time : [0f64, 1.0e+100],
        }
    }
}

mod load_meth;
