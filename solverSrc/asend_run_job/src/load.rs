use crate::constants::MAX_INT;
use crate::nd_el_set::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


#[derive(Clone)]
pub struct Load {
    pub this_type : CppStr,
    pub active_time : [f64; 2],
    pub node_set : CppStr,
    pub nd_set_ptr : usize,
    pub element_set : CppStr,
    pub el_set_ptr : usize,
    pub load : [f64; 6],
    pub normal_dir : [f64; 3],
    pub norm_tol : f64,
    pub center : [f64; 3],
    pub axis : [f64; 3],
    pub angular_vel : f64,
}

impl Load {
    pub fn new() -> Load {
        Load {
            this_type : CppStr::from("none"),
            active_time : [0f64, 1.0e+100f64],
            node_set : CppStr::from("none"),
            nd_set_ptr : MAX_INT,
            element_set : CppStr::from("none"),
            el_set_ptr : MAX_INT,
            load : [0f64; 6],
            normal_dir : [0f64; 3],
            norm_tol : 0f64,
            center : [0f64; 3],
            axis : [0f64; 3],
            angular_vel : 0f64,
        }
    }
}
