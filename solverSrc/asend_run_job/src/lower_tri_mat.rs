use crate::list_ent::*;
use crate::constraint::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


#[derive(Clone)]
pub struct LowerTriMat {
    pub mat : Vec<f64>,
    pub range : Vec<usize>,
    pub min_col : Vec<usize>,
    pub z_vec : Vec<f64>,
    pub dim : usize,
    pub size : usize,
    pub max_bandwidth : usize,
    pub allocated : bool,
}

impl LowerTriMat {
    pub fn new() -> LowerTriMat {
        LowerTriMat {
            mat : Vec::new(),
            range : Vec::new(),
            min_col : Vec::new(),
            z_vec : Vec::new(),
            dim : 0usize,
            size : 0usize,
            max_bandwidth : 0usize,
            allocated : false,
        }
    }
}

