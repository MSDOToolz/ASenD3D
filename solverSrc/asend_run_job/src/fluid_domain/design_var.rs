use crate::constants::MAX_INT;
use crate::diff_doub::*;
use crate::cpp_str::CppStr;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct DesignVariable {
    pub category : CppStr,
    pub component : usize,
    pub active_time : [f64; 2],
    pub cell_set_name : CppStr,
    pub cell_set_ptr : usize,
    pub nd_set_name : CppStr,
    pub nd_set_ptr : usize,
    pub coefs : LinkedList<f64>,
    pub val : DiffDoub0,
    pub diff_val : DiffDoub1,
    pub comp_cell_list : LinkedList<usize>,
}

impl DesignVariable {
    pub fn new() -> DesignVariable {
        DesignVariable {
            category : CppStr::new(),
            component : 0usize,
            active_time : [0f64, 1.0e+100f64],
            cell_set_name : CppStr::new(),
            cell_set_ptr : MAX_INT,
            nd_set_name : CppStr::new(),
            nd_set_ptr : MAX_INT,
            coefs : LinkedList::new(),
            value : DiffDoub0::new(),
            diff_val : DiffDoub1::new(),
            comp_cell_list : LinkedList::new(),
        }
    }
}