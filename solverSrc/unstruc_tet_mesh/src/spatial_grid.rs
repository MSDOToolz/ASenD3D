use crate::mesh_node::*;
use crate::mesh_element::*;
use crate::mesh_face::*;
use crate::cpp_str::CppStr;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct IntList {
    pub i_lst : LinkedList<usize>,
}

impl IntList {
    pub fn new() -> IntList {
        IntList {
            i_lst : LinkedList::new(),
        }
    }
}

#[derive(Clone)]
pub struct SpatialGrid {
    pub x_min : f64,
    pub x_sp : f64,
    pub x_bins : usize,
    pub y_min : f64,
    pub y_sp : f64,
    pub y_bins : usize,
    pub z_min : f64,
    pub z_sp : f64,
    pub z_bins : usize,
    pub list_ar : Vec<IntList>,
}

impl SpatialGrid {
    pub fn new() -> SpatialGrid {
        SpatialGrid {
            x_min : 0.0,
            x_sp : 0.0,
            x_bins : 0,
            y_min : 0.0,
            y_sp : 0.0,
            y_bins : 0,
            z_min : 0.0,
            z_sp : 0.0,
            z_bins : 0,
            list_ar : Vec::new(),
        }
    }
}


