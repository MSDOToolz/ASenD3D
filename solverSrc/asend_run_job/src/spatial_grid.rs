
use std::collections::LinkedList;
use std::collections::BTreeMap;
use crate::list_ent::DualInt;

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
    pub first : BTreeMap<usize, usize>,
    pub data : Vec<DualInt>,
    pub next_avail : usize,
}

impl SpatialGrid {
    pub fn new() -> SpatialGrid {
        SpatialGrid {
            x_min : 0.0,
            x_sp : 0.0,
            x_bins : 2000000,
            y_min : 0.0,
            y_sp : 0.0,
            y_bins : 2000000,
            z_min : 0.0,
            z_sp : 0.0,
            z_bins : 2000000,
            first : BTreeMap::new(),
            data : Vec::new(),
            next_avail : 0,
        }
    }
}


