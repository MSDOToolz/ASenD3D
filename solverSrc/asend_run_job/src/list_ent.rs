use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct IDCapsule {
    pub int_dat : usize,
    pub doub_dat : f64,
}

impl IDCapsule {
    pub fn new() -> IDCapsule {
        IDCapsule {
            int_dat : 0usize,
            doub_dat : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct MatrixEnt {
    pub row : usize,
    pub col : usize,
    pub value : f64,
}

impl MatrixEnt {
    pub fn new() -> MatrixEnt {
        MatrixEnt {
            row : 0usize,
            col : 0usize,
            value : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct MatrixRow {
    pub row_vec : LinkedList<MatrixEnt>,
}

impl MatrixRow {
    pub fn new() -> MatrixRow {
        MatrixRow {
            row_vec : LinkedList::new(),
        }
    }
}

#[derive(Clone)]
pub struct SparseMat {
    pub dim : usize,
    pub matrix : Vec<MatrixRow>,
}

impl SparseMat {
    pub fn new() -> SparseMat {
        SparseMat {
            dim : 0usize,
            matrix : Vec::new(),
        }
    }
}
