
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
pub struct DualInt {
    pub i1 : usize,
    pub i2 : usize,
}

impl DualInt {
    pub fn new() -> DualInt {
        DualInt {
            i1 : 0,
            i2 : 0, 
        }
    }
}

#[derive(Clone)]
pub struct DualFloat {
    pub f1 : f64,
    pub f2 : f64,
}

impl DualFloat {
    pub fn new() -> DualFloat {
        DualFloat {
            f1 : 0.0,
            f2 : 0.0,
        }
    }
}

#[derive(Clone)]
pub struct QuadFloat {
    pub f1 : f64,
    pub f2 : f64,
    pub f3 : f64,
    pub f4 : f64,
}

impl QuadFloat {
    pub fn new() -> QuadFloat {
        QuadFloat {
            f1 : 0f64,
            f2 : 0f64,
            f3 : 0f64,
            f4 : 0f64,
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

pub mod list_ent_meth;