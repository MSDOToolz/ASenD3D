use crate::diff_doub::*;

#[derive(Clone)]
pub struct FltScr {
    pub dat : Vec<f64>,
}

impl FltScr {
    pub fn new() -> FltScr {
        FltScr {
            dat : vec![0.0f64; 3600],
        }
    }
}

#[derive(Clone)]
pub struct DiffDoub0Scr {
    pub dat : Vec<DiffDoub0>,
}

impl DiffDoub0Scr {
    pub fn new() -> DiffDoub0Scr {
        DiffDoub0Scr {
            dat : vec![DiffDoub0::new(); 60],
        }
    }
}

#[derive(Clone)]
pub struct DiffDoub1Scr {
    pub dat : Vec<DiffDoub1>,
}

impl DiffDoub1Scr {
    pub fn new() -> DiffDoub1Scr {
        DiffDoub1Scr {
            dat : vec![DiffDoub1::new(); 60],
        }
    }
}