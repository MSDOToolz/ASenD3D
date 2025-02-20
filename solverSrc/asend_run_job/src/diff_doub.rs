
#[derive(Clone)]
#[derive(Copy)]
pub struct DiffDoub0 {
    pub val : f64,
    pub dval : f64,
}

impl DiffDoub0 {
    pub fn new() -> DiffDoub0 {
        DiffDoub0 {
            val : 0f64,
            dval : 0f64,
        }
    }
}

#[derive(Clone)]
#[derive(Copy)]
pub struct DiffDoub1 {
    pub val : f64,
    pub dval : f64,
    pub tmp : f64,
    pub tmp2 : f64,
}

impl DiffDoub1 {
    pub fn new() -> DiffDoub1 {
        DiffDoub1 {
            val : 0f64,
            dval : 0f64,
            tmp : 0f64,
            tmp2 : 0f64,
        }
    }
}

#[derive(Clone)]
#[derive(Copy)]
pub struct DiffDoub2 {
    pub val : f64,
    pub dv1 : f64,
    pub dv2 : f64,
    pub dv12 : f64,
    pub tmp : f64,
    pub tmp2 : f64,
    pub tmp3 : f64,
    pub tmp4 : f64,
}

impl DiffDoub2 {
    pub fn new() -> DiffDoub2 {
        DiffDoub2 {
            val : 0f64,
            dv1 : 0f64,
            dv2 : 0f64,
            dv12 : 0f64,
            tmp : 0f64,
            tmp2 : 0f64,
            tmp3 : 0f64,
            tmp4 : 0f64,
        }
    }
}
