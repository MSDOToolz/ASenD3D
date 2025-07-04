

#[derive(Clone)]
pub struct LUMat {
    pub l_mat : Vec<f64>,
    pub l_range : Vec<usize>,
    pub l_min_col : Vec<usize>,
    pub l_size : usize,
    pub u_mat : Vec<f64>,
    pub u_range : Vec<usize>,
    pub u_min_row : Vec<usize>,
    pub u_size : usize,
    pub z_vec : Vec<f64>,
    pub dim : usize,
    pub max_bandwidth : usize,
    pub allocated : bool,
}

impl LUMat {
    pub fn new() -> LUMat {
        LUMat {
            l_mat : Vec::new(),
            l_range : Vec::new(),
            l_min_col : Vec::new(),
            l_size : 0usize,
            u_mat : Vec::new(),
            u_range : Vec::new(),
            u_min_row : Vec::new(),
            u_size : 0usize,
            z_vec : Vec::new(),
            dim : 0usize,
            max_bandwidth : 0usize,
            allocated : false,
        }
    }
}
