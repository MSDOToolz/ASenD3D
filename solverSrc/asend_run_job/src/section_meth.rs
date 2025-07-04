use crate::section::*;
use crate::matrix_functions::*;
use crate::fmath::*;

impl Material {
    pub fn set_stiffness(&mut self, row : usize, col : usize, val : f64) {
        self.stiffness[6 * row + col] = val;
        self.stiffness[6 * col + row] = val;
        return;
    }

    pub fn set_damping(&mut self, row : usize, col : usize, val : f64) {
        self.damping[6 * row + col] = val;
        self.damping[6 * col + row] = val;
        return;
    }

}

impl Section {
    pub fn set_orientation(&mut self, new_ori : &mut [f64]) {
        self.orientation[0] = new_ori[0];
        self.orientation[1] = new_ori[1];
        self.orientation[2] = new_ori[2];
        self.orientation[3] = new_ori[3];
        self.orientation[4] = new_ori[4];
        self.orientation[5] = new_ori[5];
        
        let mut mag : f64 =  self.orientation[0]*self.orientation[0] + self.orientation[1]*self.orientation[1] + self.orientation[2]*self.orientation[2];
        mag = 1.0/sqrt(mag);
        self.orientation[0] = mag*self.orientation[0];
        self.orientation[1] = mag*self.orientation[1];
        self.orientation[2] = mag*self.orientation[2];
        
        let mut cp_res : [f64; 3] = [0.0, 0.0, 0.0];
        cross_prod(&mut cp_res, & self.orientation[0..], & self.orientation[3..]);
        
        mag = cp_res[0]*cp_res[0] + cp_res[1]*cp_res[1] + cp_res[2]*cp_res[2];
        mag = 1.0/sqrt(mag);
        self.orientation[6] = mag*cp_res[0];
        self.orientation[7] = mag*cp_res[1];
        self.orientation[8] = mag*cp_res[2];
        
        cross_prod(&mut cp_res, & self.orientation[6..], & self.orientation[0..]);
        for i in 0..3 {
            self.orientation[i+3] = cp_res[i];
        }
        
        return;
    }

    pub fn set_area_moment(&mut self, new_i : &mut [f64]) {
        self.area_moment[0] = new_i[0];
        self.area_moment[1] = new_i[1];
        self.area_moment[2] = new_i[2];
        self.area_moment[3] = new_i[3];
        self.area_moment[4] = new_i[4];
        return;
    }

    pub fn set_stiffness(&mut self, row : usize, col : usize, val : f64) {
        self.stiffness[6*row+col] = val;
        self.stiffness[6*col+row] = val;
        return;
    }

    pub fn set_mass(&mut self, row : usize, col : usize, val : f64) {
        self.mass[6*row+col] = val;
        self.mass[6*col+row] = val;
        return;
    }

    pub fn set_damping(&mut self, row : usize, col : usize, val : f64) {
        self.damping[6 * row + col] = val;
        self.damping[6 * col + row] = val;
        return;
    }

    pub fn set_exp_ld(&mut self, new_exp_ld : &mut [f64]) {
        self.exp_load_coef[0] = new_exp_ld[0];
        self.exp_load_coef[1] = new_exp_ld[1];
        self.exp_load_coef[2] = new_exp_ld[2];
        self.exp_load_coef[3] = new_exp_ld[3];
        self.exp_load_coef[4] = new_exp_ld[4];
        self.exp_load_coef[5] = new_exp_ld[5];
        return;
    }

    pub fn get_layer_mat_ptr(& self, layi : usize) -> usize {
        let nl = self.layers.len();
        if layi >= nl {
            panic!("Error: material requested for layer index {} of section where highest index is {}.  Remember that layer indexing begins at 0", layi, nl-1);
        }
        let mut i1 : usize =  0;
        for lay in self.layers.iter() {
            if i1 == layi {
                return lay.mat_ptr;
            }
            i1 += 1usize;
        }
        0usize
    }

}


