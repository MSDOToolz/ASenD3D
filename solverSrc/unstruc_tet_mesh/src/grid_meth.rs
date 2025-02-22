use crate::spatial_grid::*;
use crate::cpp_str::CppStr;


impl IntList {
    pub fn copy_to_vector(&mut self, in_vec : &mut Vec<usize>, st_i : usize, max_len : usize) -> usize {
        let mut i1 : usize =  st_i;
        for i2 in self.i_lst.iter_mut() {
            if (i1 >= max_len) {
                return  i1;
            }
            in_vec[i1] = *i2;
            i1 += 1usize;
        }
        return  i1;
    }

}

impl SpatialGrid {
    pub fn initialize(&mut self, x_range : &mut [f64], x_spacing : f64, y_range : &mut [f64], y_spacing : f64, z_range : &mut [f64], z_spacing : f64) {
        self.x_min = x_range[0];
        self.x_sp = x_spacing;
        self.y_min = y_range[0];
        self.y_sp = y_spacing;
        self.z_min = z_range[0];
        self.z_sp = z_spacing;
        self.x_bins = ((x_range[1] - self.x_min) / self.x_sp + 1.0) as usize;
        self.y_bins = ((y_range[1] - self.y_min) / self.y_sp + 1.0) as usize;
        self.z_bins = ((z_range[1] - self.z_min) / self.z_sp + 1.0) as usize;
        let mut tot_bins : usize =  self.x_bins * self.y_bins * self.z_bins;
        self.list_ar = vec![IntList::new(); tot_bins];
        return;
    }

    pub fn add_ent(&mut self, label : usize, crd : & [f64]) {
        let mut x_b : usize =  ((crd[0] - self.x_min) / self.x_sp) as usize;
        if (x_b < 0) {
            x_b = 0;
        }
        if (x_b >= self.x_bins) {
            x_b = self.x_bins - 1;
        }
        let mut y_b : usize =  ((crd[1] - self.y_min) / self.y_sp) as usize;
        if (y_b < 0) {
            y_b = 0;
        }
        if (y_b >= self.y_bins) {
            y_b = self.y_bins - 1;
        }
        let mut z_b : usize =  ((crd[2] - self.z_min) / self.z_sp) as usize;
        if (z_b < 0) {
            z_b = 0;
        }
        if (z_b >= self.z_bins) {
            z_b = self.z_bins - 1;
        }
        let mut ind : usize =  (z_b*self.y_bins + y_b)*self.x_bins + x_b;
        self.list_ar[ind].i_lst.push_back(label);
    }

    pub fn get_in_xyzrange(&mut self, out_lst : &mut Vec<usize>, max_len : usize, x_range : & [f64], y_range : & [f64], z_range : & [f64]) -> usize {
        let mut i_min : usize;
        let mut i_max : usize;
        let mut j_min : usize;
        let mut j_max : usize;
        let mut k_min : usize;
        let mut k_max : usize;
        
        if (x_range[0] < x_range[1]) {
            if (x_range[0] < self.x_min) {
                i_min = 0;
            }
            else {
                i_min = ((x_range[0] - self.x_min) / self.x_sp) as usize;
            }
             
            i_max = ((x_range[1] - self.x_min) / self.x_sp) as usize;
            if (i_max >= self.x_bins) {
                i_max = self.x_bins - 1;
            }
        }
        else {
            i_min = 0;
            i_max = self.x_bins - 1;
        }
        
        if (y_range[0] < y_range[1]) {
            if (y_range[0] < self.y_min) {
                j_min = 0;
            }
            else {
                j_min = ((y_range[0] - self.y_min) / self.y_sp) as usize;
            }

            j_max = ((y_range[1] - self.y_min) / self.y_sp) as usize;
            if (j_max >= self.y_bins) {
                j_max = self.y_bins - 1;
            }
        }
        else {
            j_min = 0;
            j_max = self.y_bins - 1;
        }
        
        if (z_range[0] < z_range[1]) {
            if (z_range[0] < self.z_min) {
                k_min = 0;
            }
            else {
                k_min = ((z_range[0] - self.z_min) / self.z_sp) as usize;
            }

            k_max = ((z_range[1] - self.z_min) / self.z_sp) as usize;
            if (k_max >= self.z_bins) {
                k_max = self.z_bins - 1;
            }
        }
        else {
            k_min = 0;
            k_max = self.z_bins - 1;
        }
        
        let mut ind : usize;
        let mut lst_len : usize =  0;
        for k in k_min..=k_max {
            for j in j_min..=j_max {
                for i in i_min..=i_max {
                    ind = (k * self.y_bins + j)*self.x_bins + i;
                    lst_len = self.list_ar[ind].copy_to_vector(out_lst,  lst_len,  max_len);
                }
            }
        }
        
        return  lst_len;
    }

    pub fn get_in_radius(&mut self, out_list : &mut Vec<usize>, max_len : usize, pt : &mut [f64], rad : f64) -> usize {
        let mut range : [f64; 6] = [0f64; 6];
        range[0] = pt[0] - rad;
        range[1] = pt[0] + rad;
        range[2] = pt[1] - rad;
        range[3] = pt[1] + rad;
        range[4] = pt[2] - rad;
        range[5] = pt[2] + rad;
        
        return  self.get_in_xyzrange( out_list,  max_len, & range[0..2], & range[2..4], & range[4..6]);
    }

}


