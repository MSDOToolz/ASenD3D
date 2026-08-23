use crate::constants::MAX_INT;
use crate::list_ent::DualInt;
use crate::spatial_grid::*;


impl IntList {
    pub fn copy_to_vector(&mut self, in_vec : &mut Vec<usize>, st_i : usize, max_len : usize) -> usize {
        let mut i1 : usize =  st_i;
        for i2 in self.i_lst.iter_mut() {
            if i1 >= max_len {
                return  i1;
            }
            in_vec[i1] = *i2;
            i1 += 1usize;
        }
        return  i1;
    }

}

impl SpatialGrid {
    pub fn initialize(&mut self, x_range : &mut [f64], x_spacing : f64, y_range : &mut [f64], y_spacing : f64, z_range : &mut [f64], z_spacing : f64, capacity : usize) {
        self.x_min = 0.5*(x_range[0] + x_range[1] - x_spacing*(self.x_bins as f64));
        self.x_sp = x_spacing;
        self.y_min = 0.5*(y_range[0] + y_range[1] - y_spacing*(self.y_bins as f64));
        self.y_sp = y_spacing;
        self.z_min = 0.5*(z_range[0] + z_range[1] - z_spacing*(self.z_bins as f64));
        self.z_sp = z_spacing;
        self.data = vec![DualInt {i1 : MAX_INT, i2 : MAX_INT}; capacity];
        self.next_avail = 0;
    }

    pub fn reset(&mut self) {
        self.first.clear();
        for d in self.data.iter_mut() {
            d.i1 = MAX_INT;
            d.i2 = MAX_INT;
        }
        self.next_avail = 0;
    }

    pub fn add_ent(&mut self, label : usize, crd : & [f64]) {
        let mut x_b : usize;
        if crd[0] < self.x_min {
            x_b = 0;
        }
        else {
            x_b = ((crd[0] - self.x_min) / self.x_sp) as usize;
        }
        if x_b >= self.x_bins {
            x_b = self.x_bins - 1;
        }

        let mut y_b : usize;
        if crd[1] < self.y_min {
            y_b = 0;
        }
        else {
            y_b = ((crd[1] - self.y_min) / self.y_sp) as usize;
        }
        if y_b >= self.y_bins {
            y_b = self.y_bins - 1;
        }

        let mut z_b : usize;
        if crd[2] < self.z_min {
            z_b = 0;
        }
        else {
            z_b = ((crd[2] - self.z_min) / self.z_sp) as usize;
        }
        if z_b >= self.z_bins {
            z_b = self.z_bins - 1;
        }

        let ind : usize =  (z_b*self.y_bins + y_b)*self.x_bins + x_b;
        let na = self.next_avail;
        self.data[na].i1 = label;
        let mut d_ind : usize;
        match self.first.get(&ind) {
            None => {self.first.insert(ind, na);},
            Some(x) => {d_ind = *x;
                                while self.data[d_ind].i2 < MAX_INT {
                                    d_ind = self.data[d_ind].i2;
                                }
                                self.data[d_ind].i2 = na;},
        }
        self.next_avail += 1;
    }

    pub fn get_in_xyzrange(&self, out_lst : &mut Vec<usize>, max_len : usize, x_range : & [f64], y_range : & [f64], z_range : & [f64]) -> usize {
        let i_min : usize;
        let mut i_max : usize;
        let j_min : usize;
        let mut j_max : usize;
        let k_min : usize;
        let mut k_max : usize;
        
        if x_range[0] < x_range[1] {
            if x_range[0] < self.x_min {
                i_min = 0;
            }
            else {
                i_min = ((x_range[0] - self.x_min) / self.x_sp) as usize;
            }
             
            i_max = ((x_range[1] - self.x_min) / self.x_sp) as usize;
            if i_max >= self.x_bins {
                i_max = self.x_bins - 1;
            }
        }
        else {
            i_min = 0;
            i_max = self.x_bins - 1;
        }
        
        if y_range[0] < y_range[1] {
            if y_range[0] < self.y_min {
                j_min = 0;
            }
            else {
                j_min = ((y_range[0] - self.y_min) / self.y_sp) as usize;
            }

            j_max = ((y_range[1] - self.y_min) / self.y_sp) as usize;
            if j_max >= self.y_bins {
                j_max = self.y_bins - 1;
            }
        }
        else {
            j_min = 0;
            j_max = self.y_bins - 1;
        }
        
        if z_range[0] < z_range[1] {
            if z_range[0] < self.z_min {
                k_min = 0;
            }
            else {
                k_min = ((z_range[0] - self.z_min) / self.z_sp) as usize;
            }

            k_max = ((z_range[1] - self.z_min) / self.z_sp) as usize;
            if k_max >= self.z_bins {
                k_max = self.z_bins - 1;
            }
        }
        else {
            k_min = 0;
            k_max = self.z_bins - 1;
        }
        
        let mut f_ind : usize;
        let mut lst_len : usize =  0;
        let mut d_ind : usize;
        for k in k_min..=k_max {
            for j in j_min..=j_max {
                for i in i_min..=i_max {
                    f_ind = (k*self.y_bins + j)*self.x_bins + i;
                    //d_ind = self.first[f_ind];
                    match self.first.get(&f_ind) {
                        None => {},
                        Some(x) => {d_ind = *x;
                                            while d_ind < MAX_INT && lst_len < max_len {
                                                out_lst[lst_len] = self.data[d_ind].i1;
                                                d_ind = self.data[d_ind].i2;
                                                lst_len += 1;
                                            }},
                    }
                    
                    //lst_len = self.list_ar[ind].copy_to_vector(out_lst,  lst_len,  max_len);
                }
            }
        }
        
        return  lst_len;
    }

    pub fn get_in_radius(&self, out_list : &mut Vec<usize>, max_len : usize, pt : & [f64], rad : f64) -> usize {
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


