use crate::lower_tri_mat::*;
use crate::list_ent::*;
use crate::constraint::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


impl LowerTriMat {
    pub fn set_dim(&mut self, new_dim : usize) {
        self.dim = new_dim;
        self.range = vec![0usize; new_dim+1];
        self.min_col = vec![0usize; new_dim];
        self.z_vec = vec![0f64; new_dim];
        return;
    }

    pub fn allocate_from_sparse_mat(&mut self, sp_mat : &mut SparseMat, c_list : &mut ConstraintList, block_dim : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut row : usize;
        let mut curr_block : usize;
        let mut blk_mc : usize;
        let mut const_dim : usize;
        let mut this_mat : &SparseMat;
        
        if(!self.allocated) {
            self.set_dim(sp_mat.dim);
        }
        
        for mr in sp_mat.matrix.iter() {
            row = match mr.row_vec.front() {
                None => 0usize,
                Some(x) => x.row,
            };
            self.range[row] = 1;
            curr_block = row / block_dim;
            blk_mc = block_dim * curr_block;
            for me in mr.row_vec.iter() {
                i2 = me.col;
                if(i2 <= row && i2 >= blk_mc) {
                    i3 = row - i2 + 1;
                    if(i3 > self.range[row]) {
                        self.range[row] = i3;
                    }
                }
            }
        }
        
        for cnst in c_list.const_vec.iter() {
            this_mat = &cnst.mat;
            for mr in this_mat.matrix.iter() {
                for me in mr.row_vec.iter() {
                    i2 = me.col;
                    curr_block = i2 / block_dim;
                    blk_mc = block_dim * curr_block;
                    for me2 in mr.row_vec.iter() {
                        i3 = me2.col;
                        if(i3 <= i2 && i3 >= blk_mc) {
                            i4 = i2 - i3 + 1;
                            if(i4 > self.range[i2]) {
                                self.range[i2] = i4;
                            }
                        }
                    }
                }
            }
        }
        
        self.size = 0;
        self.max_bandwidth = 0;
        for i1 in 0..self.dim {
            self.size +=  self.range[i1];
            self.min_col[i1] = i1 - self.range[i1] + 1;
            if(self.range[i1] > self.max_bandwidth) {
                self.max_bandwidth = self.range[i1];
            }
        }
        self.range[self.dim] = self.size;
        for i1 in (0..(self.dim)).rev() {
            self.range[i1] = self.range[i1+1] - self.range[i1];
        }
        
        self.mat = vec![0f64; self.size];
        
        self.allocated = true;
        
        println!("{}{}", "maxBandwidth: " , self.max_bandwidth );
        
        return;
    }

    pub fn is_allocated(&mut self) -> bool {
        return  self.allocated;
    }

    pub fn populate_from_sparse_mat(&mut self, sp_mat : &mut SparseMat, c_list : &mut ConstraintList) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut const_dim : usize;
        let mut const_sf : f64;
        let mut this_mat : &SparseMat;
        
        for me in self.mat.iter_mut() {
            *me = 0.0;
        }
        
        i1 = 0;
        for mr in sp_mat.matrix.iter() {
            for me in mr.row_vec.iter() {
                i2 = me.col;
                if(i2 <= i1 && i2 >= self.min_col[i1]) {
                    i3 = self.range[i1] + (i2 - self.min_col[i1]);
                    self.mat[i3] +=  me.value;
                }
            }
            i1 += 1usize;
        }
        
        for cnst in c_list.const_vec.iter() {
            this_mat = &cnst.mat;
            const_sf = cnst.scale_fact;
            for mr in this_mat.matrix.iter() {
                for me in mr.row_vec.iter() {
                    i2 = me.col;
                    for me2 in mr.row_vec.iter() {
                        i3 = me2.col;
                        if(i3 <= i2 && i3 >= self.min_col[i2]) {
                            i4 = self.range[i2] + (i3 - self.min_col[i2]);
                            self.mat[i4]  +=  const_sf * me.value * me2.value;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn ldl_factor(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut st_col : usize;
        let mut sum : f64;
        
        let mut ld_vec = vec![0f64; self.dim];
        
        for i1 in 0..self.dim {// i1 = column in l
            //form ld_vec
            i3 = self.range[i1];
            for i2 in self.min_col[i1]..i1 {
                i4 = self.range[i2+1] - 1;
                ld_vec[i2] = self.mat[i3]*self.mat[i4];
                i3 += 1usize;
            }
            //get d term for column i1
            st_col = self.min_col[i1];
            i2 = self.range[i1];
            sum = 0.0;
            for i3 in st_col..i1 {
                sum +=  ld_vec[i3]*self.mat[i2];
                i2 += 1usize;
            }
            i2 = self.range[i1+1] - 1;
            self.mat[i2] = self.mat[i2] - sum;
            //get l terms for column i1
            i2 = i1 + 1;
            while(i2 < (i1 + self.max_bandwidth) && i2 < self.dim) {// i2 = row in l
                if(self.min_col[i2] <= i1) {
                    st_col = self.min_col[i2];
                    if(self.min_col[i1] > st_col) {
                        st_col = self.min_col[i1];
                    }
                    i4 = self.range[i2] + (st_col - self.min_col[i2]);
                    sum = 0.0;
                    for i5 in st_col..i1 {
                        sum +=  ld_vec[i5]*self.mat[i4];
                        i4 += 1usize;
                    }
                    i3 = self.range[i2] + (i1 - self.min_col[i2]);
                    i4 = self.range[i1+1] - 1;
                    self.mat[i3] = (self.mat[i3] - sum)/self.mat[i4];
                }
                i2 += 1usize;
            }
        }
        
        return;
    }

    pub fn ldl_solve(&mut self, soln_vec : &mut Vec<f64>, rhs : &mut Vec<f64>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut stop_row : usize;
        let mut sum : f64;
        
        self.z_vec[0] = rhs[0];
        for i1 in 1..self.dim {
            i2 = self.range[i1];
            sum = 0.0;
            for i3 in self.min_col[i1]..i1 {
                sum +=  self.mat[i2]*self.z_vec[i3];
                i2 += 1usize;
            }
            self.z_vec[i1] = rhs[i1] - sum;
        }
        
        for i1 in (0..self.dim).rev() {
            sum = 0.0;
            stop_row = i1 + self.max_bandwidth + 1;
            if(stop_row > self.dim) {
                stop_row = self.dim;
            }
            for i2 in (i1+1)..stop_row {
                i3 = self.range[i2] + (i1 - self.min_col[i2]);
                if(i3 >= self.range[i2] && i3 < self.range[i2+1]) {
                    sum += self.mat[i3]*soln_vec[i2];
                }
            }
            i2 = self.range[i1+1] - 1;
            soln_vec[i1] = (self.z_vec[i1]/self.mat[i2]) - sum;
        }
        return;
    }

    pub fn pos_def(&mut self) -> bool {
        let mut i1 : usize;
        let mut d_val : f64;
        if (!self.allocated) {
            return  false;
        }
        for i1 in 0..self.dim {
            d_val = self.mat[self.range[i1 + 1] - 1];
            if (d_val < 0.0) {
                return  false;
            }
        }
        return  true;
    }

    // pub fn write_to_file(&mut self) {
    //     let mut i1 : usize;
    //     let mut i2 : usize;
    //     let mut col : usize;
        
    //     for i1 in 0..self.dim {
    //         col = self.min_col[i1];
    //         for i2 in self.range[i1]..self.range[i1 + 1] {
    //             out_file.write(format!("{}{}{}{}{}{}{}", "    - [" , i1 , ", " , col , ", " , self.mat[i2] , "]\n").as_bytes());
    //             col += 1usize;
    //         }
    //     }
        
    //     return;
    // }

}


