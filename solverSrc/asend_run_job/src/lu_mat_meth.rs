use crate::lu_mat::*;
use crate::list_ent::*;
use crate::constraint::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


impl LUMat {
    pub fn set_dim(&mut self, new_dim : usize) {
        self.dim = new_dim;
        self.l_range = vec![0usize; new_dim+1];
        self.l_min_col = vec![0usize; new_dim];
        self.u_range = vec![0usize; new_dim+1];
        self.u_min_row = vec![0usize; new_dim];
        self.z_vec = vec![0f64; new_dim];
        return;
    }

    pub fn allocate_from_sparse_mat(&mut self, sp_mat : &mut SparseMat, c_list : &mut ConstraintList, block_dim : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut curr_block : usize;
        let mut const_dim : usize;
        let mut blk_max_col : usize;
        let mut blk_min_col : usize;
        
        if(!self.allocated) {
            self.set_dim(sp_mat.dim);
        }
        
        for i1 in 0..self.dim {
            self.l_range[i1] = 0;
            self.u_range[i1] = 1;
        }
        
        for i1 in 0..self.dim {
            curr_block = i1/block_dim;
            blk_min_col = curr_block*block_dim;
            blk_max_col = blk_min_col + block_dim;
            let mut mr = &sp_mat.matrix[i1];
            for me in mr.row_vec.iter() {
                i2 = me.col;
                if(i2 >= blk_min_col && i2 < blk_max_col) {
                    if(i2 < i1) {// in l
                        i3 = i1 - i2;
                        if(i3 > self.l_range[i1]) {
                            self.l_range[i1] = i3;
                        }
                    }
                    else {
                        i3 = i2 - i1 + 1;
                        if(i3 > self.u_range[i2]) {
                            self.u_range[i2] = i3;
                        }
                    }
                }
            }
        }
        
        for cnst in c_list.const_vec.iter() {
            let mut this_mat = &cnst.mat;
            for mr in this_mat.matrix.iter() {
                for me in mr.row_vec.iter() {
                    i2 = me.col;
                    curr_block = i2/block_dim;
                    blk_min_col = curr_block*block_dim;
                    blk_max_col = blk_min_col + block_dim;
                    for me2 in mr.row_vec.iter() {
                        i3 = me2.col;
                        if (i3 >= blk_min_col && i3 < blk_max_col) {
                            if (i3 < i2) {// in lu
                                i4 = i2 - i3;
                                if (i4 > self.l_range[i2]) {
                                    self.l_range[i2] = i4;
                                }
                            }
                            else {
                                i4 = i3 - i2 + 1;
                                if(i4 > self.u_range[i3]) {
                                    self.u_range[i3] = i4;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        self.l_size = 0;
        self.u_size = 0;
        self.max_bandwidth = 0;
        for i1 in 0..self.dim {
            self.l_size  +=  self.l_range[i1];
            self.u_size  +=  self.u_range[i1];
            self.l_min_col[i1] = i1 - self.l_range[i1];
            if(self.l_range[i1] > self.max_bandwidth) {
                self.max_bandwidth = self.l_range[i1];
            }
            self.u_min_row[i1] = i1 - self.u_range[i1] + 1;
            if(self.u_range[i1] > self.max_bandwidth) {
                self.max_bandwidth = self.u_range[i1];
            }
        }
        self.l_range[self.dim] = self.l_size;
        self.u_range[self.dim] = self.u_size;
        for i1 in (0..(self.dim)).rev() {
            self.l_range[i1] = self.l_range[i1+1] - self.l_range[i1];
            self.u_range[i1] = self.u_range[i1+1] - self.u_range[i1];
        }
        
        self.l_mat = vec![0f64; self.l_size];
        self.u_mat = vec![0f64; self.u_size];
        
        self.allocated = true;
        
        println!("{}{}", "maxBandwidth: " , self.max_bandwidth );
        
        return;
    }

    pub fn populate_from_sparse_mat(&mut self, sp_mat : &mut SparseMat, c_list : &mut ConstraintList) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut const_dim : usize;
        let mut const_sf : f64;
        
        for lm in self.l_mat.iter_mut() {
            *lm = 0.0;
        }
        
        for um in self.u_mat.iter_mut() {
            *um = 0.0;
        }
        
        for i1 in 0..self.dim {
            let mut mr = &sp_mat.matrix[i1];
            for me in mr.row_vec.iter() {
                i2 = me.col;
                if(i2 >= self.l_min_col[i1] && i2 < i1) {// in l
                    i3 = self.l_range[i1] + i2 - self.l_min_col[i1];
                    self.l_mat[i3]  +=  me.value;
                }
                else if(i1 >= self.u_min_row[i2] && i1 <= i2) {
                    i3 = self.u_range[i2] + i1 - self.u_min_row[i2];
                    self.u_mat[i3]  +=  me.value;
                }
            }
        }
        
        for cnst in c_list.const_vec.iter() {
            let mut this_mat = &cnst.mat;
            const_sf = cnst.scale_fact;
            for mr in this_mat.matrix.iter() {
                for me in mr.row_vec.iter() {
                    i2 = me.col;
                    for me2 in mr.row_vec.iter() {
                        i3 = me2.col;
                        if(i3 >= self.l_min_col[i2] && i3 < i2) {// in lu
                            i4 = self.l_range[i2] + i3 - self.l_min_col[i2];
                            self.l_mat[i4]  +=  const_sf*me.value*me2.value;
                        }
                        else if(i2 >= self.u_min_row[i3] && i2 <= i3) {
                            i4 = self.u_range[i3] + i2 - self.u_min_row[i3];
                            self.u_mat[i4]  +=  const_sf*me.value*me2.value;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn lu_factor(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut this_ind : usize;
        let mut l_ind : usize;
        let mut u_ind : usize;
        let mut u_piv : usize;
        let mut max_col : usize;
        let mut max_row : usize;
        let mut tmp : f64;
        
        for i1 in 0..self.dim {
            max_col = i1 + self.max_bandwidth;
            if(max_col > self.dim) {
                max_col = self.dim;
            }
            for i2 in i1..max_col {
                if(self.u_min_row[i2] <= i1) {
                    this_ind = self.u_range[i2] + i1 - self.u_min_row[i2];
                    i3 = self.u_min_row[i2];
                    if(self.l_min_col[i1] > i3) {
                        i3 = self.l_min_col[i1];
                    }
                    l_ind = self.l_range[i1] + i3 - self.l_min_col[i1];
                    u_ind = self.u_range[i2] + i3 - self.u_min_row[i2];
                    tmp = self.u_mat[this_ind];
                    while(l_ind < self.l_range[i1+1]) {
                        tmp  -=  self.l_mat[l_ind]*self.u_mat[u_ind];
                        l_ind += 1usize;
                        u_ind += 1usize;
                    }
                    self.u_mat[this_ind] = tmp;
                }
            }
            max_row = max_col;
            for i2 in (i1+1)..max_row {
                if(self.l_min_col[i2] <= i1) {
                    this_ind = self.l_range[i2] + i1 - self.l_min_col[i2];
                    i3 = self.l_min_col[i2];
                    if(self.u_min_row[i1] > i3) {
                        i3 = self.u_min_row[i1];
                    }
                    l_ind = self.l_range[i2] + i3 - self.l_min_col[i2];
                    u_ind = self.u_range[i1] + i3 - self.u_min_row[i1];
                    tmp = self.l_mat[this_ind];
                    u_piv = self.u_range[i1+1] - 1;
                    while(u_ind < u_piv) {
                        tmp  -=  self.l_mat[l_ind]*self.u_mat[u_ind];
                        l_ind += 1usize;
                        u_ind += 1usize;
                    }
                    self.l_mat[this_ind] = tmp/self.u_mat[u_piv];
                }
            }
        }
        
        return;
    }

    pub fn lu_solve(&mut self, soln_vec : &mut Vec<f64>, rhs : &mut Vec<f64>, transpose : bool) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut max_col : usize;
        let mut max_row : usize;
        let mut u_piv : usize;
        let mut tmp : f64;
        
        if(!transpose) {
            for i1 in 0..self.dim {
                tmp = rhs[i1];
                i2 = self.l_min_col[i1];
                for i3 in self.l_range[i1]..self.l_range[i1+1] {
                    tmp  -=  self.l_mat[i3]*self.z_vec[i2];
                    i2 += 1usize;
                }
                self.z_vec[i1] = tmp;
            }
            for i1 in (0..(self.dim)).rev() {
                tmp = self.z_vec[i1];
                max_col = i1 + self.max_bandwidth;
                if(max_col > self.dim) {
                    max_col = self.dim;
                }
                for i2 in (i1+1)..max_col {
                    if(self.u_min_row[i2] <= i1) {
                        i3 = self.u_range[i2] + i1 - self.u_min_row[i2];
                        tmp  -=  self.u_mat[i3]*soln_vec[i2];
                    }
                }
                u_piv = self.u_range[i1+1] - 1;
                soln_vec[i1] = tmp/self.u_mat[u_piv];
            }
        }
        else {
            for i1 in 0..self.dim {
                tmp = rhs[i1];
                i2 = self.u_min_row[i1];
                for i3 in self.u_range[i1]..(self.u_range[i1+1]-1) {
                    tmp  -=  self.u_mat[i3]*self.z_vec[i2];
                    i2 += 1usize;
                }
                u_piv = self.u_range[i1+1] - 1;
                self.z_vec[i1] = tmp/self.u_mat[u_piv];
            }
            for i1 in (0..(self.dim)).rev() {
                tmp = self.z_vec[i1];
                max_row = i1 + self.max_bandwidth;
                if(max_row > self.dim) {
                    max_row = self.dim;
                }
                for i2 in (i1+1)..max_row {
                    if(self.l_min_col[i2] <= i1) {
                        i3 = self.l_range[i2] + i1 - self.l_min_col[i2];
                        tmp  -=  self.l_mat[i3]*soln_vec[i2];
                    }
                }
                soln_vec[i1] = tmp;
            }
        }
        
        return;
    }

}


