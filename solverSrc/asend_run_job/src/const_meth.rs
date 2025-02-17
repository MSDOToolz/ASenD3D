use crate::constraint::*;
use crate::constants::*;
use crate::list_ent::*;
use crate::nd_el_set::*;
use crate::node::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;

use std::collections::LinkedList;

impl Constraint {
    pub fn build_mat(&mut self, nd_ar : &mut Vec<Node>, set_ar : &mut Vec<Set>) {
        let mut set_len : usize =  1;
        let mut seti_len : usize;
        let mut nd_index : usize;
        let mut row : usize;
        let mut col : usize;
        let mut dof : usize;
        let mut coef : f64;
        let mut set_nm = CppStr::new();
        let mut set_found : bool;
        let mut this_set : usize;
        let mut set_labs : &LinkedList<usize>;
        
        for tm in self.terms.iter_mut() {
            this_set = tm.ns_ptr;
            seti_len = set_ar[this_set].labels.len();
            if (seti_len > set_len) {
                set_len = seti_len;
            }
        }
        
        self.mat.set_dim(set_len);
        for tm in self.terms.iter_mut() {
            coef = tm.coef;
            dof = tm.dof;
            this_set = tm.ns_ptr;
            set_labs = &set_ar[this_set].labels;
            if (set_labs.len() == 1) {
                nd_index = match set_labs.front() {
                    None => 0usize,
                    Some(x) => *x,
                };
                if (self.this_type.s == "displacement") {
                    col = nd_ar[nd_index].dof_index[dof - 1];
                }
                else {
                    col = nd_ar[nd_index].sorted_rank;
                }
                for row in 0..set_len {
                    self.mat.add_entry(row, col, coef);
                }
            }
            else {
                row = 0;
                for nd in set_labs.iter() {
                    if (self.this_type.s == "displacement") {
                        col = nd_ar[*nd].dof_index[dof - 1];
                    }
                    else {
                        col = nd_ar[*nd].sorted_rank;
                    }
                    self.mat.add_entry(row,  col,  coef);
                    row += 1usize;
                }
            }
        }
        return;
    }

    pub fn full_vec_multiply(&mut self, prod : &mut Vec<f64>, vec : &mut Vec<f64>, tmp_v : &mut Vec<f64>) {
        // compute self.scale_fact*[self.mat]^t*([self.mat]*vec + q_vec)
        let mut i1 : usize;
        for tv in tmp_v.iter_mut() {
            *tv = 0.0;
        }
        self.mat.vector_multiply(tmp_v, vec, false);
        for tv in tmp_v.iter_mut() {
            *tv  *=  self.scale_fact;
        }
        self.mat.vector_multiply(prod, tmp_v, true);
        return;
    }

    pub fn get_load(&mut self, c_ld : &mut Vec<f64>, u_vec : &mut Vec<f64>, q_vec : &mut Vec<f64>, res_dim : usize) {
        let mut i1 : usize;
        let mut dim : usize =  self.mat.dim;
        for i1 in 0..dim {
            q_vec[i1] = -self.rhs;
        }
        self.mat.vector_multiply(q_vec, u_vec, false);
        for i1 in 0..dim {
            q_vec[i1]  *=  -self.scale_fact;
        }
        self.mat.vector_multiply(c_ld, q_vec, true);
        
        return;
    }

    // pub fn write_to_file(&mut self) {
    //     out_file.write(format!("{}{}{}", "type: " , self.this_type.s , "\n").as_bytes());
    //     out_file.write(format!("{}{}{}", "scale factor: " , self.scale_fact , "\n").as_bytes());
    //     out_file.write(format!("{}{}{}", "rhs: " , self.rhs , "\n").as_bytes());
    //     out_file.write(format!("{}", "matrix: \n").as_bytes());
    //     //self.mat.write_to_file(&mut out_file);
    //     return;
    // }

}

impl ConstraintList {
    pub fn set_scale_fact(&mut self, new_sf : f64) {
        for cnst in self.const_vec.iter_mut() {
            cnst.scale_fact = new_sf;
        }
        return;
    }

    pub fn get_total_vec_mult(&mut self, prod : &mut Vec<f64>, vec : &mut Vec<f64>, tmp_v : &mut Vec<f64>) {
        for cnst in self.const_vec.iter_mut() {
            cnst.full_vec_multiply(prod, vec, tmp_v);
        }
        return;
    }

    pub fn get_total_load(&mut self, c_ld : &mut Vec<f64>, u_vec : &mut Vec<f64>, q_vec : &mut Vec<f64>, res_dim : usize) {
        for cnst in self.const_vec.iter_mut() {
            cnst.get_load(c_ld, u_vec, q_vec, res_dim);
        }
        return;
    }

    // pub fn write_all_to_file(&mut self, file_name : &mut CppStr) {
    //     out_file.open(file_name);
    //     for cnst in self.const_vec.iter_mut() {
    //         cnst.write_to_file(&mut out_file);
    //         out_file.write(format!("{}", "##------------------------------------------------------\n").as_bytes());
    //     }
    //     out_file.close();
    //     return;
    // }

}


