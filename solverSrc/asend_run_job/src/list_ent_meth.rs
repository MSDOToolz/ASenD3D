use crate::list_ent::*;

use crate::fmath::*;

impl MatrixRow {
    pub fn add_entry(&mut self, row : usize, col : usize, val : f64) {
        for me in self.row_vec.iter_mut() {
            if me.col == col {
                me.value += val;
                return;
            }
        }
        let mut new_ent = MatrixEnt::new();
        new_ent.row = row;
        new_ent.col = col;
        new_ent.value = val;
        self.row_vec.push_back(new_ent);
    }

    pub fn scale_to_sum(&mut self, new_sum : f64) {
        let mut curr_sum = 0.0f64;
        for i in self.row_vec.iter() {
            curr_sum += i.value;
        }
        let scale_fact = new_sum/curr_sum;
        for i in self.row_vec.iter_mut() {
            i.value *= scale_fact;
        }
    }

}

impl SparseMat {
    pub fn set_dim(&mut self, new_dim : usize) {
        self.dim = new_dim;
        self.matrix = vec![MatrixRow::new(); new_dim];
        return;
    }

    pub fn zero_all(&mut self) {
        for i1 in self.matrix.iter_mut() {
            for i2 in i1.row_vec.iter_mut() {
                i2.value = 0.0;
            }
        }
        
        return;
    }

    pub fn add_entry(&mut self, row : usize, col : usize, val : f64) {
        self.matrix[row].add_entry(row, col, val);
        return;
    }

    pub fn add_matrix(&mut self, inp_mat : &mut SparseMat) {
        for i1 in inp_mat.matrix.iter_mut() {
            for i2 in i1.row_vec.iter_mut() {
                self.add_entry(i2.row,   i2.col,   i2.value);
            }
        }
        return;
    }

    pub fn vector_multiply(&mut self, prod : &mut Vec<f64>, inp_vec : &mut Vec<f64>, transpose : bool) {
        if transpose {
            for i1 in self.matrix.iter_mut() {
                for i2 in i1.row_vec.iter_mut() {
                    prod[i2.col]  +=  i2.value * inp_vec[i2.row];
                }
            }
        } else {
            for i1 in self.matrix.iter_mut() {
                for i2 in i1.row_vec.iter_mut() {
                    prod[i2.row]  +=  i2.value * inp_vec[i2.col];
                }
            }
        }
        return;
    }

    pub fn get_max_abs_val(&mut self) -> f64 {
        let mut this_val : f64;
        let mut max_val : f64 =  0.0;
        for i1 in self.matrix.iter_mut() {
            for i2 in i1.row_vec.iter_mut() {
                this_val = fabs(i2.value);
                if this_val > max_val {
                    max_val = this_val;
                }
            }
        }
        
        return  max_val;
    }

    pub fn scale_row_to_sum(&mut self, row : usize, new_sum : f64) {
        self.matrix[row].scale_to_sum(new_sum);
    }

    // pub fn write_to_file(&mut self) {
    //     for i1 in self.matrix.iter_mut() {
    //         for i2 in i1.row_vec.iter_mut() {
    //             out_file.write(format!("{}{}{}{}{}{}{}", "    - [" , i2.row , ", " , i2.col , ", " , i2.value , "]\n").as_bytes());
    //         }
    //     }
    //     return;
    // }

}


