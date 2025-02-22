use crate::constants::*;
use crate::fmath::*;
use crate::cpp_str::CppStr;

use std::fs::File;
use std::io::{self, Read, BufRead};
use std::path::Path;

pub fn cross_prod(prod : &mut [f64], v1 : &mut [f64], v2 : &mut [f64]) {
    prod[0] = v1[1] * v2[2] - v1[2] * v2[1];
    prod[1] = v1[2] * v2[0] - v1[0] * v2[2];
    prod[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return;
}

pub fn q_rfactor(mat : &mut [f64], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta : f64;
    let mut sth : f64;
    let mut cth : f64;
    let mut p1 : f64;
    let mut p2 : f64;
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if (tri_diag == 0) {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if (i2_max > end_row) {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min - 1) * col_dim + i1;
            k12 = i2 * col_dim + i1;
            if (fabs(mat[k11]) < TOL) {
                mat[k11] = TOL;
            }
            theta = atan(mat[k12] / mat[k11]);
            sth = sin(theta);
            cth = cos(theta);
            i3_min = i1;
            if (tri_diag == 2) {
                i3_max = i1 + 2;
                if (i3_max > end_col) {
                    i3_max = end_col;
                }
            }
            else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min - 1) * col_dim + i3;
                k23 = i2 * col_dim + i3;
                p1 = cth * mat[k22] + sth * mat[k23];
                p2 = -sth * mat[k22] + cth * mat[k23];
                mat[k22] = p1;
                mat[k23] = p2;
            }
            mat[k12] = theta;
        }
    }
    return;
}

pub fn solveq_rx_eqb(x_vec : &mut [f64], mat : &mut [f64], b_vec : &mut [f64], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i1 : usize;
    let mut i2 : usize;
    let mut i3 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut theta : f64;
    let mut sth : f64;
    let mut cth : f64;
    let mut p1 : f64;
    let mut p2 : f64;
    let mut row_sum : f64;
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if (tri_diag == 0) {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if (i2_max > end_row) {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2 * col_dim + i1;
            theta = mat[k12];
            sth = sin(theta);
            cth = cos(theta);
            p1 = cth * b_vec[i3] + sth * b_vec[i2];
            p2 = -sth * b_vec[i3] + cth * b_vec[i2];
            b_vec[i3] = p1;
            b_vec[i2] = p2;
        }
        x_vec[i1] = 0.0;
    }
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if (tri_diag == 2) {
            i2_max = i1 + 2;
            if (i2_max > end_col) {
                i2_max = end_col;
            }
        }
        else {
            i2_max = end_col;
        }
        row_sum = 0.0;
        k11 = i3 * col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            row_sum  +=  mat[k11] * x_vec[i2];
            k11 += 1usize;
        }
        k11 = i3 * col_dim + i1;
        x_vec[i1] = (b_vec[i3] - row_sum) / mat[k11];
    }
    
    return;
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_input_line(file_line : &mut CppStr, headings : &mut Vec<CppStr>, hd_ld_space : &mut [usize], data : &mut Vec<CppStr>, data_len : &mut usize) {
    let mut i1 : usize;
    let i2 : usize;
    let mut ln_len : usize;
    let wrd_len : usize;
    i1 = file_line.find("#");
    if i1 < MAX_INT {
        *file_line = file_line.substr(0,i1);
    }
    file_line.s = file_line.s.clone() + " ";
    ln_len = file_line.len();
    i1 = file_line.find(":");
    *data_len = 0;
    if i1 < MAX_INT {
        i2 = file_line.find_first_not_of(" -\n\t");
        wrd_len = i1 - i2;
        if headings[0].s == "" || hd_ld_space[0] == i2 {
            headings[0] = file_line.substr(i2,wrd_len);
            hd_ld_space[0] = i2;
            headings[1] = CppStr::from("");
            hd_ld_space[1] = 0;
            headings[2] = CppStr::from("");
            hd_ld_space[2] = 0;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else if headings[1].s == "" || hd_ld_space[1] == i2 {
            headings[1] = file_line.substr(i2,wrd_len);
            hd_ld_space[1] = i2;
            headings[2] = CppStr::from("");
            hd_ld_space[2] = 0;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else if headings[2].s == "" || hd_ld_space[2] == i2 {
            headings[2] = file_line.substr(i2,wrd_len);
            hd_ld_space[2] = i2;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else {
            headings[3] = file_line.substr(i2,wrd_len);
            hd_ld_space[3] = i2;
        }
        i1 += 1usize;
        while i1 < ln_len {
            *file_line = file_line.substr(i1, MAX_INT);
            i1 = file_line.find_first_not_of(" ,[]\t\n");
            if i1 < MAX_INT {
                *file_line = file_line.substr(i1, MAX_INT);
                ln_len = file_line.len();
                i1 = file_line.find_first_of(" ,[]\t\n");
                if i1 < MAX_INT {
                    data[*data_len] = file_line.substr(0,i1);
                    *data_len += 1usize;
                } else {
                    i1 = ln_len;
                }
            } else {
                i1 = ln_len;
            }
        }
    } else {
        i1 = file_line.find("- ");
        if i1 < MAX_INT {
            i1 += 1usize;
            while i1 < ln_len {
                *file_line = file_line.substr(i1, MAX_INT);
                i1 = file_line.find_first_not_of(" ,[]\t\n");
                if i1 < MAX_INT {
                    *file_line = file_line.substr(i1, MAX_INT);
                    ln_len = file_line.len();
                    i1 = file_line.find_first_of(" ,[]\t\n");
                    if i1 < MAX_INT {
                        data[*data_len] = file_line.substr(0,i1);
                        *data_len += 1usize;
                    } else {
                        i1 = ln_len;
                    }
                } else {
                    i1 = ln_len;
                }
            }
        }
    }
    
    return;
}