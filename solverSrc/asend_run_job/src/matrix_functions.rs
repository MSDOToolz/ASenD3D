use crate::constants::*;
use crate::diff_doub::*;
use crate::list_ent::*;
use crate::lu_mat::*;
use crate::lower_tri_mat::*;
use crate::constraint::*;
use crate::fmath::*;

pub fn sub_vec(sub_v : &mut Vec<f64>, v_in : &mut Vec<f64>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        sub_v[i2] = v_in[i1];
        i2 += 1usize;
    }
}

pub fn return_sv(sub_v : &mut Vec<f64>, v_in : &mut Vec<f64>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        v_in[i1] = sub_v[i2];
        i2 += 1usize;
    }
}

pub fn vec_to_ar(ar : &mut [f64], vc : & Vec<f64>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        ar[i2] = vc[i1];
        i2 += 1usize;
    }
}

pub fn get_dist(p1 : & [f64], p2 : & [f64]) -> f64 {
    let dist : f64;
    let mut vec : [f64; 3] = [0f64; 3];
    vec[0] = p1[0] - p2[0];
    vec[1] = p1[1] - p2[1];
    vec[2] = p1[2] - p2[2];
    dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    return  dist;
}

pub fn cross_prod(prod : &mut [f64], v1 : & [f64], v2 : & [f64]) {
    prod[0] = v1[1]*v2[2] - v1[2]*v2[1];
    prod[1] = v1[2]*v2[0] - v1[0]*v2[2];
    prod[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return;
}

pub fn q_rfactor(mat : &mut Vec<f64>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
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
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min-1)*col_dim + i1;
            k12 = i2*col_dim + i1;
            if fabs(mat[k11]) < TOL {
                mat[k11] = TOL;
            }
            theta = atan(mat[k12]/mat[k11]);
            sth = sin(theta);
            cth = cos(theta);
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > end_col {
                    i3_max = end_col;
                }
            } else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min-1)*col_dim + i3;
                k23 = i2*col_dim + i3;
                p1 = cth*mat[k22] + sth*mat[k23];
                p2 = -sth*mat[k22] + cth*mat[k23];
                mat[k22] = p1;
                mat[k23] = p2;
            }
            mat[k12] = theta;
        }
    }
    return;
}

pub fn solveq_rx_eqb(x_vec : &mut Vec<f64>, mat : &mut Vec<f64>, b_vec : &mut Vec<f64>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
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
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2*col_dim + i1;
            theta = mat[k12];
            sth = sin(theta);
            cth = cos(theta);
            p1 = cth*b_vec[i3] + sth*b_vec[i2];
            p2 = -sth*b_vec[i3] + cth*b_vec[i2];
            b_vec[i3] = p1;
            b_vec[i2] = p2;
        }
        x_vec[i1] = 0.0;
    }
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if tri_diag == 2 {
            i2_max = i1 + 2;
            if i2_max > end_col {
                i2_max = end_col;
            }
        } else {
            i2_max = end_col;
        }
        row_sum = 0.0;
        k11 = i3*col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            row_sum +=  mat[k11]*x_vec[i2];
            k11 += 1usize;
        }
        k11 = i3*col_dim + i1;
        x_vec[i1] = (b_vec[i3] - row_sum)/mat[k11];
    }
    
    return;
}

pub fn conj_grad_sparse(soln : &mut Vec<f64>, mat : &mut SparseMat, cnst : &mut ConstraintList, pc_mat : &mut LowerTriMat, rhs : &mut Vec<f64>, conv_tol : f64, max_it : usize) {
    let mut i1 : usize;
    let mut i3 : usize;
    let dim : usize =  mat.dim;
    let status_freq : usize =  dim / 24;
    let mut res : f64;
    let mut r_next : f64;
    let mut alpha : f64;
    let mut dp : f64;
    let mut beta : f64;
    let mut h_vec = vec![0f64; dim];
    let mut g_vec = vec![0f64; dim];
    let mut z_vec = vec![0f64; dim];
    let mut w_vec = vec![0f64; dim];
    let mut t_vec = vec![0f64; dim];
    
    if !pc_mat.pos_def() {
        println!("{}", "Warning: preconditioning matrix for conjugate gradient solver is not positive definite." );
    }
    
    for i1 in 0..dim {
        soln[i1] = 0.0;
        g_vec[i1] = -rhs[i1];
    }
    
    pc_mat.ldl_solve(&mut w_vec, &mut  g_vec);
    
    res = 0.0;
    for i1 in 0..dim {
        h_vec[i1] = -w_vec[i1];
        res  +=  g_vec[i1] * w_vec[i1];
    }
    
    i1 = 0;
    i3 = 0;
    while i1 < max_it && res > conv_tol {
        for i2 in 0..dim {
            z_vec[i2] = 0.0;
        }
        mat.vector_multiply(&mut z_vec, &mut  h_vec,  false);
        cnst.get_total_vec_mult(&mut z_vec, &mut  h_vec, &mut  t_vec);
        dp = 0.0;
        for i2 in 0..dim {
            dp  +=  h_vec[i2] * z_vec[i2];
        }
        alpha = res / dp;
        for i2 in 0..dim {
            soln[i2]  +=  alpha * h_vec[i2];
            g_vec[i2]  +=  alpha * z_vec[i2];
        }
        pc_mat.ldl_solve(&mut w_vec, &mut  g_vec);
        r_next = 0.0;
        for i2 in 0..dim {
            r_next  +=  g_vec[i2] * w_vec[i2];
        }
        beta = r_next / res;
        for i2 in 0..dim {
            h_vec[i2] = -w_vec[i2] + beta * h_vec[i2];
        }
        res = r_next;
        
        if i3 == status_freq {
            println!("{}{}{}{}", "Conjugate gradient iteration: " , i1 , ", residual: " , res );
            i3 = 0;
        }
        i1 += 1usize;
        i3 += 1usize;
    }
    
    if i1 == max_it {
        println!("{}{}{}{}", "Warning: Conjugate gradient solver did not converge to specified tolerance, " , conv_tol , " within tne maximum number of iterations, " , max_it );
    }
    
    println!("{}{}", "Total CG iterations: " , i1 );
    println!("{}{}", "Final residual norm: " , res );
    
    return;
}

pub fn g_mres_sparse(soln : &mut Vec<f64>, mat : &mut SparseMat, cnst : &mut ConstraintList, pc_mat : &mut LUMat, rhs : &mut Vec<f64>, conv_tol : f64, max_it : usize, restart : usize) {
    let mut i1 : usize;
    let mut i2 : usize;
    let mut i3 : usize;
    let mut i4 : usize;
    let mut i5 : usize;
    let mut it_ct : usize;
    let mut res_nrm : f64;
    let mut tmp : f64;
    let dim : usize =  mat.dim;
    
    let mut res_vec = vec![0f64; dim];
    let mut tmp_v = vec![0f64; dim];
    let mut tmp_v2 = vec![0f64; dim];
    let mut tmp_v3 = vec![0f64; dim];
    i1 = dim * (restart + 1);
    let mut h_mat = vec![0f64; i1];
    i1 = restart * (restart + 1);
    let mut phi_mat = vec![0f64; i1];
    
    for i1 in 0..dim {
        soln[i1] = 0.0;
        tmp_v[i1] = -rhs[i1];
    }
    
    pc_mat.lu_solve(&mut res_vec, &mut  tmp_v,  false);
    //pc_mat.ldl_solve(res_vec, tmp_v);
    res_nrm = 0.0;
    for i1 in 0..dim {
        res_nrm  +=  res_vec[i1] * res_vec[i1];
    }
    res_nrm = sqrt(res_nrm);
    println!("{}{}", "Initial residual norm: " , res_nrm );
    
    it_ct = 0;
    while res_nrm > conv_tol && it_ct < max_it {
        tmp = 1.0 / res_nrm;
        for i1 in 0..dim {
            h_mat[i1] = tmp * res_vec[i1];
        }
        i2 = restart * (restart + 1);
        for i1 in 0..i2 {
            phi_mat[i1] = 0.0;
        }
        //generate basis vectors
        for i1 in 1..(restart + 1) {
            //multiply by previous vector
            for i2 in 0..dim {
                tmp_v[i2] = 0.0;
            }
            i2 = dim * (i1 - 1);
            sub_vec(&mut tmp_v3, &mut  h_mat,  i2,  i2 + dim);
            mat.vector_multiply(&mut tmp_v, &mut  tmp_v3,  false);
            cnst.get_total_vec_mult(&mut tmp_v, &mut  tmp_v3, &mut  tmp_v2);
            //pc_mat.ldl_solve(tmp_v2, tmp_v);
            pc_mat.lu_solve(&mut tmp_v2, &mut  tmp_v,  false);
            //orthogonalize with all previous vectors
            for i2 in 0..i1 {
                i4 = dim * i2;
                tmp = 0.0;
                for i3 in 0..dim {
                    tmp  +=  tmp_v2[i3] * h_mat[i4];
                    i4 += 1usize;
                }
                i5 = i2 * restart + (i1 - 1);
                phi_mat[i5] = tmp;
                i4 = dim * i2;
                for i3 in 0..dim {
                    tmp_v2[i3]  -=  tmp * h_mat[i4];
                    i4 += 1usize;
                }
            }
            //calculate magnitude and Set new unit basis std::vector
            tmp = 0.0;
            for i3 in 0..dim {
                tmp  +=  tmp_v2[i3] * tmp_v2[i3];
            }
            tmp = sqrt(tmp);
            i5 = i1 * restart + (i1 - 1);
            phi_mat[i5] = tmp;
            tmp = 1.0 / tmp;
            i4 = dim * i1;
            for i2 in 0..dim {
                h_mat[i4] = tmp * tmp_v2[i2];
                i4 += 1usize;
            }
        }
        //find the least squares solution
        i3 = 0;
        for i1 in 0..=restart {
            tmp_v[i1] = 0.0;
            for i2 in 0..dim {
                tmp_v[i1]  -=  h_mat[i3] * res_vec[i2];
                i3 += 1usize;
            }
        }
        q_rfactor(&mut phi_mat,  restart,  0,  restart,  0,  restart-1,  1);
        solveq_rx_eqb(&mut tmp_v2, &mut  phi_mat, &mut  tmp_v,  restart,  0,  restart,  0,  restart-1,  1);
        //update the solution std::vector
        i3 = 0;
        for i1 in 0..restart {
            for i2 in 0..dim {
                soln[i2]  +=  h_mat[i3] * tmp_v2[i1];
                i3 += 1usize;
            }
        }
        //update residual std::vector
        for i1 in 0..dim {
            tmp_v[i1] = -rhs[i1];
        }
        mat.vector_multiply(&mut tmp_v, soln,  false);
        cnst.get_total_vec_mult(&mut tmp_v, soln, &mut  tmp_v2);
        //pc_mat.ldl_solve(res_vec, tmp_v);
        pc_mat.lu_solve(&mut res_vec, &mut  tmp_v,  false);
        res_nrm = 0.0;
        for i1 in 0..dim {
            res_nrm  +=  res_vec[i1] * res_vec[i1];
        }
        res_nrm = sqrt(res_nrm);
        it_ct  +=  restart;
        println!("{}{}{}{}", "Iteration: " , it_ct , ",  Residual Norm: " , res_nrm );
    }
    
    if res_nrm > conv_tol {
        println!("{}{}", "Warning: GMRES solver did not converge to the requested tolerance of " , conv_tol );
        println!("{}{}{}{}", "Residual norm after " , it_ct , " iterations: " , res_nrm );
    }
    
    return;
}

pub fn sym_factor(mat : &mut Vec<f64>, q_mat : &mut Vec<f64>, mat_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    let mut i5 : usize;
    let mut theta : f64;
    let mut sth : f64;
    let mut cth : f64;
    let mut a1 : f64;
    let mut a2 : f64;
    
    i3 = 0;
    for i1 in 0..mat_dim {
        for i2 in 0..mat_dim {
            if i2 == i1 {
                q_mat[i3] = 1.0;
            } else {
                q_mat[i3] = 0.0;
            }
            i3 += 1usize;
        }
    }
    
    for i1 in 0..mat_dim {
        for i2 in i1+2..mat_dim {
            i3 = (i1+1)*mat_dim + i1;
            if fabs(mat[i3]) < TOL {
                theta = PI_2;
            } else {
                i4 = i2*mat_dim + i1;
                theta = atan(mat[i4]/mat[i3]);
            }
            sth = sin(theta);
            cth = cos(theta);
            i4 = (i1+1)*mat_dim + i1;
            i5 = i2*mat_dim + i1;
            for _i3 in i1..mat_dim {
                a1 = cth*mat[i4] + sth*mat[i5];
                a2 = -sth*mat[i4] + cth*mat[i5];
                mat[i4] = a1;
                mat[i5] = a2;
                i4 += 1usize;
                i5 += 1usize;
            }
            i4 = i1 + 1;
            i5 = i2;
            for _i3 in 0..mat_dim {
                a1 = cth*mat[i4] + sth*mat[i5];
                a2 = -sth*mat[i4] + cth*mat[i5];
                mat[i4] = a1;
                mat[i5] = a2;
                a1 = cth*q_mat[i4] + sth*q_mat[i5];
                a2 = -sth*q_mat[i4] + cth*q_mat[i5];
                q_mat[i4] = a1;
                q_mat[i5] = a2;
                i4 +=  mat_dim;
                i5 +=  mat_dim;
            }
        }
    }
    
    return;
}

pub fn get_char_fun(c_fun : &mut DiffDoub1, mat : &mut Vec<DiffDoub1>, mat_dim : usize, e_vals : &mut Vec<f64>, lam : f64, tri_diag : usize) {
    let mut i2 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta = DiffDoub1::new();
    let mut sth = DiffDoub1::new();
    let mut cth = DiffDoub1::new();
    let mut p1 = DiffDoub1::new();
    let mut p2 = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    
    for i1 in 0..=(mat_dim - 1) {
        i2_min = i1 + 1;
        if tri_diag == 0 {
            i2_max = mat_dim - 1;
        } else {
            i2_max = i1 + 1;
            if i2_max > (mat_dim - 1) {
                i2_max = mat_dim - 1;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min-1)*mat_dim + i1;
            k12 = i2*mat_dim + i1;
            if fabs(mat[k11].val) < TOL {
                mat[k11].val = TOL;
            }
            theta.set_val_dfd1(&mut mat[k12]);
            tmp.set_val_dfd1(&mut mat[k11]);
            theta.dvd(& tmp);
            theta.atn();
            sth.set_val_dfd1(&mut theta);
            sth.sn();
            cth.set_val_dfd1(&mut theta);
            cth.cs();
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > (mat_dim - 1) {
                    i3_max = mat_dim - 1;
                }
            } else {
                i3_max = mat_dim - 1;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min-1)*mat_dim + i3;
                k23 = i2*mat_dim + i3;
                p1.set_val_dfd1(&mut cth);
                p1.mult(& mat[k22]);
                tmp.set_val_dfd1(&mut sth);
                tmp.mult(& mat[k23]);
                p1.add(& tmp);
                p2.set_val_dfd1(&mut sth);
                p2.neg();
                p2.mult(& mat[k22]);
                tmp.set_val_dfd1(&mut cth);
                tmp.mult(& mat[k23]);
                p2.add(& tmp);
                mat[k22].set_val_dfd1(&mut p1);
                mat[k23].set_val_dfd1(&mut p2);
            }
            mat[k12].set_val_dfd1(&mut theta);
        }
    }
    
    c_fun.set_val_dfd1(&mut mat[0]);
    i2 = mat_dim + 1;
    for _i1 in 1..mat_dim {
        c_fun.mult(& mat[i2]);
        while fabs(c_fun.val) > MAX_MAG {
            c_fun.val *=  MIN_MAG;
            c_fun.dval *=  MIN_MAG;
        }
        while fabs(c_fun.val) < MIN_MAG {
            c_fun.val *=  MAX_MAG;
            c_fun.dval *=  MAX_MAG;
        }
        i2 +=  mat_dim + 1;
    }
    
    let mut term : f64;
    for ev in e_vals.iter_mut() {
        term = *ev - lam;
        tmp.set_val_2(term,   -1.0);
        c_fun.dvd(& tmp);
        while fabs(c_fun.val) > MAX_MAG {
            c_fun.val *=  MIN_MAG;
            c_fun.dval *=  MIN_MAG;
        }
        while fabs(c_fun.val) < MIN_MAG {
            c_fun.val *=  MAX_MAG;
            c_fun.dval *=  MAX_MAG;
        }
    }
    
    return;
}

pub fn get_evals(e_vals : &mut Vec<f64>, mat : &mut Vec<f64>, mat_dim : usize, lam_init : f64, conv_tol : f64, tri_diag : usize) {
    let mut i1 : usize;
    let mut i3 : usize;
    let mat_size : usize =  mat_dim*mat_dim;
    let mut e_val_list : Vec<f64> = Vec::new();
    let mut mat_copy = vec![DiffDoub1::new(); mat_size];
    let mut c_fun = DiffDoub1::new();
    let mut lam : f64 =  lam_init;
    let mut d_lam : f64;
    let mut d_trm : f64;
    let back_step : f64 =  1.0e+6*conv_tol;
    let mut num_found : usize =  0;
    let max_it : usize =  100*mat_dim;
    let mut it : usize =  0;
    while num_found < mat_dim && it < max_it {
        i3 = 0;
        for i1 in 0..mat_dim {
            for i2 in 0..mat_dim {
                if i1 == i2 {
                    d_trm = mat[i3] - lam;
                    mat_copy[i3].set_val_2(d_trm,   -1.0);
                } else {
                    mat_copy[i3].set_val(mat[i3]);
                }
                i3 += 1usize;
            }
        }
        get_char_fun(&mut c_fun, &mut mat_copy, mat_dim, &mut e_val_list, lam, tri_diag);
        d_lam = -c_fun.val/c_fun.dval;
        lam +=  d_lam;
        if fabs(d_lam) < conv_tol {
            //e_val_list.add_entry(lam);
            e_val_list.push(lam);
            num_found += 1usize;
            lam -=  back_step;
        }
        it += 1usize;
    }
    //doub_list_ent *this_ent = e_val_list.get_first();
    i1 = 0;
    for ev in e_val_list.iter_mut() {
        e_vals[i1] = *ev;
        i1 += 1usize;
    }
    return;
}

pub fn get_low_eval(mat : &mut Vec<f64>, mat_dim : usize) -> f64 {
    let mut vec = vec![0f64; mat_dim];
    let mut prev_vec = vec![0f64; mat_dim];
    let high : f64;
    let mut tmp : f64;
    let mut mag : f64;
    let mut dp : f64 = 0.0;
    let mut i3 : usize;
    let mut it : usize;
    let max_it : usize;
    
    mag = 0.0;
    for i1 in 0..mat_dim {
        tmp = sin(1.0*(i1 as f64));
        prev_vec[i1] = tmp;
        mag +=  tmp*tmp;
    }
    mag = 1.0/sqrt(mag);
    for i1 in 0..mat_dim {
        prev_vec[i1] *=  mag;
    }
    
    let mut conv : bool =  false;
    max_it = 100*mat_dim;
    it = 0;
    while !conv && it < max_it {
        i3 = 0;
        mag = 0.0;
        for i1 in 0..mat_dim {
            vec[i1] = 0.0;
            for i2 in 0..mat_dim {
                vec[i1] +=  mat[i3]*prev_vec[i2];
                i3 += 1usize;
            }
            mag +=  vec[i1]*vec[i1];
        }
        mag = sqrt(mag);
        dp = 0.0;
        for i1 in 0..mat_dim {
            dp +=  vec[i1]*prev_vec[i1];
        }
        if fabs(dp) > 0.99999999*mag {
            conv = true;
        } else {
            mag = 1.0/mag;
            for i1 in 0..mat_dim {
                prev_vec[i1] = mag*vec[i1];
            }
        }
        it += 1usize;
    }
    
    if dp < 0.0 {
        return  dp;
    }
    
    high = dp;
    
    i3 = 0;
    for _i1 in 0..mat_dim {
        mat[i3] -=  high;
        i3 +=  mat_dim + 1;
    }
    
    mag = 0.0;
    for i1 in 0..mat_dim {
        tmp = sin(1.0*(i1 as f64));
        prev_vec[i1] = tmp;
        mag +=  tmp*tmp;
    }
    mag = 1.0/sqrt(mag);
    for i1 in 0..mat_dim {
        prev_vec[i1] *=  mag;
    }
    
    conv = false;
    it = 0;
    while !conv && it < max_it {
        i3 = 0;
        mag = 0.0;
        for i1 in 0..mat_dim {
            vec[i1] = 0.0;
            for i2 in 0..mat_dim {
                vec[i1] +=  mat[i3]*prev_vec[i2];
                i3 += 1usize;
            }
            mag +=  vec[i1]*vec[i1];
        }
        mag = sqrt(mag);
        dp = 0.0;
        for i1 in 0..mat_dim {
            dp +=  vec[i1]*prev_vec[i1];
        }
        if fabs(dp) > 0.99999999*mag {
            conv = true;
        } else {
            mag = 1.0/mag;
            for i1 in 0..mat_dim {
                prev_vec[i1] = mag*vec[i1];
            }
        }
        it += 1usize;
    }
    
    i3 = 0;
    for _i1 in 0..mat_dim {
        mat[i3] +=  high;
        i3 +=  mat_dim + 1;
    }
    
    dp +=  high;
    return  dp;
}

pub fn eigen_solve(e_vals : &mut Vec<f64>, e_vecs : &mut Vec<f64>, mat : &mut Vec<f64>, mat_dim : usize, tri_diag : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    let mut tmp : f64;
    let mat_size : usize =  mat_dim*mat_dim;
    let mut new_td : usize =  tri_diag;
    let mut q_mat = vec![0f64; 1];
    if tri_diag == 0 {
        q_mat = vec![0f64; mat_size];
        sym_factor(mat, &mut q_mat, mat_dim);
        new_td = 1;
    }
    
    let mut low : f64 =  get_low_eval(mat, mat_dim);
    
    let mut mat_mag : f64 =  0.0;
    for i1 in 0..mat_size {
        tmp = fabs(mat[i1]);
        if tmp > mat_mag {
            mat_mag = tmp;
        }
    }
    
    low -=  0.001*mat_mag;
    let conv_tol : f64 =  1.0e-12*mat_mag;
    get_evals(e_vals, mat, mat_dim, low, conv_tol, new_td);
    
    //doub_list_ent *this_val = e_vals.get_first();
    let mut mat_copy = vec![0f64; mat_size];
    let mut vec = vec![0f64; mat_dim];
    let mut prev_vec = vec![0f64; mat_dim];
    let mut b_vec = vec![0f64; mat_dim];
    let mut shift : f64;
    let mut conv : bool;
    let mut mag : f64;
    let mut dp : f64;
    let mut it : usize;
    let max_it : usize =  100*mat_dim;
    for i5 in 0..mat_dim {
        shift = e_vals[i5] + 1.0e-8*mat_mag;
        i3 = 0;
        for i1 in 0..mat_dim {
            for i2 in 0..mat_dim {
                if i1 == i2 {
                    mat_copy[i3] = mat[i3] - shift;
                }
                else {
                    mat_copy[i3] = mat[i3];
                }
                i3 += 1usize;
            }
        }
        q_rfactor(&mut mat_copy,  mat_dim,  0,  mat_dim - 1,  0,  mat_dim - 1,  new_td);
        mag = 0.0;
        for i1 in 0..mat_dim {
            tmp = sin((i1 * (i5 + 1)) as f64);
            prev_vec[i1] = tmp;
            mag  +=  tmp * tmp;
        }
        mag = 1.0 / sqrt(mag);
        for i1 in 0..mat_dim {
            prev_vec[i1] = mag * prev_vec[i1];
        }
        conv = false;
        it = 0;
        while !conv && it < max_it {
            for i1 in 0..mat_dim {
                b_vec[i1] = prev_vec[i1];
            }
            solveq_rx_eqb(&mut vec, &mut  mat_copy, &mut  b_vec,  mat_dim,  0,  mat_dim - 1,  0,  mat_dim - 1,  new_td);
            mag = 0.0;
            dp = 0.0;
            for i1 in 0..mat_dim {
                mag  +=  vec[i1] * vec[i1];
                dp  +=  vec[i1] * prev_vec[i1];
            }
            mag = sqrt(mag);
            if fabs(dp) > 0.999999 * mag {
                conv = true;
                if tri_diag == 0 {
                    i3 = 0;
                    //i4 = i5;
                    i4 = i5 * mat_dim;
                    for _i1 in 0..mat_dim {
                        //i4 = i5*mat_dim + i1
                        e_vecs[i4] = 0.0;
                        for i2 in 0..mat_dim {
                            e_vecs[i4]  +=  q_mat[i3] * prev_vec[i2];
                            i3 += 1usize;
                        }
                        //i4 += mat_dim;
                        i4 += 1usize;
                    }
                }
                else {
                    //i4 = i5;
                    i4 = i5 * mat_dim;
                    for i1 in 0..mat_dim {
                        //i4 = i5*mat_dim + i1;
                        e_vecs[i4] = prev_vec[i1];
                        //i4 += mat_dim;
                        i4 += 1usize;
                    }
                }
            }
            else {
                mag = 1.0 / mag;
                for i1 in 0..mat_dim {
                    prev_vec[i1] = mag * vec[i1];
                }
            }
            it += 1usize;
        }
    }
    return;
}

pub fn sym_eigen_solve(e_vals : &mut Vec<f64>, e_vecs : &mut Vec<f64>, mat : &mut Vec<f64>, mat_dim : usize, tri_diag : usize) {
    let mut i2 : usize;
    let mut i3 : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut i4 : usize;
    let mut i5 : usize;
    let mut loop_ct : usize;
    let mut tmp : f64;
    let mut mag : f64;
    let mut dp : f64;
    
    let mut temp_v1 = vec![0f64; mat_dim];
    let mut temp_v2 = vec![0f64; mat_dim];
    let mut e_vtmp = vec![0f64; mat_dim*mat_dim];
    let mut q_mat = vec![0f64; mat_dim];
    let mut srt_order = vec![0usize; mat_dim];
    
    if tri_diag == 0 {
        sym_factor(mat, &mut q_mat, mat_dim);
    }
    
    for i1 in 0..mat_dim {
        for i2 in 0..mat_dim {
            tmp = sin((i2 * (i1+1)) as f64);
            temp_v1[i2] = tmp;
        }
        for i2 in 0..i1 {
            dp = 0.0;
            i4 = i2 * mat_dim;
            for i3 in 0..mat_dim {
                dp  +=  temp_v1[i3] * e_vtmp[i4];
                i4 += 1usize;
            }
            i4 = i2 * mat_dim;
            for i3 in 0..mat_dim {
                temp_v1[i3]  -=  dp*e_vtmp[i4];
                i4 += 1usize;
            }
        }
        mag = 0.0;
        for i2 in 0..mat_dim {
            mag  +=  temp_v1[i2] * temp_v1[i2];
        }
        mag = 1.0 / sqrt(mag);
        for i2 in 0..mat_dim {
            temp_v1[i2]  *=  mag;
        }
        dp = 0.0;
        loop_ct = 0;
        while fabs(dp) < 0.99999999*mag && loop_ct < 10000 {
            // multiply by mat
            for i2 in 0..mat_dim {
                temp_v2[i2] = 0.0;
                if tri_diag == 1 {
                    if i2 == 0 {
                        i3_min = 0;
                    }
                    else {
                        i3_min = i2 - 1;
                    }
                    i3_max = i2 + 2;
                    if i3_max > mat_dim {
                        i3_max = mat_dim;
                    }
                }
                else {
                    i3_min = 0;
                    i3_max = mat_dim;
                }
                i4 = i2 * mat_dim + i3_min;
                for i3 in i3_min..i3_max {
                    //i4 = i2 * mat_dim + i3;
                    temp_v2[i2]  +=  mat[i4] * temp_v1[i3];
                    i4 += 1usize;
                }
            }
            // orthogonalize with previous std::vectors
            for i2 in 0..i1 {
                dp = 0.0;
                i4 = i2 * mat_dim;
                for i3 in 0..mat_dim {
                    dp  +=  temp_v2[i3] * e_vtmp[i4];
                    i4 += 1usize;
                }
                i4 = i2 * mat_dim;
                for i3 in 0..mat_dim {
                    temp_v2[i3]  -=  dp * e_vtmp[i4];
                    i4 += 1usize;
                }
            }
            // get mag, dp
            mag = 0.0;
            dp = 0.0;
            for i2 in 0..mat_dim {
                mag  +=  temp_v2[i2] * temp_v2[i2];
                dp  +=  temp_v1[i2] * temp_v2[i2];
            }
            mag = sqrt(mag);
            // update temp_v1
            tmp = 1.0 / mag;
            for i2 in 0..mat_dim {
                temp_v1[i2] = tmp * temp_v2[i2];
            }
            loop_ct += 1usize;
        }
        e_vals[i1] = dp;
        if tri_diag == 0 {
            i4 = 0;
            for i2 in 0..mat_dim {
                i5 = i1 * mat_dim + i2;
                e_vtmp[i5] = 0.0;
                for i3 in 0..mat_dim {
                    e_vtmp[i5]  +=  q_mat[i4] * temp_v1[i3];
                    i4 += 1usize;
                }
            }
        }
        else {
            i3 = i1 * mat_dim;
            for i2 in 0..mat_dim {
                e_vtmp[i3] = temp_v1[i2];
                i3 += 1usize;
            }
        }
    }
    
    // sort pairs
    
    for i1 in 0..mat_dim {
        srt_order[i1] = i1;
    }
    
    for _i1 in 0..mat_dim {
        for i2 in 0..(mat_dim - 1) {
            if e_vals[i2 + 1] < e_vals[i2] {
                tmp = e_vals[i2];
                e_vals[i2] = e_vals[i2 + 1];
                e_vals[i2 + 1] = tmp;
                i3 = srt_order[i2];
                srt_order[i2] = srt_order[i2 + 1];
                srt_order[i2 + 1] = i3;
            }
        }
    }
    
    for i1 in 0..mat_dim {
        i2 = mat_dim * i1;
        i3 = mat_dim * srt_order[i1];
        for _i4 in 0..mat_dim {
            e_vecs[i2] = e_vtmp[i3];
            i2 += 1usize;
            i3 += 1usize;
        }
    }
    
    return;
}

pub fn eigen_sparse_direct(e_vals : &mut Vec<f64>, e_vecs : &mut Vec<f64>, num_pairs : usize, mat : &mut LowerTriMat, mass_mat : &mut Vec<f64>, mat_dim : usize) {
    let mut i1 : usize;
    let mut i2 : usize;
    let mut i3 : usize;
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    
    let h_cols : usize =  5*num_pairs;
    i1 = h_cols*mat_dim;
    let mut h_mat = vec![0f64; i1];
    let mut t_vec1 = vec![0f64; mat_dim];
    let mut t_vec2 = vec![0f64; mat_dim];
    i1 = h_cols*(h_cols - 1);
    let mut coef_mat = vec![0f64; i1];
    
    for i2 in 0..i1 {
        coef_mat[i2] = 0.0;
    }
    
    let mut mag : f64;
    let mut tmp : f64;
    for i2 in 0..mat_dim {
        mass_mat[i2] = sqrt(mass_mat[i2]);
    }
    
    for i2 in mat_dim..(2 * mat_dim) {
        h_mat[i2] = sin(i2 as f64);
    }
    
    sub_vec(&mut t_vec1, &mut  h_mat,  0,  mat_dim);
    sub_vec(&mut t_vec2, &mut  h_mat,  mat_dim,  mat_dim + mat_dim);
    mat.ldl_solve(&mut t_vec1, &mut  t_vec2);
    return_sv(&mut t_vec1, &mut h_mat, 0, mat_dim);
    
    mag = 0.0;
    for i2 in 0..mat_dim {
        mag  +=  h_mat[i2] * h_mat[i2];
    }
    mag = 1.0 / sqrt(mag);
    
    for i2 in 0..mat_dim {
        h_mat[i2]  *=  mag;
    }
    
    let mut dp : f64;
    let mut st_vec : usize;
    for i1 in 1..h_cols {
        i3 = mat_dim*(i1-1);
        for i2 in 0..mat_dim {
            t_vec1[i2] = mass_mat[i2]*h_mat[i3];
            i3 += 1usize;
        }
        mat.ldl_solve(&mut t_vec2, &mut t_vec1);
        for i2 in 0..mat_dim {
            t_vec1[i2] = t_vec2[i2]*mass_mat[i2];
        }
        if i1 == 1 {
            st_vec = 0;
        } else {
            st_vec = i1 - 2;
        }
        for i2 in st_vec..i1 {
            dp = 0.0;
            i3 = mat_dim*i2;
            for i4 in 0..mat_dim {
                dp +=  t_vec1[i4]*h_mat[i3];
                i3 += 1usize;
            }
            i3 = mat_dim*i2;
            for i4 in 0..mat_dim {
                t_vec1[i4] -=  dp*h_mat[i3];
                i3 += 1usize;
            }
            i3 = (h_cols-1)*i2 + (i1-1);
            coef_mat[i3] = dp;
        }
        mag = 0.0;
        for i2 in 0..mat_dim {
            mag +=  t_vec1[i2]*t_vec1[i2];
        }
        mag = sqrt(mag);
        i3 = h_cols*i1 - 1;// (h_cols-1)*i1 + (i1 - 1)
        coef_mat[i3] = mag;
        mag = 1.0/mag;
        i3 = mat_dim*i1;
        for i2 in 0..mat_dim {
            h_mat[i3] = mag*t_vec1[i2];
            i3 += 1usize;
        }
    }
    
    i1 = h_cols - 1;
    let mut coef_vals = vec![0f64; i1];
    let mut coef_vecs = vec![0f64; i1*i1];
    //eigen_solve(coef_vals, coef_vecs, coef_mat, i1, 2);
    sym_eigen_solve(&mut coef_vals,  &mut  coef_vecs,  &mut  coef_mat,   i1,   1);
    
    for i3 in 0..mat_dim {
        mass_mat[i3] = 1.0/mass_mat[i3];
    }
    
    i1 = 0;
    i2 = h_cols - 2;
    while i1 < num_pairs && i2 > 0 {
        for i3 in 0..mat_dim {
            t_vec1[i3] = 0.0;
        }
        i5 = 0;
        for i3 in 0..(h_cols-1) {
            //i6 = (h_cols-1)*i3 + i2;
            i6 = i2 * (h_cols - 1) + i3;
            for i4 in 0..mat_dim {
                t_vec1[i4]  +=  h_mat[i5] * coef_vecs[i6];
                i5 += 1usize;
            }
        }
        mag = 0.0;
        for i3 in 0..mat_dim {
            tmp = t_vec1[i3]*mass_mat[i3];
            t_vec2[i3] = tmp;
            mag +=  tmp*tmp;
        }
        mag = 1.0/sqrt(mag);
        for i3 in 0..mat_dim {
            t_vec2[i3]  *=  mag;
        }
        if i1 > 0 {
            tmp = 0.0;
            i4 = mat_dim * (i1 - 1);
            for i3 in 0..mat_dim {
                tmp  +=  e_vecs[i4] * t_vec2[i3];
                i4 += 1usize;
            }
            if fabs(tmp) < 0.9 {
                for i3 in 0..mat_dim {
                    e_vecs[i4] = t_vec2[i3];
                    i4 += 1usize;
                }
                e_vals[i1] = 1.0 / coef_vals[i2];
                i1 += 1usize;
                i2 -= 1usize;
            }
            else {
                i2 -= 1usize;
            }
        }
        else {
            i4 = 0;
            for i3 in 0..mat_dim {
                //i4 = i1 * mat_dim + i3;
                e_vecs[i4] = t_vec2[i3];
                i4 += 1usize;
            }
            e_vals[i1] = 1.0 / coef_vals[i2];
            i1 += 1usize;
            i2 -= 1usize;
        }
    }
    
    for i1 in 0..mat_dim {
        tmp = 1.0/mass_mat[i1];
        mass_mat[i1] = tmp*tmp;
    }
    
    return;
}

pub fn ray_quot(grad : &mut Vec<f64>, kv : &mut Vec<f64>, mv : &mut Vec<f64>, mat : &mut SparseMat, cnst : &mut ConstraintList, mass_mat : &mut Vec<f64>, in_vec : &mut Vec<f64>) -> f64 {
    let r_c : f64;
    let dim : usize =  mat.dim;
    let mut v_kv : f64;
    let mut v_mv : f64 =  0.0;
    for i1 in 0..dim {
        kv[i1] = 0.0;
        mv[i1] = mass_mat[i1] * in_vec[i1];
        v_mv  +=  in_vec[i1] * mv[i1];
    }
    mat.vector_multiply(kv, in_vec, false);
    cnst.get_total_vec_mult(kv, in_vec, grad);
    v_kv = 0.0;
    for i1 in 0..dim {
        v_kv  +=  in_vec[i1] * kv[i1];
    }
    r_c = v_kv / v_mv;
    for i1 in 0..dim {
        grad[i1] = 2.0 * ((kv[i1] / v_mv) - r_c * (mv[i1] / v_mv));
    }
    return  r_c;
}

pub fn unitize_vec(vec : &mut Vec<f64>, dim : usize) -> f64 {
    let mut mag : f64;
    let minv : f64;
    mag = 0.0;
    for i1 in 0..dim {
        mag  +=  vec[i1] * vec[i1];
    }
    mag = sqrt(mag);
    minv = 1.0 / mag;
    for i1 in 0..dim {
        vec[i1]  *=  minv;
    }
    return  mag;
}

pub fn get_nearest_evec_rq(mat : &mut SparseMat, cnst : &mut ConstraintList, mass_mat : &mut Vec<f64>, in_vecs : &mut Vec<f64>, e_vals : &mut Vec<f64>, num_vecs : usize, max_it : usize) {
    let mut i2 : usize;
    let mut i3 : usize;
    let dim : usize;
    let mut dp : f64;
    let mut res : f64;
    let mut r_q0 : f64;
    let mut _r_q1 : f64;
    
    dim = mat.dim;
    
    let mut grad0 = vec![0f64; dim];
    let mut grad1 = vec![0f64; dim];
    let mut d2_rq = vec![0f64; dim];
    let mut v_step = vec![0f64; dim];
    let mut kv = vec![0f64; dim];
    let mut mv = vec![0f64; dim];
    let mut t_vec1 = vec![0f64; dim];
    
    for i1 in 0..num_vecs {
        i2 = i1 * dim;
        sub_vec(&mut t_vec1, in_vecs,  i2,  i2 + dim);
        unitize_vec(&mut t_vec1,  dim);
        r_q0 = ray_quot(&mut grad0, &mut  kv, &mut  mv, mat, cnst,  mass_mat, &mut  t_vec1);
        return_sv(&mut t_vec1, in_vecs, i2, i2 + dim);
        res = 1.0;
        i3 = 0;
        while i3 < max_it && res > 1.0e-6 {
            for i4 in 0..dim {
                v_step[i4] = -grad0[i4];
                if v_step[i4] == 0.0 {
                    v_step[i4] = 1.0e-12;
                }
            }
            unitize_vec(&mut v_step,  dim);
            for i4 in 0..dim {
                in_vecs[i2 + i4]  +=  0.01*v_step[i4];
            }
            sub_vec(&mut t_vec1, in_vecs,  i2,  i2 + dim);
            _r_q1 = ray_quot(&mut grad1, &mut  kv, &mut  mv, mat, cnst,  mass_mat, &mut  t_vec1);
            for i4 in 0..dim {
                in_vecs[i2 + i4]  -=  0.01 * v_step[i4];
            }
            for i4 in 0..dim {
                d2_rq[i4] = (grad1[i4] - grad0[i4]) / (0.01 * v_step[i4]);
                if d2_rq[i4] != 0.0 {
                    v_step[i4] = -grad0[i4] / d2_rq[i4];
                }
                else {
                    v_step[i4] = 0.0;
                }
                in_vecs[i2 + i4]  +=  v_step[i4];
            }
            sub_vec(&mut t_vec1, in_vecs,  i2,  i2 + dim);
            unitize_vec(&mut t_vec1,  dim);
            r_q0 = ray_quot(&mut grad0, &mut  kv, &mut  mv, mat,  cnst, mass_mat, &mut  t_vec1);
            return_sv(&mut t_vec1, in_vecs, i2, i2 + dim);
            unitize_vec(&mut kv,  dim);
            unitize_vec(&mut mv,  dim);
            dp = 0.0;
            for i4 in 0..dim {
                dp  +=  kv[i4] * mv[i4];
            }
            res = 1.0 - fabs(dp);
            i3 += 1usize;
        }
        println!("{}{}{}", "Warning: eigenstd::vector " , i1 , " did not converge within the max iterations." );
        e_vals[i1] = r_q0;
    }
    
}

//dup1
pub fn q_rfactor_dfd0(mat : &mut Vec<DiffDoub0>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta = DiffDoub0::new();
    let mut sth = DiffDoub0::new();
    let mut cth = DiffDoub0::new();
    let mut p1 = DiffDoub0::new();
    let mut p2 = DiffDoub0::new();
    let mut tmp = DiffDoub0::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min-1)*col_dim + i1;
            k12 = i2*col_dim + i1;
            if fabs(mat[k11].val) < TOL {
                mat[k11].val = TOL;
            }
            theta.set_val_dfd0(& mat[k12]);
            tmp.set_val_dfd0(& mat[k11]);
            theta.dvd(& tmp);
            theta.atn();
            sth.set_val_dfd0(& theta);
            sth.sn();
            cth.set_val_dfd0(& theta);
            cth.cs();
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > end_col {
                    i3_max = end_col;
                }
            } else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min-1)*col_dim + i3;
                k23 = i2*col_dim + i3;
                p1.set_val_dfd0(& cth);
                p1.mult(& mat[k22]);
                tmp.set_val_dfd0(& sth);
                tmp.mult(& mat[k23]);
                p1.add(& tmp);
                p2.set_val_dfd0(& sth);
                p2.neg();
                p2.mult(& mat[k22]);
                tmp.set_val_dfd0(& cth);
                tmp.mult(& mat[k23]);
                p2.add(& tmp);
                mat[k22].set_val_dfd0(& p1);
                mat[k23].set_val_dfd0(& p2);
            }
            mat[k12].set_val_dfd0(& theta);
        }
    }
    return;
}

pub fn q_rfactor_ar_dfd0(mat : &mut [DiffDoub0], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta = DiffDoub0::new();
    let mut sth = DiffDoub0::new();
    let mut cth = DiffDoub0::new();
    let mut p1 = DiffDoub0::new();
    let mut p2 = DiffDoub0::new();
    let mut tmp = DiffDoub0::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min - 1) * col_dim + i1;
            k12 = i2 * col_dim + i1;
            if fabs(mat[k11].val) < TOL {
                mat[k11].val = TOL;
            }
            theta.set_val_dfd0(& mat[k12]);
            tmp.set_val_dfd0(& mat[k11]);
            theta.dvd(& tmp);
            theta.atn();
            sth.set_val_dfd0(& theta);
            sth.sn();
            cth.set_val_dfd0(& theta);
            cth.cs();
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > end_col {
                    i3_max = end_col;
                }
            }
            else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min - 1) * col_dim + i3;
                k23 = i2 * col_dim + i3;
                p1.set_val_dfd0(& cth);
                p1.mult(& mat[k22]);
                tmp.set_val_dfd0(& sth);
                tmp.mult(& mat[k23]);
                p1.add(& tmp);
                p2.set_val_dfd0(& sth);
                p2.neg();
                p2.mult(& mat[k22]);
                tmp.set_val_dfd0(& cth);
                tmp.mult(& mat[k23]);
                p2.add(& tmp);
                mat[k22].set_val_dfd0(& p1);
                mat[k23].set_val_dfd0(& p2);
            }
            mat[k12].set_val_dfd0(& theta);
        }
    }
    return;
}

pub fn solveq_rx_eqb_dfd0(x_vec : &mut Vec<DiffDoub0>, mat : &mut Vec<DiffDoub0>, b_vec : &mut Vec<DiffDoub0>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i3 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut theta = DiffDoub0::new();
    let mut sth = DiffDoub0::new();
    let mut cth = DiffDoub0::new();
    let mut p1 = DiffDoub0::new();
    let mut p2 = DiffDoub0::new();
    let mut tmp = DiffDoub0::new();
    let mut row_sum = DiffDoub0::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2*col_dim + i1;
            theta.set_val_dfd0(& mat[k12]);
            sth.set_val_dfd0(& theta);
            sth.sn();
            cth.set_val_dfd0(& theta);
            cth.cs();
            p1.set_val_dfd0(& cth);
            p1.mult(& b_vec[i3]);
            tmp.set_val_dfd0(& sth);
            tmp.mult(& b_vec[i2]);
            p1.add(& tmp);
            p2.set_val_dfd0(& sth);
            p2.neg();
            p2.mult(& b_vec[i3]);
            tmp.set_val_dfd0(& cth);
            tmp.mult(& b_vec[i2]);
            p2.add(& tmp);
            b_vec[i3].set_val_dfd0(& p1);
            b_vec[i2].set_val_dfd0(& p2);
        }
        x_vec[i1].set_val(0.0);
    }
    
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if tri_diag == 2 {
            i2_max = i1 + 2;
            if i2_max > end_col {
                i2_max = end_col;
            }
        } else {
            i2_max = end_col;
        }
        row_sum.set_val(0.0);
        k11 = i3*col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            tmp.set_val_dfd0(& mat[k11]);
            tmp.mult(& x_vec[i2]);
            row_sum.add(& tmp);
            k11 += 1usize;
        }
        tmp.set_val_dfd0(& b_vec[i3]);
        tmp.sub(& row_sum);
        k11 = i3*col_dim + i1;
        tmp.dvd(& mat[k11]);
        x_vec[i1].set_val_dfd0(& tmp);
    }
    
    return;
}

pub fn solveq_rx_eqb_ar_dfd0(x_vec : &mut [DiffDoub0], mat : &mut [DiffDoub0], b_vec : &mut [DiffDoub0], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i3 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut theta = DiffDoub0::new();
    let mut sth = DiffDoub0::new();
    let mut cth = DiffDoub0::new();
    let mut p1 = DiffDoub0::new();
    let mut p2 = DiffDoub0::new();
    let mut tmp = DiffDoub0::new();
    let mut row_sum = DiffDoub0::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2 * col_dim + i1;
            theta.set_val_dfd0(& mat[k12]);
            sth.set_val_dfd0(& theta);
            sth.sn();
            cth.set_val_dfd0(& theta);
            cth.cs();
            p1.set_val_dfd0(& cth);
            p1.mult(& b_vec[i3]);
            tmp.set_val_dfd0(& sth);
            tmp.mult(& b_vec[i2]);
            p1.add(& tmp);
            p2.set_val_dfd0(& sth);
            p2.neg();
            p2.mult(& b_vec[i3]);
            tmp.set_val_dfd0(& cth);
            tmp.mult(& b_vec[i2]);
            p2.add(& tmp);
            b_vec[i3].set_val_dfd0(& p1);
            b_vec[i2].set_val_dfd0(& p2);
        }
        x_vec[i1].set_val(0.0);
    }
    
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if tri_diag == 2 {
            i2_max = i1 + 2;
            if i2_max > end_col {
                i2_max = end_col;
            }
        }
        else {
            i2_max = end_col;
        }
        row_sum.set_val(0.0);
        k11 = i3 * col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            tmp.set_val_dfd0(& mat[k11]);
            tmp.mult(& x_vec[i2]);
            row_sum.add(& tmp);
            k11 += 1usize;
        }
        tmp.set_val_dfd0(& b_vec[i3]);
        tmp.sub(& row_sum);
        k11 = i3 * col_dim + i1;
        tmp.dvd(& mat[k11]);
        x_vec[i1].set_val_dfd0(& tmp);
    }
    
    return;
}

pub fn get_det_inv_dfd0(det : &mut DiffDoub0, inv : &mut Vec<DiffDoub0>, mat : &mut Vec<DiffDoub0>, col_dim : usize, tri_diag : usize, x_vec : &mut Vec<DiffDoub0>, b_vec : &mut Vec<DiffDoub0>) {
    q_rfactor_dfd0(mat,  col_dim,  0,  col_dim-1,  0,  col_dim-1,  tri_diag);
    let mut i3 : usize;
    det.set_val(1.0);
    for i1 in 0..col_dim {
        for i2 in 0..col_dim {
            if i1 == i2 {
                b_vec[i2].set_val(1.0);
            } else {
                b_vec[i2].set_val(0.0);
            }
        }
        solveq_rx_eqb_dfd0(x_vec, mat, b_vec, col_dim, 0, col_dim-1, 0, col_dim-1, tri_diag);
        for i2 in 0..col_dim {
            i3 = i2*col_dim + i1;
            inv[i3].set_val_dfd0(& x_vec[i2]);
        }
        i3 = i1*col_dim + i1;
        det.mult(& mat[i3]);
    }
    return;
}

pub fn get_det_inv_ar_dfd0(det : &mut DiffDoub0, inv : &mut [DiffDoub0], mat : &mut [DiffDoub0], col_dim : usize, tri_diag : usize, x_vec : &mut [DiffDoub0], b_vec : &mut [DiffDoub0]) {
    q_rfactor_ar_dfd0(mat,  col_dim,  0,  col_dim - 1,  0,  col_dim - 1,  tri_diag);
    let mut i3 : usize;
    det.set_val(1.0);
    for i1 in 0..col_dim {
        for i2 in 0..col_dim {
            if i1 == i2 {
                b_vec[i2].set_val(1.0);
            }
            else {
                b_vec[i2].set_val(0.0);
            }
        }
        solveq_rx_eqb_ar_dfd0(x_vec, mat, b_vec, col_dim, 0, col_dim - 1, 0, col_dim - 1, tri_diag);
        for i2 in 0..col_dim {
            i3 = i2 * col_dim + i1;
            inv[i3].set_val_dfd0(& x_vec[i2]);
        }
        i3 = i1 * col_dim + i1;
        det.mult(& mat[i3]);
    }
    return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
pub fn q_rfactor_dfd1(mat : &mut Vec<DiffDoub1>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta = DiffDoub1::new();
    let mut sth = DiffDoub1::new();
    let mut cth = DiffDoub1::new();
    let mut p1 = DiffDoub1::new();
    let mut p2 = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min-1)*col_dim + i1;
            k12 = i2*col_dim + i1;
            if fabs(mat[k11].val) < TOL {
                mat[k11].val = TOL;
            }
            theta.set_val_dfd1(& mat[k12]);
            tmp.set_val_dfd1(& mat[k11]);
            theta.dvd(& tmp);
            theta.atn();
            sth.set_val_dfd1(& theta);
            sth.sn();
            cth.set_val_dfd1(& theta);
            cth.cs();
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > end_col {
                    i3_max = end_col;
                }
            } else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min-1)*col_dim + i3;
                k23 = i2*col_dim + i3;
                p1.set_val_dfd1(& cth);
                p1.mult(& mat[k22]);
                tmp.set_val_dfd1(& sth);
                tmp.mult(& mat[k23]);
                p1.add(& tmp);
                p2.set_val_dfd1(& sth);
                p2.neg();
                p2.mult(& mat[k22]);
                tmp.set_val_dfd1(& cth);
                tmp.mult(& mat[k23]);
                p2.add(& tmp);
                mat[k22].set_val_dfd1(& p1);
                mat[k23].set_val_dfd1(& p2);
            }
            mat[k12].set_val_dfd1(& theta);
        }
    }
    return;
}

pub fn q_rfactor_ar_dfd1(mat : &mut [DiffDoub1], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut i3_min : usize;
    let mut i3_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut k22 : usize;
    let mut k23 : usize;
    let mut theta = DiffDoub1::new();
    let mut sth = DiffDoub1::new();
    let mut cth = DiffDoub1::new();
    let mut p1 = DiffDoub1::new();
    let mut p2 = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        for i2 in i2_min..=i2_max {
            k11 = (i2_min - 1) * col_dim + i1;
            k12 = i2 * col_dim + i1;
            if fabs(mat[k11].val) < TOL {
                mat[k11].val = TOL;
            }
            theta.set_val_dfd1(& mat[k12]);
            tmp.set_val_dfd1(& mat[k11]);
            theta.dvd(& tmp);
            theta.atn();
            sth.set_val_dfd1(& theta);
            sth.sn();
            cth.set_val_dfd1(& theta);
            cth.cs();
            i3_min = i1;
            if tri_diag == 2 {
                i3_max = i1 + 2;
                if i3_max > end_col {
                    i3_max = end_col;
                }
            }
            else {
                i3_max = end_col;
            }
            for i3 in i3_min..=i3_max {
                k22 = (i2_min - 1) * col_dim + i3;
                k23 = i2 * col_dim + i3;
                p1.set_val_dfd1(& cth);
                p1.mult(& mat[k22]);
                tmp.set_val_dfd1(& sth);
                tmp.mult(& mat[k23]);
                p1.add(& tmp);
                p2.set_val_dfd1(& sth);
                p2.neg();
                p2.mult(& mat[k22]);
                tmp.set_val_dfd1(& cth);
                tmp.mult(& mat[k23]);
                p2.add(& tmp);
                mat[k22].set_val_dfd1(& p1);
                mat[k23].set_val_dfd1(& p2);
            }
            mat[k12].set_val_dfd1(& theta);
        }
    }
    return;
}

pub fn solveq_rx_eqb_dfd1(x_vec : &mut Vec<DiffDoub1>, mat : &mut Vec<DiffDoub1>, b_vec : &mut Vec<DiffDoub1>, col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i3 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut theta = DiffDoub1::new();
    let mut sth = DiffDoub1::new();
    let mut cth = DiffDoub1::new();
    let mut p1 = DiffDoub1::new();
    let mut p2 = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    let mut row_sum = DiffDoub1::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        } else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2*col_dim + i1;
            theta.set_val_dfd1(& mat[k12]);
            sth.set_val_dfd1(& theta);
            sth.sn();
            cth.set_val_dfd1(& theta);
            cth.cs();
            p1.set_val_dfd1(& cth);
            p1.mult(& b_vec[i3]);
            tmp.set_val_dfd1(& sth);
            tmp.mult(& b_vec[i2]);
            p1.add(& tmp);
            p2.set_val_dfd1(& sth);
            p2.neg();
            p2.mult(& b_vec[i3]);
            tmp.set_val_dfd1(& cth);
            tmp.mult(& b_vec[i2]);
            p2.add(& tmp);
            b_vec[i3].set_val_dfd1(& p1);
            b_vec[i2].set_val_dfd1(& p2);
        }
        x_vec[i1].set_val(0.0);
    }
    
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if tri_diag == 2 {
            i2_max = i1 + 2;
            if i2_max > end_col {
                i2_max = end_col;
            }
        } else {
            i2_max = end_col;
        }
        row_sum.set_val(0.0);
        k11 = i3*col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            tmp.set_val_dfd1(& mat[k11]);
            tmp.mult(& x_vec[i2]);
            row_sum.add(& tmp);
            k11 += 1usize;
        }
        tmp.set_val_dfd1(& b_vec[i3]);
        tmp.sub(& row_sum);
        k11 = i3*col_dim + i1;
        tmp.dvd(& mat[k11]);
        x_vec[i1].set_val_dfd1(& tmp);
    }
    
    return;
}

pub fn solveq_rx_eqb_ar_dfd1(x_vec : &mut [DiffDoub1], mat : &mut [DiffDoub1], b_vec : &mut [DiffDoub1], col_dim : usize, st_row : usize, end_row : usize, st_col : usize, end_col : usize, tri_diag : usize) {
    let mut i3 : usize;
    let mut i2_min : usize;
    let mut i2_max : usize;
    let mut k11 : usize;
    let mut k12 : usize;
    let mut theta = DiffDoub1::new();
    let mut sth = DiffDoub1::new();
    let mut cth = DiffDoub1::new();
    let mut p1 = DiffDoub1::new();
    let mut p2 = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    let mut row_sum = DiffDoub1::new();
    
    for i1 in st_col..=end_col {
        i2_min = st_row + (i1 - st_col) + 1;
        if tri_diag == 0 {
            i2_max = end_row;
        }
        else {
            i2_max = st_row + (i1 - st_col) + 1;
            if i2_max > end_row {
                i2_max = end_row;
            }
        }
        i3 = i2_min - 1;
        for i2 in i2_min..=i2_max {
            k12 = i2 * col_dim + i1;
            theta.set_val_dfd1(& mat[k12]);
            sth.set_val_dfd1(& theta);
            sth.sn();
            cth.set_val_dfd1(& theta);
            cth.cs();
            p1.set_val_dfd1(& cth);
            p1.mult(& b_vec[i3]);
            tmp.set_val_dfd1(& sth);
            tmp.mult(& b_vec[i2]);
            p1.add(& tmp);
            p2.set_val_dfd1(& sth);
            p2.neg();
            p2.mult(& b_vec[i3]);
            tmp.set_val_dfd1(& cth);
            tmp.mult(& b_vec[i2]);
            p2.add(& tmp);
            b_vec[i3].set_val_dfd1(& p1);
            b_vec[i2].set_val_dfd1(& p2);
        }
        x_vec[i1].set_val(0.0);
    }
    
    
    for i1 in (st_col..(end_col+1)).rev() {
        i3 = st_row + (i1 - st_col);
        i2_min = i1 + 1;
        if tri_diag == 2 {
            i2_max = i1 + 2;
            if i2_max > end_col {
                i2_max = end_col;
            }
        }
        else {
            i2_max = end_col;
        }
        row_sum.set_val(0.0);
        k11 = i3 * col_dim + i2_min;
        for i2 in i2_min..=i2_max {
            tmp.set_val_dfd1(& mat[k11]);
            tmp.mult(& x_vec[i2]);
            row_sum.add(& tmp);
            k11 += 1usize;
        }
        tmp.set_val_dfd1(& b_vec[i3]);
        tmp.sub(& row_sum);
        k11 = i3 * col_dim + i1;
        tmp.dvd(& mat[k11]);
        x_vec[i1].set_val_dfd1(& tmp);
    }
    
    return;
}

pub fn get_det_inv_dfd1(det : &mut DiffDoub1, inv : &mut Vec<DiffDoub1>, mat : &mut Vec<DiffDoub1>, col_dim : usize, tri_diag : usize, x_vec : &mut Vec<DiffDoub1>, b_vec : &mut Vec<DiffDoub1>) {
    q_rfactor_dfd1(mat,  col_dim,  0,  col_dim-1,  0,  col_dim-1,  tri_diag);
    let mut i3 : usize;
    det.set_val(1.0);
    for i1 in 0..col_dim {
        for i2 in 0..col_dim {
            if i1 == i2 {
                b_vec[i2].set_val(1.0);
            } else {
                b_vec[i2].set_val(0.0);
            }
        }
        solveq_rx_eqb_dfd1(x_vec, mat, b_vec, col_dim, 0, col_dim-1, 0, col_dim-1, tri_diag);
        for i2 in 0..col_dim {
            i3 = i2*col_dim + i1;
            inv[i3].set_val_dfd1(& x_vec[i2]);
        }
        i3 = i1*col_dim + i1;
        det.mult(& mat[i3]);
    }
    return;
}

pub fn get_det_inv_ar_dfd1(det : &mut DiffDoub1, inv : &mut [DiffDoub1], mat : &mut [DiffDoub1], col_dim : usize, tri_diag : usize, x_vec : &mut [DiffDoub1], b_vec : &mut [DiffDoub1]) {
    q_rfactor_ar_dfd1(mat,  col_dim,  0,  col_dim - 1,  0,  col_dim - 1,  tri_diag);
    let mut i3 : usize;
    det.set_val(1.0);
    for i1 in 0..col_dim {
        for i2 in 0..col_dim {
            if i1 == i2 {
                b_vec[i2].set_val(1.0);
            }
            else {
                b_vec[i2].set_val(0.0);
            }
        }
        solveq_rx_eqb_ar_dfd1(x_vec, mat, b_vec, col_dim, 0, col_dim - 1, 0, col_dim - 1, tri_diag);
        for i2 in 0..col_dim {
            i3 = i2 * col_dim + i1;
            inv[i3].set_val_dfd1(& x_vec[i2]);
        }
        i3 = i1 * col_dim + i1;
        det.mult(& mat[i3]);
    }
    return;
}

//end dup
 
//end skip 
 
 
//dup2

pub fn sub_vec_dfd0(sub_v : &mut Vec<DiffDoub0>, in_v : &mut Vec<DiffDoub0>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        sub_v[i2].set_val_dfd0(& in_v[i1]);
        i2 += 1usize;
    }
}

pub fn return_sv_dfd0(sub_v : &mut Vec<DiffDoub0>, in_v : &mut Vec<DiffDoub0>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        in_v[i1].set_val_dfd0(& sub_v[i2]);
        i2 += 1usize;
    }
}

pub fn vec_to_ar_dfd0(ar : &mut [DiffDoub0], vc : &mut Vec<DiffDoub0>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        ar[i2].set_val_dfd0(& vc[i1]);
        i2 += 1usize;
    }
}

pub fn ar_to_vec_dfd0(ar : &mut [DiffDoub0], vc : &mut Vec<DiffDoub0>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        vc[i1].set_val_dfd0(& ar[i2]);
        i2 += 1usize;
    }
}

pub fn mat_mul_dfd0(prod : &mut Vec<DiffDoub0>, mat1 : & Vec<DiffDoub0>, mat2 : & Vec<DiffDoub0>, m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub0::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd0(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6 +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn mat_mul_ar_dfd0(prod : &mut [DiffDoub0], mat1 : & [DiffDoub0], mat2 : & [DiffDoub0], m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub0::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd0(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6  +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn transpose_dfd0(mat_t : &mut Vec<DiffDoub0>, mat : &mut Vec<DiffDoub0>, row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd0(& mat[i3]);
            i3 += 1usize;
            i4 +=  row_dim;
        }
    }
    return;
}

pub fn transpose_ar_dfd0(mat_t : &mut [DiffDoub0], mat : &mut [DiffDoub0], row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd0(& mat[i3]);
            i3 += 1usize;
            i4  +=  row_dim;
        }
    }
    return;
}

pub fn cross_prod_dfd0(prod : &mut [DiffDoub0], v1 : & [DiffDoub0], v2 : & [DiffDoub0]) {
    let mut tmp = DiffDoub0::new();
    
    prod[0].set_val_dfd0(& v1[1]);
    prod[0].mult(& v2[2]);
    tmp.set_val_dfd0(& v1[2]);
    tmp.mult(& v2[1]);
    prod[0].sub(& tmp);
    
    prod[1].set_val_dfd0(& v1[2]);
    prod[1].mult(& v2[0]);
    tmp.set_val_dfd0(& v1[0]);
    tmp.mult(& v2[2]);
    prod[1].sub(& tmp);
    
    prod[2].set_val_dfd0(& v1[0]);
    prod[2].mult(& v2[1]);
    tmp.set_val_dfd0(& v1[1]);
    tmp.mult(& v2[0]);
    prod[2].sub(& tmp);
    return;
}

pub fn rotate_orient_dfd0(inst_ori : &mut [DiffDoub0], loc_ori : &mut [DiffDoub0], rot : &mut [DiffDoub0]) {
    let mut mag = DiffDoub0::new();
    let mut tmp = DiffDoub0::new();
    let mut tmp2 = DiffDoub0::new();
    let mut one_half = DiffDoub0::new();
    let mut a1 = [DiffDoub0::new(); 9];
    let mut a2 = [DiffDoub0::new(); 9];
    let mut a3 = [DiffDoub0::new(); 9];
    let mut tmp_v = [DiffDoub0::new(); 3];
    let mut i1 : usize;
    let mut i3 : usize;
    
    one_half.set_val(0.5);
    
    mag.set_val_dfd0(& rot[0]);
    mag.sqr();
    tmp.set_val_dfd0(& rot[1]);
    tmp.sqr();
    tmp2.set_val_dfd0(& rot[2]);
    tmp2.sqr();
    mag.add(& tmp);
    mag.add(& tmp2);
    mag.sqt();
    if mag.val < MAG_TOL {
        let mut loc_rot = [DiffDoub0::new(); 3];
        i3 = 0;
        for i1 in 0..3 {
            loc_rot[i1].set_val(0.0);
            for i2 in 0..3 {
                tmp.set_val_dfd0(& loc_ori[i3]);
                tmp.mult(& rot[i2]);
                loc_rot[i1].add(& tmp);
                i3 += 1usize;
            }
        }
        a1[0].set_val(1.0);
        tmp.set_val_dfd0(& loc_rot[1]);
        tmp.sqr();
        tmp2.set_val_dfd0(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[0].sub(& tmp);
        
        a1[1].set_val(0.5);
        a1[1].mult(& loc_rot[0]);
        a1[1].mult(& loc_rot[1]);
        a1[1].add(& loc_rot[2]);
        
        a1[2].set_val(0.5);
        a1[2].mult(& loc_rot[0]);
        a1[2].mult(& loc_rot[2]);
        a1[2].sub(& loc_rot[1]);
        
        a1[3].set_val(0.5);
        a1[3].mult(& loc_rot[0]);
        a1[3].mult(& loc_rot[1]);
        a1[3].sub(& loc_rot[2]);
        
        a1[4].set_val(1.0);
        tmp.set_val_dfd0(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd0(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[4].sub(& tmp);
        
        a1[5].set_val(0.5);
        a1[5].mult(& loc_rot[1]);
        a1[5].mult(& loc_rot[2]);
        a1[5].add(& loc_rot[0]);
        
        a1[6].set_val(0.5);
        a1[6].mult(& loc_rot[0]);
        a1[6].mult(& loc_rot[2]);
        a1[6].add(& loc_rot[1]);
        
        a1[7].set_val(0.5);
        a1[7].mult(& loc_rot[1]);
        a1[7].mult(& loc_rot[2]);
        a1[7].sub(& loc_rot[0]);
        
        a1[8].set_val(1.0);
        tmp.set_val_dfd0(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd0(& loc_rot[1]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[8].sub(& tmp);
        
        mat_mul_ar_dfd0(inst_ori, &mut  a1, loc_ori,  3,  3,  3);
    } else {
        let mut sth = DiffDoub0::new();
        let mut cth = DiffDoub0::new();
        
        sth.set_val_dfd0(& mag);
        sth.sn();
        cth.set_val_dfd0(& mag);
        cth.cs();
        
        let mut unit_rot = [DiffDoub0::new(); 3];
        tmp.set_val(1.0);
        tmp.dvd(& mag);
        unit_rot[0].set_val_dfd0(& rot[0]);
        unit_rot[0].mult(& tmp);
        unit_rot[1].set_val_dfd0(& rot[1]);
        unit_rot[1].mult(& tmp);
        unit_rot[2].set_val_dfd0(& rot[2]);
        unit_rot[2].mult(& tmp);
        
        a1[0].set_val_dfd0(& unit_rot[0]);
        a1[1].set_val_dfd0(& unit_rot[1]);
        a1[2].set_val_dfd0(& unit_rot[2]);
        
        i1 = 0;
        if fabs(unit_rot[1].val) < fabs(unit_rot[0].val) {
            i1 = 1;
        }
        if fabs(unit_rot[2].val) < fabs(unit_rot[i1].val) {
            i1 = 2;
        }
        tmp.set_val(1.0);
        tmp2.set_val_dfd0(& unit_rot[i1]);
        tmp2.sqr();
        tmp.sub(& tmp2);
        tmp.sqt();// tmp = sqrt(1 - unit_rot[i1]^2)
        a1[3+i1].set_val_dfd0(& tmp);
        for i2 in 0..3 {
            if i2 != i1 {
                i3 = 3 + i2;
                tmp2.set_val_dfd0(& a1[i1]);
                a1[i3].set_val_dfd0(& tmp2);
                a1[i3].neg();
                tmp2.set_val_dfd0(& a1[i2]);
                a1[i3].mult(& tmp2);
                a1[i3].dvd(& tmp);
            }
        }
        
        cross_prod_dfd0(&mut tmp_v, & a1[0..], & a1[3..]);
        for i in 0..3 {
            a1[i+6].set_val_dfd0(&tmp_v[i]);
        }
        
        a2[0].set_val_dfd0(& a1[0]);
        a2[1].set_val_dfd0(& a1[1]);
        a2[2].set_val_dfd0(& a1[2]);
        
        a2[3].set_val_dfd0(& cth);
        a2[3].mult(& a1[3]);
        tmp.set_val_dfd0(& sth);
        tmp.mult(& a1[6]);
        a2[3].sub(& tmp);
        
        a2[4].set_val_dfd0(& cth);
        a2[4].mult(& a1[4]);
        tmp.set_val_dfd0(& sth);
        tmp.mult(& a1[7]);
        a2[4].sub(& tmp);
        
        a2[5].set_val_dfd0(& cth);
        a2[5].mult(& a1[5]);
        tmp.set_val_dfd0(& sth);
        tmp.mult(& a1[8]);
        a2[5].sub(& tmp);
        
        a2[6].set_val_dfd0(& sth);
        a2[6].mult(& a1[3]);
        tmp.set_val_dfd0(& cth);
        tmp.mult(& a1[6]);
        a2[6].add(& tmp);
        
        a2[7].set_val_dfd0(& sth);
        a2[7].mult(& a1[4]);
        tmp.set_val_dfd0(& cth);
        tmp.mult(& a1[7]);
        a2[7].add(& tmp);
        
        a2[8].set_val_dfd0(& sth);
        a2[8].mult(& a1[5]);
        tmp.set_val_dfd0(& cth);
        tmp.mult(& a1[8]);
        a2[8].add(& tmp);
        
        transpose_ar_dfd0(&mut a3, &mut a2, 3, 3);// a3 = a2^t
        mat_mul_ar_dfd0(&mut a2, &mut a3, &mut a1, 3, 3, 3);//a2 = a2^t*a1
        mat_mul_ar_dfd0(inst_ori, loc_ori, &mut a2, 3, 3, 3);
    }
    return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup2

pub fn sub_vec_dfd1(sub_v : &mut Vec<DiffDoub1>, in_v : &mut Vec<DiffDoub1>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        sub_v[i2].set_val_dfd1(& in_v[i1]);
        i2 += 1usize;
    }
}

pub fn return_sv_dfd1(sub_v : &mut Vec<DiffDoub1>, in_v : &mut Vec<DiffDoub1>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        in_v[i1].set_val_dfd1(& sub_v[i2]);
        i2 += 1usize;
    }
}

pub fn vec_to_ar_dfd1(ar : &mut [DiffDoub1], vc : &mut Vec<DiffDoub1>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        ar[i2].set_val_dfd1(& vc[i1]);
        i2 += 1usize;
    }
}

pub fn ar_to_vec_dfd1(ar : &mut [DiffDoub1], vc : &mut Vec<DiffDoub1>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        vc[i1].set_val_dfd1(& ar[i2]);
        i2 += 1usize;
    }
}

pub fn mat_mul_dfd1(prod : &mut Vec<DiffDoub1>, mat1 : & Vec<DiffDoub1>, mat2 : & Vec<DiffDoub1>, m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub1::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd1(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6 +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn mat_mul_ar_dfd1(prod : &mut [DiffDoub1], mat1 : & [DiffDoub1], mat2 : & [DiffDoub1], m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub1::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd1(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6  +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn transpose_dfd1(mat_t : &mut Vec<DiffDoub1>, mat : &mut Vec<DiffDoub1>, row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd1(& mat[i3]);
            i3 += 1usize;
            i4 +=  row_dim;
        }
    }
    return;
}

pub fn transpose_ar_dfd1(mat_t : &mut [DiffDoub1], mat : &mut [DiffDoub1], row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd1(& mat[i3]);
            i3 += 1usize;
            i4  +=  row_dim;
        }
    }
    return;
}

pub fn cross_prod_dfd1(prod : &mut [DiffDoub1], v1 : & [DiffDoub1], v2 : & [DiffDoub1]) {
    let mut tmp = DiffDoub1::new();
    
    prod[0].set_val_dfd1(& v1[1]);
    prod[0].mult(& v2[2]);
    tmp.set_val_dfd1(& v1[2]);
    tmp.mult(& v2[1]);
    prod[0].sub(& tmp);
    
    prod[1].set_val_dfd1(& v1[2]);
    prod[1].mult(& v2[0]);
    tmp.set_val_dfd1(& v1[0]);
    tmp.mult(& v2[2]);
    prod[1].sub(& tmp);
    
    prod[2].set_val_dfd1(& v1[0]);
    prod[2].mult(& v2[1]);
    tmp.set_val_dfd1(& v1[1]);
    tmp.mult(& v2[0]);
    prod[2].sub(& tmp);
    return;
}

pub fn rotate_orient_dfd1(inst_ori : &mut [DiffDoub1], loc_ori : &mut [DiffDoub1], rot : &mut [DiffDoub1]) {
    let mut mag = DiffDoub1::new();
    let mut tmp = DiffDoub1::new();
    let mut tmp2 = DiffDoub1::new();
    let mut one_half = DiffDoub1::new();
    let mut a1 = [DiffDoub1::new(); 9];
    let mut a2 = [DiffDoub1::new(); 9];
    let mut a3 = [DiffDoub1::new(); 9];
    let mut tmp_v = [DiffDoub1::new(); 3];
    let mut i1 : usize;
    let mut i3 : usize;
    
    one_half.set_val(0.5);
    
    mag.set_val_dfd1(& rot[0]);
    mag.sqr();
    tmp.set_val_dfd1(& rot[1]);
    tmp.sqr();
    tmp2.set_val_dfd1(& rot[2]);
    tmp2.sqr();
    mag.add(& tmp);
    mag.add(& tmp2);
    mag.sqt();
    if mag.val < MAG_TOL {
        let mut loc_rot = [DiffDoub1::new(); 3];
        i3 = 0;
        for i1 in 0..3 {
            loc_rot[i1].set_val(0.0);
            for i2 in 0..3 {
                tmp.set_val_dfd1(& loc_ori[i3]);
                tmp.mult(& rot[i2]);
                loc_rot[i1].add(& tmp);
                i3 += 1usize;
            }
        }
        a1[0].set_val(1.0);
        tmp.set_val_dfd1(& loc_rot[1]);
        tmp.sqr();
        tmp2.set_val_dfd1(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[0].sub(& tmp);
        
        a1[1].set_val(0.5);
        a1[1].mult(& loc_rot[0]);
        a1[1].mult(& loc_rot[1]);
        a1[1].add(& loc_rot[2]);
        
        a1[2].set_val(0.5);
        a1[2].mult(& loc_rot[0]);
        a1[2].mult(& loc_rot[2]);
        a1[2].sub(& loc_rot[1]);
        
        a1[3].set_val(0.5);
        a1[3].mult(& loc_rot[0]);
        a1[3].mult(& loc_rot[1]);
        a1[3].sub(& loc_rot[2]);
        
        a1[4].set_val(1.0);
        tmp.set_val_dfd1(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd1(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[4].sub(& tmp);
        
        a1[5].set_val(0.5);
        a1[5].mult(& loc_rot[1]);
        a1[5].mult(& loc_rot[2]);
        a1[5].add(& loc_rot[0]);
        
        a1[6].set_val(0.5);
        a1[6].mult(& loc_rot[0]);
        a1[6].mult(& loc_rot[2]);
        a1[6].add(& loc_rot[1]);
        
        a1[7].set_val(0.5);
        a1[7].mult(& loc_rot[1]);
        a1[7].mult(& loc_rot[2]);
        a1[7].sub(& loc_rot[0]);
        
        a1[8].set_val(1.0);
        tmp.set_val_dfd1(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd1(& loc_rot[1]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[8].sub(& tmp);
        
        mat_mul_ar_dfd1(inst_ori, &mut  a1, loc_ori,  3,  3,  3);
    } else {
        let mut sth = DiffDoub1::new();
        let mut cth = DiffDoub1::new();
        
        sth.set_val_dfd1(& mag);
        sth.sn();
        cth.set_val_dfd1(& mag);
        cth.cs();
        
        let mut unit_rot = [DiffDoub1::new(); 3];
        tmp.set_val(1.0);
        tmp.dvd(& mag);
        unit_rot[0].set_val_dfd1(& rot[0]);
        unit_rot[0].mult(& tmp);
        unit_rot[1].set_val_dfd1(& rot[1]);
        unit_rot[1].mult(& tmp);
        unit_rot[2].set_val_dfd1(& rot[2]);
        unit_rot[2].mult(& tmp);
        
        a1[0].set_val_dfd1(& unit_rot[0]);
        a1[1].set_val_dfd1(& unit_rot[1]);
        a1[2].set_val_dfd1(& unit_rot[2]);
        
        i1 = 0;
        if fabs(unit_rot[1].val) < fabs(unit_rot[0].val) {
            i1 = 1;
        }
        if fabs(unit_rot[2].val) < fabs(unit_rot[i1].val) {
            i1 = 2;
        }
        tmp.set_val(1.0);
        tmp2.set_val_dfd1(& unit_rot[i1]);
        tmp2.sqr();
        tmp.sub(& tmp2);
        tmp.sqt();// tmp = sqrt(1 - unit_rot[i1]^2)
        a1[3+i1].set_val_dfd1(& tmp);
        for i2 in 0..3 {
            if i2 != i1 {
                i3 = 3 + i2;
                tmp2.set_val_dfd1(& a1[i1]);
                a1[i3].set_val_dfd1(& tmp2);
                a1[i3].neg();
                tmp2.set_val_dfd1(& a1[i2]);
                a1[i3].mult(& tmp2);
                a1[i3].dvd(& tmp);
            }
        }
        
        cross_prod_dfd1(&mut tmp_v, & a1[0..], & a1[3..]);
        for i in 0..3 {
            a1[i+6].set_val_dfd1(&tmp_v[i]);
        }
        
        a2[0].set_val_dfd1(& a1[0]);
        a2[1].set_val_dfd1(& a1[1]);
        a2[2].set_val_dfd1(& a1[2]);
        
        a2[3].set_val_dfd1(& cth);
        a2[3].mult(& a1[3]);
        tmp.set_val_dfd1(& sth);
        tmp.mult(& a1[6]);
        a2[3].sub(& tmp);
        
        a2[4].set_val_dfd1(& cth);
        a2[4].mult(& a1[4]);
        tmp.set_val_dfd1(& sth);
        tmp.mult(& a1[7]);
        a2[4].sub(& tmp);
        
        a2[5].set_val_dfd1(& cth);
        a2[5].mult(& a1[5]);
        tmp.set_val_dfd1(& sth);
        tmp.mult(& a1[8]);
        a2[5].sub(& tmp);
        
        a2[6].set_val_dfd1(& sth);
        a2[6].mult(& a1[3]);
        tmp.set_val_dfd1(& cth);
        tmp.mult(& a1[6]);
        a2[6].add(& tmp);
        
        a2[7].set_val_dfd1(& sth);
        a2[7].mult(& a1[4]);
        tmp.set_val_dfd1(& cth);
        tmp.mult(& a1[7]);
        a2[7].add(& tmp);
        
        a2[8].set_val_dfd1(& sth);
        a2[8].mult(& a1[5]);
        tmp.set_val_dfd1(& cth);
        tmp.mult(& a1[8]);
        a2[8].add(& tmp);
        
        transpose_ar_dfd1(&mut a3, &mut a2, 3, 3);// a3 = a2^t
        mat_mul_ar_dfd1(&mut a2, &mut a3, &mut a1, 3, 3, 3);//a2 = a2^t*a1
        mat_mul_ar_dfd1(inst_ori, loc_ori, &mut a2, 3, 3, 3);
    }
    return;
}

//end dup
 
//DiffDoub2 versions: 
//dup2

pub fn sub_vec_dfd2(sub_v : &mut Vec<DiffDoub2>, in_v : &mut Vec<DiffDoub2>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        sub_v[i2].set_val_dfd2(& in_v[i1]);
        i2 += 1usize;
    }
}

pub fn return_sv_dfd2(sub_v : &mut Vec<DiffDoub2>, in_v : &mut Vec<DiffDoub2>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        in_v[i1].set_val_dfd2(& sub_v[i2]);
        i2 += 1usize;
    }
}

pub fn vec_to_ar_dfd2(ar : &mut [DiffDoub2], vc : &mut Vec<DiffDoub2>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        ar[i2].set_val_dfd2(& vc[i1]);
        i2 += 1usize;
    }
}

pub fn ar_to_vec_dfd2(ar : &mut [DiffDoub2], vc : &mut Vec<DiffDoub2>, st : usize, end : usize) {
    let mut i2 : usize =  0;
    for i1 in st..end {
        vc[i1].set_val_dfd2(& ar[i2]);
        i2 += 1usize;
    }
}

pub fn mat_mul_dfd2(prod : &mut Vec<DiffDoub2>, mat1 : & Vec<DiffDoub2>, mat2 : & Vec<DiffDoub2>, m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub2::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd2(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6 +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn mat_mul_ar_dfd2(prod : &mut [DiffDoub2], mat1 : & [DiffDoub2], mat2 : & [DiffDoub2], m1_rows : usize, m1_cols : usize, m2_cols : usize) {
    let mut i4 : usize;
    let mut i5 : usize;
    let mut i6 : usize;
    let mut tmp = DiffDoub2::new();
    
    i4 = 0;
    i5 = 0;
    for _i1 in 0..m1_rows {
        for i2 in 0..m2_cols {
            i6 = i2;
            prod[i4].set_val(0.0);
            for _i3 in 0..m1_cols {
                tmp.set_val_dfd2(& mat1[i5]);
                tmp.mult(& mat2[i6]);
                prod[i4].add(& tmp);
                i5 += 1usize;
                i6  +=  m2_cols;
            }
            i5  -=  m1_cols;
            i4 += 1usize;
        }
        i5  +=  m1_cols;
    }
    return;
}

pub fn transpose_dfd2(mat_t : &mut Vec<DiffDoub2>, mat : &mut Vec<DiffDoub2>, row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd2(& mat[i3]);
            i3 += 1usize;
            i4 +=  row_dim;
        }
    }
    return;
}

pub fn transpose_ar_dfd2(mat_t : &mut [DiffDoub2], mat : &mut [DiffDoub2], row_dim : usize, col_dim : usize) {
    let mut i3 : usize;
    let mut i4 : usize;
    
    i3 = 0;
    for i1 in 0..row_dim {
        i4 = i1;
        for _i2 in 0..col_dim {
            mat_t[i4].set_val_dfd2(& mat[i3]);
            i3 += 1usize;
            i4  +=  row_dim;
        }
    }
    return;
}

pub fn cross_prod_dfd2(prod : &mut [DiffDoub2], v1 : & [DiffDoub2], v2 : & [DiffDoub2]) {
    let mut tmp = DiffDoub2::new();
    
    prod[0].set_val_dfd2(& v1[1]);
    prod[0].mult(& v2[2]);
    tmp.set_val_dfd2(& v1[2]);
    tmp.mult(& v2[1]);
    prod[0].sub(& tmp);
    
    prod[1].set_val_dfd2(& v1[2]);
    prod[1].mult(& v2[0]);
    tmp.set_val_dfd2(& v1[0]);
    tmp.mult(& v2[2]);
    prod[1].sub(& tmp);
    
    prod[2].set_val_dfd2(& v1[0]);
    prod[2].mult(& v2[1]);
    tmp.set_val_dfd2(& v1[1]);
    tmp.mult(& v2[0]);
    prod[2].sub(& tmp);
    return;
}

pub fn rotate_orient_dfd2(inst_ori : &mut [DiffDoub2], loc_ori : &mut [DiffDoub2], rot : &mut [DiffDoub2]) {
    let mut mag = DiffDoub2::new();
    let mut tmp = DiffDoub2::new();
    let mut tmp2 = DiffDoub2::new();
    let mut one_half = DiffDoub2::new();
    let mut a1 = [DiffDoub2::new(); 9];
    let mut a2 = [DiffDoub2::new(); 9];
    let mut a3 = [DiffDoub2::new(); 9];
    let mut tmp_v = [DiffDoub2::new(); 3];
    let mut i1 : usize;
    let mut i3 : usize;
    
    one_half.set_val(0.5);
    
    mag.set_val_dfd2(& rot[0]);
    mag.sqr();
    tmp.set_val_dfd2(& rot[1]);
    tmp.sqr();
    tmp2.set_val_dfd2(& rot[2]);
    tmp2.sqr();
    mag.add(& tmp);
    mag.add(& tmp2);
    mag.sqt();
    if mag.val < MAG_TOL {
        let mut loc_rot = [DiffDoub2::new(); 3];
        i3 = 0;
        for i1 in 0..3 {
            loc_rot[i1].set_val(0.0);
            for i2 in 0..3 {
                tmp.set_val_dfd2(& loc_ori[i3]);
                tmp.mult(& rot[i2]);
                loc_rot[i1].add(& tmp);
                i3 += 1usize;
            }
        }
        a1[0].set_val(1.0);
        tmp.set_val_dfd2(& loc_rot[1]);
        tmp.sqr();
        tmp2.set_val_dfd2(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[0].sub(& tmp);
        
        a1[1].set_val(0.5);
        a1[1].mult(& loc_rot[0]);
        a1[1].mult(& loc_rot[1]);
        a1[1].add(& loc_rot[2]);
        
        a1[2].set_val(0.5);
        a1[2].mult(& loc_rot[0]);
        a1[2].mult(& loc_rot[2]);
        a1[2].sub(& loc_rot[1]);
        
        a1[3].set_val(0.5);
        a1[3].mult(& loc_rot[0]);
        a1[3].mult(& loc_rot[1]);
        a1[3].sub(& loc_rot[2]);
        
        a1[4].set_val(1.0);
        tmp.set_val_dfd2(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd2(& loc_rot[2]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[4].sub(& tmp);
        
        a1[5].set_val(0.5);
        a1[5].mult(& loc_rot[1]);
        a1[5].mult(& loc_rot[2]);
        a1[5].add(& loc_rot[0]);
        
        a1[6].set_val(0.5);
        a1[6].mult(& loc_rot[0]);
        a1[6].mult(& loc_rot[2]);
        a1[6].add(& loc_rot[1]);
        
        a1[7].set_val(0.5);
        a1[7].mult(& loc_rot[1]);
        a1[7].mult(& loc_rot[2]);
        a1[7].sub(& loc_rot[0]);
        
        a1[8].set_val(1.0);
        tmp.set_val_dfd2(& loc_rot[0]);
        tmp.sqr();
        tmp2.set_val_dfd2(& loc_rot[1]);
        tmp2.sqr();
        tmp.add(& tmp2);
        tmp.mult(& one_half);
        a1[8].sub(& tmp);
        
        mat_mul_ar_dfd2(inst_ori, &mut  a1, loc_ori,  3,  3,  3);
    } else {
        let mut sth = DiffDoub2::new();
        let mut cth = DiffDoub2::new();
        
        sth.set_val_dfd2(& mag);
        sth.sn();
        cth.set_val_dfd2(& mag);
        cth.cs();
        
        let mut unit_rot = [DiffDoub2::new(); 3];
        tmp.set_val(1.0);
        tmp.dvd(& mag);
        unit_rot[0].set_val_dfd2(& rot[0]);
        unit_rot[0].mult(& tmp);
        unit_rot[1].set_val_dfd2(& rot[1]);
        unit_rot[1].mult(& tmp);
        unit_rot[2].set_val_dfd2(& rot[2]);
        unit_rot[2].mult(& tmp);
        
        a1[0].set_val_dfd2(& unit_rot[0]);
        a1[1].set_val_dfd2(& unit_rot[1]);
        a1[2].set_val_dfd2(& unit_rot[2]);
        
        i1 = 0;
        if fabs(unit_rot[1].val) < fabs(unit_rot[0].val) {
            i1 = 1;
        }
        if fabs(unit_rot[2].val) < fabs(unit_rot[i1].val) {
            i1 = 2;
        }
        tmp.set_val(1.0);
        tmp2.set_val_dfd2(& unit_rot[i1]);
        tmp2.sqr();
        tmp.sub(& tmp2);
        tmp.sqt();// tmp = sqrt(1 - unit_rot[i1]^2)
        a1[3+i1].set_val_dfd2(& tmp);
        for i2 in 0..3 {
            if i2 != i1 {
                i3 = 3 + i2;
                tmp2.set_val_dfd2(& a1[i1]);
                a1[i3].set_val_dfd2(& tmp2);
                a1[i3].neg();
                tmp2.set_val_dfd2(& a1[i2]);
                a1[i3].mult(& tmp2);
                a1[i3].dvd(& tmp);
            }
        }
        
        cross_prod_dfd2(&mut tmp_v, & a1[0..], & a1[3..]);
        for i in 0..3 {
            a1[i+6].set_val_dfd2(&tmp_v[i]);
        }
        
        a2[0].set_val_dfd2(& a1[0]);
        a2[1].set_val_dfd2(& a1[1]);
        a2[2].set_val_dfd2(& a1[2]);
        
        a2[3].set_val_dfd2(& cth);
        a2[3].mult(& a1[3]);
        tmp.set_val_dfd2(& sth);
        tmp.mult(& a1[6]);
        a2[3].sub(& tmp);
        
        a2[4].set_val_dfd2(& cth);
        a2[4].mult(& a1[4]);
        tmp.set_val_dfd2(& sth);
        tmp.mult(& a1[7]);
        a2[4].sub(& tmp);
        
        a2[5].set_val_dfd2(& cth);
        a2[5].mult(& a1[5]);
        tmp.set_val_dfd2(& sth);
        tmp.mult(& a1[8]);
        a2[5].sub(& tmp);
        
        a2[6].set_val_dfd2(& sth);
        a2[6].mult(& a1[3]);
        tmp.set_val_dfd2(& cth);
        tmp.mult(& a1[6]);
        a2[6].add(& tmp);
        
        a2[7].set_val_dfd2(& sth);
        a2[7].mult(& a1[4]);
        tmp.set_val_dfd2(& cth);
        tmp.mult(& a1[7]);
        a2[7].add(& tmp);
        
        a2[8].set_val_dfd2(& sth);
        a2[8].mult(& a1[5]);
        tmp.set_val_dfd2(& cth);
        tmp.mult(& a1[8]);
        a2[8].add(& tmp);
        
        transpose_ar_dfd2(&mut a3, &mut a2, 3, 3);// a3 = a2^t
        mat_mul_ar_dfd2(&mut a2, &mut a3, &mut a1, 3, 3, 3);//a2 = a2^t*a1
        mat_mul_ar_dfd2(inst_ori, loc_ori, &mut a2, 3, 3, 3);
    }
    return;
}

//end dup
 
//end skip 
 
 
//dup1

pub fn d_orid_thet_dfd0(inst_ori : &mut [DiffDoub0], loc_ori : &mut [DiffDoub0], rot : &mut [DiffDoub0], v1 : usize, v2 : usize) {
    if v1 + v2 == 0 {
        rotate_orient_dfd0(inst_ori, loc_ori, rot);
        return;
    }
    //preserve
    let mut d_rot = [DiffDoub1::new(); 3];
    let mut d2_rot = [DiffDoub2::new(); 3];
    let mut d_loc_ori = [DiffDoub1::new(); 9];
    let mut d2_loc_ori = [DiffDoub2::new(); 9];
    let mut d_inst_ori = [DiffDoub1::new(); 9];
    let mut d2_inst_ori = [DiffDoub2::new(); 9];
    //end preserve
    let is_diff : bool =  loc_ori[0].diff_type();
    let mut param : f64;
    let mut param2 : f64;
    if !is_diff {
        if v1*v2 == 0 {
            for i1 in 0..9 {
                param = loc_ori[i1].val;
                d_loc_ori[i1].set_val(param);
            }
            let vnz = v1 + v2;
            //v1 = v1 + v2;
            for i1 in 1..4 {
                param = rot[i1-1].val;
                if i1 == vnz {
                    d_rot[i1-1].set_val_2(param,   1.0);
                } else {
                    d_rot[i1-1].set_val_2(param,   0.0);
                }
            }
            rotate_orient_dfd1(&mut d_inst_ori, &mut d_loc_ori, &mut d_rot);
            for i1 in 0..9 {
                param = d_inst_ori[i1].dval;
                inst_ori[i1].set_val(param);
            }
        } else {
            for i1 in 0..9 {
                param = loc_ori[i1].val;
                param2 = loc_ori[i1].dval;
                d2_loc_ori[i1].set_val_2(param,   param2);
            }
            if v1 == v2 {
                for i1 in 1..4 {
                    param = rot[i1-1].val;
                    if i1 == v1 {
                        d2_rot[i1-1].set_val_3(param, 1.0, 1.0);
                    } else {
                        d2_rot[i1-1].set_val(param);
                    }
                }
            } else {
                for i1 in 1..4 {
                    param = rot[i1-1].val;
                    if i1 == v1 {
                        d2_rot[i1-1].set_val_2(param,   1.0);
                    } else if i1 == v2 {
                        d2_rot[i1-1].set_val_3(param, 0.0, 1.0);
                    } else {
                        d2_rot[i1-1].set_val_2(param,   0.0);
                    }
                }
            }
            rotate_orient_dfd2(&mut d2_inst_ori, &mut d2_loc_ori, &mut d2_rot);
            for i1 in 0..9 {
                param = d2_inst_ori[i1].dv12;
                inst_ori[i1].set_val(param);
            }
        }
    } else {
        for i1 in 0..9 {
            param = loc_ori[i1].val;
            param2 = loc_ori[i1].dval;
            d2_loc_ori[i1].set_val_2(param,   param2);
        }
        let vnz = v1 + v2;
        //v1 = v1 + v2;
        for i1 in 1..4 {
            param = rot[i1-1].val;
            if i1 == vnz {
                d2_rot[i1-1].set_val_3(param, 0.0, 1.0);
            } else {
                d2_rot[i1-1].set_val_2(param,   0.0);
            }
        }
        rotate_orient_dfd2(&mut d2_inst_ori, &mut d2_loc_ori, &mut d2_rot);
        for i1 in 0..9 {
            param = d2_inst_ori[i1].dv2;
            param2 = d2_inst_ori[i1].dv12;
            inst_ori[i1].set_val_2(param, param2);
        }
    }
    return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

pub fn d_orid_thet_dfd1(inst_ori : &mut [DiffDoub1], loc_ori : &mut [DiffDoub1], rot : &mut [DiffDoub1], v1 : usize, v2 : usize) {
    if v1 + v2 == 0 {
        rotate_orient_dfd1(inst_ori, loc_ori, rot);
        return;
    }
    let mut d_rot = [DiffDoub1::new(); 3];
    let mut d2_rot = [DiffDoub2::new(); 3];
    let mut d_loc_ori = [DiffDoub1::new(); 9];
    let mut d2_loc_ori = [DiffDoub2::new(); 9];
    let mut d_inst_ori = [DiffDoub1::new(); 9];
    let mut d2_inst_ori = [DiffDoub2::new(); 9];
    let is_diff : bool =  loc_ori[0].diff_type();
    let mut param : f64;
    let mut param2 : f64;
    if !is_diff {
        if v1*v2 == 0 {
            for i1 in 0..9 {
                param = loc_ori[i1].val;
                d_loc_ori[i1].set_val(param);
            }
            let vnz = v1 + v2;
            //v1 = v1 + v2;
            for i1 in 1..4 {
                param = rot[i1-1].val;
                if i1 == vnz {
                    d_rot[i1-1].set_val_2(param,   1.0);
                } else {
                    d_rot[i1-1].set_val_2(param,   0.0);
                }
            }
            rotate_orient_dfd1(&mut d_inst_ori, &mut d_loc_ori, &mut d_rot);
            for i1 in 0..9 {
                param = d_inst_ori[i1].dval;
                inst_ori[i1].set_val(param);
            }
        } else {
            for i1 in 0..9 {
                param = loc_ori[i1].val;
                param2 = loc_ori[i1].dval;
                d2_loc_ori[i1].set_val_2(param,   param2);
            }
            if v1 == v2 {
                for i1 in 1..4 {
                    param = rot[i1-1].val;
                    if i1 == v1 {
                        d2_rot[i1-1].set_val_3(param, 1.0, 1.0);
                    } else {
                        d2_rot[i1-1].set_val(param);
                    }
                }
            } else {
                for i1 in 1..4 {
                    param = rot[i1-1].val;
                    if i1 == v1 {
                        d2_rot[i1-1].set_val_2(param,   1.0);
                    } else if i1 == v2 {
                        d2_rot[i1-1].set_val_3(param, 0.0, 1.0);
                    } else {
                        d2_rot[i1-1].set_val_2(param,   0.0);
                    }
                }
            }
            rotate_orient_dfd2(&mut d2_inst_ori, &mut d2_loc_ori, &mut d2_rot);
            for i1 in 0..9 {
                param = d2_inst_ori[i1].dv12;
                inst_ori[i1].set_val(param);
            }
        }
    } else {
        for i1 in 0..9 {
            param = loc_ori[i1].val;
            param2 = loc_ori[i1].dval;
            d2_loc_ori[i1].set_val_2(param,   param2);
        }
        let vnz = v1 + v2;
        //v1 = v1 + v2;
        for i1 in 1..4 {
            param = rot[i1-1].val;
            if i1 == vnz {
                d2_rot[i1-1].set_val_3(param, 0.0, 1.0);
            } else {
                d2_rot[i1-1].set_val_2(param,   0.0);
            }
        }
        rotate_orient_dfd2(&mut d2_inst_ori, &mut d2_loc_ori, &mut d2_rot);
        for i1 in 0..9 {
            param = d2_inst_ori[i1].dv2;
            param2 = d2_inst_ori[i1].dv12;
            inst_ori[i1].set_val_2(param, param2);
        }
    }
    return;
}

//end dup
 
//end skip 
 
 
