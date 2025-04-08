use crate::element::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::list_ent::*;
use crate::design_var::*;
use crate::node::*;
use crate::face::*;
use crate::section::*;
use crate::load::*;
use crate::job::*;
use crate::matrix_functions::*;
use crate::scratch::*;
use crate::cpp_str::CppStr;

use std::collections::linked_list::IterMut;


impl Element {
    pub fn condense_mat(&mut self, mat : &mut Vec<f64>, scr1 : &mut Vec<f64>, scr2 : &mut Vec<f64>) {
        let st_row : usize;
        let end_row : usize;
        let end_col : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        
        st_row = self.num_nds*self.dof_per_nd;
        end_row = st_row + self.num_int_dof - 1;
        end_col = self.num_int_dof - 1;
        q_rfactor(&mut self.internal_mat,  self.num_int_dof,  st_row,  end_row,  0,  end_col,  0);
        
        for i1 in 0..st_row {
            scr1[i1] = 0.0;
        }
        
        for i1 in 0..st_row {
            i3 = st_row;
            i4 = i1*self.num_int_dof;
            for _i2 in 0..self.num_int_dof {
                scr1[i3] = self.internal_mat[i4];
                i3 += 1usize;
                i4 += 1usize;
            }
            solveq_rx_eqb(scr2, &mut self.internal_mat, scr1, self.num_int_dof, st_row, end_row, 0, end_col, 0);
            i4 = i1;
            for i2 in 0..st_row {
                i5 = i2*self.num_int_dof;
                for i3 in 0..self.num_int_dof {
                    mat[i4] -=  self.internal_mat[i5]*scr2[i3];
                    i5 += 1usize;
                }
                i4 +=  end_row+1;
            }
        }
        
        return;
    }

    pub fn update_external(&mut self, ext_vec : &mut Vec<f64>, for_soln : usize, nd_ar : &mut Vec<Node>, sc_it : &mut IterMut<'_,FltScr>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let st_row : usize;
        let end_row : usize;
        let end_col : usize;
        let mut nd : usize;
        let mut dof : usize;
        let mut glb_ind : usize;
        
        if self.num_int_dof == 0 {
            return;
        }

        let scr1 = match sc_it.next() {
            None => panic!("Error: ran out of scratch matrices."),
            Some(x) => &mut x.dat,
        };

        let scr2 = match sc_it.next() {
            None => panic!("Error: ran out of scratch matrices."),
            Some(x) => &mut x.dat,
        };
        
        st_row = self.num_nds*self.dof_per_nd;
        end_row = st_row + self.num_int_dof - 1;
        end_col = self.num_int_dof - 1;
        for i1 in 0..st_row {
            scr1[i1] = 0.0;
        }
        
        if for_soln == 1 {
            i2 = st_row;
            for i1 in 0..self.num_int_dof {
                scr1[i2] = self.internal_ru[i1].val;
                i2 += 1usize;
            }
        } else {
            i2 = st_row;
            for i1 in 0..self.num_int_dof {
                scr1[i2] = -self.internald_ldu[i1];
                i2 += 1usize;
            }
        }
        
        solveq_rx_eqb(scr2, &mut self.internal_mat, scr1, self.num_int_dof, st_row, end_row, 0, end_col, 0);
        i3 = 0;
        i4 = 0;
        for _i1 in 0..st_row {
            nd = self.nodes[self.dof_table[i4]];
            dof = self.dof_table[i4+1];
            glb_ind = nd_ar[nd].dof_index[dof];
            for i2 in 0..self.num_int_dof {
                ext_vec[glb_ind] +=  self.internal_mat[i3]*scr2[i2];
                i3 += 1usize;
            }
            i4  +=  2;
        }
        
        return;
    }

    pub fn update_internal(&mut self, ext_vec : &mut Vec<f64>, for_soln : usize, nd_ar : &mut Vec<Node>, sc_it : &mut IterMut<'_,FltScr>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let st_row : usize;
        let end_row : usize;
        let end_col : usize;
        let mut nd : usize;
        let mut dof : usize;
        let mut glb_ind : usize;
        
        if self.num_int_dof == 0 {
            return;
        }

        let scr1 = match sc_it.next() {
            None => panic!("Error: ran out of scratch matrices."),
            Some(x) => &mut x.dat,
        };

        let scr2 = match sc_it.next() {
            None => panic!("Error: ran out of scratch matrices."),
            Some(x) => &mut x.dat,
        };
        
        st_row = self.num_nds*self.dof_per_nd;
        end_row = st_row + self.num_int_dof - 1;
        end_col = self.num_int_dof - 1;
        for i1 in 0..st_row {
            scr1[i1] = 0.0;
        }
        
        if for_soln == 1 {
            i2 = st_row;
            for i1 in 0..self.num_int_dof {
                scr1[i2] = -self.internal_ru[i1].val;
                i2 += 1usize;
            }
        } else {
            i2 = st_row;
            for i1 in 0..self.num_int_dof {
                scr1[i2] = self.internald_ldu[i1];
                i2 += 1usize;
            }
        }
        
        for i1 in 0..st_row {
            nd = self.nodes[self.dof_table[2*i1]];
            dof = self.dof_table[2*i1+1];
            glb_ind = nd_ar[nd].dof_index[dof];
            i3 = i1*self.num_int_dof;
            i4 = st_row;
            for _i2 in 0..self.num_int_dof {
                scr1[i4]  -=  self.internal_mat[i3]*ext_vec[glb_ind];
                i3 += 1usize;
                i4 += 1usize;
            }
        }
        
        solveq_rx_eqb(scr2, &mut self.internal_mat, scr1, self.num_int_dof, st_row, end_row, 0, end_col, 0);
        
        if for_soln == 1 {
            for i1 in 0..self.num_int_dof {
                self.internal_disp[i1]  +=  scr2[i1];
            }
        } else {
            for i1 in 0..self.num_int_dof {
                self.internal_adj[i1] = scr2[i1];
            }
        }
        
        return;
    }

    pub fn get_int_adjd_rd_d(&mut self) -> f64 {
        let mut prod : f64 =  0.0;
        for i1 in 0..self.num_int_dof {
            prod  +=  self.internal_adj[i1] * self.internal_ru[i1].dval;
        }
        return  prod;
    }

    //dup1

    pub fn get_ruk_dfd0(&mut self, rvec : &mut Vec<DiffDoub0>, d_rdu : &mut Vec<f64>, d_rd_t : &mut Vec<f64>, get_matrix : bool, n_lgeom : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let tot_dof : usize;
        
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut d_jwt = DiffDoub0::new();
        
        let mut strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut thrm_stn = [DiffDoub0::new(); 6];
        let mut ip_temp = DiffDoub0::new();
        let mut cte = [DiffDoub0::new(); 6];
        let mut cten = [DiffDoub0::new(); 90];
        let mut ux = [DiffDoub0::new(); 9];
        let mut sec_def = [DiffDoub0::new(); 9];
        let mut sec_fc_mom = [DiffDoub0::new(); 9];
        
        let mut tmp = DiffDoub0::new();
        let mut tmp_c = [DiffDoub0::new(); 81];
        let mut tmp_gd = [DiffDoub0::new(); 60];
        let mut tmp_te = [DiffDoub0::new(); 6];
        
        vec_to_ar_dfd0(&mut tmp_c, &mut  pre.cmat,  0,  81);
        vec_to_ar_dfd0(&mut tmp_gd, &mut  pre.glob_disp,  0,  60);
        vec_to_ar_dfd0(&mut tmp_te, &mut  pre.therm_exp,  0,  6);
        
        tot_dof = self.num_nds*self.dof_per_nd + self.num_int_dof;
        
        i3 = 0;
        i4 = 0;
        for i1 in 0..tot_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..tot_dof {
                    d_rdu[i3] = 0.0;
                    i3 += 1usize;
                }
                for _i2 in 0..self.num_nds {
                    d_rd_t[i4] = 0.0;
                    i4 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        if self.dof_per_nd == 6 {
            if n_lgeom {
                self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
        }
        
        for i1 in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.loc_nds, & self.int_pts[(3*i1)..]);
            d_jwt.set_val_dfd0(& det_j);
            tmp.set_val(self.ip_wt[i1]);
            d_jwt.mult(& tmp);
            ip_temp.set_val(0.0);
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd0(& pre.glob_temp[i2]);
                tmp.mult(& n_vec[i2]);
                ip_temp.add(& tmp);
            }
            if self.dof_per_nd == 3 {
                mat_mul_ar_dfd0(&mut ux, &mut tmp_gd, &mut d_ndx, 3, self.n_dim, 3);
                self.get_solid_strain_dfd0(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, MAX_INT, MAX_INT, n_lgeom);
                for i2 in 0..6 {
                    thrm_stn[i2].set_val_dfd0(& pre.therm_exp[i2]);
                    thrm_stn[i2].mult(& ip_temp);
                    strain[i2].sub(& thrm_stn[i2]);
                    strain[i2].sub(& pre.einit[i2]);
                }
                mat_mul_ar_dfd0(&mut stress, &mut tmp_c, &mut strain, 6, 6, 1);
                if get_matrix {
                    mat_mul_ar_dfd0(&mut cte, &mut  tmp_c, &mut  tmp_te,  6,  6,  1);
                    mat_mul_ar_dfd0(&mut cten, &mut  cte, &mut  n_vec,  6,  1,  self.num_nds);
                }
                for i2 in 0..tot_dof {
                    self.get_solid_strain_dfd0(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, i2, MAX_INT, n_lgeom);
                    i4 = i2;
                    for i3 in 0..6 {
                        tmp.set_val_dfd0(& stress[i3]);
                        tmp.mult(& strain[i3]);
                        tmp.mult(& d_jwt);
                        rvec[i2].add(& tmp);
                        if get_matrix {
                            pre.bmat[i4].set_val_dfd0(& strain[i3]);
                            i4 +=  tot_dof;
                        }
                    }
                    if n_lgeom && get_matrix {
                        i5 = (tot_dof + 1)*i2;
                        i6 = i5;
                        for i3 in i2..tot_dof {
                            self.get_solid_strain_dfd0(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, i2, i3, n_lgeom);
                            for i4 in 0..6 {
                                d_rdu[i5] +=  stress[i4].val*strain[i4].val*d_jwt.val;
                            }
                            d_rdu[i6] = d_rdu[i5];
                            i5 += 1usize;
                            i6 +=  tot_dof;
                        }
                    }
                }
            } else {
                self.get_section_def_dfd0(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, MAX_INT, MAX_INT);
                mat_mul_ar_dfd0(&mut sec_fc_mom, &mut tmp_c, &mut sec_def, self.def_dim, self.def_dim, 1);
                for i2 in 0..6 {
                    tmp.set_val_dfd0(& pre.therm_exp[i2]);
                    tmp.mult(& ip_temp);
                    tmp.add(& pre.einit[i2]);
                    sec_fc_mom[i2].sub(& tmp);
                }
                if get_matrix {
                    mat_mul_ar_dfd0(&mut cten, &mut  tmp_te, &mut  n_vec,  6,  1,  self.num_nds);
                }
                for i2 in 0..tot_dof {
                    self.get_section_def_dfd0(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, i2, MAX_INT);
                    i4 = i2;
                    for i3 in 0..self.def_dim {
                        tmp.set_val_dfd0(& sec_fc_mom[i3]);
                        tmp.mult(& sec_def[i3]);
                        tmp.mult(& d_jwt);
                        rvec[i2].add(& tmp);
                        if get_matrix {
                            pre.bmat[i4].set_val_dfd0(& sec_def[i3]);
                            i4 +=  tot_dof;
                        }
                    }
                    if n_lgeom && get_matrix {
                        i5 = (tot_dof + 1)*i2;
                        i6 = i5;
                        for i3 in i2..tot_dof {
                            self.get_section_def_dfd0(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, i2, i3);
                            for i4 in 0..self.def_dim {
                                d_rdu[i5] +=  sec_fc_mom[i4].val*sec_def[i4].val*d_jwt.val;
                            }
                            d_rdu[i6] = d_rdu[i5];
                            i5 += 1usize;
                            i6 +=  tot_dof;
                        }
                    }
                }
            }
            if get_matrix {
                mat_mul_dfd0(&mut pre.cbmat, &mut pre.cmat, &mut pre.bmat, self.def_dim, self.def_dim, tot_dof);
                i3 = self.def_dim*tot_dof;
                for i2 in 0..i3 {
                    pre.cbmat[i2].mult(& d_jwt);
                }
                i5 = 0;
                for i2 in 0..tot_dof {
                    for i3 in 0..tot_dof {
                        i6 = i2;
                        i7 = i3;
                        for _i4 in 0..self.def_dim {
                            d_rdu[i5] +=  pre.bmat[i6].val*pre.cbmat[i7].val;
                            i6 +=  tot_dof;
                            i7 +=  tot_dof;
                        }
                        i5 += 1usize;
                    }
                }
                for i2 in 0..6 * self.num_nds {
                    cten[i2].mult(& d_jwt);
                }
                i5 = 0;
                for i2 in 0..tot_dof {
                    for i3 in 0..self.num_nds {
                        i6 = i2;
                        i7 = i3;
                        for _i4 in 0..6 {
                            d_rd_t[i5]  -=  pre.bmat[i6].val * cten[i7].val;
                            i6  +=  tot_dof;
                            i7  +=  self.num_nds;
                        }
                        i5 += 1usize;
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_rum_dfd0(& self, rvec : &mut Vec<DiffDoub0>, d_rd_a : &mut Vec<f64>, get_matrix : bool, actual_props : bool, n_lgeom : bool, 
        pre : &mut DiffDoub0StressPrereq, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd1 : usize;
        let mut dof1 : usize;
        let mut nd2 : usize;
        let mut dof2 : usize;
        let nd_dof : usize =  self.num_nds * self.dof_per_nd;
        
        let mut inst_disp = [DiffDoub0::new(); 60];
        
        let mut tmp_s : [f64; 3] = [0f64; 3];
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        
        let mut tmp = DiffDoub0::new();
        let mut tmp61 = [DiffDoub0::new(); 6];
        let mut det_jwt = DiffDoub0::new();
        
        let mut save_m = [DiffDoub0::new(); 36];
        
        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v2 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        
        i3 = 0;
        for i1 in 0..nd_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..nd_dof {
                    d_rd_a[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 {
            i2 = 0;
            for i1 in 0..3 {
                rvec[i1].set_val_dfd0(& pre.glob_acc[i1]);
                rvec[i1].mult(& pre.mass_per_el);
                if get_matrix {
                    d_rd_a[i2] = pre.mass_per_el.val;
                    i2  +=  4;
                }
            }
            return;
        }
        else if self.this_type == 21 {
            return;
        }
        
        if self.dof_per_nd == 6 {
            if n_lgeom {
                self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  1);
            }
            if !actual_props {
                for i1 in 0..36 {
                    save_m[i1].set_val_dfd0(& pre.mmat[i1]);
                    pre.mmat[i1].set_val(0.0);
                }
                i1 = 0usize;
                while i1 < 36 {
                    pre.mmat[i1].set_val(1.0);
                    i1 += 7;
                }
            }
        }
        else {
            if !actual_props {
                save_m[0].set_val_dfd0(& pre.mmat[0]);
                pre.mmat[0].set_val(1.0);
            }
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            det_jwt.set_val(self.ip_wt[i1]);
            det_jwt.mult(& det_j);
            if self.dof_per_nd == 6 {
                // build matrix [d_u^i/d_u^g] * {n}
                i7 = 0;
                for i2 in 0..nd_dof {
                    //nd1 = self.dof_table[i7];
                    //dof1 = self.dof_table[i7 + 1];
                    self.get_inst_disp_dfd0(&mut inst_disp, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds,  n_lgeom,  i2,  MAX_INT);
                    i6 = 0;
                    i5 = i2;
                    for _i3 in 0..6 {
                        pre.bmat[i5].set_val(0.0);
                        for i4 in 0..self.num_nds {
                            tmp.set_val_dfd0(& n_vec[i4]);
                            tmp.mult(& inst_disp[i6]);
                            pre.bmat[i5].add(& tmp);
                            i6 += 1usize;
                        }
                        i5  +=  nd_dof;
                        i6  +=  self.n_dim - self.num_nds;
                    }
                    i7  +=  2;
                }
                // tmp61 = bmat*{acc}
                i5 = 0;
                for i2 in 0..6 {
                    i4 = 0;
                    tmp61[i2].set_val(0.0);
                    for _i3 in 0..nd_dof {
                        nd1 = self.dof_table[i4];
                        dof1 = self.dof_table[i4 + 1];
                        i6 = dof1 * self.num_nds + nd1;
                        tmp.set_val_dfd0(& pre.bmat[i5]);
                        tmp.mult(& pre.glob_acc[i6]);
                        tmp61[i2].add(& tmp);
                        i4  +=  2;
                        i5 += 1usize;
                    }
                }
                // pre.scr_vec5 = [m][b]{a}, originally tmp62
                ar_to_vec_dfd0(&mut tmp61, &mut scr_v1, 0, 6);
                mat_mul_dfd0(&mut scr_v2, &mut pre.mmat, &mut  scr_v1,  6,  6,  1);
                mat_mul_dfd0(&mut scr_v1, &mut  scr_v2, &mut  pre.bmat,  1,  6,  nd_dof);
                // update rvec
                for i2 in 0..nd_dof {
                    tmp.set_val_dfd0(& scr_v1[i2]);
                    tmp.mult(& det_jwt);
                    rvec[i2].add(& tmp);
                }
                if get_matrix {
                    mat_mul_dfd0(&mut pre.cbmat, &mut  pre.mmat, &mut  pre.bmat,  6,  6,  nd_dof);
                    i4 = 6 * nd_dof;
                    for i2 in 0..i4 {
                        pre.cbmat[i2].mult(& det_jwt);
                    }
                    i5 = 0;
                    for i2 in 0..nd_dof {
                        for i3 in 0..nd_dof {
                            i6 = i2;
                            i7 = i3;
                            for _i4 in 0..6 {
                                d_rd_a[i5]  +=  pre.cbmat[i6].val * pre.bmat[i7].val;
                                i6  +=  nd_dof;
                                i7  +=  nd_dof;
                            }
                            i5 += 1usize;
                        }
                    }
                }
            }
            else {
                ar_to_vec_dfd0(&mut n_vec, &mut  scr_v1,  0,  11);
                mat_mul_dfd0(&mut pre.bmat, &mut  pre.glob_acc, &mut  scr_v1,  3,  self.num_nds,  1);
                for i2 in 0..nd_dof {
                    nd1 = self.dof_table[2 * i2];
                    dof1 = self.dof_table[2 * i2 + 1];
                    tmp.set_val_dfd0(& n_vec[nd1]);
                    tmp.mult(& pre.mmat[0]);
                    tmp.mult(& pre.bmat[dof1]);
                    tmp.mult(& det_jwt);
                    rvec[i2].add(& tmp);
                    if get_matrix {
                        for i3 in 0..nd_dof {
                            nd2 = self.dof_table[2 * i3];
                            dof2 = self.dof_table[2 * i3 + 1];
                            if dof2 == dof1 {
                                tmp.set_val_dfd0(& n_vec[nd1]);
                                tmp.mult(& n_vec[nd2]);
                                tmp.mult(& pre.mmat[0]);
                                tmp.mult(& det_jwt);
                                i4 = i2 * nd_dof + i3;
                                d_rd_a[i4]  +=  tmp.val;
                            }
                        }
                    }
                }
            }
        }
        
        if self.dof_per_nd == 6 {
            if !actual_props {
                for i1 in 0..36 {
                    pre.mmat[i1].set_val_dfd0(& save_m[i1]);
                }
            }
        }
        else {
            if !actual_props {
                pre.mmat[0].set_val_dfd0(& save_m[0]);
            }
        }
        
        return;
    }

    pub fn get_rud_dfd0(&mut self, rvec : &mut Vec<DiffDoub0>, d_rd_v : &mut Vec<f64>, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub0StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let mut dof : usize;
        let nd_dof : usize =  self.num_nds * self.dof_per_nd;
        let mut tmp = DiffDoub0::new();
        //let mut rtmp = &mut scr.scr_v3;
        let rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v2 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rtmp[1089];
        //let mut d_rtmp = &mut scr.scr_m4;
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rd_t[330];
        //let mut d_rd_t = &mut scr.scr_m5;
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_non_zero : bool;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut wt_det_j = DiffDoub0::new();
        let mut strain = [DiffDoub0::new(); 6];
        let mut sec_def = [DiffDoub0::new(); 9];
        let mut ux = [DiffDoub0::new(); 9];
        let mut bvel = [DiffDoub0::new(); 9];
        
        let mut tmp_s : [f64; 3] = [0f64; 3];
        
        i3 = 0;
        for i1 in 0..nd_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..nd_dof {
                    d_rd_v[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        if cmd.ray_damp_cm > 0.0 {
            tmp.set_val(cmd.ray_damp_cm);
            for i1 in 0..36 {
                pre.mmat[i1].mult(& tmp);
            }
            for i1 in 0..nd_dof {
                pre.glob_acc[i1].set_val_dfd0(& pre.glob_vel[i1]);
            }
            self.get_rum_dfd0(rtmp, d_rtmp,  get_matrix,  true,  cmd.nonlinear_geom, pre, scr_dfd);
            for i1 in 0..nd_dof {
                rvec[i1].add(& rtmp[i1]);
            }
            if get_matrix {
                i3 = 0;
                for _i1 in 0..nd_dof {
                    for _i2 in 0..nd_dof {
                        d_rd_v[i3]  +=  d_rtmp[i3];
                        i3 += 1usize;
                    }
                }
            }
        }
        if cmd.ray_damp_ck > 0.0 {
            tmp.set_val(cmd.ray_damp_ck);
            for i1 in 0..81 {
                pre.cmat[i1].mult(& tmp);
            }
            i3 = 0;
            i4 = 0;
            for _i1 in 0..self.dof_per_nd {
                for i2 in 0..self.n_dim {
                    if i2 >= self.num_nds {
                        pre.glob_disp[i3].set_val(0.0);
                    }
                    else {
                        pre.glob_disp[i3].set_val_dfd0(& pre.glob_vel[i4]);
                        i4 += 1usize;
                    }
                    i3 += 1usize;
                }
            }
            self.get_ruk_dfd0(rtmp, &mut  d_rtmp, &mut  d_rd_t,  get_matrix,  cmd.nonlinear_geom,  pre);
            for i1 in 0..nd_dof {
                rvec[i1].add(& rtmp[i1]);
            }
            if get_matrix {
                i3 = 0;
                i4 = 0;
                for _i1 in 0..nd_dof {
                    for _i2 in 0..nd_dof {
                        d_rd_v[i3]  +=  d_rtmp[i4];
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                    i3  +=  self.num_int_dof;
                }
            }
        }
        
        // Material damping
        d_non_zero = false;
        i1 = 0;
        while !d_non_zero && i1 < 36 {
            if pre.dmat[i1].val > 0.0 {
                d_non_zero = true;
            }
            i1 += 1usize;
        }
        
        if d_non_zero {
            if cmd.nonlinear_geom {
                self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            for i1 in 0..self.num_ip {
                i2 = i1 * 3;
                vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
                self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
                wt_det_j.set_val(self.ip_wt[i1]);
                wt_det_j.mult(& det_j);
                if self.dof_per_nd == 3 {
                    ar_to_vec_dfd0(&mut d_ndx, &mut  scr_v1,  0,  33);
                    mat_mul_dfd0(&mut scr_v2, &mut  pre.glob_disp, &mut  scr_v1,  3,  self.n_dim,  3);
                    vec_to_ar_dfd0(&mut ux, &mut  scr_v2,  0,  9);
                    for i2 in 0..nd_dof {
                        self.get_solid_strain_dfd0(&mut strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  i2,  MAX_INT,  cmd.nonlinear_geom);
                        i4 = i2;
                        for i3 in 0..6 {
                            //i4 = i3 * nd_dof + i2;
                            pre.bmat[i4].set_val_dfd0(& strain[i3]);
                            i4  +=  nd_dof;
                        }
                    }
                }
                else {
                    for i2 in 0..nd_dof {
                        self.get_section_def_dfd0(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  cmd.nonlinear_geom,  i2,  MAX_INT);
                        i4 = i2;
                        for i3 in 0..6 {
                            pre.bmat[i4].set_val_dfd0(& sec_def[i3]);
                            i4  +=  nd_dof;
                        }
                    }
                }
                i5 = 0;
                for i2 in 0..6 {
                    bvel[i2].set_val(0.0);
                    i4 = 0;
                    for _i3 in 0..nd_dof {
                        nd = self.dof_table[i4];
                        dof = self.dof_table[i4 + 1];
                        tmp.set_val_dfd0(& pre.bmat[i5]);
                        tmp.mult(& pre.glob_vel[dof * self.num_nds + nd]);
                        bvel[i2].add(& tmp);
                        i4  +=  2;
                        i5 += 1usize;
                    }
                }
                ar_to_vec_dfd0(&mut bvel, &mut  scr_v1,  0,  9);
                mat_mul_dfd0(&mut scr_v2, &mut  pre.dmat, &mut  scr_v1,  self.def_dim,  self.def_dim,  1);
                mat_mul_dfd0(&mut scr_v1, &mut  scr_v2, &mut  pre.bmat,  1,  self.def_dim,  nd_dof);
                for i2 in 0..nd_dof {
                    scr_v1[i2].mult(& wt_det_j);
                    rvec[i2].add(& scr_v1[i2]);
                }
                if get_matrix {
                    mat_mul_dfd0(&mut pre.cbmat, &mut  pre.dmat, &mut  pre.bmat,  self.def_dim,  self.def_dim,  nd_dof);
                    i3 = nd_dof * self.def_dim;
                    for i2 in 0..i3 {
                        pre.cbmat[i2].mult(& wt_det_j);
                    }
                    i5 = 0;
                    for i2 in 0..nd_dof {
                        for i3 in 0..nd_dof {
                            i6 = i2;
                            i7 = i3;
                            for _i4 in 0..6 {
                                //i5 = i2 * nd_dof + i3;
                                //i6 = i4 * nd_dof + i2;
                                //i7 = i4 * nd_dof + i3;
                                d_rd_v[i5]  +=  pre.bmat[i6].val * pre.cbmat[i7].val;
                                i6  +=  nd_dof;
                                i7  +=  nd_dof;
                            }
                            i5 += 1usize;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_ru_dfd0(&mut self, glob_r : &mut Vec<DiffDoub0>, globd_rdu : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub0StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut dof : usize;
        let mut dof2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let nd_dof : usize;
        let tot_dof : usize;
        let mut temp_acc = [DiffDoub0::new(); 30];
        //let mut rvec = &mut scr.scr_v1;
        let mut rvec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //let mut rtmp = &mut scr.scr_v2;
        let mut rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rdu[1089];
        //let mut d_rdu = &mut scr.scr_m1;
        let mut d_rdu = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rd_t[330];
        //let mut d_rd_t = &mut scr.scr_m2;
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rtmp[1089];
        //let mut d_rtmp = &mut scr.scr_m3;
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut c1 : f64;
        let c2 : f64;
        
        nd_dof = self.num_nds*self.dof_per_nd;
        tot_dof = nd_dof + self.num_int_dof;
        
        for i1 in 0..tot_dof {
            rvec[i1].set_val(0.0);
            if get_matrix && i1 < nd_dof {
                i3 = i1 * tot_dof;
                for _i2 in 0..tot_dof {
                    d_rdu[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 21 {
            self.get_ru_frc_fld_dfd0(glob_r, globd_rdu,  get_matrix, cmd, pre, nd_ar);
            return;
        }
        
        if self.this_type != 1 {
            self.get_ruk_dfd0(&mut rvec, &mut  d_rdu, &mut  d_rd_t,  get_matrix,  cmd.nonlinear_geom, pre);
        }
        
        if self.num_int_dof > 0 {
            i2 = nd_dof;
            for i1 in 0..self.num_int_dof {
                self.internal_ru[i1].set_val_dfd0(& rvec[i2]);
                i2 += 1usize;
            }
            if get_matrix {
                i4 = 0;
                for i1 in 0..tot_dof {
                    i3 = i1 * tot_dof + nd_dof;
                    for _i2 in nd_dof..tot_dof {
                        self.internal_mat[i4] = d_rdu[i3];
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
                //condense matrix, d_rd_t and d_rtmp used for scratch
                self.condense_mat(&mut d_rdu, &mut d_rd_t, &mut d_rtmp);
            }
        }
        
        if cmd.dynamic {
            if get_matrix {
                if cmd.lump_mass {
                    i2 = 0;
                    for i1 in 0..nd_dof {
                        nd = self.dof_table[i2];
                        dof = self.dof_table[i2 + 1];
                        i3 = dof * self.num_nds + nd;
                        temp_acc[i1].set_val_dfd0(& pre.glob_acc[i3]);
                        pre.glob_acc[i3].set_val(1.0);
                        i2  +=  2;
                    }
                    self.get_rum_dfd0(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    c1 = cmd.time_step;
                    c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
                    i2 = 0;
                    i3 = tot_dof + 1;
                    for i1 in 0..nd_dof {
                        temp_acc[i1].mult(& rtmp[i1]);
                        rvec[i1].add(& temp_acc[i1]);
                        d_rdu[i2]  +=  c1 * rtmp[i1].val;
                        i2  +=  i3;
                    }
                }
                else {
                    self.get_rum_dfd0(&mut rtmp, &mut d_rtmp,  true,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                    c1 = cmd.time_step;
                    c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
                    i3 = 0;
                    i4 = 0;
                    for _i1 in 0..nd_dof {
                        for _i2 in 0..nd_dof {
                            d_rdu[i3]  +=  c1 * d_rtmp[i4];
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                        i3  +=  self.num_int_dof;
                    }
                }
                if self.this_type != 1 {
                    self.get_rud_dfd0(&mut rtmp, &mut  d_rtmp,  true, cmd, pre, scr, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                    c2 = cmd.time_step * cmd.newmark_gamma * c1;
                    i3 = 0;
                    i4 = 0;
                    for _i1 in 0..nd_dof {
                        for _i2 in 0..nd_dof {
                            d_rdu[i3]  +=  c2 * d_rtmp[i4];
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                        i3  +=  self.num_int_dof;
                    }
                }
            }
            else {
                if cmd.lump_mass {
                    i2 = 0;
                    for i1 in 0..nd_dof {
                        nd = self.dof_table[i2];
                        dof = self.dof_table[i2 + 1];
                        i3 = dof * self.num_nds + nd;
                        temp_acc[i1].set_val_dfd0(& pre.glob_acc[i3]);
                        pre.glob_acc[i3].set_val(1.0);
                        i2  +=  2;
                    }
                    self.get_rum_dfd0(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        temp_acc[i1].mult(& rtmp[i1]);
                        rvec[i1].add(& temp_acc[i1]);
                    }
                }
                else {
                    self.get_rum_dfd0(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                }
                if self.this_type != 1 {
                    self.get_rud_dfd0(&mut rtmp, &mut  d_rtmp,  false,  cmd,  pre, scr, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                }
            }
        }
        
        i4 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i4]];
            dof = self.dof_table[i4+1];
            glob_ind = nd_ar[nd].dof_index[dof];
            glob_r[glob_ind].add(& rvec[i1]);
            if get_matrix {
                i3 = i1*tot_dof;
                i5 = 0;
                for _i2 in 0..nd_dof {
                    nd2 = self.nodes[self.dof_table[i5]];
                    dof2 = self.dof_table[i5+1];
                    glob_ind2 = nd_ar[nd2].dof_index[dof2];
                    globd_rdu.add_entry(glob_ind,   glob_ind2,   d_rdu[i3]);
                    i3 += 1usize;
                    i5  +=  2;
                }
            }
            i4 +=  2;
        }
        
        return;
    }

    pub fn get_rtk_dfd0(&mut self, rvec : &mut Vec<DiffDoub0>, d_rd_t : &mut Vec<f64>, get_matrix : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut d_nt = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut grad_t = [DiffDoub0::new(); 3];
        let mut q_vec = [DiffDoub0::new(); 3];
        let mut rtmp = [DiffDoub0::new(); 10];
        let mut d_rtmp = [DiffDoub0::new(); 100];
        let mut tmp_s : [f64; 3] = [0f64; 3];
        let mut tmp_t = [DiffDoub0::new(); 10];
        let mut tmp_tc = [DiffDoub0::new(); 9];
        let mut tmp_ar = [DiffDoub0::new(); 30];
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..self.num_nds {
                    d_rd_t[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            tmp.set_val(self.ip_wt[i1]);
            det_j.mult(& tmp);
            vec_to_ar_dfd0(&mut tmp_t, &mut  pre.glob_temp,  0,  10);
            mat_mul_ar_dfd0(&mut grad_t, &mut  tmp_t, &mut  d_ndx,  1,  self.num_nds,  3);
            vec_to_ar_dfd0(&mut tmp_tc, &mut  pre.tcmat,  0,  9);
            mat_mul_ar_dfd0(&mut q_vec, &mut  tmp_tc, &mut  grad_t,  3,  3,  1);
            mat_mul_ar_dfd0(&mut rtmp, &mut  d_ndx, &mut  q_vec,  self.num_nds,  3,  1);
            for i2 in 0..self.num_nds {
                rtmp[i2].mult(& det_j);
                rvec[i2].add(& rtmp[i2]);
            }
            if get_matrix {
                transpose_ar_dfd0(&mut d_nt, &mut  d_ndx,  self.num_nds, 3);
                mat_mul_ar_dfd0(&mut tmp_ar, &mut  tmp_tc, &mut  d_nt,  3,  3,  self.num_nds);
                mat_mul_ar_dfd0(&mut d_rtmp, &mut  d_ndx, &mut  tmp_ar,  self.num_nds,  3,  self.num_nds);
                i3 = self.num_nds * self.num_nds;
                for i2 in 0..i3 {
                    d_rd_t[i2]  +=  d_rtmp[i2].val*det_j.val;
                }
            }
        }
        
        return;
    }

    pub fn get_rtm_dfd0(& self, rvec : &mut Vec<DiffDoub0>, d_rd_tdot : &mut Vec<f64>, get_matrix : bool, actual_props : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut pt_tdot = DiffDoub0::new();
        let mut d_rtmp = [DiffDoub0::new(); 100];
        let mut save_cp = DiffDoub0::new();
        let mut tmp_s : [f64; 3] = [0f64; 3];
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..self.num_nds {
                    d_rd_tdot[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 21 {
            return;
        }
        else if self.this_type == 1 {
            tmp.set_val_dfd0(& pre.mass_per_el);
            tmp.mult(& pre.spec_heat);
            d_rd_tdot[0] = tmp.val;
            tmp.mult(& pre.glob_temp[0]);
            rvec[0].set_val_dfd0(& tmp);
            return;
        }
        
        if !actual_props {
            save_cp.set_val_dfd0(& pre.spec_heat);
            pre.spec_heat.set_val(1.0);
        }
        else if self.dof_per_nd == 3 {
            pre.spec_heat.mult(& pre.mmat[0]);
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            tmp.set_val(self.ip_wt[i1]);
            det_j.mult(& tmp);
            pt_tdot.set_val(0.0);
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd0(& n_vec[i2]);
                tmp.mult(& pre.glob_tdot[i2]);
                pt_tdot.add(& tmp);
            }
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd0(& n_vec[i2]);
                tmp.mult(& pre.spec_heat);
                tmp.mult(& pt_tdot);
                tmp.mult(& det_j);
                rvec[i2].add(& tmp);
            }
            if get_matrix {
                mat_mul_ar_dfd0(&mut d_rtmp, & n_vec, & n_vec, self.num_nds, 1, self.num_nds);
                i3 = self.num_nds * self.num_nds;
                for i2 in 0..i3 {
                    d_rd_tdot[i2]  +=  d_rtmp[i2].val * pre.spec_heat.val * det_j.val;
                }
            }
        }
        
        if !actual_props {
            pre.spec_heat.set_val_dfd0(& save_cp);
        }
        else if self.dof_per_nd == 3 {
            pre.spec_heat.dvd(& pre.mmat[0]);
        }
        
        return;
    }

    pub fn get_rt_dfd0(&mut self, glob_r : &mut Vec<DiffDoub0>, globd_rd_t : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub0StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut glob_ind1 : usize;
        let mut glob_ind2 : usize;
        let c1 : f64 =  1.0/(cmd.time_step*cmd.newmark_gamma);
        //let mut rvec = &mut scr.scr_v1;
        let mut rvec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //let mut d_rd_t = &mut scr.scr_m1;
        //let mut rtmp = &mut scr.scr_v2;
        //let mut d_rtmp = &mut scr.scr_m2;
        
        if self.this_type == 21 {
            //self.get_rt_frc_fld_dfd0(glob_r, globd_rd_t, &mut  scr.scr_m1, &mut  scr.scr_m2,  get_matrix, cmd, &mut  pre, &mut  nd_ar);
            self.get_rt_frc_fld_dfd0(glob_r, globd_rd_t, &mut  d_rd_t, &mut  d_rtmp,  get_matrix, pre, nd_ar);
            return;
        }
        
        if self.this_type != 1 {
            self.get_rtk_dfd0(&mut rvec, &mut  d_rd_t,  get_matrix,  pre);
        }
        
        if cmd.dynamic {
            self.get_rtm_dfd0(&mut rtmp, &mut  d_rtmp,  get_matrix,  true,  pre);
            i3 = 0;
            for i1 in 0..self.num_nds {
                rvec[i1].add(& rtmp[i1]);
                if get_matrix {
                    for _i2 in 0..self.num_nds {
                        d_rd_t[i3]  +=  c1 * d_rtmp[i3];
                        i3 += 1usize;
                    }
                }
            }
        }
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            glob_ind1 = nd_ar[self.nodes[i1]].sorted_rank;
            glob_r[glob_ind1].add(& rvec[i1]);
            if get_matrix {
                for i2 in 0..self.num_nds {
                    glob_ind2 = nd_ar[self.nodes[i2]].sorted_rank;
                    globd_rd_t.add_entry(glob_ind1,   glob_ind2,   d_rd_t[i3]);
                    i3 += 1usize;
                }
            }
        }
        
        return;
    }

    pub fn get_ru_frc_fld_dfd0(&mut self, glob_r : &mut Vec<DiffDoub0>, globd_rdu : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub0StressPrereq, nd_ar : &mut Vec<Node>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let nd_dof : usize =  6;
        let tot_dof : usize =  6;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut dof : usize;
        let mut dof2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let mut rvec = [DiffDoub0::new(); 6];
        let mut d_rd_u : [f64; 36] = [0f64; 36];
        let mut d_vec = [DiffDoub0::new(); 3];
        let mut dist = DiffDoub0::new();
        let mut dv_vec = [DiffDoub0::new(); 3];
        let mut d_dvecd_u = [DiffDoub0::new(); 18];
        let mut d_distd_u = [DiffDoub0::new(); 6];
        let mut f_n1 = [DiffDoub0::new(); 3];
        let mut df_n1d_u = [DiffDoub0::new(); 18];
        let mut dto_p = DiffDoub0::new();
        let mut dto_p1 = DiffDoub0::new();
        let mut dto_p2 = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        
        // potential force
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd0(& pre.glob_nds[i2]);
            tmp.add(& pre.glob_disp[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_nds[i2]);
            tmp.sub(& pre.glob_disp[i2]);
            d_vec[i1].set_val_dfd0(& tmp);
        }
        
        dist.set_val_dfd0(& d_vec[0]);
        dist.sqr();
        tmp.set_val_dfd0(& d_vec[1]);
        tmp.sqr();
        dist.add(& tmp);
        tmp.set_val_dfd0(& d_vec[2]);
        tmp.sqr();
        dist.add(& tmp);
        dist.sqt();
        
        d_dvecd_u[0].set_val(-1.0);
        d_dvecd_u[3].set_val(1.0);
        d_dvecd_u[7].set_val(-1.0);
        d_dvecd_u[10].set_val(1.0);
        d_dvecd_u[14].set_val(-1.0);
        d_dvecd_u[17].set_val(1.0);
        
        mat_mul_ar_dfd0(&mut d_distd_u, &mut  d_vec, &mut  d_dvecd_u,  1,  3,  6);
        tmp.set_val(1.0);
        tmp.dvd(& dist);
        for i1 in 0..6 {
            d_distd_u[i1].mult(& tmp);
        }
        
        dto_p.set_val_dfd0(& dist);
        i1 = 1;
        while i1 < (pre.frc_fld_exp[0].val as usize) {
            dto_p.mult(& dist);
            i1 += 1usize;
        }
        dto_p1.set_val_dfd0(& dto_p);
        dto_p1.mult(& dist);
        dto_p2.set_val_dfd0(& dto_p1);
        dto_p2.mult(& dist);
        
        tmp.set_val_dfd0(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd0(& d_vec[i1]);
            f_n1[i1].mult(& tmp);
        }
        
        mat_mul_ar_dfd0(&mut df_n1d_u, &mut  d_vec, &mut  d_distd_u,  3,  1,  6);
        tmp.set_val(1.0);
        tmp.add(& pre.frc_fld_exp[0]);
        tmp.neg();
        tmp.mult(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p2);
        for i1 in 0..18 {
            df_n1d_u[i1].mult(& tmp);
        }
        
        tmp.set_val_dfd0(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..18 {
            tmp2.set_val_dfd0(& tmp);
            tmp2.mult(& d_dvecd_u[i1]);
            df_n1d_u[i1].add(& tmp2);
        }
        
        for i1 in 0..3 {
            rvec[i1].set_val_dfd0(& f_n1[i1]);
            rvec[i1].neg();
            rvec[i1 + 3].set_val_dfd0(& f_n1[i1]);
        }
        
        for i1 in 0..18 {
            d_rd_u[i1] = -df_n1d_u[i1].val;
            d_rd_u[i1 + 18] = df_n1d_u[i1].val;
        }
        
        // damping force
        
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd0(& pre.glob_vel[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_vel[i2]);
            dv_vec[i1].set_val_dfd0(& tmp);
        }
        
        tmp.set_val_dfd0(& pre.frc_fld_coef[1]);
        tmp.dvd(& dto_p);
        
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd0(& tmp);
            f_n1[i1].mult(& dv_vec[i1]);
        }
        
        
        tmp2.set_val(-cmd.newmark_gamma / (cmd.time_step * (cmd.newmark_beta - cmd.newmark_gamma)));
        tmp.mult(& tmp2);
        
        for i1 in 0..18 {
            tmp2.set_val_dfd0(& tmp);
            tmp2.mult(& d_dvecd_u[i1]);
            //df_n1d_u[i1].add(tmp2);
            df_n1d_u[i1].set_val_dfd0(& tmp2);
        }
        
        for i1 in 0..3 {
            rvec[i1].sub(& f_n1[i1]);
            rvec[i1 + 3].add(& f_n1[i1]);
        }
        
        for i1 in 0..18 {
            d_rd_u[i1]  -=  df_n1d_u[i1].val;
            d_rd_u[i1 + 18]  +=  df_n1d_u[i1].val;
        }
        
        i4 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i4]];
            dof = self.dof_table[i4 + 1];
            glob_ind = nd_ar[nd].dof_index[dof];
            glob_r[glob_ind].add(& rvec[i1]);
            if get_matrix {
                i3 = i1 * tot_dof;
                i5 = 0;
                for _i2 in 0..nd_dof {
                    nd2 = self.nodes[self.dof_table[i5]];
                    dof2 = self.dof_table[i5 + 1];
                    glob_ind2 = nd_ar[nd2].dof_index[dof2];
                    globd_rdu.add_entry(glob_ind,   glob_ind2,   d_rd_u[i3]);
                    i3 += 1usize;
                    i5  +=  2;
                }
            }
            i4  +=  2;
        }
        
        return;
    }

    pub fn get_rt_frc_fld_dfd0(&mut self, glob_r : &mut Vec<DiffDoub0>, globd_rd_t : &mut SparseMat, d_rd_u : &mut Vec<f64>, d_rd_v : &mut Vec<f64>, get_matrix : bool, pre : &mut DiffDoub0StressPrereq, nd_ar : &mut Vec<Node>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let mut rvec = [DiffDoub0::new(); 2];
        let mut d_rd_t : [f64; 4] = [0f64; 4];
        let mut d_vec = [DiffDoub0::new(); 3];
        let mut dist = DiffDoub0::new();
        let mut dv_vec = [DiffDoub0::new(); 3];
        let mut d_dvecd_u = [DiffDoub0::new(); 18];
        let mut d_distd_u = [DiffDoub0::new(); 6];
        let mut f_n1 = [DiffDoub0::new(); 3];
        let mut df_n1d_u = [DiffDoub0::new(); 18];
        let mut df_n1d_v = [DiffDoub0::new(); 18];
        let mut dto_p = DiffDoub0::new();
        let mut dto_p1 = DiffDoub0::new();
        let mut dto_p2 = DiffDoub0::new();
        let mut t1to3 = DiffDoub0::new();
        let mut t1to4 = DiffDoub0::new();
        let mut t2to3 = DiffDoub0::new();
        let mut t2to4 = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        let mut tmp_mat = [DiffDoub0::new(); 18];
        let mut tmp_mat2 = [DiffDoub0::new(); 18];
        
        // potential force
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd0(& pre.glob_nds[i2]);
            tmp.add(& pre.glob_disp[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_nds[i2]);
            tmp.sub(& pre.glob_disp[i2]);
            d_vec[i1].set_val_dfd0(& tmp);
        }
        
        dist.set_val_dfd0(& d_vec[0]);
        dist.sqr();
        tmp.set_val_dfd0(& d_vec[1]);
        tmp.sqr();
        dist.add(& tmp);
        tmp.set_val_dfd0(& d_vec[2]);
        tmp.sqr();
        dist.add(& tmp);
        dist.sqt();
        
        dto_p.set_val_dfd0(& dist);
        i1 = 1;
        while i1 < (pre.frc_fld_exp[0].val as usize) {
            dto_p.mult(& dist);
            i1 += 1usize;
        }
        dto_p1.set_val_dfd0(& dto_p);
        dto_p1.mult(& dist);// ||d||^(p+1)
        dto_p2.set_val_dfd0(& dto_p1);
        dto_p2.mult(& dist);// ||d||^(p+2)
        
        tmp.set_val_dfd0(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd0(& d_vec[i1]);
            f_n1[i1].mult(& tmp);
        }
        
        // damping force
        
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd0(& pre.glob_vel[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_vel[i2]);
            dv_vec[i1].set_val_dfd0(& tmp);
        }
        
        tmp.set_val_dfd0(& pre.frc_fld_coef[1]);
        tmp.dvd(& dto_p);// tmp = c_d/||d||^p
        
        for i1 in 0..3 {
            tmp2.set_val_dfd0(& tmp);
            tmp2.mult(& dv_vec[i1]);
            f_n1[i1].add(& tmp2);
        }
        
        // calculate absolute temps
        tmp.set_val_dfd0(& pre.glob_temp[0]);
        tmp.add(& pre.ref_temp);
        t1to4.set_val_dfd0(& tmp);
        t1to4.sqr();
        t1to4.sqr();
        t1to3.set_val_dfd0(& t1to4);
        t1to3.dvd(& tmp);
        
        tmp.set_val_dfd0(& pre.glob_temp[1]);
        tmp.add(& pre.ref_temp);
        t2to4.set_val_dfd0(& tmp);
        t2to4.sqr();
        t2to4.sqr();
        t2to3.set_val_dfd0(& t2to4);
        t2to3.dvd(& tmp);
        
        // conduction and radiation terms
        tmp.set_val_dfd0(& pre.glob_temp[0]);
        tmp.sub(& pre.glob_temp[1]);
        tmp.mult(& pre.thrm_fld_coef[0]);
        tmp.dvd(& dist);
        rvec[0].set_val_dfd0(& tmp);
        
        tmp.set_val_dfd0(& t1to4);
        tmp.sub(& t2to4);
        tmp.mult(& pre.thrm_fld_coef[1]);
        tmp.dvd(& dist);
        tmp.dvd(& dist);
        rvec[0].add(& tmp);
        
        tmp.set_val_dfd0(& rvec[0]);
        rvec[1].set_val_dfd0(& tmp);
        rvec[1].neg();
        
        // work dissipation term
        
        tmp.set_val_dfd0(& f_n1[0]);
        tmp.mult(& dv_vec[0]);
        tmp2.set_val_dfd0(& f_n1[1]);
        tmp2.mult(& dv_vec[1]);
        tmp.add(& tmp2);
        tmp2.set_val_dfd0(& f_n1[2]);
        tmp2.mult(& dv_vec[2]);
        tmp.add(& tmp2);
        tmp2.set_val(0.5);
        tmp.mult(& tmp2);// tmp = 0.5*dot(fn1,d_v)
        
        rvec[0].sub(& tmp);
        rvec[1].sub(& tmp);
        for i1 in 0..2 {
            nd = self.nodes[i1];
            glob_ind = nd_ar[nd].sorted_rank;
            glob_r[glob_ind].add(& rvec[i1]);
        }
        
        if get_matrix {
            d_dvecd_u[0].set_val(-1.0);
            d_dvecd_u[3].set_val(1.0);
            d_dvecd_u[7].set_val(-1.0);
            d_dvecd_u[10].set_val(1.0);
            d_dvecd_u[14].set_val(-1.0);
            d_dvecd_u[17].set_val(1.0);
            
            mat_mul_ar_dfd0(&mut d_distd_u, &mut  d_vec, &mut  d_dvecd_u,  1,  3,  6);
            tmp.set_val(1.0);
            tmp.dvd(& dist);
            for i1 in 0..6 {
                d_distd_u[i1].mult(& tmp);
            }
            
            mat_mul_ar_dfd0(&mut df_n1d_u, &mut  d_vec, &mut  d_distd_u,  3,  1,  6);
            tmp.set_val(1.0);
            tmp.add(& pre.frc_fld_exp[0]);
            tmp.neg();
            tmp.mult(& pre.frc_fld_coef[0]);
            tmp.dvd(& dto_p2);
            for i1 in 0..18 {
                df_n1d_u[i1].mult(& tmp);
            }
            
            tmp.set_val_dfd0(& pre.frc_fld_coef[0]);
            tmp.dvd(& dto_p1);
            for i1 in 0..18 {
                tmp2.set_val_dfd0(& tmp);
                tmp2.mult(& d_dvecd_u[i1]);
                df_n1d_u[i1].add(& tmp2);
            }
            
            mat_mul_ar_dfd0(&mut tmp_mat, &mut  dv_vec, &mut  d_distd_u,  3,  1,  6);
            tmp.set_val_dfd0(& pre.frc_fld_coef[1]);
            tmp.mult(& pre.frc_fld_exp[1]);
            tmp.dvd(& dto_p1);
            tmp.neg();
            for i1 in 0..18 {
                tmp2.set_val_dfd0(& tmp);
                tmp2.mult(& tmp_mat[i1]);
                df_n1d_u[i1].add(& tmp2);
            }
            
            tmp.set_val_dfd0(& pre.frc_fld_coef[1]);
            tmp.dvd(& dto_p);
            for i1 in 0..18 {
                df_n1d_v[i1].set_val_dfd0(& d_dvecd_u[i1]);
                df_n1d_v[i1].mult(& tmp);
            }
            
            // d_rd_t
            tmp.set_val_dfd0(& pre.thrm_fld_coef[0]);
            tmp.dvd(& dist);
            tmp2.set_val_dfd0(& pre.thrm_fld_coef[1]);
            tmp2.dvd(& dist);
            tmp2.dvd(& dist);
            d_rd_t[0] = tmp.val + tmp2.val * t1to3.val;
            d_rd_t[1] = -tmp.val;
            d_rd_t[2] = -tmp.val;
            d_rd_t[3] = tmp.val + tmp2.val * t2to3.val;
            
            i3 = 0;
            for i1 in 0..2 {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                for i2 in 0..2 {
                    nd2 = self.nodes[i2];
                    glob_ind2 = nd_ar[nd2].sorted_rank;
                    globd_rd_t.add_entry(glob_ind,   glob_ind2,   d_rd_t[i3]);
                    i3 += 1usize;
                }
            }
            
            // d_rd_u
            for i1 in 0..12 {
                d_rd_u[i1] = 0.0;
            }
            
            // conduction term
            tmp.set_val_dfd0(& pre.glob_temp[0]);
            tmp.sub(& pre.glob_temp[1]);
            tmp.mult(& pre.thrm_fld_coef[0]);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.neg();
            tmp_mat[0].set_val_dfd0(& tmp);
            tmp_mat[1].set_val_dfd0(& tmp);
            tmp_mat[1].neg();
            
            mat_mul_ar_dfd0(&mut tmp_mat2, &mut  tmp_mat, &mut  d_distd_u,  2,  1,  6);
            
            for i1 in 0..12 {
                d_rd_u[i1]  +=  tmp_mat2[i1].val;
            }
            
            // radiation term
            tmp.set_val_dfd0(& t1to4);
            tmp.sub(& t2to4);
            tmp.mult(& pre.thrm_fld_coef[1]);
            tmp2.set_val(2.0);
            tmp.mult(& tmp2);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.neg();
            tmp_mat[0].set_val_dfd0(& tmp);
            tmp_mat[1].set_val_dfd0(& tmp);
            tmp_mat[1].neg();
            
            mat_mul_ar_dfd0(&mut tmp_mat2, &mut  tmp_mat, &mut  d_distd_u,  2,  1,  6);
            
            for i1 in 0..12 {
                d_rd_u[i1]  +=  tmp_mat2[i1].val;
            }
            
            // work dissipation term
            mat_mul_ar_dfd0(&mut tmp_mat, &mut  dv_vec, &mut  df_n1d_u,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_u[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
            // d_rd_v
            
            for i1 in 0..12 {
                d_rd_v[i1] = 0.0;
            }
            
            mat_mul_ar_dfd0(&mut tmp_mat, &mut  dv_vec, &mut  df_n1d_v,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_v[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
            mat_mul_ar_dfd0(&mut tmp_mat, &mut  f_n1, &mut  d_dvecd_u,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_v[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
        }
        
        return;
    }

    pub fn get_app_load_dfd0(& self, app_ld : &mut Vec<DiffDoub0>, ld_pt : & Load, n_lgeom : bool, pre : &mut DiffDoub0StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>, sec_ar : &mut Vec<Section>, fc_ar : &mut Vec<Face>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut glob_ind : usize;
        let mut nd : usize;
        let mut dof : usize;
        let nd_dof : usize =  self.num_nds*self.dof_per_nd;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let ld_type : CppStr = ld_pt.this_type.clone();
        //let mut d_rd_a = &mut scr.scr_m1;
        let mut d_rd_a = match scr.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        let mut fc_num_nds : usize;
        let mut nd_in_face : bool;
        //let mut el_app_ld = &mut scr.scr_v1;
        let mut el_app_ld = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut tot_nd_f = [DiffDoub0::new(); 6];
        let mut inp_mag = DiffDoub0::new();
        let mut vec_mag = DiffDoub0::new();
        let mut el_vol = DiffDoub0::new();
        let mut sc_fact = DiffDoub0::new();
        let mut el_cent = [DiffDoub0::new(); 3];
        let mut cent_to_el = [DiffDoub0::new(); 3];
        let mut ax_to_el = [DiffDoub0::new(); 3];
        let mut ang_vel2 = DiffDoub0::new();
        let mut nn_inv = DiffDoub0::new();
        let mut fc_area = DiffDoub0::new();
        let mut fc_norm = [DiffDoub0::new(); 3];
        let mut trac = [DiffDoub0::new(); 3];
        let mut dp = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        
        for i1 in 0..nd_dof {
            pre.glob_acc[i1].set_val(0.0);
        }
        
        if ld_type.s == "bodyForce" {
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                pre.glob_acc[i3].set_val(ld_pt.load[dof]);
                i2  +=  2;
            }
            self.get_rum_dfd0(&mut el_app_ld, &mut  d_rd_a,  false,  false,  n_lgeom, pre, scr_dfd);
            i2 = 0;
            for i1 in 0..nd_dof {
                dof = self.dof_table[i2 + 1];
                tot_nd_f[dof].add(& el_app_ld[i1]);
                i2  +=  2;
            }
            for i1 in 0..self.dof_per_nd {
                tmp.set_val(ld_pt.load[i1]);
                tmp.sqr();
                inp_mag.add(& tmp);
                tmp.set_val_dfd0(& tot_nd_f[i1]);
                tmp.sqr();
                vec_mag.add(& tmp);
            }
            inp_mag.sqt();
            vec_mag.sqt();
            if self.this_type == 41 || self.this_type == 3 {
                el_vol.set_val(0.0);
                for i1 in 0..num_lay {
                    self.get_volume_dfd0(&mut tmp, pre,  i1, sec_ar, dv_ar);
                    el_vol.add(& tmp);
                }
            }
            else {
                self.get_volume_dfd0(&mut el_vol, pre,  0, sec_ar, dv_ar);
            }
            sc_fact.set_val_dfd0(& inp_mag);
            sc_fact.mult(& el_vol);
            sc_fact.dvd(& vec_mag);
            for i1 in 0..nd_dof {
                el_app_ld[i1].mult(& sc_fact);
            }
        }
        else if ld_type.s == "gravitational" {
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                if dof < 3 {
                    pre.glob_acc[i3].set_val(ld_pt.load[dof]);
                }
                i2  +=  2;
            }
            self.get_rum_dfd0(&mut el_app_ld, &mut  d_rd_a,  false,  true,  n_lgeom,  pre, scr_dfd);
        }
        else if ld_type.s == "centrifugal" {
            ang_vel2.set_val(ld_pt.angular_vel);
            ang_vel2.sqr();
            i3 = 0;
            nn_inv.set_val(1.0 / (self.num_nds as f64));
            dp.set_val(0.0);
            for i1 in 0..3 {
                for _i2 in 0..self.num_nds {
                    el_cent[i1].add(& pre.glob_nds[i3]);
                    i3 += 1usize;
                }
                el_cent[i1].mult(& nn_inv);
                cent_to_el[i1].set_val_dfd0(& el_cent[i1]);
                tmp.set_val(ld_pt.center[i1]);
                cent_to_el[i1].sub(& tmp);
                tmp.set_val(ld_pt.axis[i1]);
                tmp.mult(& cent_to_el[i1]);
                dp.add(& tmp);
            }
            for i1 in 0..3 {
                ax_to_el[i1].set_val_dfd0(& cent_to_el[i1]);
                tmp.set_val(ld_pt.axis[i1]);
                tmp.mult(& dp);
                ax_to_el[i1].sub(& tmp);
                ax_to_el[i1].mult(& ang_vel2);
            }
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                if dof < 3 {
                    pre.glob_acc[i3].set_val_dfd0(& ax_to_el[dof]);
                }
                i2  +=  2;
            }
            self.get_rum_dfd0(&mut el_app_ld, &mut  d_rd_a,  false,  true,  n_lgeom, pre, scr_dfd);
        }
        else if ld_type.s == "surfaceTraction" || ld_type.s == "surfacePressure" {
            let mut this_fc : &Face;
            for fi in self.faces.iter() {
                this_fc = &fc_ar[*fi];
                if this_fc.on_surf {
                    this_fc.get_area_normal_dfd0(&mut fc_area, &mut  fc_norm, &pre.glob_nds, self.num_nds);
                    dp.set_val(0.0);
                    for i1 in 0..3 {
                        tmp.set_val(ld_pt.normal_dir[i1]);
                        tmp.mult(& fc_norm[i1]);
                        dp.add(& tmp);
                    }
                    tmp.set_val(R_PI_0180* ld_pt.norm_tol);
                    tmp.cs();
                    if dp.val > tmp.val {
                        if ld_type.s == "surfacePressure" {
                            for i1 in 0..3 {
                                trac[i1].set_val(ld_pt.load[0]);
                                trac[i1].mult(& fc_norm[i1]);
                                trac[i1].neg();
                            }
                        }
                        else {
                            for i1 in 0..3 {
                                trac[i1].set_val(ld_pt.load[i1]);
                            }
                        }
                        fc_num_nds = this_fc.num_nds;
                        for i1 in 0..fc_num_nds {
                            i4 = this_fc.loc_nodes[i1];
                            for i3 in 0..3 {
                                //i4 = i3 * self.num_nds + fc_loc_nd[i1];
                                pre.glob_acc[i4].set_val_dfd0(& trac[i3]);
                                i4  +=  self.num_nds;
                            }
                        }
                        self.get_rum_dfd0(&mut el_app_ld, &mut  d_rd_a,  false,  false,  n_lgeom, pre, scr_dfd);
                        i2 = 0;
                        for i1 in 0..nd_dof {
                            nd = self.dof_table[i2];
                            dof = self.dof_table[i2 + 1];
                            nd_in_face = false;
                            for i3 in 0..fc_num_nds {
                                if this_fc.loc_nodes[i3] == nd {
                                    nd_in_face = true;
                                }
                            }
                            if !nd_in_face {
                                el_app_ld[i1].set_val(0.0);
                            }
                            tot_nd_f[dof].add(& el_app_ld[i1]);
                            i2  +=  2;
                        }
                        inp_mag.set_val(0.0);
                        vec_mag.set_val(0.0);
                        for i1 in 0..3 {
                            tmp.set_val_dfd0(& tot_nd_f[i1]);
                            tmp.sqr();
                            vec_mag.add(& tmp);
                            tmp.set_val_dfd0(& trac[i1]);
                            tmp.sqr();
                            inp_mag.add(& tmp);
                        }
                        inp_mag.sqt();
                        vec_mag.sqt();
                        tmp.set_val_dfd0(& fc_area);
                        tmp.mult(& inp_mag);
                        tmp.dvd(& vec_mag);
                        for i1 in 0..nd_dof {
                            el_app_ld[i1].mult(& tmp);
                        }
                    }
                }
            }
        }
        
        i2 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i2]];
            dof = self.dof_table[i2 + 1];
            glob_ind = nd_ar[nd].dof_index[dof];
            app_ld[glob_ind].add(& el_app_ld[i1]);
            i2  +=  2;
        }
        
        return;
    }

    pub fn get_app_therm_load_dfd0(& self, app_ld : &mut Vec<DiffDoub0>, ld_pt : & Load, pre : &mut DiffDoub0StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>, sec_ar : &mut Vec<Section>, fc_ar : &mut Vec<Face>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut glob_ind : usize;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let ld_type : CppStr = ld_pt.this_type.clone();
        //let mut el_app_ld = &mut scr.scr_v1;
        let el_app_ld = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        //let mut d_rd_t = &mut scr.scr_m1;
        let d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        let mut fc_num_nds : usize;
        let mut nd_in_face : bool;
        let mut fc_area = DiffDoub0::new();
        let mut fc_norm = [DiffDoub0::new(); 3];
        let mut tot_hg = DiffDoub0::new();
        let mut el_vol = DiffDoub0::new();
        let mut dp = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        
        for i1 in 0..self.num_nds {
            pre.glob_tdot[i1].set_val(0.0);
        }
        
        if ld_type.s == "bodyHeatGen" {
            for i1 in 0..self.num_nds {
                pre.glob_tdot[i1].set_val(ld_pt.load[0]);
            }
            self.get_rtm_dfd0(el_app_ld, d_rd_t, false, false, pre);
            tot_hg.set_val(0.0);
            for i1 in 0..self.num_nds {
                tot_hg.add(& el_app_ld[i1]);
            }
            if num_lay > 0 {
                el_vol.set_val(0.0);
                for i1 in 0..num_lay {
                    self.get_volume_dfd0(&mut tmp, pre,  i1, sec_ar, dv_ar);
                    el_vol.add(& tmp);
                }
            }
            else {
                self.get_volume_dfd0(&mut el_vol, pre,  0, sec_ar, dv_ar);
            }
            tmp.set_val(ld_pt.load[0]);
            tmp.mult(& el_vol);
            tmp.dvd(& tot_hg);
            for i1 in 0..self.num_nds {
                el_app_ld[i1].mult(& tmp);
            }
        }
        else if ld_type.s == "surfaceFlux" {
            let mut this_fc : &Face;
            for fi in self.faces.iter() {
                this_fc = &fc_ar[*fi];
                if this_fc.on_surf {
                    this_fc.get_area_normal_dfd0(&mut fc_area, &mut  fc_norm, &pre.glob_nds, self.num_nds);
                    for i1 in 0..3 {
                        tmp.set_val(ld_pt.normal_dir[i1]);
                        tmp.mult(& fc_norm[i1]);
                        dp.add(& tmp);
                    }
                    tmp.set_val(R_PI_0180 * ld_pt.norm_tol);
                    tmp.cs();
                    if dp.val > tmp.val {
                        fc_num_nds = this_fc.num_nds;
                        for i1 in 0..fc_num_nds {
                            i2 = this_fc.loc_nodes[i1];
                            pre.glob_tdot[i2].set_val(ld_pt.load[0]);
                        }
                        self.get_rtm_dfd0(el_app_ld, d_rd_t,  false,  false, pre);
                        tot_hg.set_val(0.0);
                        for i1 in 0..self.num_nds {
                            nd_in_face = false;
                            for i2 in 0..fc_num_nds {
                                if this_fc.loc_nodes[i2] == i1 {
                                    nd_in_face = true;
                                }
                            }
                            if !nd_in_face {
                                el_app_ld[i1].set_val(0.0);
                            }
                            tot_hg.add(& el_app_ld[i1]);
                        }
                        tmp.set_val(ld_pt.load[0]);
                        tmp.mult(& fc_area);
                        tmp.dvd(& tot_hg);
                        for i1 in 0..self.num_nds {
                            el_app_ld[i1].mult(& tmp);
                        }
                    }
                }
            }
        }
        
        for i1 in 0..self.num_nds {
            glob_ind = nd_ar[self.nodes[i1]].sorted_rank;
            app_ld[glob_ind].add(& el_app_ld[i1]);
        }
        
        return;
    }

    pub fn get_rf_unsteady_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq) {
        let mut rho = DiffDoub0::new();
        let mut rho_dot = DiffDoub0::new();
        let mut fl_vel = [DiffDoub0::new(); 3];
        let mut vel_dot = [DiffDoub0::new(); 3];
        let mut temp = DiffDoub0::new();
        let mut t_dot = DiffDoub0::new();
        let mut turb_e = DiffDoub0::new();
        let mut turb_e_dot = DiffDoub0::new();
        let mut enth = DiffDoub0::new();
        let mut enth_dot = DiffDoub0::new();

        let mut si = 0usize;
        let mut vi : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut r_sum = [DiffDoub0::new(); 6];
        let mut ipw_j = DiffDoub0::new();
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        let mut td3 = DiffDoub0::new();

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            rho_dot.set_val(0.0);
            temp.set_val(0.0);
            t_dot.set_val(0.0);
            turb_e.set_val(0.0);
            turb_e_dot.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den_dot[j]);
                rho_dot.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_temp[j]);
                temp.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_tdot[j]);
                t_dot.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_edot[j]);
                turb_e_dot.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                vel_dot[k].set_val(0.0);
                vi = self.num_nds*k;
                for j in 0..self.num_nds {
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel[vi]);
                    fl_vel[k].add(&td);
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel_dot[vi]);
                    vel_dot[k].add(&td);
                    vi += 1;
                }
            }
            enth.set_val_dfd0(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);
            enth_dot.set_val_dfd0(&pre.spec_heat);
            enth_dot.mult(&t_dot);

            //mass
            td.set_val_dfd0(&rho_dot);
            td.mult(&ipw_j);
            r_sum[0].add(&td);

            //momentum
            for j in 0..3 {
                td.set_val_dfd0(&rho_dot);
                td.mult(&fl_vel[j]);
                td2.set_val_dfd0(&rho);
                td2.mult(&vel_dot[j]);
                td.add(&td2);
                td.mult(&ipw_j);
                r_sum[j+1].add(&td);
            }
            
            //energy
            td.set_val(0.5);
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd0(&fl_vel[j]);
                td3.sqr();
                td2.add(&td3);
            }
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_dot); // td = rho_dot(0.5v*v + h + tau)
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd0(&fl_vel[j]);
                td3.mult(&vel_dot[j]);
                td2.add(&td3);
            }
            td2.add(&enth_dot);
            td2.add(&turb_e_dot);
            td2.mult(&rho);  //td2 = rho(v*v_dot + h_dot + tau_dot)
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[4].add(&td);

            //turbulence
            td.set_val_dfd0(&rho_dot);
            td.mult(&turb_e);
            td2.set_val_dfd0(&rho);
            td2.mult(&turb_e_dot);
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[5].add(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_convective_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_convective_dfd0(r_vec,pre);
        }
    }

    pub fn get_rf_viscous_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub0::new();
        let mut turb_e = DiffDoub0::new();
        let mut grad_v = [DiffDoub0::new(); 9];
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();

        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ipw_j = DiffDoub0::new();
        let mut r_sum = DiffDoub0::new();
        let mut si = 0usize;
        let mut gvi : usize;
        let mut nvi : usize;
        let mut gni : usize;

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
            }
            gvi = 0;
            for j in 0..3 {
                for k in 0..3 {
                    grad_v[gvi].set_val(0.0);
                    nvi = j*self.num_nds;
                    gni = k;
                    for _m  in 0..self.num_nds {
                        //nvi = j*self.num_nds + m;
                        //gni = m*3 + k;
                        td.set_val_dfd0(&pre.fl_vel[nvi]);
                        td.mult(&d_ndx[gni]);
                        grad_v[gvi].add(&td);
                        nvi += 1;
                        gni += 3;
                    }
                    gvi += 1;
                }
            }
            td.set_val(0.0);
            for j in 0..9 {
                td2.set_val_dfd0(&grad_v[j]);
                td2.sqr();
                td.add(&td2);
            }
            td.sqt(); // td = mag(grad_v)
            
            td2.set_val_dfd0(&pre.grad_turb_coef);
            td2.mult(&rho);
            td2.mult(&td);
            td2.mult(&ipw_j);
            r_sum.sub(&td2);

            td2.set_val_dfd0(&pre.diss_turb_coef);
            td2.mult(&rho);
            td2.mult(&turb_e);
            td2.mult(&ipw_j);
            r_sum.add(&td2);

            si += 3;
        }

        let mut ri = 5usize;
        for j in 0..self.num_nds {
            r_vec[ri].add(&r_sum);
            ri += 6;
        }

        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_viscous_dfd0(r_vec, pre, fc_vec, ec_vec);
        }

    }

    pub fn get_rf_pressure_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_pressure_dfd0(r_vec,pre, fc_vec, ec_vec);
        }
    }

    pub fn get_rf_load_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq) {
        let mut rho = DiffDoub0::new();
        let mut vel = [DiffDoub0::new(); 3];

        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ipw_j = DiffDoub0::new();
        let mut si = 0usize;
        let mut nvi : usize;
        let mut r_sum = [DiffDoub0::new(); 6];
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        
        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
            }
            nvi = 0;
            for k in 0..3 {
                vel[k].set_val(0.0);
                for j in 0..self.num_nds {
                    //nvi = k*self.num_nds + j;
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel[nvi]);
                    vel[k].add(&td);
                    nvi += 1;
                }
            }

            //momentum
            for j in 1..4 {
                td.set_val_dfd0(&pre.f_per_mass[j-1]);
                td.mult(&rho);
                td.mult(&ipw_j);
                r_sum[j].sub(&td);
            }

            //energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&vel[j]);
                td2.mult(&pre.f_per_mass[j]);
                td.add(&td2);
            }
            td.mult(&rho);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            td.set_val_dfd0(&pre.hg_per_vol);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_dfd0(&self, r_glob : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq, dynamic : bool, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_dfd0"),
            Some(x) => &mut x.dat,
        };
        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd0(r_vec, pre);
        }
        self.get_rf_convective_dfd0( r_vec, pre, fc_vec);
        self.get_rf_viscous_dfd0(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_pressure_dfd0(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_load_dfd0(r_vec,pre);

        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                r_glob[rg+j].add(&r_vec[ri]);
                ri += 1;
            }
        }
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn get_ruk_dfd1(&mut self, rvec : &mut Vec<DiffDoub1>, d_rdu : &mut Vec<f64>, d_rd_t : &mut Vec<f64>, get_matrix : bool, n_lgeom : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let tot_dof : usize;
        
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut d_jwt = DiffDoub1::new();
        
        let mut strain = [DiffDoub1::new(); 6];
        let mut stress = [DiffDoub1::new(); 6];
        let mut thrm_stn = [DiffDoub1::new(); 6];
        let mut ip_temp = DiffDoub1::new();
        let mut cte = [DiffDoub1::new(); 6];
        let mut cten = [DiffDoub1::new(); 90];
        let mut ux = [DiffDoub1::new(); 9];
        let mut sec_def = [DiffDoub1::new(); 9];
        let mut sec_fc_mom = [DiffDoub1::new(); 9];
        
        let mut tmp = DiffDoub1::new();
        let mut tmp_c = [DiffDoub1::new(); 81];
        let mut tmp_gd = [DiffDoub1::new(); 60];
        let mut tmp_te = [DiffDoub1::new(); 6];
        
        vec_to_ar_dfd1(&mut tmp_c, &mut  pre.cmat,  0,  81);
        vec_to_ar_dfd1(&mut tmp_gd, &mut  pre.glob_disp,  0,  60);
        vec_to_ar_dfd1(&mut tmp_te, &mut  pre.therm_exp,  0,  6);
        
        tot_dof = self.num_nds*self.dof_per_nd + self.num_int_dof;
        
        i3 = 0;
        i4 = 0;
        for i1 in 0..tot_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..tot_dof {
                    d_rdu[i3] = 0.0;
                    i3 += 1usize;
                }
                for _i2 in 0..self.num_nds {
                    d_rd_t[i4] = 0.0;
                    i4 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        if self.dof_per_nd == 6 {
            if n_lgeom {
                self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
        }
        
        for i1 in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.loc_nds, & self.int_pts[(3*i1)..]);
            d_jwt.set_val_dfd1(& det_j);
            tmp.set_val(self.ip_wt[i1]);
            d_jwt.mult(& tmp);
            ip_temp.set_val(0.0);
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd1(& pre.glob_temp[i2]);
                tmp.mult(& n_vec[i2]);
                ip_temp.add(& tmp);
            }
            if self.dof_per_nd == 3 {
                mat_mul_ar_dfd1(&mut ux, &mut tmp_gd, &mut d_ndx, 3, self.n_dim, 3);
                self.get_solid_strain_dfd1(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, MAX_INT, MAX_INT, n_lgeom);
                for i2 in 0..6 {
                    thrm_stn[i2].set_val_dfd1(& pre.therm_exp[i2]);
                    thrm_stn[i2].mult(& ip_temp);
                    strain[i2].sub(& thrm_stn[i2]);
                    strain[i2].sub(& pre.einit[i2]);
                }
                mat_mul_ar_dfd1(&mut stress, &mut tmp_c, &mut strain, 6, 6, 1);
                if get_matrix {
                    mat_mul_ar_dfd1(&mut cte, &mut  tmp_c, &mut  tmp_te,  6,  6,  1);
                    mat_mul_ar_dfd1(&mut cten, &mut  cte, &mut  n_vec,  6,  1,  self.num_nds);
                }
                for i2 in 0..tot_dof {
                    self.get_solid_strain_dfd1(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, i2, MAX_INT, n_lgeom);
                    i4 = i2;
                    for i3 in 0..6 {
                        tmp.set_val_dfd1(& stress[i3]);
                        tmp.mult(& strain[i3]);
                        tmp.mult(& d_jwt);
                        rvec[i2].add(& tmp);
                        if get_matrix {
                            pre.bmat[i4].set_val_dfd1(& strain[i3]);
                            i4 +=  tot_dof;
                        }
                    }
                    if n_lgeom && get_matrix {
                        i5 = (tot_dof + 1)*i2;
                        i6 = i5;
                        for i3 in i2..tot_dof {
                            self.get_solid_strain_dfd1(&mut strain, &mut ux, &mut d_ndx, &mut pre.loc_ori, i2, i3, n_lgeom);
                            for i4 in 0..6 {
                                d_rdu[i5] +=  stress[i4].val*strain[i4].val*d_jwt.val;
                            }
                            d_rdu[i6] = d_rdu[i5];
                            i5 += 1usize;
                            i6 +=  tot_dof;
                        }
                    }
                }
            } else {
                self.get_section_def_dfd1(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, MAX_INT, MAX_INT);
                mat_mul_ar_dfd1(&mut sec_fc_mom, &mut tmp_c, &mut sec_def, self.def_dim, self.def_dim, 1);
                for i2 in 0..6 {
                    tmp.set_val_dfd1(& pre.therm_exp[i2]);
                    tmp.mult(& ip_temp);
                    tmp.add(& pre.einit[i2]);
                    sec_fc_mom[i2].sub(& tmp);
                }
                if get_matrix {
                    mat_mul_ar_dfd1(&mut cten, &mut  tmp_te, &mut  n_vec,  6,  1,  self.num_nds);
                }
                for i2 in 0..tot_dof {
                    self.get_section_def_dfd1(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, i2, MAX_INT);
                    i4 = i2;
                    for i3 in 0..self.def_dim {
                        tmp.set_val_dfd1(& sec_fc_mom[i3]);
                        tmp.mult(& sec_def[i3]);
                        tmp.mult(& d_jwt);
                        rvec[i2].add(& tmp);
                        if get_matrix {
                            pre.bmat[i4].set_val_dfd1(& sec_def[i3]);
                            i4 +=  tot_dof;
                        }
                    }
                    if n_lgeom && get_matrix {
                        i5 = (tot_dof + 1)*i2;
                        i6 = i5;
                        for i3 in i2..tot_dof {
                            self.get_section_def_dfd1(&mut sec_def, &mut pre.glob_disp, &mut pre.inst_ori, &mut pre.loc_ori, &mut pre.glob_nds, &mut d_ndx, &mut n_vec, n_lgeom, i2, i3);
                            for i4 in 0..self.def_dim {
                                d_rdu[i5] +=  sec_fc_mom[i4].val*sec_def[i4].val*d_jwt.val;
                            }
                            d_rdu[i6] = d_rdu[i5];
                            i5 += 1usize;
                            i6 +=  tot_dof;
                        }
                    }
                }
            }
            if get_matrix {
                mat_mul_dfd1(&mut pre.cbmat, &mut pre.cmat, &mut pre.bmat, self.def_dim, self.def_dim, tot_dof);
                i3 = self.def_dim*tot_dof;
                for i2 in 0..i3 {
                    pre.cbmat[i2].mult(& d_jwt);
                }
                i5 = 0;
                for i2 in 0..tot_dof {
                    for i3 in 0..tot_dof {
                        i6 = i2;
                        i7 = i3;
                        for _i4 in 0..self.def_dim {
                            d_rdu[i5] +=  pre.bmat[i6].val*pre.cbmat[i7].val;
                            i6 +=  tot_dof;
                            i7 +=  tot_dof;
                        }
                        i5 += 1usize;
                    }
                }
                for i2 in 0..6 * self.num_nds {
                    cten[i2].mult(& d_jwt);
                }
                i5 = 0;
                for i2 in 0..tot_dof {
                    for i3 in 0..self.num_nds {
                        i6 = i2;
                        i7 = i3;
                        for _i4 in 0..6 {
                            d_rd_t[i5]  -=  pre.bmat[i6].val * cten[i7].val;
                            i6  +=  tot_dof;
                            i7  +=  self.num_nds;
                        }
                        i5 += 1usize;
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_rum_dfd1(& self, rvec : &mut Vec<DiffDoub1>, d_rd_a : &mut Vec<f64>, get_matrix : bool, actual_props : bool, n_lgeom : bool, 
        pre : &mut DiffDoub1StressPrereq, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {
        
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd1 : usize;
        let mut dof1 : usize;
        let mut nd2 : usize;
        let mut dof2 : usize;
        let nd_dof : usize =  self.num_nds * self.dof_per_nd;
        
        let mut inst_disp = [DiffDoub1::new(); 60];
        
        let mut tmp_s : [f64; 3] = [0f64; 3];
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        
        let mut tmp = DiffDoub1::new();
        let mut tmp61 = [DiffDoub1::new(); 6];
        let mut det_jwt = DiffDoub1::new();
        
        let mut save_m = [DiffDoub1::new(); 36];
        
        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v2 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        
        i3 = 0;
        for i1 in 0..nd_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..nd_dof {
                    d_rd_a[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 {
            i2 = 0;
            for i1 in 0..3 {
                rvec[i1].set_val_dfd1(& pre.glob_acc[i1]);
                rvec[i1].mult(& pre.mass_per_el);
                if get_matrix {
                    d_rd_a[i2] = pre.mass_per_el.val;
                    i2  +=  4;
                }
            }
            return;
        }
        else if self.this_type == 21 {
            return;
        }
        
        if self.dof_per_nd == 6 {
            if n_lgeom {
                self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  1);
            }
            if !actual_props {
                for i1 in 0..36 {
                    save_m[i1].set_val_dfd1(& pre.mmat[i1]);
                    pre.mmat[i1].set_val(0.0);
                }
                i1 = 0usize;
                while i1 < 36 {
                    pre.mmat[i1].set_val(1.0);
                    i1 += 7;
                }
            }
        }
        else {
            if !actual_props {
                save_m[0].set_val_dfd1(& pre.mmat[0]);
                pre.mmat[0].set_val(1.0);
            }
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            det_jwt.set_val(self.ip_wt[i1]);
            det_jwt.mult(& det_j);
            if self.dof_per_nd == 6 {
                // build matrix [d_u^i/d_u^g] * {n}
                i7 = 0;
                for i2 in 0..nd_dof {
                    //nd1 = self.dof_table[i7];
                    //dof1 = self.dof_table[i7 + 1];
                    self.get_inst_disp_dfd1(&mut inst_disp, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds,  n_lgeom,  i2,  MAX_INT);
                    i6 = 0;
                    i5 = i2;
                    for _i3 in 0..6 {
                        pre.bmat[i5].set_val(0.0);
                        for i4 in 0..self.num_nds {
                            tmp.set_val_dfd1(& n_vec[i4]);
                            tmp.mult(& inst_disp[i6]);
                            pre.bmat[i5].add(& tmp);
                            i6 += 1usize;
                        }
                        i5  +=  nd_dof;
                        i6  +=  self.n_dim - self.num_nds;
                    }
                    i7  +=  2;
                }
                // tmp61 = bmat*{acc}
                i5 = 0;
                for i2 in 0..6 {
                    i4 = 0;
                    tmp61[i2].set_val(0.0);
                    for _i3 in 0..nd_dof {
                        nd1 = self.dof_table[i4];
                        dof1 = self.dof_table[i4 + 1];
                        i6 = dof1 * self.num_nds + nd1;
                        tmp.set_val_dfd1(& pre.bmat[i5]);
                        tmp.mult(& pre.glob_acc[i6]);
                        tmp61[i2].add(& tmp);
                        i4  +=  2;
                        i5 += 1usize;
                    }
                }
                // pre.scr_vec5 = [m][b]{a}, originally tmp62
                ar_to_vec_dfd1(&mut tmp61, &mut scr_v1, 0, 6);
                mat_mul_dfd1(&mut scr_v2, &mut pre.mmat, &mut  scr_v1,  6,  6,  1);
                mat_mul_dfd1(&mut scr_v1, &mut  scr_v2, &mut  pre.bmat,  1,  6,  nd_dof);
                // update rvec
                for i2 in 0..nd_dof {
                    tmp.set_val_dfd1(& scr_v1[i2]);
                    tmp.mult(& det_jwt);
                    rvec[i2].add(& tmp);
                }
                if get_matrix {
                    mat_mul_dfd1(&mut pre.cbmat, &mut  pre.mmat, &mut  pre.bmat,  6,  6,  nd_dof);
                    i4 = 6 * nd_dof;
                    for i2 in 0..i4 {
                        pre.cbmat[i2].mult(& det_jwt);
                    }
                    i5 = 0;
                    for i2 in 0..nd_dof {
                        for i3 in 0..nd_dof {
                            i6 = i2;
                            i7 = i3;
                            for _i4 in 0..6 {
                                d_rd_a[i5]  +=  pre.cbmat[i6].val * pre.bmat[i7].val;
                                i6  +=  nd_dof;
                                i7  +=  nd_dof;
                            }
                            i5 += 1usize;
                        }
                    }
                }
            }
            else {
                ar_to_vec_dfd1(&mut n_vec, &mut  scr_v1,  0,  11);
                mat_mul_dfd1(&mut pre.bmat, &mut  pre.glob_acc, &mut  scr_v1,  3,  self.num_nds,  1);
                for i2 in 0..nd_dof {
                    nd1 = self.dof_table[2 * i2];
                    dof1 = self.dof_table[2 * i2 + 1];
                    tmp.set_val_dfd1(& n_vec[nd1]);
                    tmp.mult(& pre.mmat[0]);
                    tmp.mult(& pre.bmat[dof1]);
                    tmp.mult(& det_jwt);
                    rvec[i2].add(& tmp);
                    if get_matrix {
                        for i3 in 0..nd_dof {
                            nd2 = self.dof_table[2 * i3];
                            dof2 = self.dof_table[2 * i3 + 1];
                            if dof2 == dof1 {
                                tmp.set_val_dfd1(& n_vec[nd1]);
                                tmp.mult(& n_vec[nd2]);
                                tmp.mult(& pre.mmat[0]);
                                tmp.mult(& det_jwt);
                                i4 = i2 * nd_dof + i3;
                                d_rd_a[i4]  +=  tmp.val;
                            }
                        }
                    }
                }
            }
        }
        
        if self.dof_per_nd == 6 {
            if !actual_props {
                for i1 in 0..36 {
                    pre.mmat[i1].set_val_dfd1(& save_m[i1]);
                }
            }
        }
        else {
            if !actual_props {
                pre.mmat[0].set_val_dfd1(& save_m[0]);
            }
        }
        
        return;
    }

    pub fn get_rud_dfd1(&mut self, rvec : &mut Vec<DiffDoub1>, d_rd_v : &mut Vec<f64>, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub1StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let mut dof : usize;
        let nd_dof : usize =  self.num_nds * self.dof_per_nd;
        let mut tmp = DiffDoub1::new();
        //let mut rtmp = &mut scr.scr_v3;
        let rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut scr_v2 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rtmp[1089];
        //let mut d_rtmp = &mut scr.scr_m4;
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rd_t[330];
        //let mut d_rd_t = &mut scr.scr_m5;
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_non_zero : bool;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut wt_det_j = DiffDoub1::new();
        let mut strain = [DiffDoub1::new(); 6];
        let mut sec_def = [DiffDoub1::new(); 9];
        let mut ux = [DiffDoub1::new(); 9];
        let mut bvel = [DiffDoub1::new(); 9];
        
        let mut tmp_s : [f64; 3] = [0f64; 3];
        
        i3 = 0;
        for i1 in 0..nd_dof {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..nd_dof {
                    d_rd_v[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        if cmd.ray_damp_cm > 0.0 {
            tmp.set_val(cmd.ray_damp_cm);
            for i1 in 0..36 {
                pre.mmat[i1].mult(& tmp);
            }
            for i1 in 0..nd_dof {
                pre.glob_acc[i1].set_val_dfd1(& pre.glob_vel[i1]);
            }
            self.get_rum_dfd1(rtmp, d_rtmp,  get_matrix,  true,  cmd.nonlinear_geom, pre, scr_dfd);
            for i1 in 0..nd_dof {
                rvec[i1].add(& rtmp[i1]);
            }
            if get_matrix {
                i3 = 0;
                for _i1 in 0..nd_dof {
                    for _i2 in 0..nd_dof {
                        d_rd_v[i3]  +=  d_rtmp[i3];
                        i3 += 1usize;
                    }
                }
            }
        }
        if cmd.ray_damp_ck > 0.0 {
            tmp.set_val(cmd.ray_damp_ck);
            for i1 in 0..81 {
                pre.cmat[i1].mult(& tmp);
            }
            i3 = 0;
            i4 = 0;
            for _i1 in 0..self.dof_per_nd {
                for i2 in 0..self.n_dim {
                    if i2 >= self.num_nds {
                        pre.glob_disp[i3].set_val(0.0);
                    }
                    else {
                        pre.glob_disp[i3].set_val_dfd1(& pre.glob_vel[i4]);
                        i4 += 1usize;
                    }
                    i3 += 1usize;
                }
            }
            self.get_ruk_dfd1(rtmp, &mut  d_rtmp, &mut  d_rd_t,  get_matrix,  cmd.nonlinear_geom,  pre);
            for i1 in 0..nd_dof {
                rvec[i1].add(& rtmp[i1]);
            }
            if get_matrix {
                i3 = 0;
                i4 = 0;
                for _i1 in 0..nd_dof {
                    for _i2 in 0..nd_dof {
                        d_rd_v[i3]  +=  d_rtmp[i4];
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                    i3  +=  self.num_int_dof;
                }
            }
        }
        
        // Material damping
        d_non_zero = false;
        i1 = 0;
        while !d_non_zero && i1 < 36 {
            if pre.dmat[i1].val > 0.0 {
                d_non_zero = true;
            }
            i1 += 1usize;
        }
        
        if d_non_zero {
            if cmd.nonlinear_geom {
                self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            for i1 in 0..self.num_ip {
                i2 = i1 * 3;
                vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
                self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
                wt_det_j.set_val(self.ip_wt[i1]);
                wt_det_j.mult(& det_j);
                if self.dof_per_nd == 3 {
                    ar_to_vec_dfd1(&mut d_ndx, &mut  scr_v1,  0,  33);
                    mat_mul_dfd1(&mut scr_v2, &mut  pre.glob_disp, &mut  scr_v1,  3,  self.n_dim,  3);
                    vec_to_ar_dfd1(&mut ux, &mut  scr_v2,  0,  9);
                    for i2 in 0..nd_dof {
                        self.get_solid_strain_dfd1(&mut strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  i2,  MAX_INT,  cmd.nonlinear_geom);
                        i4 = i2;
                        for i3 in 0..6 {
                            //i4 = i3 * nd_dof + i2;
                            pre.bmat[i4].set_val_dfd1(& strain[i3]);
                            i4  +=  nd_dof;
                        }
                    }
                }
                else {
                    for i2 in 0..nd_dof {
                        self.get_section_def_dfd1(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  cmd.nonlinear_geom,  i2,  MAX_INT);
                        i4 = i2;
                        for i3 in 0..6 {
                            pre.bmat[i4].set_val_dfd1(& sec_def[i3]);
                            i4  +=  nd_dof;
                        }
                    }
                }
                i5 = 0;
                for i2 in 0..6 {
                    bvel[i2].set_val(0.0);
                    i4 = 0;
                    for _i3 in 0..nd_dof {
                        nd = self.dof_table[i4];
                        dof = self.dof_table[i4 + 1];
                        tmp.set_val_dfd1(& pre.bmat[i5]);
                        tmp.mult(& pre.glob_vel[dof * self.num_nds + nd]);
                        bvel[i2].add(& tmp);
                        i4  +=  2;
                        i5 += 1usize;
                    }
                }
                ar_to_vec_dfd1(&mut bvel, &mut  scr_v1,  0,  9);
                mat_mul_dfd1(&mut scr_v2, &mut  pre.dmat, &mut  scr_v1,  self.def_dim,  self.def_dim,  1);
                mat_mul_dfd1(&mut scr_v1, &mut  scr_v2, &mut  pre.bmat,  1,  self.def_dim,  nd_dof);
                for i2 in 0..nd_dof {
                    scr_v1[i2].mult(& wt_det_j);
                    rvec[i2].add(& scr_v1[i2]);
                }
                if get_matrix {
                    mat_mul_dfd1(&mut pre.cbmat, &mut  pre.dmat, &mut  pre.bmat,  self.def_dim,  self.def_dim,  nd_dof);
                    i3 = nd_dof * self.def_dim;
                    for i2 in 0..i3 {
                        pre.cbmat[i2].mult(& wt_det_j);
                    }
                    i5 = 0;
                    for i2 in 0..nd_dof {
                        for i3 in 0..nd_dof {
                            i6 = i2;
                            i7 = i3;
                            for _i4 in 0..6 {
                                //i5 = i2 * nd_dof + i3;
                                //i6 = i4 * nd_dof + i2;
                                //i7 = i4 * nd_dof + i3;
                                d_rd_v[i5]  +=  pre.bmat[i6].val * pre.cbmat[i7].val;
                                i6  +=  nd_dof;
                                i7  +=  nd_dof;
                            }
                            i5 += 1usize;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_ru_dfd1(&mut self, glob_r : &mut Vec<DiffDoub1>, globd_rdu : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub1StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut dof : usize;
        let mut dof2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let nd_dof : usize;
        let tot_dof : usize;
        let mut temp_acc = [DiffDoub1::new(); 30];
        //let mut rvec = &mut scr.scr_v1;
        let mut rvec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //let mut rtmp = &mut scr.scr_v2;
        let mut rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rdu[1089];
        //let mut d_rdu = &mut scr.scr_m1;
        let mut d_rdu = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rd_t[330];
        //let mut d_rd_t = &mut scr.scr_m2;
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //double d_rtmp[1089];
        //let mut d_rtmp = &mut scr.scr_m3;
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut c1 : f64;
        let c2 : f64;
        
        nd_dof = self.num_nds*self.dof_per_nd;
        tot_dof = nd_dof + self.num_int_dof;
        
        for i1 in 0..tot_dof {
            rvec[i1].set_val(0.0);
            if get_matrix && i1 < nd_dof {
                i3 = i1 * tot_dof;
                for _i2 in 0..tot_dof {
                    d_rdu[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 21 {
            self.get_ru_frc_fld_dfd1(glob_r, globd_rdu,  get_matrix, cmd, pre, nd_ar);
            return;
        }
        
        if self.this_type != 1 {
            self.get_ruk_dfd1(&mut rvec, &mut  d_rdu, &mut  d_rd_t,  get_matrix,  cmd.nonlinear_geom, pre);
        }
        
        if self.num_int_dof > 0 {
            i2 = nd_dof;
            for i1 in 0..self.num_int_dof {
                self.internal_ru[i1].set_val_dfd1(& rvec[i2]);
                i2 += 1usize;
            }
            if get_matrix {
                i4 = 0;
                for i1 in 0..tot_dof {
                    i3 = i1 * tot_dof + nd_dof;
                    for _i2 in nd_dof..tot_dof {
                        self.internal_mat[i4] = d_rdu[i3];
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
                //condense matrix, d_rd_t and d_rtmp used for scratch
                self.condense_mat(&mut d_rdu, &mut d_rd_t, &mut d_rtmp);
            }
        }
        
        if cmd.dynamic {
            if get_matrix {
                if cmd.lump_mass {
                    i2 = 0;
                    for i1 in 0..nd_dof {
                        nd = self.dof_table[i2];
                        dof = self.dof_table[i2 + 1];
                        i3 = dof * self.num_nds + nd;
                        temp_acc[i1].set_val_dfd1(& pre.glob_acc[i3]);
                        pre.glob_acc[i3].set_val(1.0);
                        i2  +=  2;
                    }
                    self.get_rum_dfd1(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    c1 = cmd.time_step;
                    c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
                    i2 = 0;
                    i3 = tot_dof + 1;
                    for i1 in 0..nd_dof {
                        temp_acc[i1].mult(& rtmp[i1]);
                        rvec[i1].add(& temp_acc[i1]);
                        d_rdu[i2]  +=  c1 * rtmp[i1].val;
                        i2  +=  i3;
                    }
                }
                else {
                    self.get_rum_dfd1(&mut rtmp, &mut d_rtmp,  true,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                    c1 = cmd.time_step;
                    c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
                    i3 = 0;
                    i4 = 0;
                    for _i1 in 0..nd_dof {
                        for _i2 in 0..nd_dof {
                            d_rdu[i3]  +=  c1 * d_rtmp[i4];
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                        i3  +=  self.num_int_dof;
                    }
                }
                if self.this_type != 1 {
                    self.get_rud_dfd1(&mut rtmp, &mut  d_rtmp,  true, cmd, pre, scr, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                    c2 = cmd.time_step * cmd.newmark_gamma * c1;
                    i3 = 0;
                    i4 = 0;
                    for _i1 in 0..nd_dof {
                        for _i2 in 0..nd_dof {
                            d_rdu[i3]  +=  c2 * d_rtmp[i4];
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                        i3  +=  self.num_int_dof;
                    }
                }
            }
            else {
                if cmd.lump_mass {
                    i2 = 0;
                    for i1 in 0..nd_dof {
                        nd = self.dof_table[i2];
                        dof = self.dof_table[i2 + 1];
                        i3 = dof * self.num_nds + nd;
                        temp_acc[i1].set_val_dfd1(& pre.glob_acc[i3]);
                        pre.glob_acc[i3].set_val(1.0);
                        i2  +=  2;
                    }
                    self.get_rum_dfd1(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        temp_acc[i1].mult(& rtmp[i1]);
                        rvec[i1].add(& temp_acc[i1]);
                    }
                }
                else {
                    self.get_rum_dfd1(&mut rtmp, &mut  d_rtmp,  false,  true,  cmd.nonlinear_geom, pre, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                }
                if self.this_type != 1 {
                    self.get_rud_dfd1(&mut rtmp, &mut  d_rtmp,  false,  cmd,  pre, scr, scr_dfd);
                    for i1 in 0..nd_dof {
                        rvec[i1].add(& rtmp[i1]);
                    }
                }
            }
        }
        
        i4 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i4]];
            dof = self.dof_table[i4+1];
            glob_ind = nd_ar[nd].dof_index[dof];
            glob_r[glob_ind].add(& rvec[i1]);
            if get_matrix {
                i3 = i1*tot_dof;
                i5 = 0;
                for _i2 in 0..nd_dof {
                    nd2 = self.nodes[self.dof_table[i5]];
                    dof2 = self.dof_table[i5+1];
                    glob_ind2 = nd_ar[nd2].dof_index[dof2];
                    globd_rdu.add_entry(glob_ind,   glob_ind2,   d_rdu[i3]);
                    i3 += 1usize;
                    i5  +=  2;
                }
            }
            i4 +=  2;
        }
        
        return;
    }

    pub fn get_rtk_dfd1(&mut self, rvec : &mut Vec<DiffDoub1>, d_rd_t : &mut Vec<f64>, get_matrix : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut d_nt = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut grad_t = [DiffDoub1::new(); 3];
        let mut q_vec = [DiffDoub1::new(); 3];
        let mut rtmp = [DiffDoub1::new(); 10];
        let mut d_rtmp = [DiffDoub1::new(); 100];
        let mut tmp_s : [f64; 3] = [0f64; 3];
        let mut tmp_t = [DiffDoub1::new(); 10];
        let mut tmp_tc = [DiffDoub1::new(); 9];
        let mut tmp_ar = [DiffDoub1::new(); 30];
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..self.num_nds {
                    d_rd_t[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 1 || self.this_type == 21 {
            return;
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            tmp.set_val(self.ip_wt[i1]);
            det_j.mult(& tmp);
            vec_to_ar_dfd1(&mut tmp_t, &mut  pre.glob_temp,  0,  10);
            mat_mul_ar_dfd1(&mut grad_t, &mut  tmp_t, &mut  d_ndx,  1,  self.num_nds,  3);
            vec_to_ar_dfd1(&mut tmp_tc, &mut  pre.tcmat,  0,  9);
            mat_mul_ar_dfd1(&mut q_vec, &mut  tmp_tc, &mut  grad_t,  3,  3,  1);
            mat_mul_ar_dfd1(&mut rtmp, &mut  d_ndx, &mut  q_vec,  self.num_nds,  3,  1);
            for i2 in 0..self.num_nds {
                rtmp[i2].mult(& det_j);
                rvec[i2].add(& rtmp[i2]);
            }
            if get_matrix {
                transpose_ar_dfd1(&mut d_nt, &mut  d_ndx,  self.num_nds, 3);
                mat_mul_ar_dfd1(&mut tmp_ar, &mut  tmp_tc, &mut  d_nt,  3,  3,  self.num_nds);
                mat_mul_ar_dfd1(&mut d_rtmp, &mut  d_ndx, &mut  tmp_ar,  self.num_nds,  3,  self.num_nds);
                i3 = self.num_nds * self.num_nds;
                for i2 in 0..i3 {
                    d_rd_t[i2]  +=  d_rtmp[i2].val*det_j.val;
                }
            }
        }
        
        return;
    }

    pub fn get_rtm_dfd1(& self, rvec : &mut Vec<DiffDoub1>, d_rd_tdot : &mut Vec<f64>, get_matrix : bool, actual_props : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut pt_tdot = DiffDoub1::new();
        let mut d_rtmp = [DiffDoub1::new(); 100];
        let mut save_cp = DiffDoub1::new();
        let mut tmp_s : [f64; 3] = [0f64; 3];
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            rvec[i1].set_val(0.0);
            if get_matrix {
                for _i2 in 0..self.num_nds {
                    d_rd_tdot[i3] = 0.0;
                    i3 += 1usize;
                }
            }
        }
        
        if self.this_type == 21 {
            return;
        }
        else if self.this_type == 1 {
            tmp.set_val_dfd1(& pre.mass_per_el);
            tmp.mult(& pre.spec_heat);
            d_rd_tdot[0] = tmp.val;
            tmp.mult(& pre.glob_temp[0]);
            rvec[0].set_val_dfd1(& tmp);
            return;
        }
        
        if !actual_props {
            save_cp.set_val_dfd1(& pre.spec_heat);
            pre.spec_heat.set_val(1.0);
        }
        else if self.dof_per_nd == 3 {
            pre.spec_heat.mult(& pre.mmat[0]);
        }
        
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut tmp_s, & self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  tmp_s);
            tmp.set_val(self.ip_wt[i1]);
            det_j.mult(& tmp);
            pt_tdot.set_val(0.0);
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd1(& n_vec[i2]);
                tmp.mult(& pre.glob_tdot[i2]);
                pt_tdot.add(& tmp);
            }
            for i2 in 0..self.num_nds {
                tmp.set_val_dfd1(& n_vec[i2]);
                tmp.mult(& pre.spec_heat);
                tmp.mult(& pt_tdot);
                tmp.mult(& det_j);
                rvec[i2].add(& tmp);
            }
            if get_matrix {
                mat_mul_ar_dfd1(&mut d_rtmp, & n_vec, & n_vec, self.num_nds, 1, self.num_nds);
                i3 = self.num_nds * self.num_nds;
                for i2 in 0..i3 {
                    d_rd_tdot[i2]  +=  d_rtmp[i2].val * pre.spec_heat.val * det_j.val;
                }
            }
        }
        
        if !actual_props {
            pre.spec_heat.set_val_dfd1(& save_cp);
        }
        else if self.dof_per_nd == 3 {
            pre.spec_heat.dvd(& pre.mmat[0]);
        }
        
        return;
    }

    pub fn get_rt_dfd1(&mut self, glob_r : &mut Vec<DiffDoub1>, globd_rd_t : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub1StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut glob_ind1 : usize;
        let mut glob_ind2 : usize;
        let c1 : f64 =  1.0/(cmd.time_step*cmd.newmark_gamma);
        //let mut rvec = &mut scr.scr_v1;
        let mut rvec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut rtmp = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut d_rtmp = match scr.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        //let mut d_rd_t = &mut scr.scr_m1;
        //let mut rtmp = &mut scr.scr_v2;
        //let mut d_rtmp = &mut scr.scr_m2;
        
        if self.this_type == 21 {
            //self.get_rt_frc_fld_dfd1(glob_r, globd_rd_t, &mut  scr.scr_m1, &mut  scr.scr_m2,  get_matrix, cmd, &mut  pre, &mut  nd_ar);
            self.get_rt_frc_fld_dfd1(glob_r, globd_rd_t, &mut  d_rd_t, &mut  d_rtmp,  get_matrix, pre, nd_ar);
            return;
        }
        
        if self.this_type != 1 {
            self.get_rtk_dfd1(&mut rvec, &mut  d_rd_t,  get_matrix,  pre);
        }
        
        if cmd.dynamic {
            self.get_rtm_dfd1(&mut rtmp, &mut  d_rtmp,  get_matrix,  true,  pre);
            i3 = 0;
            for i1 in 0..self.num_nds {
                rvec[i1].add(& rtmp[i1]);
                if get_matrix {
                    for _i2 in 0..self.num_nds {
                        d_rd_t[i3]  +=  c1 * d_rtmp[i3];
                        i3 += 1usize;
                    }
                }
            }
        }
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            glob_ind1 = nd_ar[self.nodes[i1]].sorted_rank;
            glob_r[glob_ind1].add(& rvec[i1]);
            if get_matrix {
                for i2 in 0..self.num_nds {
                    glob_ind2 = nd_ar[self.nodes[i2]].sorted_rank;
                    globd_rd_t.add_entry(glob_ind1,   glob_ind2,   d_rd_t[i3]);
                    i3 += 1usize;
                }
            }
        }
        
        return;
    }

    pub fn get_ru_frc_fld_dfd1(&mut self, glob_r : &mut Vec<DiffDoub1>, globd_rdu : &mut SparseMat, get_matrix : bool, cmd : & JobCommand, pre : &mut DiffDoub1StressPrereq, nd_ar : &mut Vec<Node>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let nd_dof : usize =  6;
        let tot_dof : usize =  6;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut dof : usize;
        let mut dof2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let mut rvec = [DiffDoub1::new(); 6];
        let mut d_rd_u : [f64; 36] = [0f64; 36];
        let mut d_vec = [DiffDoub1::new(); 3];
        let mut dist = DiffDoub1::new();
        let mut dv_vec = [DiffDoub1::new(); 3];
        let mut d_dvecd_u = [DiffDoub1::new(); 18];
        let mut d_distd_u = [DiffDoub1::new(); 6];
        let mut f_n1 = [DiffDoub1::new(); 3];
        let mut df_n1d_u = [DiffDoub1::new(); 18];
        let mut dto_p = DiffDoub1::new();
        let mut dto_p1 = DiffDoub1::new();
        let mut dto_p2 = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        
        // potential force
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd1(& pre.glob_nds[i2]);
            tmp.add(& pre.glob_disp[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_nds[i2]);
            tmp.sub(& pre.glob_disp[i2]);
            d_vec[i1].set_val_dfd1(& tmp);
        }
        
        dist.set_val_dfd1(& d_vec[0]);
        dist.sqr();
        tmp.set_val_dfd1(& d_vec[1]);
        tmp.sqr();
        dist.add(& tmp);
        tmp.set_val_dfd1(& d_vec[2]);
        tmp.sqr();
        dist.add(& tmp);
        dist.sqt();
        
        d_dvecd_u[0].set_val(-1.0);
        d_dvecd_u[3].set_val(1.0);
        d_dvecd_u[7].set_val(-1.0);
        d_dvecd_u[10].set_val(1.0);
        d_dvecd_u[14].set_val(-1.0);
        d_dvecd_u[17].set_val(1.0);
        
        mat_mul_ar_dfd1(&mut d_distd_u, &mut  d_vec, &mut  d_dvecd_u,  1,  3,  6);
        tmp.set_val(1.0);
        tmp.dvd(& dist);
        for i1 in 0..6 {
            d_distd_u[i1].mult(& tmp);
        }
        
        dto_p.set_val_dfd1(& dist);
        i1 = 1;
        while i1 < (pre.frc_fld_exp[0].val as usize) {
            dto_p.mult(& dist);
            i1 += 1usize;
        }
        dto_p1.set_val_dfd1(& dto_p);
        dto_p1.mult(& dist);
        dto_p2.set_val_dfd1(& dto_p1);
        dto_p2.mult(& dist);
        
        tmp.set_val_dfd1(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd1(& d_vec[i1]);
            f_n1[i1].mult(& tmp);
        }
        
        mat_mul_ar_dfd1(&mut df_n1d_u, &mut  d_vec, &mut  d_distd_u,  3,  1,  6);
        tmp.set_val(1.0);
        tmp.add(& pre.frc_fld_exp[0]);
        tmp.neg();
        tmp.mult(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p2);
        for i1 in 0..18 {
            df_n1d_u[i1].mult(& tmp);
        }
        
        tmp.set_val_dfd1(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..18 {
            tmp2.set_val_dfd1(& tmp);
            tmp2.mult(& d_dvecd_u[i1]);
            df_n1d_u[i1].add(& tmp2);
        }
        
        for i1 in 0..3 {
            rvec[i1].set_val_dfd1(& f_n1[i1]);
            rvec[i1].neg();
            rvec[i1 + 3].set_val_dfd1(& f_n1[i1]);
        }
        
        for i1 in 0..18 {
            d_rd_u[i1] = -df_n1d_u[i1].val;
            d_rd_u[i1 + 18] = df_n1d_u[i1].val;
        }
        
        // damping force
        
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd1(& pre.glob_vel[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_vel[i2]);
            dv_vec[i1].set_val_dfd1(& tmp);
        }
        
        tmp.set_val_dfd1(& pre.frc_fld_coef[1]);
        tmp.dvd(& dto_p);
        
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd1(& tmp);
            f_n1[i1].mult(& dv_vec[i1]);
        }
        
        
        tmp2.set_val(-cmd.newmark_gamma / (cmd.time_step * (cmd.newmark_beta - cmd.newmark_gamma)));
        tmp.mult(& tmp2);
        
        for i1 in 0..18 {
            tmp2.set_val_dfd1(& tmp);
            tmp2.mult(& d_dvecd_u[i1]);
            //df_n1d_u[i1].add(tmp2);
            df_n1d_u[i1].set_val_dfd1(& tmp2);
        }
        
        for i1 in 0..3 {
            rvec[i1].sub(& f_n1[i1]);
            rvec[i1 + 3].add(& f_n1[i1]);
        }
        
        for i1 in 0..18 {
            d_rd_u[i1]  -=  df_n1d_u[i1].val;
            d_rd_u[i1 + 18]  +=  df_n1d_u[i1].val;
        }
        
        i4 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i4]];
            dof = self.dof_table[i4 + 1];
            glob_ind = nd_ar[nd].dof_index[dof];
            glob_r[glob_ind].add(& rvec[i1]);
            if get_matrix {
                i3 = i1 * tot_dof;
                i5 = 0;
                for _i2 in 0..nd_dof {
                    nd2 = self.nodes[self.dof_table[i5]];
                    dof2 = self.dof_table[i5 + 1];
                    glob_ind2 = nd_ar[nd2].dof_index[dof2];
                    globd_rdu.add_entry(glob_ind,   glob_ind2,   d_rd_u[i3]);
                    i3 += 1usize;
                    i5  +=  2;
                }
            }
            i4  +=  2;
        }
        
        return;
    }

    pub fn get_rt_frc_fld_dfd1(&mut self, glob_r : &mut Vec<DiffDoub1>, globd_rd_t : &mut SparseMat, d_rd_u : &mut Vec<f64>, d_rd_v : &mut Vec<f64>, get_matrix : bool, pre : &mut DiffDoub1StressPrereq, nd_ar : &mut Vec<Node>) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut nd : usize;
        let mut nd2 : usize;
        let mut glob_ind : usize;
        let mut glob_ind2 : usize;
        let mut rvec = [DiffDoub1::new(); 2];
        let mut d_rd_t : [f64; 4] = [0f64; 4];
        let mut d_vec = [DiffDoub1::new(); 3];
        let mut dist = DiffDoub1::new();
        let mut dv_vec = [DiffDoub1::new(); 3];
        let mut d_dvecd_u = [DiffDoub1::new(); 18];
        let mut d_distd_u = [DiffDoub1::new(); 6];
        let mut f_n1 = [DiffDoub1::new(); 3];
        let mut df_n1d_u = [DiffDoub1::new(); 18];
        let mut df_n1d_v = [DiffDoub1::new(); 18];
        let mut dto_p = DiffDoub1::new();
        let mut dto_p1 = DiffDoub1::new();
        let mut dto_p2 = DiffDoub1::new();
        let mut t1to3 = DiffDoub1::new();
        let mut t1to4 = DiffDoub1::new();
        let mut t2to3 = DiffDoub1::new();
        let mut t2to4 = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut tmp_mat = [DiffDoub1::new(); 18];
        let mut tmp_mat2 = [DiffDoub1::new(); 18];
        
        // potential force
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd1(& pre.glob_nds[i2]);
            tmp.add(& pre.glob_disp[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_nds[i2]);
            tmp.sub(& pre.glob_disp[i2]);
            d_vec[i1].set_val_dfd1(& tmp);
        }
        
        dist.set_val_dfd1(& d_vec[0]);
        dist.sqr();
        tmp.set_val_dfd1(& d_vec[1]);
        tmp.sqr();
        dist.add(& tmp);
        tmp.set_val_dfd1(& d_vec[2]);
        tmp.sqr();
        dist.add(& tmp);
        dist.sqt();
        
        dto_p.set_val_dfd1(& dist);
        i1 = 1;
        while i1 < (pre.frc_fld_exp[0].val as usize) {
            dto_p.mult(& dist);
            i1 += 1usize;
        }
        dto_p1.set_val_dfd1(& dto_p);
        dto_p1.mult(& dist);// ||d||^(p+1)
        dto_p2.set_val_dfd1(& dto_p1);
        dto_p2.mult(& dist);// ||d||^(p+2)
        
        tmp.set_val_dfd1(& pre.frc_fld_coef[0]);
        tmp.dvd(& dto_p1);
        for i1 in 0..3 {
            f_n1[i1].set_val_dfd1(& d_vec[i1]);
            f_n1[i1].mult(& tmp);
        }
        
        // damping force
        
        for i1 in 0..3 {
            i2 = i1 * 2 + 1;
            tmp.set_val_dfd1(& pre.glob_vel[i2]);
            i2  -=  1;
            tmp.sub(& pre.glob_vel[i2]);
            dv_vec[i1].set_val_dfd1(& tmp);
        }
        
        tmp.set_val_dfd1(& pre.frc_fld_coef[1]);
        tmp.dvd(& dto_p);// tmp = c_d/||d||^p
        
        for i1 in 0..3 {
            tmp2.set_val_dfd1(& tmp);
            tmp2.mult(& dv_vec[i1]);
            f_n1[i1].add(& tmp2);
        }
        
        // calculate absolute temps
        tmp.set_val_dfd1(& pre.glob_temp[0]);
        tmp.add(& pre.ref_temp);
        t1to4.set_val_dfd1(& tmp);
        t1to4.sqr();
        t1to4.sqr();
        t1to3.set_val_dfd1(& t1to4);
        t1to3.dvd(& tmp);
        
        tmp.set_val_dfd1(& pre.glob_temp[1]);
        tmp.add(& pre.ref_temp);
        t2to4.set_val_dfd1(& tmp);
        t2to4.sqr();
        t2to4.sqr();
        t2to3.set_val_dfd1(& t2to4);
        t2to3.dvd(& tmp);
        
        // conduction and radiation terms
        tmp.set_val_dfd1(& pre.glob_temp[0]);
        tmp.sub(& pre.glob_temp[1]);
        tmp.mult(& pre.thrm_fld_coef[0]);
        tmp.dvd(& dist);
        rvec[0].set_val_dfd1(& tmp);
        
        tmp.set_val_dfd1(& t1to4);
        tmp.sub(& t2to4);
        tmp.mult(& pre.thrm_fld_coef[1]);
        tmp.dvd(& dist);
        tmp.dvd(& dist);
        rvec[0].add(& tmp);
        
        tmp.set_val_dfd1(& rvec[0]);
        rvec[1].set_val_dfd1(& tmp);
        rvec[1].neg();
        
        // work dissipation term
        
        tmp.set_val_dfd1(& f_n1[0]);
        tmp.mult(& dv_vec[0]);
        tmp2.set_val_dfd1(& f_n1[1]);
        tmp2.mult(& dv_vec[1]);
        tmp.add(& tmp2);
        tmp2.set_val_dfd1(& f_n1[2]);
        tmp2.mult(& dv_vec[2]);
        tmp.add(& tmp2);
        tmp2.set_val(0.5);
        tmp.mult(& tmp2);// tmp = 0.5*dot(fn1,d_v)
        
        rvec[0].sub(& tmp);
        rvec[1].sub(& tmp);
        for i1 in 0..2 {
            nd = self.nodes[i1];
            glob_ind = nd_ar[nd].sorted_rank;
            glob_r[glob_ind].add(& rvec[i1]);
        }
        
        if get_matrix {
            d_dvecd_u[0].set_val(-1.0);
            d_dvecd_u[3].set_val(1.0);
            d_dvecd_u[7].set_val(-1.0);
            d_dvecd_u[10].set_val(1.0);
            d_dvecd_u[14].set_val(-1.0);
            d_dvecd_u[17].set_val(1.0);
            
            mat_mul_ar_dfd1(&mut d_distd_u, &mut  d_vec, &mut  d_dvecd_u,  1,  3,  6);
            tmp.set_val(1.0);
            tmp.dvd(& dist);
            for i1 in 0..6 {
                d_distd_u[i1].mult(& tmp);
            }
            
            mat_mul_ar_dfd1(&mut df_n1d_u, &mut  d_vec, &mut  d_distd_u,  3,  1,  6);
            tmp.set_val(1.0);
            tmp.add(& pre.frc_fld_exp[0]);
            tmp.neg();
            tmp.mult(& pre.frc_fld_coef[0]);
            tmp.dvd(& dto_p2);
            for i1 in 0..18 {
                df_n1d_u[i1].mult(& tmp);
            }
            
            tmp.set_val_dfd1(& pre.frc_fld_coef[0]);
            tmp.dvd(& dto_p1);
            for i1 in 0..18 {
                tmp2.set_val_dfd1(& tmp);
                tmp2.mult(& d_dvecd_u[i1]);
                df_n1d_u[i1].add(& tmp2);
            }
            
            mat_mul_ar_dfd1(&mut tmp_mat, &mut  dv_vec, &mut  d_distd_u,  3,  1,  6);
            tmp.set_val_dfd1(& pre.frc_fld_coef[1]);
            tmp.mult(& pre.frc_fld_exp[1]);
            tmp.dvd(& dto_p1);
            tmp.neg();
            for i1 in 0..18 {
                tmp2.set_val_dfd1(& tmp);
                tmp2.mult(& tmp_mat[i1]);
                df_n1d_u[i1].add(& tmp2);
            }
            
            tmp.set_val_dfd1(& pre.frc_fld_coef[1]);
            tmp.dvd(& dto_p);
            for i1 in 0..18 {
                df_n1d_v[i1].set_val_dfd1(& d_dvecd_u[i1]);
                df_n1d_v[i1].mult(& tmp);
            }
            
            // d_rd_t
            tmp.set_val_dfd1(& pre.thrm_fld_coef[0]);
            tmp.dvd(& dist);
            tmp2.set_val_dfd1(& pre.thrm_fld_coef[1]);
            tmp2.dvd(& dist);
            tmp2.dvd(& dist);
            d_rd_t[0] = tmp.val + tmp2.val * t1to3.val;
            d_rd_t[1] = -tmp.val;
            d_rd_t[2] = -tmp.val;
            d_rd_t[3] = tmp.val + tmp2.val * t2to3.val;
            
            i3 = 0;
            for i1 in 0..2 {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                for i2 in 0..2 {
                    nd2 = self.nodes[i2];
                    glob_ind2 = nd_ar[nd2].sorted_rank;
                    globd_rd_t.add_entry(glob_ind,   glob_ind2,   d_rd_t[i3]);
                    i3 += 1usize;
                }
            }
            
            // d_rd_u
            for i1 in 0..12 {
                d_rd_u[i1] = 0.0;
            }
            
            // conduction term
            tmp.set_val_dfd1(& pre.glob_temp[0]);
            tmp.sub(& pre.glob_temp[1]);
            tmp.mult(& pre.thrm_fld_coef[0]);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.neg();
            tmp_mat[0].set_val_dfd1(& tmp);
            tmp_mat[1].set_val_dfd1(& tmp);
            tmp_mat[1].neg();
            
            mat_mul_ar_dfd1(&mut tmp_mat2, &mut  tmp_mat, &mut  d_distd_u,  2,  1,  6);
            
            for i1 in 0..12 {
                d_rd_u[i1]  +=  tmp_mat2[i1].val;
            }
            
            // radiation term
            tmp.set_val_dfd1(& t1to4);
            tmp.sub(& t2to4);
            tmp.mult(& pre.thrm_fld_coef[1]);
            tmp2.set_val(2.0);
            tmp.mult(& tmp2);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.dvd(& dist);
            tmp.neg();
            tmp_mat[0].set_val_dfd1(& tmp);
            tmp_mat[1].set_val_dfd1(& tmp);
            tmp_mat[1].neg();
            
            mat_mul_ar_dfd1(&mut tmp_mat2, &mut  tmp_mat, &mut  d_distd_u,  2,  1,  6);
            
            for i1 in 0..12 {
                d_rd_u[i1]  +=  tmp_mat2[i1].val;
            }
            
            // work dissipation term
            mat_mul_ar_dfd1(&mut tmp_mat, &mut  dv_vec, &mut  df_n1d_u,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_u[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
            // d_rd_v
            
            for i1 in 0..12 {
                d_rd_v[i1] = 0.0;
            }
            
            mat_mul_ar_dfd1(&mut tmp_mat, &mut  dv_vec, &mut  df_n1d_v,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_v[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
            mat_mul_ar_dfd1(&mut tmp_mat, &mut  f_n1, &mut  d_dvecd_u,  1,  3,  6);
            i3 = 0;
            for _i1 in 0..2 {
                for i2 in 0..6 {
                    d_rd_v[i3]  -=  0.5 * tmp_mat[i2].val;
                    i3 += 1usize;
                }
            }
            
        }
        
        return;
    }

    pub fn get_app_load_dfd1(& self, app_ld : &mut Vec<DiffDoub1>, ld_pt : & Load, n_lgeom : bool, pre : &mut DiffDoub1StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>, sec_ar : &mut Vec<Section>, fc_ar : &mut Vec<Face>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut glob_ind : usize;
        let mut nd : usize;
        let mut dof : usize;
        let nd_dof : usize =  self.num_nds*self.dof_per_nd;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let ld_type : CppStr = ld_pt.this_type.clone();
        //let mut d_rd_a = &mut scr.scr_m1;
        let mut d_rd_a = match scr.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        let mut fc_num_nds : usize;
        let mut nd_in_face : bool;
        //let mut el_app_ld = &mut scr.scr_v1;
        let mut el_app_ld = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        let mut tot_nd_f = [DiffDoub1::new(); 6];
        let mut inp_mag = DiffDoub1::new();
        let mut vec_mag = DiffDoub1::new();
        let mut el_vol = DiffDoub1::new();
        let mut sc_fact = DiffDoub1::new();
        let mut el_cent = [DiffDoub1::new(); 3];
        let mut cent_to_el = [DiffDoub1::new(); 3];
        let mut ax_to_el = [DiffDoub1::new(); 3];
        let mut ang_vel2 = DiffDoub1::new();
        let mut nn_inv = DiffDoub1::new();
        let mut fc_area = DiffDoub1::new();
        let mut fc_norm = [DiffDoub1::new(); 3];
        let mut trac = [DiffDoub1::new(); 3];
        let mut dp = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        
        for i1 in 0..nd_dof {
            pre.glob_acc[i1].set_val(0.0);
        }
        
        if ld_type.s == "bodyForce" {
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                pre.glob_acc[i3].set_val(ld_pt.load[dof]);
                i2  +=  2;
            }
            self.get_rum_dfd1(&mut el_app_ld, &mut  d_rd_a,  false,  false,  n_lgeom, pre, scr_dfd);
            i2 = 0;
            for i1 in 0..nd_dof {
                dof = self.dof_table[i2 + 1];
                tot_nd_f[dof].add(& el_app_ld[i1]);
                i2  +=  2;
            }
            for i1 in 0..self.dof_per_nd {
                tmp.set_val(ld_pt.load[i1]);
                tmp.sqr();
                inp_mag.add(& tmp);
                tmp.set_val_dfd1(& tot_nd_f[i1]);
                tmp.sqr();
                vec_mag.add(& tmp);
            }
            inp_mag.sqt();
            vec_mag.sqt();
            if self.this_type == 41 || self.this_type == 3 {
                el_vol.set_val(0.0);
                for i1 in 0..num_lay {
                    self.get_volume_dfd1(&mut tmp, pre,  i1, sec_ar, dv_ar);
                    el_vol.add(& tmp);
                }
            }
            else {
                self.get_volume_dfd1(&mut el_vol, pre,  0, sec_ar, dv_ar);
            }
            sc_fact.set_val_dfd1(& inp_mag);
            sc_fact.mult(& el_vol);
            sc_fact.dvd(& vec_mag);
            for i1 in 0..nd_dof {
                el_app_ld[i1].mult(& sc_fact);
            }
        }
        else if ld_type.s == "gravitational" {
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                if dof < 3 {
                    pre.glob_acc[i3].set_val(ld_pt.load[dof]);
                }
                i2  +=  2;
            }
            self.get_rum_dfd1(&mut el_app_ld, &mut  d_rd_a,  false,  true,  n_lgeom,  pre, scr_dfd);
        }
        else if ld_type.s == "centrifugal" {
            ang_vel2.set_val(ld_pt.angular_vel);
            ang_vel2.sqr();
            i3 = 0;
            nn_inv.set_val(1.0 / (self.num_nds as f64));
            dp.set_val(0.0);
            for i1 in 0..3 {
                for _i2 in 0..self.num_nds {
                    el_cent[i1].add(& pre.glob_nds[i3]);
                    i3 += 1usize;
                }
                el_cent[i1].mult(& nn_inv);
                cent_to_el[i1].set_val_dfd1(& el_cent[i1]);
                tmp.set_val(ld_pt.center[i1]);
                cent_to_el[i1].sub(& tmp);
                tmp.set_val(ld_pt.axis[i1]);
                tmp.mult(& cent_to_el[i1]);
                dp.add(& tmp);
            }
            for i1 in 0..3 {
                ax_to_el[i1].set_val_dfd1(& cent_to_el[i1]);
                tmp.set_val(ld_pt.axis[i1]);
                tmp.mult(& dp);
                ax_to_el[i1].sub(& tmp);
                ax_to_el[i1].mult(& ang_vel2);
            }
            i2 = 0;
            for _i1 in 0..nd_dof {
                nd = self.dof_table[i2];
                dof = self.dof_table[i2 + 1];
                i3 = dof * self.num_nds + nd;
                if dof < 3 {
                    pre.glob_acc[i3].set_val_dfd1(& ax_to_el[dof]);
                }
                i2  +=  2;
            }
            self.get_rum_dfd1(&mut el_app_ld, &mut  d_rd_a,  false,  true,  n_lgeom, pre, scr_dfd);
        }
        else if ld_type.s == "surfaceTraction" || ld_type.s == "surfacePressure" {
            let mut this_fc : &Face;
            for fi in self.faces.iter() {
                this_fc = &fc_ar[*fi];
                if this_fc.on_surf {
                    this_fc.get_area_normal_dfd1(&mut fc_area, &mut  fc_norm, &pre.glob_nds, self.num_nds);
                    dp.set_val(0.0);
                    for i1 in 0..3 {
                        tmp.set_val(ld_pt.normal_dir[i1]);
                        tmp.mult(& fc_norm[i1]);
                        dp.add(& tmp);
                    }
                    tmp.set_val(R_PI_0180* ld_pt.norm_tol);
                    tmp.cs();
                    if dp.val > tmp.val {
                        if ld_type.s == "surfacePressure" {
                            for i1 in 0..3 {
                                trac[i1].set_val(ld_pt.load[0]);
                                trac[i1].mult(& fc_norm[i1]);
                                trac[i1].neg();
                            }
                        }
                        else {
                            for i1 in 0..3 {
                                trac[i1].set_val(ld_pt.load[i1]);
                            }
                        }
                        fc_num_nds = this_fc.num_nds;
                        for i1 in 0..fc_num_nds {
                            i4 = this_fc.loc_nodes[i1];
                            for i3 in 0..3 {
                                //i4 = i3 * self.num_nds + fc_loc_nd[i1];
                                pre.glob_acc[i4].set_val_dfd1(& trac[i3]);
                                i4  +=  self.num_nds;
                            }
                        }
                        self.get_rum_dfd1(&mut el_app_ld, &mut  d_rd_a,  false,  false,  n_lgeom, pre, scr_dfd);
                        i2 = 0;
                        for i1 in 0..nd_dof {
                            nd = self.dof_table[i2];
                            dof = self.dof_table[i2 + 1];
                            nd_in_face = false;
                            for i3 in 0..fc_num_nds {
                                if this_fc.loc_nodes[i3] == nd {
                                    nd_in_face = true;
                                }
                            }
                            if !nd_in_face {
                                el_app_ld[i1].set_val(0.0);
                            }
                            tot_nd_f[dof].add(& el_app_ld[i1]);
                            i2  +=  2;
                        }
                        inp_mag.set_val(0.0);
                        vec_mag.set_val(0.0);
                        for i1 in 0..3 {
                            tmp.set_val_dfd1(& tot_nd_f[i1]);
                            tmp.sqr();
                            vec_mag.add(& tmp);
                            tmp.set_val_dfd1(& trac[i1]);
                            tmp.sqr();
                            inp_mag.add(& tmp);
                        }
                        inp_mag.sqt();
                        vec_mag.sqt();
                        tmp.set_val_dfd1(& fc_area);
                        tmp.mult(& inp_mag);
                        tmp.dvd(& vec_mag);
                        for i1 in 0..nd_dof {
                            el_app_ld[i1].mult(& tmp);
                        }
                    }
                }
            }
        }
        
        i2 = 0;
        for i1 in 0..nd_dof {
            nd = self.nodes[self.dof_table[i2]];
            dof = self.dof_table[i2 + 1];
            glob_ind = nd_ar[nd].dof_index[dof];
            app_ld[glob_ind].add(& el_app_ld[i1]);
            i2  +=  2;
        }
        
        return;
    }

    pub fn get_app_therm_load_dfd1(& self, app_ld : &mut Vec<DiffDoub1>, ld_pt : & Load, pre : &mut DiffDoub1StressPrereq, 
        scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>, sec_ar : &mut Vec<Section>, fc_ar : &mut Vec<Face>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut glob_ind : usize;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let ld_type : CppStr = ld_pt.this_type.clone();
        //let mut el_app_ld = &mut scr.scr_v1;
        let el_app_ld = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        //let mut d_rd_t = &mut scr.scr_m1;
        let d_rd_t = match scr.next() {
            None => panic!("Error: ran out of scratch matrices"),
            Some(x) => &mut x.dat,
        };
        let mut fc_num_nds : usize;
        let mut nd_in_face : bool;
        let mut fc_area = DiffDoub1::new();
        let mut fc_norm = [DiffDoub1::new(); 3];
        let mut tot_hg = DiffDoub1::new();
        let mut el_vol = DiffDoub1::new();
        let mut dp = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        
        for i1 in 0..self.num_nds {
            pre.glob_tdot[i1].set_val(0.0);
        }
        
        if ld_type.s == "bodyHeatGen" {
            for i1 in 0..self.num_nds {
                pre.glob_tdot[i1].set_val(ld_pt.load[0]);
            }
            self.get_rtm_dfd1(el_app_ld, d_rd_t, false, false, pre);
            tot_hg.set_val(0.0);
            for i1 in 0..self.num_nds {
                tot_hg.add(& el_app_ld[i1]);
            }
            if num_lay > 0 {
                el_vol.set_val(0.0);
                for i1 in 0..num_lay {
                    self.get_volume_dfd1(&mut tmp, pre,  i1, sec_ar, dv_ar);
                    el_vol.add(& tmp);
                }
            }
            else {
                self.get_volume_dfd1(&mut el_vol, pre,  0, sec_ar, dv_ar);
            }
            tmp.set_val(ld_pt.load[0]);
            tmp.mult(& el_vol);
            tmp.dvd(& tot_hg);
            for i1 in 0..self.num_nds {
                el_app_ld[i1].mult(& tmp);
            }
        }
        else if ld_type.s == "surfaceFlux" {
            let mut this_fc : &Face;
            for fi in self.faces.iter() {
                this_fc = &fc_ar[*fi];
                if this_fc.on_surf {
                    this_fc.get_area_normal_dfd1(&mut fc_area, &mut  fc_norm, &pre.glob_nds, self.num_nds);
                    for i1 in 0..3 {
                        tmp.set_val(ld_pt.normal_dir[i1]);
                        tmp.mult(& fc_norm[i1]);
                        dp.add(& tmp);
                    }
                    tmp.set_val(R_PI_0180 * ld_pt.norm_tol);
                    tmp.cs();
                    if dp.val > tmp.val {
                        fc_num_nds = this_fc.num_nds;
                        for i1 in 0..fc_num_nds {
                            i2 = this_fc.loc_nodes[i1];
                            pre.glob_tdot[i2].set_val(ld_pt.load[0]);
                        }
                        self.get_rtm_dfd1(el_app_ld, d_rd_t,  false,  false, pre);
                        tot_hg.set_val(0.0);
                        for i1 in 0..self.num_nds {
                            nd_in_face = false;
                            for i2 in 0..fc_num_nds {
                                if this_fc.loc_nodes[i2] == i1 {
                                    nd_in_face = true;
                                }
                            }
                            if !nd_in_face {
                                el_app_ld[i1].set_val(0.0);
                            }
                            tot_hg.add(& el_app_ld[i1]);
                        }
                        tmp.set_val(ld_pt.load[0]);
                        tmp.mult(& fc_area);
                        tmp.dvd(& tot_hg);
                        for i1 in 0..self.num_nds {
                            el_app_ld[i1].mult(& tmp);
                        }
                    }
                }
            }
        }
        
        for i1 in 0..self.num_nds {
            glob_ind = nd_ar[self.nodes[i1]].sorted_rank;
            app_ld[glob_ind].add(& el_app_ld[i1]);
        }
        
        return;
    }

    pub fn get_rf_unsteady_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq) {
        let mut rho = DiffDoub1::new();
        let mut rho_dot = DiffDoub1::new();
        let mut fl_vel = [DiffDoub1::new(); 3];
        let mut vel_dot = [DiffDoub1::new(); 3];
        let mut temp = DiffDoub1::new();
        let mut t_dot = DiffDoub1::new();
        let mut turb_e = DiffDoub1::new();
        let mut turb_e_dot = DiffDoub1::new();
        let mut enth = DiffDoub1::new();
        let mut enth_dot = DiffDoub1::new();

        let mut si = 0usize;
        let mut vi : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut r_sum = [DiffDoub1::new(); 6];
        let mut ipw_j = DiffDoub1::new();
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        let mut td3 = DiffDoub1::new();

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            rho_dot.set_val(0.0);
            temp.set_val(0.0);
            t_dot.set_val(0.0);
            turb_e.set_val(0.0);
            turb_e_dot.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den_dot[j]);
                rho_dot.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_temp[j]);
                temp.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_tdot[j]);
                t_dot.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_edot[j]);
                turb_e_dot.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                vel_dot[k].set_val(0.0);
                vi = self.num_nds*k;
                for j in 0..self.num_nds {
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel[vi]);
                    fl_vel[k].add(&td);
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel_dot[vi]);
                    vel_dot[k].add(&td);
                    vi += 1;
                }
            }
            enth.set_val_dfd1(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);
            enth_dot.set_val_dfd1(&pre.spec_heat);
            enth_dot.mult(&t_dot);

            //mass
            td.set_val_dfd1(&rho_dot);
            td.mult(&ipw_j);
            r_sum[0].add(&td);

            //momentum
            for j in 0..3 {
                td.set_val_dfd1(&rho_dot);
                td.mult(&fl_vel[j]);
                td2.set_val_dfd1(&rho);
                td2.mult(&vel_dot[j]);
                td.add(&td2);
                td.mult(&ipw_j);
                r_sum[j+1].add(&td);
            }
            
            //energy
            td.set_val(0.5);
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd1(&fl_vel[j]);
                td3.sqr();
                td2.add(&td3);
            }
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_dot); // td = rho_dot(0.5v*v + h + tau)
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd1(&fl_vel[j]);
                td3.mult(&vel_dot[j]);
                td2.add(&td3);
            }
            td2.add(&enth_dot);
            td2.add(&turb_e_dot);
            td2.mult(&rho);  //td2 = rho(v*v_dot + h_dot + tau_dot)
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[4].add(&td);

            //turbulence
            td.set_val_dfd1(&rho_dot);
            td.mult(&turb_e);
            td2.set_val_dfd1(&rho);
            td2.mult(&turb_e_dot);
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[5].add(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_convective_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_convective_dfd1(r_vec,pre);
        }
    }

    pub fn get_rf_viscous_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub1::new();
        let mut turb_e = DiffDoub1::new();
        let mut grad_v = [DiffDoub1::new(); 9];
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();

        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ipw_j = DiffDoub1::new();
        let mut r_sum = DiffDoub1::new();
        let mut si = 0usize;
        let mut gvi : usize;
        let mut nvi : usize;
        let mut gni : usize;

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
            }
            gvi = 0;
            for j in 0..3 {
                for k in 0..3 {
                    grad_v[gvi].set_val(0.0);
                    nvi = j*self.num_nds;
                    gni = k;
                    for _m  in 0..self.num_nds {
                        //nvi = j*self.num_nds + m;
                        //gni = m*3 + k;
                        td.set_val_dfd1(&pre.fl_vel[nvi]);
                        td.mult(&d_ndx[gni]);
                        grad_v[gvi].add(&td);
                        nvi += 1;
                        gni += 3;
                    }
                    gvi += 1;
                }
            }
            td.set_val(0.0);
            for j in 0..9 {
                td2.set_val_dfd1(&grad_v[j]);
                td2.sqr();
                td.add(&td2);
            }
            td.sqt(); // td = mag(grad_v)
            
            td2.set_val_dfd1(&pre.grad_turb_coef);
            td2.mult(&rho);
            td2.mult(&td);
            td2.mult(&ipw_j);
            r_sum.sub(&td2);

            td2.set_val_dfd1(&pre.diss_turb_coef);
            td2.mult(&rho);
            td2.mult(&turb_e);
            td2.mult(&ipw_j);
            r_sum.add(&td2);

            si += 3;
        }

        let mut ri = 5usize;
        for j in 0..self.num_nds {
            r_vec[ri].add(&r_sum);
            ri += 6;
        }

        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
        }

    }

    pub fn get_rf_pressure_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_pressure_dfd1(r_vec,pre, fc_vec, ec_vec);
        }
    }

    pub fn get_rf_load_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq) {
        let mut rho = DiffDoub1::new();
        let mut vel = [DiffDoub1::new(); 3];

        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ipw_j = DiffDoub1::new();
        let mut si = 0usize;
        let mut nvi : usize;
        let mut r_sum = [DiffDoub1::new(); 6];
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        
        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
            }
            nvi = 0;
            for k in 0..3 {
                vel[k].set_val(0.0);
                for j in 0..self.num_nds {
                    //nvi = k*self.num_nds + j;
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel[nvi]);
                    vel[k].add(&td);
                    nvi += 1;
                }
            }

            //momentum
            for j in 1..4 {
                td.set_val_dfd1(&pre.f_per_mass[j-1]);
                td.mult(&rho);
                td.mult(&ipw_j);
                r_sum[j].sub(&td);
            }

            //energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&vel[j]);
                td2.mult(&pre.f_per_mass[j]);
                td.add(&td2);
            }
            td.mult(&rho);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            td.set_val_dfd1(&pre.hg_per_vol);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_dfd1(&self, r_glob : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq, dynamic : bool, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {
        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_dfd1"),
            Some(x) => &mut x.dat,
        };
        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd1(r_vec, pre);
        }
        self.get_rf_convective_dfd1( r_vec, pre, fc_vec);
        self.get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_pressure_dfd1(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_load_dfd1(r_vec,pre);

        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                r_glob[rg+j].add(&r_vec[ri]);
                ri += 1;
            }
        }
    }

    //end dup
 
//end skip 
 
    pub fn get_rf_mat_col(&self, glob_mat : &mut SparseMat, pre : &mut DiffDoub1FlPrereq, dynamic : bool, glob_col : usize, prim_var : bool, mult_fact : f64,
        fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {

        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_mat()"),
            Some(x) => &mut x.dat,
        };

        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd1(r_vec, pre);
        }
        if prim_var {
            self.get_rf_convective_dfd1(r_vec, pre, fc_vec);
            self.get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
            self.get_rf_pressure_dfd1(r_vec, pre, fc_vec, ec_vec);
            self.get_rf_load_dfd1(r_vec,pre);
        }
        
        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                glob_mat.add_entry(rg+j, glob_col, mult_fact*r_vec[ri].dval);
                ri += 1;
            }
        }
    }
 
}


