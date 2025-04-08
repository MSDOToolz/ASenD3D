use crate::element::*;
use crate::face::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::list_ent::*;
use crate::design_var::*;
use crate::node::*;
use crate::section::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;


impl Element {

    //dup1

    pub fn get_nd_disp_dfd0(&self, glob_disp : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_disp[i3].set_val(this_nd.displacement[i2]);
                i3 +=  self.n_dim;
            }
        }
        
        if self.num_int_dof > 0 {
            i2 = 2*self.num_nds*self.dof_per_nd;
            for i3 in 0..self.num_int_dof {
                i4 = self.dof_table[i2];
                i5 = self.dof_table[i2+1];
                i6 = self.n_dim*i5 + i4;
                glob_disp[i6].set_val(self.internal_disp[i3]);
                i2 +=  2;
            }
        }
        
        return;
    }

    pub fn get_nd_vel_dfd0(& self, glob_vel : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_vel[i3].set_val(this_nd.velocity[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_acc_dfd0(& self, glob_acc : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_acc[i3].set_val(this_nd.acceleration[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_fl_vel_dfd0(& self, fl_vel : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..3 {
                fl_vel[i3].set_val(this_nd.fl_vel[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_fl_vdot_dfd0(& self, fl_vdot : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..3 {
                fl_vdot[i3].set_val(this_nd.fl_vel_dot[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_temp_dfd0(& self, glob_temp : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            glob_temp[i1].set_val(this_nd.temperature);
        }
        return;
    }

    pub fn get_nd_tdot_dfd0(& self, glob_tdot : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            glob_tdot[i1].set_val(this_nd.temp_change_rate);
        }
        return;
    }

    pub fn get_nd_fl_den_dfd0(& self, fl_den : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            fl_den[i1].set_val(this_nd.fl_den);
        }
        
        return;
    }

    pub fn get_nd_fl_den_dot_dfd0(& self, fl_den_dot : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            fl_den_dot[i1].set_val(this_nd.fl_den_dot);
        }
        
        return;
    }

    pub fn get_nd_turb_e_dfd0(& self, turb_e : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            turb_e[i1].set_val(this_nd.turb_e);
        }
        
        return;
    }

    pub fn get_nd_turb_edot_dfd0(& self, turb_edot : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            turb_edot[i1].set_val(this_nd.turb_edot);
        }
        
        return;
    }

    pub fn eval_n_dfd0(&self, n_vec : &mut [DiffDoub0], d_nds : &mut [DiffDoub0], spt : & [f64]) {
        if self.this_type == 4 || self.this_type == 400 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]-spt[2]);
            d_nds[0].set_val(-1.0);
            d_nds[1].set_val(-1.0);
            d_nds[2].set_val(-1.0);
            
            n_vec[1].set_val(spt[0]);
            d_nds[3].set_val(1.0);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(spt[1]);
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(1.0);
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(spt[2]);
            d_nds[9].set_val(0.0);
            d_nds[10].set_val(0.0);
            d_nds[11].set_val(1.0);
        } else if self.this_type == 6 || self.this_type == 600 {
            n_vec[0].set_val(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
            d_nds[0].set_val(-0.5*(1.0-spt[2]));
            d_nds[1].set_val(-0.5*(1.0-spt[2]));
            d_nds[2].set_val(-0.5*(1.0-spt[0]-spt[1]));
            
            n_vec[1].set_val(0.5*spt[0]*(1.0-spt[2]));
            d_nds[3].set_val(0.5*(1.0-spt[2]));
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(-0.5*spt[0]);
            
            n_vec[2].set_val(0.5*spt[1]*(1.0-spt[2]));
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(0.5*(1.0-spt[2]));
            d_nds[8].set_val(-0.5*spt[1]);
            
            n_vec[3].set_val(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
            d_nds[9].set_val(-0.5*(1.0+spt[2]));
            d_nds[10].set_val(-0.5*(1.0+spt[2]));
            d_nds[11].set_val(0.5*(1.0-spt[0]-spt[1]));
            
            n_vec[4].set_val(0.5*spt[0]*(1.0+spt[2]));
            d_nds[12].set_val(0.5*(1.0+spt[2]));
            d_nds[13].set_val(0.0);
            d_nds[14].set_val(0.5*spt[0]);
            
            n_vec[5].set_val(0.5*spt[1]*(1.0+spt[2]));
            d_nds[15].set_val(0.0);
            d_nds[16].set_val(0.5*(1.0+spt[2]));
            d_nds[17].set_val(0.5*spt[1]);
        } else if self.this_type == 8 || self.this_type == 81 || self.this_type == 800 {
            n_vec[0].set_val(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
            d_nds[0].set_val(-0.125*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[1].set_val(-0.125*(1.0-spt[0])*(1.0-spt[2]));
            d_nds[2].set_val(-0.125*(1.0-spt[0])*(1.0-spt[1]));
            
            n_vec[1].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[3].set_val(0.125*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[4].set_val(-0.125*(1.0+spt[0])*(1.0-spt[2]));
            d_nds[5].set_val(-0.125*(1.0+spt[0])*(1.0-spt[1]));
            
            n_vec[2].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[6].set_val(0.125*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[7].set_val(0.125*(1.0+spt[0])*(1.0-spt[2]));
            d_nds[8].set_val(-0.125*(1.0+spt[0])*(1.0+spt[1]));
            
            n_vec[3].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[9].set_val(-0.125*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[10].set_val(0.125*(1.0-spt[0])*(1.0-spt[2]));
            d_nds[11].set_val(-0.125*(1.0-spt[0])*(1.0+spt[1]));
            
            n_vec[4].set_val(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[12].set_val(-0.125*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[13].set_val(-0.125*(1.0-spt[0])*(1.0+spt[2]));
            d_nds[14].set_val(0.125*(1.0-spt[0])*(1.0-spt[1]));
            
            n_vec[5].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[15].set_val(0.125*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[16].set_val(-0.125*(1.0+spt[0])*(1.0+spt[2]));
            d_nds[17].set_val(0.125*(1.0+spt[0])*(1.0-spt[1]));
            
            n_vec[6].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[18].set_val(0.125*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[19].set_val(0.125*(1.0+spt[0])*(1.0+spt[2]));
            d_nds[20].set_val(0.125*(1.0+spt[0])*(1.0+spt[1]));
            
            n_vec[7].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[21].set_val(-0.125*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[22].set_val(0.125*(1.0-spt[0])*(1.0+spt[2]));
            d_nds[23].set_val(0.125*(1.0-spt[0])*(1.0+spt[1]));
            
            if self.this_type == 81 {
                n_vec[8].set_val(1.0-spt[0]*spt[0]);
                d_nds[24].set_val(-2.0*spt[0]);
                d_nds[25].set_val(0.0);
                d_nds[26].set_val(0.0);
                
                n_vec[9].set_val(1.0-spt[1]*spt[1]);
                d_nds[27].set_val(0.0);
                d_nds[28].set_val(-2.0*spt[1]);
                d_nds[29].set_val(0.0);
                
                n_vec[10].set_val(1.0-spt[2]*spt[2]);
                d_nds[30].set_val(0.0);
                d_nds[31].set_val(0.0);
                d_nds[32].set_val(-2.0*spt[2]);
            }
        }
        else if self.this_type == 10 || self.this_type == 1000 {
            let p1 : f64 =  1.0 - spt[0] - spt[1] - spt[2];
            let p2 : f64 =  p1 - 0.5;
            n_vec[0].set_val(2.0 * p1 * p2);
            d_nds[0].set_val(2.0 * (-p2 - p1));
            d_nds[1].set_val(2.0 * (-p2 - p1));
            d_nds[2].set_val(2.0 * (-p2 - p1));
            
            n_vec[1].set_val(-2.0 * spt[0] * (0.5 - spt[0]));
            d_nds[3].set_val(-2.0 * (0.5 - spt[0] - spt[0]));
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(-2.0 * spt[1] * (0.5 - spt[1]));
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(-2.0 * (0.5 - spt[1] - spt[1]));
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(-2.0 * spt[2] * (0.5 - spt[2]));
            d_nds[9].set_val(0.0);
            d_nds[10].set_val(0.0);
            d_nds[11].set_val(-2.0 * (0.5 - spt[2] - spt[2]));
            
            n_vec[4].set_val(4.0 * spt[0] * p1);
            d_nds[12].set_val(4.0 * (p1 - spt[0]));
            d_nds[13].set_val(-4.0 * spt[0]);
            d_nds[14].set_val(-4.0 * spt[0]);
            
            n_vec[5].set_val(4.0 * spt[0] * spt[1]);
            d_nds[15].set_val(4.0 * spt[1]);
            d_nds[16].set_val(4.0 * spt[0]);
            d_nds[17].set_val(0.0);
            
            n_vec[6].set_val(4.0 * spt[1] * p1);
            d_nds[18].set_val(-4.0 * spt[1]);
            d_nds[19].set_val(4.0 * (p1 - spt[1]));
            d_nds[20].set_val(-4.0 * spt[1]);
            
            n_vec[7].set_val(4.0 * spt[2] * p1);
            d_nds[21].set_val(-4.0 * spt[2]);
            d_nds[22].set_val(-4.0 * spt[2]);
            d_nds[23].set_val(4.0 * (p1 - spt[2]));
            
            n_vec[8].set_val(4.0 * spt[0] * spt[2]);
            d_nds[24].set_val(4.0 * spt[2]);
            d_nds[25].set_val(0.0);
            d_nds[26].set_val(4.0 * spt[0]);
            
            n_vec[9].set_val(4.0 * spt[1] * spt[2]);
            d_nds[27].set_val(0.0);
            d_nds[28].set_val(4.0 * spt[2]);
            d_nds[29].set_val(4.0 * spt[1]);
        }
        else if self.this_type == 3 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]);
            d_nds[0].set_val(-1.0);
            d_nds[1].set_val(-1.0);
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(spt[0]);
            d_nds[3].set_val(1.0);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(spt[1]);
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(1.0);
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(spt[0]*(1.0-spt[0]-spt[1]));
            d_nds[9].set_val(1.0-spt[0]-spt[1] - spt[0]);
            d_nds[10].set_val(-spt[0]);
            d_nds[11].set_val(0.0);
            
            n_vec[4].set_val(spt[0]*spt[1]);
            d_nds[12].set_val(spt[1]);
            d_nds[13].set_val(spt[0]);
            d_nds[14].set_val(0.0);
            
            n_vec[5].set_val(spt[1]*(1.0-spt[0]-spt[1]));
            d_nds[15].set_val(-spt[1]);
            d_nds[16].set_val(1.0-spt[0]-spt[1]-spt[1]);
            d_nds[17].set_val(0.0);
            
        } else if self.this_type == 41 {
            n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));
            d_nds[0].set_val(-0.25*(1.0-spt[1]));
            d_nds[1].set_val(-0.25*(1.0-spt[0]));
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
            d_nds[3].set_val(0.25*(1.0-spt[1]));
            d_nds[4].set_val(-0.25*(1.0+spt[0]));
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
            d_nds[6].set_val(0.25*(1.0+spt[1]));
            d_nds[7].set_val(0.25*(1.0+spt[0]));
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
            d_nds[9].set_val(-0.25*(1.0+spt[1]));
            d_nds[10].set_val(0.25*(1.0-spt[0]));
            d_nds[11].set_val(0.0);
            
            n_vec[4].set_val(1.0-spt[0]*spt[0]);
            d_nds[12].set_val(-2.0*spt[0]);
            d_nds[13].set_val(0.0);
            d_nds[14].set_val(0.0);
            
            n_vec[5].set_val(1.0-spt[1]*spt[1]);
            d_nds[15].set_val(0.0);
            d_nds[16].set_val(-2.0*spt[1]);
            d_nds[17].set_val(0.0);
            
            n_vec[6].set_val((1.0-spt[0]*spt[0])*(1.0-spt[1]));
            d_nds[18].set_val(-2.0*spt[0]*(1.0-spt[1]));
            d_nds[19].set_val(-1.0+spt[0]*spt[0]);
            d_nds[20].set_val(0.0);
            
            n_vec[7].set_val((1.0-spt[1]*spt[1])*(1.0+spt[0]));
            d_nds[21].set_val(1.0-spt[1]*spt[1]);
            d_nds[22].set_val(-2.0*spt[1]*(1.0+spt[0]));
            d_nds[23].set_val(0.0);
            
            n_vec[8].set_val((1.0-spt[0]*spt[0])*(1.0+spt[1]));
            d_nds[24].set_val(-2.0*spt[0]*(1.0+spt[1]));
            d_nds[25].set_val(1.0-spt[0]*spt[0]);
            d_nds[26].set_val(0.0);
            
            n_vec[9].set_val((1.0-spt[1]*spt[1])*(1.0-spt[0]));
            d_nds[27].set_val(-1.0+spt[1]*spt[1]);
            d_nds[28].set_val(-2.0*spt[1]*(1.0-spt[0]));
            d_nds[29].set_val(0.0);
            
        } else if self.this_type == 2 {
            n_vec[0].set_val(0.5*(1.0 - spt[0]));
            d_nds[0].set_val(-0.5);
            d_nds[1].set_val(0.0);
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(0.5*(1.0 + spt[0]));
            d_nds[3].set_val(0.5);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(1.0 - spt[0]*spt[0]);
            d_nds[0].set_val(-2.0*spt[0]);
            d_nds[1].set_val(0.0);
            d_nds[2].set_val(0.0);
        }
        return;
    }

    pub fn get_ip_data_dfd0(&self, n_vec : &mut [DiffDoub0], d_ndx : &mut [DiffDoub0], det_j : &mut DiffDoub0, loc_nds : &mut Vec<DiffDoub0>, spt : & [f64]) {
        let i1 : usize;
        let i2 : usize;
        let mut n_cent = [DiffDoub0::new(); 11];
        let mut d_nds = [DiffDoub0::new(); 33];
        let mut d_nds_cent = [DiffDoub0::new(); 33];
        let mut j_mat = [DiffDoub0::new(); 9];
        let mut tmp_nds = [DiffDoub0::new(); 30];
        let mut j_cent = [DiffDoub0::new(); 9];
        let mut det_cent = DiffDoub0::new();
        let mut j_inv = [DiffDoub0::new(); 9];
        let mut j_inv_cent = [DiffDoub0::new(); 9];
        let mut x_vec = [DiffDoub0::new(); 3];
        let mut b_vec = [DiffDoub0::new(); 3];
        let mut z_dir = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        
        self.eval_n_dfd0(n_vec, &mut d_nds, spt);
        vec_to_ar_dfd0(&mut tmp_nds, loc_nds,  0,  30);
        mat_mul_ar_dfd0(&mut j_mat, &mut tmp_nds, &mut d_nds, 3, self.num_nds, 3);
        
        if self.this_type == 41 || self.this_type == 3 {
            z_dir.set_val_dfd0(& j_mat[0]);
            z_dir.mult(& j_mat[4]);
            tmp.set_val_dfd0(& j_mat[3]);
            tmp.mult(& j_mat[1]);
            z_dir.sub(& tmp);
            if z_dir.val > 0.0 {
                j_mat[8].set_val(1.0);
            } else {
                j_mat[8].set_val(-1.0);
            }
        } else if self.this_type == 2 {
            j_mat[4].set_val(1.0);
            if j_mat[0].val > 0.0 {
                j_mat[8].set_val(1.0);
            } else {
                j_mat[8].set_val(-1.0);
            }
        }
        
        get_det_inv_ar_dfd0(det_j, &mut j_inv, &mut j_mat, 3, 0, &mut x_vec, &mut b_vec);
        
        // mat_mul_dfd0(d_nds,d_nds,j_inv,self.n_dim,3,3);
        mat_mul_ar_dfd0(d_ndx, &mut d_nds, &mut j_inv, self.num_nds, 3, 3);
        
        if self.n_dim > self.num_nds {
            self.eval_n_dfd0(&mut n_cent, &mut  d_nds_cent, & self.s_cent);
            mat_mul_ar_dfd0(&mut j_cent, &mut  tmp_nds, &mut  d_nds_cent,  3,  self.num_nds,  3);
            
            if self.this_type == 41 || self.this_type == 3 {
                z_dir.set_val_dfd0(& j_cent[0]);
                z_dir.mult(& j_cent[4]);
                tmp.set_val_dfd0(& j_cent[3]);
                tmp.mult(& j_cent[1]);
                z_dir.sub(& tmp);
                if z_dir.val > 0.0 {
                    j_cent[8].set_val(1.0);
                }
                else {
                    j_cent[8].set_val(-1.0);
                }
            }
            else if self.this_type == 2 {
                j_cent[4].set_val(1.0);
                if j_cent[0].val > 0.0 {
                    j_cent[8].set_val(1.0);
                }
                else {
                    j_cent[8].set_val(-1.0);
                }
            }
            
            get_det_inv_ar_dfd0(&mut det_cent, &mut  j_inv_cent, &mut  j_cent,  3,  0, &mut  x_vec, &mut  b_vec);
            
            i1 = 3 * self.num_nds;
            i2 = self.n_dim - self.num_nds;
            mat_mul_ar_dfd0(&mut d_ndx[i1..], &mut d_nds[i1..], &mut  j_inv_cent,  i2,  3,  3);
        }
        return;
    }

    pub fn get_inst_ori_dfd0(& self, inst_ori_mat : &mut Vec<DiffDoub0>, loc_ori : &mut Vec<DiffDoub0>, glob_disp : &mut Vec<DiffDoub0>, stat : usize) {
        // stat = 1: nonlinear geometry, Element ori from nodal theta, no derivatives, 1st order diff for self.nodes
        // stat = 2: nonlinear geometry, 2nd order diff for DiffDoub0 version, 1st order for DiffDoub1 version
        let mut rot = [DiffDoub0::new(); 3];
        let mut nnds = DiffDoub0::new();
        let mut one = DiffDoub0::new();
        let mut tmp_ori = [DiffDoub0::new(); 9];
        let mut tmp_inst = [DiffDoub0::new(); 9];
        let mut tmp = DiffDoub0::new();
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut st_index : usize;
        let i_ori_size : usize =  (self.num_nds+1)*144;
        for i1 in 0..i_ori_size {
            inst_ori_mat[i1].set_val(0.0);
        }
        let is_diff : bool =  loc_ori[0].diff_type();
        
        nnds.set_val(0.0);
        one.set_val(1.0);
        rot[0].set_val(0.0);
        rot[1].set_val(0.0);
        rot[2].set_val(0.0);
        i2 = 3 * self.n_dim;
        for _i1 in 0..self.num_nds {
            rot[0].add(& glob_disp[i2]);
            rot[1].add(& glob_disp[i2 + self.n_dim]);
            rot[2].add(& glob_disp[i2 + 2 * self.n_dim]);
            nnds.add(& one);
            i2 += 1usize;
        }
        rot[0].dvd(& nnds);
        rot[1].dvd(& nnds);
        rot[2].dvd(& nnds);
        
        vec_to_ar_dfd0(&mut tmp_ori, loc_ori,  0,  9);
        
        if stat == 1 {
            d_orid_thet_dfd0(&mut tmp_inst, &mut  tmp_ori, &mut  rot,  0,  0);
            ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  0,  9);
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd0(& glob_disp[i2]);
                rot[1].set_val_dfd0(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd0(& glob_disp[i2+2*self.n_dim]);
                st_index = 144 * i1;
                d_orid_thet_dfd0(&mut tmp_inst, &mut  tmp_ori, &mut  rot,  0,  0);
                ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                for i2 in 1..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd0(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd0(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
        } else if is_diff {
            for i2 in 0..4 {
                st_index = 36*i2;
                d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                i3 = st_index;
                i4 = 9*i2;
                for _i5 in 0..9 {
                    tmp.set_val_dfd0(& inst_ori_mat[i3]);
                    inst_ori_mat[i4].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4 += 1usize;
                }
            }
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd0(& glob_disp[i2]);
                rot[1].set_val_dfd0(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd0(& glob_disp[i2+2*self.n_dim]);
                for i2 in 0..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd0(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd0(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
        } else {
            for i2 in 0..4 {
                st_index = 36*i2;
                d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                i3 = st_index;
                i4 = 9*i2;
                for _i5 in 0..9 {
                    tmp.set_val_dfd0(& inst_ori_mat[i3]);
                    inst_ori_mat[i4].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4 += 1usize;
                }
                for i6 in i2..4 {
                    st_index = 36*i2 + 9*i6;
                    d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, i6);
                    ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 36*i6 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd0(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd0(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd0(& glob_disp[i2]);
                rot[1].set_val_dfd0(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd0(& glob_disp[i2+2*self.n_dim]);
                for i2 in 0..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd0(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd0(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                    for i6 in i2..4 {
                        st_index = 144*i1 + 36*i2 + 9*i6;
                        d_orid_thet_dfd0(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, i6);
                        ar_to_vec_dfd0(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                        i3 = st_index;
                        i4 = 144*i1 + 36*i6 + 9*i2;
                        for _i5 in 0..9 {
                            tmp.set_val_dfd0(& inst_ori_mat[i3]);
                            inst_ori_mat[i4].set_val_dfd0(& tmp);
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_inst_disp_dfd0(& self, inst_disp : &mut [DiffDoub0], glob_disp : &mut Vec<DiffDoub0>, inst_ori_mat : &mut Vec<DiffDoub0>, loc_ori : &mut Vec<DiffDoub0>, x_glob : &mut Vec<DiffDoub0>, n_lgeom : bool, dv1 : usize, dv2 : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let mut nd_oind : usize;
        let mut dof : usize;
        let mut dof_oind : usize;
        let nd2 : usize;
        let mut dof2 : usize;
        let dof2_oind : usize;
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        let mut nn_inv = DiffDoub0::new();
        let mut nn_inv2 = DiffDoub0::new();
        
        i2 = 6*self.n_dim;
        for i1 in 0..i2 {
            inst_disp[i1].set_val(0.0);
        }
        
        if !n_lgeom {
            if dv1 == MAX_INT && dv2 == MAX_INT {
                i7 = 3 * self.n_dim;
                for i1 in 0..3 {
                    i4 = i1 * self.n_dim;
                    for i2 in 0..self.num_nds {
                        i5 = i1 * 3;
                        i6 = i2;
                        for _i3 in 0..3 {
                            //i4 = i1 * self.n_dim + i2;
                            //i5 = i1 * 3 + i3;
                            //i6 = i3 * self.n_dim + i2;
                            tmp.set_val_dfd0(& loc_ori[i5]);
                            tmp.mult(& glob_disp[i6]);
                            inst_disp[i4].add(& tmp);
                            tmp.set_val_dfd0(& loc_ori[i5]);
                            tmp.mult(& glob_disp[i6 + i7]);
                            inst_disp[i4 + i7].add(& tmp);
                            i5 += 1usize;
                            i6  +=  self.n_dim;
                        }
                        i4 += 1usize;
                    }
                }
                i2 = 2 * self.num_nds * self.dof_per_nd;
                for _i1 in 0..self.num_int_dof {
                    nd = self.dof_table[i2];
                    dof = self.dof_table[i2 + 1];
                    i3 = dof * self.n_dim + nd;
                    inst_disp[i3].set_val_dfd0(& glob_disp[i3]);
                    i2  +=  2;
                }
            }
            else if (dv1 + dv2) >= MAX_INT {
                let dnz = dv1 + dv2 - MAX_INT;
                //dv1 = dv1 + dv2 - MAX_INT;
                nd = self.dof_table[2 * dnz];
                dof = self.dof_table[2 * dnz + 1];
                if dof < 3 {
                    if nd < self.num_nds {
                        i2 = nd;
                        i3 = dof;
                        for _i1 in 0..3 {
                            //i2 = i1 * self.n_dim + nd;
                            //i3 = i1 * 3 + dof;
                            inst_disp[i2].set_val_dfd0(& loc_ori[i3]);
                            i2  +=  self.n_dim;
                            i3  +=  3;
                        }
                    }
                    else {
                        i1 = dof * self.n_dim + nd;
                        inst_disp[i1].set_val(1.0);
                    }
                }
                else {
                    i2 = 3 * self.n_dim + nd;
                    i3 = dof - 3;
                    for _i1 in 0..3 {
                        //i2 = (i1 + 3) * self.n_dim + nd;
                        //i3 = i1 * 3 + (dof - 3);
                        inst_disp[i2].set_val_dfd0(& loc_ori[i3]);
                        i2  +=  self.n_dim;
                        i3  +=  3;
                    }
                }
            }
        }
        else {
            if dv1 == MAX_INT && dv2 == MAX_INT {
                for i1 in 0..3 {
                    for i2 in 0..self.num_nds {
                        i4 = i1 * self.n_dim + i2;
                        i5 = i2;
                        i6 = i2;
                        i7 = i1 * 3;
                        for _i3 in 0..3 {
                            tmp.set_val_dfd0(& glob_disp[i5]);
                            tmp.add(& x_glob[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp2.set_val_dfd0(& loc_ori[i7]);
                            tmp2.mult(& x_glob[i6]);
                            tmp.sub(& tmp2);
                            inst_disp[i4].add(& tmp);
                            i5  +=  self.n_dim;
                            i6  +=  self.num_nds;
                            i7 += 1usize;
                        }
                    }
                }
                
                i2 = self.num_nds * self.dof_per_nd;
                i3 = i2 + self.num_int_dof;
                for i1 in i2..i3 {
                    i4 = self.dof_table[2 * i1];
                    i5 = self.dof_table[2 * i1 + 1];
                    i6 = i5 * self.n_dim + i4;
                    inst_disp[i6].set_val_dfd0(& glob_disp[i6]);
                }
                
                for i1 in 0..self.num_nds {
                    i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                    i4 = 4 * self.n_dim + i1;
                    i5 = 5 * self.n_dim + i1;
                    nd_oind = 144 * (i1 + 1);
                    for i2 in 0..3 {
                        i6 = nd_oind + 3 + i2;
                        i7 = 6 + i2;
                        tmp.set_val_dfd0(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i3].add(& tmp);
                        i6 = nd_oind + 6 + i2;
                        i7 = i2;
                        tmp.set_val_dfd0(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i4].add(& tmp);
                        i6 = nd_oind + i2;
                        i7 = 3 + i2;
                        tmp.set_val_dfd0(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i5].add(& tmp);
                    }
                }
            }
            else if (dv1 + dv2) >= MAX_INT {
                let dnz = dv1 + dv2 - MAX_INT;
                //dv1 = dv1 + dv2 - MAX_INT;
                nd = self.dof_table[2 * dnz];
                dof = self.dof_table[2 * dnz + 1];
                if dof < 3 {
                    if nd < self.num_nds {
                        i1 = nd;// index in inst_disp
                        i2 = dof;//index in inst_ori_mat
                        inst_disp[i1].set_val_dfd0(& inst_ori_mat[i2]);
                        i1 = self.n_dim + nd;
                        i2 = 3 + dof;
                        inst_disp[i1].set_val_dfd0(& inst_ori_mat[i2]);
                        i1 = 2 * self.n_dim + nd;
                        i2 = 6 + dof;
                        inst_disp[i1].set_val_dfd0(& inst_ori_mat[i2]);
                    }
                    else {
                        i1 = dof * self.n_dim + nd;
                        inst_disp[i1].set_val(1.0);
                    }
                }
                else {// dof is rotation
                    nn_inv.set_val(1.0 / (self.num_nds as f64));
                    dof_oind = 36 * (dof - 2);
                    for i1 in 0..3 {
                        for i2 in 0..self.num_nds {
                            i4 = i1 * self.n_dim + i2;
                            i5 = i2;
                            i6 = i2;
                            i7 = dof_oind + i1 * 3;
                            for _i3 in 0..3 {
                                tmp.set_val_dfd0(& glob_disp[i5]);
                                tmp.add(& x_glob[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i5  +=  self.n_dim;
                                i6  +=  self.num_nds;
                                i7 += 1usize;
                            }
                        }
                    }
                    
                    for i1 in 0..self.num_nds {
                        i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                        i4 = 4 * self.n_dim + i1;
                        i5 = 5 * self.n_dim + i1;
                        nd_oind = 144 * (i1 + 1);
                        for i2 in 0..3 {
                            i6 = nd_oind + 3 + i2;
                            i7 = dof_oind + 6 + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i3].add(& tmp);
                            i6 = nd_oind + 6 + i2;
                            i7 = dof_oind + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i4].add(& tmp);
                            i6 = nd_oind + i2;
                            i7 = dof_oind + 3 + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i5].add(& tmp);
                            if i1 == nd {
                                i6 = nd_oind + dof_oind + 3 + i2;
                                i7 = 6 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + 6 + i2;
                                i7 = i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + i2;
                                i7 = 3 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i5].add(& tmp);
                            }
                        }
                    }
                }
            }
            else {
                nd = self.dof_table[2 * dv1];
                dof = self.dof_table[2 * dv1 + 1];
                nd2 = self.dof_table[2 * dv2];
                dof2 = self.dof_table[2 * dv2 + 1];
                nn_inv.set_val(1.0 / (self.num_nds as f64));
                nn_inv2.set_val_dfd0(& nn_inv);
                nn_inv2.sqr();
                if dof > 2 && dof2 > 2 {
                    for i1 in 0..3 {
                        dof_oind = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
                        for i2 in 0..self.num_nds {
                            i4 = self.n_dim * i1 + i2;
                            i5 = i2;
                            i6 = i2;
                            i7 = dof_oind;
                            for _i3 in 0..3 {
                                tmp.set_val_dfd0(& glob_disp[i5]);
                                tmp.add(& x_glob[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv2);
                                inst_disp[i4].add(& tmp);
                                i5  +=  self.n_dim;
                                i6  +=  self.num_nds;
                                i7 += 1usize;
                            }
                        }
                    }
                    
                    dof_oind = 36 * (dof - 2);
                    dof2_oind = 9 * (dof2 - 2);
                    for i1 in 0..self.num_nds {
                        i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                        i4 = 4 * self.n_dim + i1;
                        i5 = 5 * self.n_dim + i1;
                        nd_oind = 144 * (i1 + 1);
                        for i2 in 0..3 {
                            i6 = nd_oind + 3 + i2;
                            i7 = dof_oind + dof2_oind + 6 + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i3].add(& tmp);
                            i6 = nd_oind + 6 + i2;
                            i7 = dof_oind + dof2_oind + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i4].add(& tmp);
                            i6 = nd_oind + i2;
                            i7 = dof_oind + dof2_oind + 3 + i2;
                            tmp.set_val_dfd0(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i5].add(& tmp);
                            if i1 == nd {
                                i6 = nd_oind + dof_oind + 3 + i2;
                                i7 = dof2_oind + 6 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + 6 + i2;
                                i7 = dof2_oind + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + i2;
                                i7 = dof2_oind + 3 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i5].add(& tmp);
                            }
                            if i1 == nd2 {
                                i6 = nd_oind + dof2_oind + 3 + i2;
                                i7 = dof_oind + 6 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof2_oind + 6 + i2;
                                i7 = dof_oind + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof2_oind + i2;
                                i7 = dof_oind + 3 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i5].add(& tmp);
                            }
                            if i1 == nd && i1 == nd2 {
                                i6 = nd_oind + dof_oind + dof2_oind + 3 + i2;
                                i7 = 6 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + dof2_oind + 6 + i2;
                                i7 = i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + dof2_oind + i2;
                                i7 = 3 + i2;
                                tmp.set_val_dfd0(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i5].add(& tmp);
                            }
                        }
                    }
                }
                else if dof < 3 && dof2 < 3 {
                    return;
                }
                else {
                    if dof > 2 {
                        i1 = dof;
                        dof = dof2;
                        dof2 = i1;
                        //i1 = nd;
                        nd = nd2;
                        //nd2 = i1;
                    }
                    i1 = 36 * (dof2 - 2) + dof;
                    tmp.set_val_dfd0(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = nd;
                    inst_disp[i2].add(& tmp);
                    i1 = 36 * (dof2 - 2) + 3 + dof;
                    tmp.set_val_dfd0(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = self.n_dim + nd;
                    inst_disp[i2].add(& tmp);
                    i1 = 36 * (dof2 - 2) + 6 + dof;
                    tmp.set_val_dfd0(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = 2 * self.n_dim + nd;
                    inst_disp[i2].add(& tmp);
                }
            }
        }

        return;
    }

    pub fn get_stress_prereq_dfd0(& self, pre : &mut DiffDoub0StressPrereq, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut offset = DiffDoub0::new();
        self.get_nd_crds_dfd0(&mut pre.glob_nds, nd_ar, dv_ar);
        self.get_loc_ori_dfd0(&mut pre.loc_ori, sec_ar, dv_ar);
        self.get_nd_disp_dfd0(&mut pre.glob_disp, nd_ar);
        self.get_nd_vel_dfd0(&mut pre.glob_vel, nd_ar);
        self.get_nd_acc_dfd0(&mut pre.glob_acc, nd_ar);
        self.get_nd_temp_dfd0(&mut pre.glob_temp, nd_ar);
        self.get_nd_tdot_dfd0(&mut pre.glob_tdot, nd_ar);
        if self.dof_per_nd == 6 {
            self.correct_orient_dfd0(&mut pre.loc_ori, &mut  pre.glob_nds);
            if self.this_type != 2 {
                self.get_layer_thk_z_dfd0(&mut pre.layer_thk, &mut  pre.layer_z, &mut  offset, sec_ar, dv_ar);
                self.get_layer_angle_dfd0(&mut pre.layer_ang, sec_ar, dv_ar);
                self.get_layer_q_dfd0(&mut pre.layer_q, sec_ar, mat_ar, dv_ar);
                self.get_layer_d_dfd0(&mut pre.layer_d, sec_ar, mat_ar, dv_ar);
                self.get_layer_th_exp_dfd0(&mut pre.layer_te, sec_ar, mat_ar, dv_ar);
                self.get_layer_einit_dfd0(&mut pre.layer_e0, sec_ar, dv_ar);
                self.get_layer_den_dfd0(&mut pre.layer_den, sec_ar, mat_ar, dv_ar);
                self.get_layer_cond_dfd0(&mut pre.layer_tc, sec_ar, mat_ar, dv_ar);
                self.get_layer_spec_heat_dfd0(&mut pre.layer_sh, sec_ar, mat_ar, dv_ar);
                self.get_abd_dfd0(&mut pre.cmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_q, &mut  pre.layer_ang, sec_ar);
                self.get_shell_damp_dfd0(&mut pre.dmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_d, &mut  pre.layer_ang, sec_ar);
                self.get_shell_exp_load_dfd0(&mut pre.therm_exp, &mut  pre.einit, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_q, &mut  pre.layer_te, &mut  pre.layer_e0, &mut  pre.layer_ang, sec_ar);
                self.get_shell_mass_dfd0(&mut pre.mmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_den, sec_ar);
                self.get_shell_cond_dfd0(&mut pre.tcmat, &mut  pre.layer_thk, &mut  pre.layer_ang, &mut  pre.layer_tc, sec_ar);
                self.get_shell_spec_heat_dfd0(&mut pre.spec_heat, &mut  pre.layer_thk, &mut  pre.layer_sh, &mut  pre.layer_den, sec_ar);
            }
            else {
                self.get_beam_stiff_dfd0(&mut pre.cmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_damp_dfd0(&mut pre.dmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_exp_load_dfd0(&mut pre.therm_exp, &mut  pre.einit, sec_ar, mat_ar, dv_ar);
                self.get_beam_mass_dfd0(&mut pre.mmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_cond_dfd0(&mut pre.tcmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_spec_heat_dfd0(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
            }
        }
        else if self.this_type == 21 {
            self.get_frc_fld_const_dfd0(&mut pre.frc_fld_coef, &mut  pre.frc_fld_exp, sec_ar, dv_ar);
            self.get_thrm_fld_const_dfd0(&mut pre.thrm_fld_coef, &mut  pre.ref_temp, sec_ar, dv_ar);
        }
        else if self.this_type == 1 {
            self.get_mass_per_el_dfd0(&mut pre.mass_per_el, sec_ar, dv_ar);
            self.get_specific_heat_dfd0(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
        }
        else {
            self.get_solid_stiff_dfd0(&mut pre.cmat, sec_ar, mat_ar, dv_ar);
            self.get_solid_damp_dfd0(&mut pre.dmat, sec_ar, mat_ar, dv_ar);
            self.get_thermal_exp_dfd0(&mut pre.therm_exp, &mut  pre.einit, sec_ar, mat_ar, dv_ar);
            self.get_density_dfd0(&mut pre.mmat[0],  0, sec_ar, mat_ar, dv_ar);
            self.get_conductivity_dfd0(&mut pre.tcmat, sec_ar, mat_ar, dv_ar);
            self.get_specific_heat_dfd0(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
        }
        mat_mul_dfd0(&mut pre.loc_nds, &mut  pre.loc_ori, &mut  pre.glob_nds,  3,  3,  self.num_nds);
        
        
        return;
    }

    pub fn get_fluid_prereq_dfd0(&mut self, pre : &mut DiffDoub0FlPrereq, sec_ar : &mut Vec<Section>, fl_ar : &mut Vec<Fluid>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>, el_cent : & Vec<CentData>) {
        let i2 : usize;
        let mut sec_prop : f64;
        
        pre.num_nds = self.num_nds;
        self.get_nd_crds_dfd0(&mut pre.glob_nds, nd_ar, dv_ar);
        self.get_nd_disp_dfd0(&mut pre.glob_disp, nd_ar);
        i2 = 3 * self.num_nds;
        for i1 in 0..i2 {
            pre.glob_nds[i1].add(& pre.glob_disp[i1]);
        }
        self.get_nd_vel_dfd0(&mut pre.glob_vel, nd_ar);
        self.get_nd_fl_den_dfd0(&mut pre.fl_den, nd_ar);
        self.get_nd_fl_vel_dfd0(&mut pre.fl_vel, nd_ar);
        self.get_nd_temp_dfd0(&mut pre.fl_temp, nd_ar);
        self.get_nd_turb_e_dfd0(&mut pre.fl_turb_e, nd_ar);
        self.get_nd_fl_den_dot_dfd0(&mut pre.fl_den_dot, nd_ar);
        self.get_nd_fl_vdot_dfd0(&mut pre.fl_vel_dot, nd_ar);
        self.get_nd_tdot_dfd0(&mut pre.fl_tdot, nd_ar);
        self.get_nd_turb_edot_dfd0(&mut pre.fl_turb_edot, nd_ar);

        let ec : &CentData = &el_cent[self.label];
        for i in 0..3 {
            pre.f_per_mass[i].set_val_dfd0(&ec.f_dfd0[i]);
        }
        self.get_f_per_mass_dfd0(&mut pre.f_per_mass, dv_ar);
        pre.hg_per_vol.set_val_dfd0(&ec.hg_dfd0);
        self.get_gen_prop_dfd0(&mut pre.hg_per_vol, &mut CppStr::from("bodyHeatGen"), dv_ar);
        
        let this_sec = &sec_ar[self.sect_ptr];
        let this_fl = &fl_ar[this_sec.fl_ptr];

        pre.compressible = this_fl.compressible;
        
        sec_prop = this_fl.grad_turb_coef;
        pre.grad_turb_coef.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.grad_turb_coef, &mut CppStr::from("gradVTurbCoef"), dv_ar);
        
        sec_prop = this_fl.diss_turb_coef;
        pre.diss_turb_coef.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.diss_turb_coef, &mut CppStr::from("dissTurbCoef"), dv_ar);
        
        sec_prop = this_fl.viscosity;
        pre.ref_visc.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.ref_visc, &mut CppStr::from("viscosity"), dv_ar);
        
        sec_prop = this_fl.temp_vis_coef;
        pre.temp_vis_coef.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.temp_vis_coef, &mut CppStr::from("tempVisCoef"), dv_ar);
        
        sec_prop = this_fl.turb_vis_coef;
        pre.turb_vis_coef.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.turb_vis_coef, &mut CppStr::from("turbVisCoef"), dv_ar);
        
        sec_prop = this_fl.ref_enth;
        pre.ref_enth.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.ref_enth, &mut CppStr::from("refEnth"), dv_ar);
        
        sec_prop = this_fl.ref_den;
        pre.ref_den.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.ref_den, &mut CppStr::from("refDen"), dv_ar);
        
        sec_prop = this_fl.ref_temp;
        pre.ref_temp.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.ref_temp, &mut CppStr::from("refTemp"), dv_ar);
        
        sec_prop = this_fl.therm_cond;
        pre.therm_cond.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.therm_cond, &mut CppStr::from("thermCond"), dv_ar);

        sec_prop = this_fl.expansion;
        pre.expansion.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.expansion, &mut CppStr::from("thermExp"), dv_ar);
        
        sec_prop = this_fl.spec_heat;
        pre.spec_heat.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.spec_heat, &mut CppStr::from("specHeat"), dv_ar);
        
        sec_prop = this_fl.ideal_gas;
        pre.i_gconst.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.i_gconst, &mut CppStr::from("iGConst"), dv_ar);

        sec_prop = this_fl.bulk_modulus;
        pre.bulk_mod.set_val(sec_prop);
        self.get_gen_prop_dfd0(&mut pre.bulk_mod, &mut CppStr::from("bulkModulus"), dv_ar);

        if pre.compressible {
            pre.den_vis_coef.set_val_dfd0(&pre.ref_visc);
            pre.den_vis_coef.dvd(&pre.ref_den);
            pre.den_vis_coef.dvd(&pre.ref_temp);
            
            pre.den_pres_coef.set_val_dfd0(&pre.bulk_mod);
            pre.den_pres_coef.dvd(&pre.ref_den);
            
            pre.temp_pres_coef.set_val(3.0);
            pre.temp_pres_coef.mult(&pre.bulk_mod);
            pre.temp_pres_coef.mult(&pre.expansion);
        }
        else {
            pre.den_vis_coef.set_val(0.0);
            pre.den_pres_coef.set_val(0.0);
            pre.temp_pres_coef.set_val(0.0);
        }
        
        
        return;
    }

    pub fn get_volume_dfd0(& self, vol : &mut DiffDoub0, pre : &mut DiffDoub0StressPrereq, layer : usize, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut thk = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut s_tmp : [f64; 3] = [0f64; 3];
        let mut this_dv : &DesignVariable;
        
        if self.this_type == 2 {
            thk.set_val(sec_ar[self.sect_ptr].area);//rem: update to factor in dvars;
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    thk.add(& dv_val);
                }
            }
        } else if self.this_type == 3 || self.this_type == 41 {
            thk.set_val_dfd0(& pre.layer_thk[layer]);
        } else {
            thk.set_val(1.0);
        }
        
        vol.set_val(0.0);
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut s_tmp, &self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  s_tmp);
            tmp.set_val(self.ip_wt[i1]);
            tmp.mult(& det_j);
            tmp.mult(& thk);
            vol.add(& tmp);
        }
        
        return;
    }

    pub fn get_section_def_dfd0(&mut self, sec_def : &mut [DiffDoub0], glob_disp : &mut Vec<DiffDoub0>, inst_ori_mat : &mut Vec<DiffDoub0>, loc_ori : &mut Vec<DiffDoub0>, x_glob : &mut Vec<DiffDoub0>, d_ndx : &mut [DiffDoub0], n_vec : &mut [DiffDoub0], n_lgeom : bool, dv1 : usize, dv2 : usize) {
        let i1 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let dof : usize;
        let nd2 : usize;
        let dof2 : usize;
        let mut inst_disp = [DiffDoub0::new(); 60];
        let mut ux = [DiffDoub0::new(); 9];
        let mut rx = [DiffDoub0::new(); 9];
        let mut rot = [DiffDoub0::new(); 3];
        let mut tmp = DiffDoub0::new();
        
        if self.dof_per_nd != 6 {
            for i1 in 0..9 {
                sec_def[i1].set_val(0.0);
            }
            return;
        }
        
        self.get_inst_disp_dfd0(&mut inst_disp, glob_disp, inst_ori_mat, loc_ori, x_glob,  n_lgeom,  dv1,  dv2);
        
        if dv1 == MAX_INT && dv2 == MAX_INT {
            mat_mul_ar_dfd0(&mut ux, &mut inst_disp, d_ndx, 3, self.n_dim, 3);
            i1 = 3*self.n_dim;
            mat_mul_ar_dfd0(&mut rx, &mut inst_disp[i1..], d_ndx, 3, self.n_dim, 3);
            mat_mul_ar_dfd0(&mut rot, &mut inst_disp[i1..], n_vec, 3, self.n_dim, 1);
        } else if (dv1 + dv2) >= MAX_INT {
            let dnz = dv1 + dv2 - MAX_INT;
            //dv1 = dv1 + dv2 - MAX_INT;
            nd = self.dof_table[2*dnz];
            dof = self.dof_table[2*dnz+1];
            if dof < 3 {
                i3 = 0;
                i4 = nd;
                for i1 in 0..3 {
                    i5 = 3*nd;
                    for _i2 in 0..3 {
                        ux[i3].set_val_dfd0(& inst_disp[i4]);
                        ux[i3].mult(& d_ndx[i5]);
                        rx[i3].set_val(0.0);
                        i3 += 1usize;
                        i5 += 1usize;
                    }
                    rot[i1].set_val(0.0);
                    i4 +=  self.n_dim;
                }
            } else {
                i4 = 0;
                for i1 in 0..3 {
                    for i2 in 0..3 {
                        i5 = self.n_dim*i1;
                        i6 = i2;
                        i7 = self.n_dim*(i1+3);
                        ux[i4].set_val(0.0);
                        rx[i4].set_val(0.0);
                        for _i3 in 0..self.num_nds {
                            tmp.set_val_dfd0(& inst_disp[i5]);
                            tmp.mult(& d_ndx[i6]);
                            ux[i4].add(& tmp);
                            tmp.set_val_dfd0(& inst_disp[i7]);
                            tmp.mult(& d_ndx[i6]);
                            rx[i4].add(& tmp);
                            i5 += 1usize;
                            i6 +=  3;
                            i7 += 1usize;
                        }
                        i4 += 1usize;
                    }
                }
                for i1 in 0..3 {
                    rot[i1].set_val(0.0);
                    i3 = self.n_dim*(i1+3);
                    for i2 in 0..self.num_nds {
                        tmp.set_val_dfd0(& inst_disp[i3]);
                        tmp.mult(& n_vec[i2]);
                        rot[i1].add(& tmp);
                        i3 += 1usize;
                    }
                }
            }
        } else {
            nd = self.dof_table[2*dv1];
            dof = self.dof_table[2*dv1+1];
            nd2 = self.dof_table[2*dv2];
            dof2 = self.dof_table[2*dv2+1];
            if dof < 3 && dof2 < 3 {
                for i1 in 0..9 {
                    sec_def[i1].set_val(0.0);
                }
                return;
            } else if dof > 2 && dof2 > 2 {
                i4 = 0;
                for i1 in 0..3 {
                    for i2 in 0..3 {
                        i5 = self.n_dim*i1;
                        i6 = i2;
                        i7 = self.n_dim*(i1+3);
                        ux[i4].set_val(0.0);
                        for _i3 in 0..self.num_nds {
                            tmp.set_val_dfd0(& inst_disp[i5]);
                            tmp.mult(& d_ndx[i6]);
                            ux[i4].add(& tmp);
                            tmp.set_val_dfd0(& inst_disp[i7]);
                            tmp.mult(& d_ndx[i6]);
                            rx[i4].add(& tmp);
                            i5 += 1usize;
                            i6 +=  3;
                            i7 += 1usize;
                        }
                        i4 += 1usize;
                    }
                }
                for i1 in 0..3 {
                    rot[i1].set_val(0.0);
                    i3 = self.n_dim*(i1+3);
                    for i2 in 0..self.num_nds {
                        tmp.set_val_dfd0(& inst_disp[i3]);
                        tmp.mult(& n_vec[i2]);
                        rot[i1].add(& tmp);
                        i3 += 1usize;
                    }
                }
            } else {
                if dof > dof2 {
                    //i1 = dof;
                    //dof = dof2;
                    //dof2 = i1;
                    //i1 = nd;
                    nd = nd2;
                    //nd2 = i1;
                }
                i3 = 0;
                i4 = nd;
                for i1 in 0..3 {
                    i5 = 3*nd;
                    for _i2 in 0..3 {
                        ux[i3].set_val_dfd0(& inst_disp[i4]);
                        ux[i3].mult(& d_ndx[i5]);
                        rx[i3].set_val(0.0);
                        i3 += 1usize;
                        i5 += 1usize;
                    }
                    rot[i1].set_val(0.0);
                    i4 +=  self.n_dim;
                }
            }
        }
        
        if self.this_type == 2 {
            sec_def[0].set_val_dfd0(& ux[0]);
            sec_def[1].set_val_dfd0(& ux[3]);
            sec_def[1].sub(& rot[2]);
            sec_def[2].set_val_dfd0(& ux[6]);
            sec_def[2].add(& rot[1]);
            sec_def[3].set_val_dfd0(& rx[0]);
            sec_def[4].set_val_dfd0(& rx[3]);
            sec_def[5].set_val_dfd0(& rx[6]);
        } else {
            sec_def[0].set_val_dfd0(& ux[0]);
            sec_def[1].set_val_dfd0(& ux[4]);
            sec_def[2].set_val_dfd0(& ux[1]);
            sec_def[2].add(& ux[3]);
            sec_def[3].set_val_dfd0(& rx[3]);
            sec_def[4].set_val_dfd0(& rx[1]);
            sec_def[4].neg();
            sec_def[5].set_val_dfd0(& rx[4]);
            sec_def[5].sub(& rx[0]);
            sec_def[6].set_val_dfd0(& ux[6]);
            sec_def[6].add(& rot[1]);
            sec_def[7].set_val_dfd0(& ux[7]);
            sec_def[7].sub(& rot[0]);
            sec_def[8].set_val(2.0);
            sec_def[8].mult(& rot[2]);
            sec_def[8].sub(& ux[3]);
            sec_def[8].add(& ux[1]);
        }
        
        return;
    }

    pub fn get_solid_strain_dfd0(&mut self, strain : &mut [DiffDoub0], ux : &mut [DiffDoub0], d_ndx : &mut [DiffDoub0], loc_ori : &mut Vec<DiffDoub0>, dv1 : usize, dv2 : usize, n_lgeom : bool) {
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let nd : usize;
        let dof : usize;
        let nd2 : usize;
        let dof2 : usize;
        let mut ux_l = [DiffDoub0::new(); 9];
        let mut strn_mat = [DiffDoub0::new(); 9];
        let mut tmp_ori = [DiffDoub0::new(); 9];
        let mut tmp = DiffDoub0::new();
        
        for i1 in 0..9 {
            strn_mat[i1].set_val(0.0);
        }
        if dv1 == MAX_INT && dv2 == MAX_INT {
            vec_to_ar_dfd0(&mut tmp_ori, loc_ori,  0,  9);
            mat_mul_ar_dfd0(&mut ux_l, &mut tmp_ori, ux, 3, 3, 3);
            for i1 in 0..3 {
                i4 = 4*i1;
                i5 = 4*i1;
                for i2 in i1..3 {
                    strn_mat[i4].add(& ux_l[i4]);
                    strn_mat[i4].add(& ux_l[i5]);
                    if n_lgeom {
                        i6 = i1;
                        i7 = i2;
                        for _i3 in 0..3 {
                            tmp.set_val_dfd0(& ux[i6]);
                            tmp.mult(& ux[i7]);
                            strn_mat[i4].add(& tmp);
                            i6 +=  3;
                            i7 +=  3;
                        }
                    }
                    i4 += 1usize;
                    i5 +=  3;
                }
            }
        } else if (dv1 + dv2) >= MAX_INT {
            if dv1 < MAX_INT {
                nd = self.dof_table[2*dv1];
                dof = self.dof_table[2*dv1+1];
            } else {
                nd = self.dof_table[2*dv2];
                dof = self.dof_table[2*dv2+1];
            }
            for i1 in 0..3 {
                i4 = 4*i1;
                for i2 in i1..3 {
                    i5 = 3*i1 + dof;
                    i6 = 3*nd + i2;
                    tmp.set_val_dfd0(& loc_ori[i5]);
                    tmp.mult(& d_ndx[i6]);
                    strn_mat[i4].add(& tmp);
                    i5 = 3*i2 + dof;
                    i6 = 3*nd + i1;
                    tmp.set_val_dfd0(& loc_ori[i5]);
                    tmp.mult(& d_ndx[i6]);
                    strn_mat[i4].add(& tmp);
                    if n_lgeom {
                        i5 = 3*nd + i1;
                        i6 = 3*dof + i2;
                        tmp.set_val_dfd0(& d_ndx[i5]);
                        tmp.mult(& ux[i6]);
                        strn_mat[i4].add(& tmp);
                        i5 = 3*nd + i2;
                        i6 = 3*dof + i1;
                        tmp.set_val_dfd0(& d_ndx[i5]);
                        tmp.mult(& ux[i6]);
                        strn_mat[i4].add(& tmp);
                    }
                    i4 += 1usize;
                }
            }
        } else {
            if n_lgeom {
                nd = self.dof_table[2*dv1];
                dof = self.dof_table[2*dv1+1];
                nd2 = self.dof_table[2*dv2];
                dof2 = self.dof_table[2*dv2+1];
                if dof == dof2 {
                    for i1 in 0..3 {
                        i4 = 4*i1;
                        for i2 in i1..3 {
                            i5 = 3*nd + i1;
                            i6 = 3*nd2 + i2;
                            tmp.set_val_dfd0(& d_ndx[i5]);
                            tmp.mult(& d_ndx[i6]);
                            strn_mat[i4].add(& tmp);
                            i5 = 3*nd2 + i1;
                            i6 = 3*nd + i2;
                            tmp.set_val_dfd0(& d_ndx[i5]);
                            tmp.mult(& d_ndx[i6]);
                            strn_mat[i4].add(& tmp);
                            i4 += 1usize;
                        }
                    }
                }
            }
        }
        
        tmp.set_val(0.5);
        strain[0].set_val_dfd0(& strn_mat[0]);
        strain[0].mult(& tmp);
        strain[1].set_val_dfd0(& strn_mat[4]);
        strain[1].mult(& tmp);
        strain[2].set_val_dfd0(& strn_mat[8]);
        strain[2].mult(& tmp);
        strain[3].set_val_dfd0(& strn_mat[1]);
        strain[4].set_val_dfd0(& strn_mat[2]);
        strain[5].set_val_dfd0(& strn_mat[5]);
        
        return;
    }

    pub fn get_stress_strain_dfd0(&mut self, stress : &mut [DiffDoub0], strain : &mut [DiffDoub0], t_strain : &mut [DiffDoub0], spt : &mut [f64], layer : usize, n_lgeom : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i2 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ux = [DiffDoub0::new(); 9];
        let mut sec_def = [DiffDoub0::new(); 9];
        let mut sect_strn = [DiffDoub0::new(); 3];
        let mut adj_stn = [DiffDoub0::new(); 6];
        let mut tmp = DiffDoub0::new();
        let mut ip_temp = DiffDoub0::new();
        let mut tmp_ar = [DiffDoub0::new(); 60];

        for i in 0..6 {
            stress[i].set_val(0f64);
            strain[i].set_val(0f64);
            t_strain[i].set_val(0f64);
        }
        
        self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        
        ip_temp.set_val(0.0);
        for i1 in 0..self.num_nds {
            tmp.set_val_dfd0(& pre.glob_temp[i1]);
            tmp.mult(& n_vec[i1]);
            ip_temp.add(& tmp);
        }
        
        if self.this_type == 41 || self.this_type == 3 {
            if n_lgeom {
                self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            
            self.get_section_def_dfd0(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  MAX_INT,  MAX_INT);
            
            sect_strn[0].set_val_dfd0(& sec_def[0]);
            tmp.set_val_dfd0(& pre.layer_z[layer]);
            tmp.mult(& sec_def[3]);
            sect_strn[0].add(& tmp);
            
            sect_strn[1].set_val_dfd0(& sec_def[1]);
            tmp.set_val_dfd0(& pre.layer_z[layer]);
            tmp.mult(& sec_def[4]);
            sect_strn[1].sub(& tmp);
            
            sect_strn[2].set_val_dfd0(& sec_def[2]);
            tmp.set_val_dfd0(& pre.layer_z[layer]);
            tmp.mult(& sec_def[5]);
            sect_strn[2].add(& tmp);
            
            tmp.set_val_dfd0(& pre.layer_ang[layer]);
            tmp.neg();
            self.transform_strain_dfd0(strain, &mut  sect_strn, &mut  tmp);
            i2 = 3 * layer;
            for i1 in 0..3 {
                tmp.set_val_dfd0(& pre.layer_te[i2]);
                tmp.mult(& ip_temp);
                t_strain[i1].set_val_dfd0(&tmp);
                adj_stn[i1].set_val_dfd0(& strain[i1]);
                adj_stn[i1].sub(& tmp);
                adj_stn[i1].sub(& pre.layer_e0[i2]);
                i2 += 1usize;
            }
            mat_mul_ar_dfd0(stress, &mut  pre.layer_q[(9*layer)..], &mut  adj_stn,  3,  3,  1);
            
            tmp.set_val_dfd0(& strain[2]);
            strain[3].set_val_dfd0(& tmp);
            strain[2].set_val(0.0);

            tmp.set_val_dfd0(& t_strain[2]);
            t_strain[3].set_val_dfd0(&tmp);
            t_strain[2].set_val(0.0);
            
            tmp.set_val_dfd0(& stress[2]);
            stress[3].set_val_dfd0(& tmp);
            stress[2].set_val(0.0);
            
        } else if self.this_type != 2 {
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.glob_disp,  0,  60);
            mat_mul_ar_dfd0(&mut ux, &mut  tmp_ar, &mut  d_ndx,  3,  self.n_dim,  3);
            self.get_solid_strain_dfd0(strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  MAX_INT,  MAX_INT,  n_lgeom);
            for i1 in 0..6 {
                tmp.set_val_dfd0(& pre.therm_exp[i1]);
                tmp.mult(& ip_temp);
                t_strain[i1].set_val_dfd0(&tmp);
                adj_stn[i1].set_val_dfd0(& strain[i1]);
                adj_stn[i1].sub(& tmp);
                adj_stn[i1].sub(& pre.einit[i1]);
            }
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.cmat,  0,  36);
            mat_mul_ar_dfd0(stress, &mut tmp_ar, &mut adj_stn,  6,  6,  1);
        }
        return;
    }

    pub fn d_stress_straind_u_dfd0(&mut self, dsd_u : &mut Vec<DiffDoub0>, ded_u : &mut Vec<DiffDoub0>, dsd_t : &mut Vec<DiffDoub0>, spt : &mut [f64], layer : usize, n_lgeom : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ux = [DiffDoub0::new(); 9];
        let mut sec_def = [DiffDoub0::new(); 9];
        let mut sect_strn = [DiffDoub0::new(); 3];
        let mut d_strain = [DiffDoub0::new(); 6];
        let mut d_stress = [DiffDoub0::new(); 6];
        let mut cte = [DiffDoub0::new(); 6];
        let mut cten = [DiffDoub0::new(); 60];
        let mut tmp = DiffDoub0::new();
        let mut tmp_ar = [DiffDoub0::new(); 60];
        let mut tmp_ar2 = [DiffDoub0::new(); 6];
        
        let tot_dof : usize =  self.num_nds * self.dof_per_nd + self.num_int_dof;
        i2 = 6 * tot_dof;
        for i1 in 0..i2 {
            dsd_u[i1].set_val(0.0);
            ded_u[i1].set_val(0.0);
        }
        i2 = 6 * self.num_nds;
        for i1 in 0..i2 {
            dsd_t[i1].set_val(0.0);
        }
        
        if self.this_type == 41 || self.this_type == 3 {
            if n_lgeom {
                self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
            for i1 in 0..tot_dof {
                self.get_section_def_dfd0(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  i1,  MAX_INT);
                
                sect_strn[0].set_val_dfd0(& sec_def[0]);
                tmp.set_val_dfd0(& pre.layer_z[layer]);
                tmp.mult(& sec_def[3]);
                sect_strn[0].add(& tmp);
                
                sect_strn[1].set_val_dfd0(& sec_def[1]);
                tmp.set_val_dfd0(& pre.layer_z[layer]);
                tmp.mult(& sec_def[4]);
                sect_strn[1].sub(& tmp);
                
                sect_strn[2].set_val_dfd0(& sec_def[2]);
                tmp.set_val_dfd0(& pre.layer_z[layer]);
                tmp.mult(& sec_def[5]);
                sect_strn[2].add(& tmp);
                
                tmp.set_val_dfd0(& pre.layer_ang[layer]);
                tmp.neg();
                self.transform_strain_dfd0(&mut d_strain, &mut  sect_strn, &mut  tmp);
                mat_mul_ar_dfd0(&mut d_stress, &mut pre.layer_q[(9*layer)..], &mut  d_strain,  3,  3,  1);
                
                ded_u[i1].set_val_dfd0(& d_strain[0]);
                ded_u[i1 + tot_dof].set_val_dfd0(& d_strain[1]);
                ded_u[i1 + 3 * tot_dof].set_val_dfd0(& d_strain[2]);
                
                dsd_u[i1].set_val_dfd0(& d_stress[0]);
                dsd_u[tot_dof + i1].set_val_dfd0(& d_stress[1]);
                dsd_u[3 * tot_dof + i1].set_val_dfd0(& d_stress[2]);
            }
            mat_mul_ar_dfd0(&mut cte, &mut pre.layer_q[(9*layer)..], &mut pre.layer_te[(3*layer)..],  3,  3,  1);
            cte[0].neg();
            cte[1].neg();
            cte[2].neg();
            mat_mul_ar_dfd0(&mut cten, &mut  cte, &mut  n_vec,  3,  1,  self.num_nds);
            for i1 in 0..self.num_nds {
                dsd_t[i1].set_val_dfd0(& cten[i1]);
                dsd_t[i1 + self.num_nds].set_val_dfd0(& cten[i1 + self.num_nds]);
                dsd_t[i1 + 3 * self.num_nds].set_val_dfd0(& cten[i1 + 2 * self.num_nds]);
            }
        }
        else if self.this_type != 2 {
            self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.glob_disp,  0,  60);
            mat_mul_ar_dfd0(&mut ux, &mut  tmp_ar, &mut  d_ndx,  3,  self.n_dim,  3);
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.cmat,  0,  36);
            for i1 in 0..tot_dof {
                self.get_solid_strain_dfd0(&mut d_strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  i1,  MAX_INT,  n_lgeom);
                mat_mul_ar_dfd0(&mut d_stress, &mut  tmp_ar, &mut  d_strain,  6,  6,  1);
                i3 = i1;
                for i2 in 0..6 {
                    ded_u[i3].set_val_dfd0(& d_strain[i2]);
                    dsd_u[i3].set_val_dfd0(& d_stress[i2]);
                    i3  +=  tot_dof;
                }
            }
            vec_to_ar_dfd0(&mut tmp_ar2, &mut  pre.therm_exp,  0,  6);
            mat_mul_ar_dfd0(&mut cte, &mut  tmp_ar, &mut  tmp_ar2,  6,  6,  1);
            for i1 in 0..6 {
                cte[i1].neg();
            }
            mat_mul_ar_dfd0(&mut tmp_ar, &mut  cte, &mut  n_vec,  6,  1,  self.num_nds);
            ar_to_vec_dfd0(&mut tmp_ar, dsd_t,  0,  60);
        }
        
        return;
    }

    pub fn get_def_frc_mom_dfd0(&mut self, def : &mut [DiffDoub0], frc_mom : &mut [DiffDoub0], spt : &mut [f64], n_lgeom : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut pt_temp = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut tmp_ar = [DiffDoub0::new(); 36];
        
        if self.dof_per_nd != 6 {
            for i1 in 0..9 {
                def[i1].set_val(0.0);
            }
            return;
        }
        
        if n_lgeom {
            self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
        }
        
        self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        self.get_section_def_dfd0(def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  MAX_INT,  MAX_INT);
        
        pt_temp.set_val(0.0);
        for i1 in 0..self.num_nds {
            tmp.set_val_dfd0(& pre.glob_temp[i1]);
            tmp.mult(& n_vec[i1]);
            pt_temp.add(& tmp);
        }
        
        vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.cmat,  0,  36);
        mat_mul_ar_dfd0(frc_mom, &mut  tmp_ar, def,  self.def_dim,  self.def_dim,  1);
        
        for i1 in 0..6 {
            frc_mom[i1].sub(& pre.einit[i1]);
            tmp.set_val_dfd0(& pt_temp);
            tmp.mult(& pre.therm_exp[i1]);
            frc_mom[i1].sub(& tmp);
        }
        
        return;
    }

    pub fn d_def_frc_momd_u_dfd0(&mut self, d_defd_u : &mut Vec<DiffDoub0>, d_frc_momd_u : &mut Vec<DiffDoub0>, d_frc_momd_t : &mut Vec<DiffDoub0>, spt : &mut [f64], n_lgeom : bool, pre : &mut DiffDoub0StressPrereq) {
        let mut i2 : usize;
        let tot_dof : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut def = [DiffDoub0::new(); 9];
        let mut tmp_ar = [DiffDoub0::new(); 60];
        let mut tmp_ar2 = [DiffDoub0::new(); 6];
        
        if n_lgeom {
            self.get_inst_ori_dfd0(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
        }
        
        self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        tot_dof = self.num_nds * self.dof_per_nd + self.num_int_dof;
        for i1 in 0..tot_dof {
            self.get_section_def_dfd0(&mut def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  i1,  MAX_INT);
            i2 = i1;
            for i3 in 0..self.def_dim {
                d_defd_u[i2].set_val_dfd0(& def[i3]);
                i2  +=  tot_dof;
            }
        }
        
        mat_mul_dfd0(d_frc_momd_u, &mut  pre.cmat, d_defd_u,  self.def_dim,  self.def_dim,  tot_dof);
        
        vec_to_ar_dfd0(&mut tmp_ar2, &mut  pre.therm_exp,  0,  6);
        mat_mul_ar_dfd0(&mut tmp_ar, &mut  tmp_ar2, &mut  n_vec,  6,  1,  self.num_nds);
        ar_to_vec_dfd0(&mut tmp_ar, d_frc_momd_t,  0,  60);
        
        return;
    }

    pub fn get_flux_tgrad_dfd0(&mut self, flux : &mut [DiffDoub0], t_grad : &mut [DiffDoub0], spt : &mut [f64], layer : usize, pre : &mut DiffDoub0StressPrereq) {
        let i1 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut tmp_ar = [DiffDoub0::new(); 10];
        
        self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.glob_temp,  0,  10);
        mat_mul_ar_dfd0(t_grad, &mut  tmp_ar, &mut  d_ndx,  1,  self.num_nds,  3);
        if self.this_type == 41 || self.this_type == 3 {
            i1 = 9 * layer;
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.layer_tc,  i1,  i1 + 9);
            mat_mul_ar_dfd0(flux, &mut  tmp_ar, t_grad,  3,  3,  1);
        }
        else {
            vec_to_ar_dfd0(&mut tmp_ar, &mut  pre.tcmat,  0,  9);
            mat_mul_ar_dfd0(flux, &mut  tmp_ar, t_grad,  3,  3,  1);
        }
        flux[0].neg();
        flux[1].neg();
        flux[2].neg();
        
        return;
    }

    pub fn d_flux_tgradd_t_dfd0(&mut self, d_fd_t : &mut Vec<DiffDoub0>, d_tg : &mut Vec<DiffDoub0>, spt : &mut [f64], layer : usize, pre : &mut DiffDoub0StressPrereq) {
        let i1 : usize;
        let i2 : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut tmp_ar1 = [DiffDoub0::new(); 33];
        let mut tmp_ar2 = [DiffDoub0::new(); 9];
        let mut tmp_ar3 = [DiffDoub0::new(); 30];
        
        self.get_ip_data_dfd0(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        transpose_ar_dfd0(&mut tmp_ar1, &mut  d_ndx,  self.num_nds,  3);
        ar_to_vec_dfd0(&mut tmp_ar1, d_tg,  0,  33);
        if self.this_type == 41 || self.this_type == 3 {
            i1 = 9 * layer;
            vec_to_ar_dfd0(&mut tmp_ar2, &mut  pre.layer_tc,  i1,  i1 + 9);
            mat_mul_ar_dfd0(&mut tmp_ar3, &mut  tmp_ar2, &mut  tmp_ar1,  3,  3,  self.num_nds);
            ar_to_vec_dfd0(&mut tmp_ar3, d_fd_t,  0,  30);
        }
        else {
            mat_mul_dfd0(d_fd_t, &mut  pre.tcmat, d_tg,  3,  3,  self.num_nds);
        }
        i2 = 3 * self.num_nds;
        for i1 in 0..i2 {
            d_fd_t[i1].neg();
        }
        
        return;
    }

    pub fn get_cent_v_grad_dfd0(&mut self, v_grad : &mut [DiffDoub0], t_grad : &mut [DiffDoub0], pre : &mut DiffDoub0FlPrereq, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable> ) {
        // get fluid prerequisites before calling
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut spt = [0f64; 3];
       
        for i in 0..3 {
            spt[i] = self.s_cent[i];
        }

        // self.get_nd_crds_dfd0(&mut pre.glob_nds, nd_ar, dv_ar);
        // self.get_nd_disp_dfd0(&mut pre.glob_disp, nd_ar);
        // let i2 = 3 * self.num_nds;
        // for i1 in 0..i2 {
        //     pre.glob_nds[i1].add(& pre.glob_disp[i1]);
        // }
        // self.get_nd_fl_vel_dfd0(&mut pre.fl_vel, nd_ar);
        // self.get_nd_temp_dfd0(&mut pre.fl_temp, nd_ar);

        self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &mut spt);
        mat_mul_ar_dfd0(v_grad, &pre.fl_vel, &d_ndx, 3, self.num_nds, 3);
        mat_mul_ar_dfd0(t_grad, &pre.fl_temp, &d_ndx, 1, self.num_nds, 3);
        
    }

    pub fn put_vec_to_glob_mat_dfd0(&mut self, q_mat : &mut SparseMat, el_qvec : &mut Vec<DiffDoub0>, for_therm : bool, mat_row : usize, nd_ar : &mut Vec<Node>) {
        let nd_dof : usize =  self.num_nds*self.dof_per_nd;
        let tot_dof : usize =  nd_dof + self.num_int_dof;
        let mut nd : usize;
        let mut dof : usize;
        let mut glob_ind : usize;
        
        if for_therm {
            for i1 in 0..self.num_nds {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                q_mat.add_entry(mat_row,   glob_ind,   el_qvec[i1].val);
            }
        }
        else {
            for i1 in 0..tot_dof {
                if i1 < nd_dof {
                    nd = self.nodes[self.dof_table[2 * i1]];
                    dof = self.dof_table[2 * i1 + 1];
                    glob_ind = nd_ar[nd].dof_index[dof];
                    q_mat.add_entry(mat_row,  glob_ind,  el_qvec[i1].val);
                }
                else {
                    dof = i1 - nd_dof;
                    glob_ind = self.int_dof_index + dof;
                    q_mat.add_entry(mat_row,   glob_ind,   el_qvec[i1].val);
                }
            }
        }
        
        return;
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn get_nd_disp_dfd1(&self, glob_disp : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_disp[i3].set_val(this_nd.displacement[i2]);
                i3 +=  self.n_dim;
            }
        }
        
        if self.num_int_dof > 0 {
            i2 = 2*self.num_nds*self.dof_per_nd;
            for i3 in 0..self.num_int_dof {
                i4 = self.dof_table[i2];
                i5 = self.dof_table[i2+1];
                i6 = self.n_dim*i5 + i4;
                glob_disp[i6].set_val(self.internal_disp[i3]);
                i2 +=  2;
            }
        }
        
        return;
    }

    pub fn get_nd_vel_dfd1(& self, glob_vel : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_vel[i3].set_val(this_nd.velocity[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_acc_dfd1(& self, glob_acc : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..self.dof_per_nd {
                glob_acc[i3].set_val(this_nd.acceleration[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_fl_vel_dfd1(& self, fl_vel : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..3 {
                fl_vel[i3].set_val(this_nd.fl_vel[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_fl_vdot_dfd1(& self, fl_vdot : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut i3 : usize;
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            i3 = i1;
            for i2 in 0..3 {
                fl_vdot[i3].set_val(this_nd.fl_vel_dot[i2]);
                i3  +=  self.num_nds;
            }
        }
        
        return;
    }

    pub fn get_nd_temp_dfd1(& self, glob_temp : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            glob_temp[i1].set_val(this_nd.temperature);
        }
        return;
    }

    pub fn get_nd_tdot_dfd1(& self, glob_tdot : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            glob_tdot[i1].set_val(this_nd.temp_change_rate);
        }
        return;
    }

    pub fn get_nd_fl_den_dfd1(& self, fl_den : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            fl_den[i1].set_val(this_nd.fl_den);
        }
        
        return;
    }

    pub fn get_nd_fl_den_dot_dfd1(& self, fl_den_dot : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            fl_den_dot[i1].set_val(this_nd.fl_den_dot);
        }
        
        return;
    }

    pub fn get_nd_turb_e_dfd1(& self, turb_e : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            turb_e[i1].set_val(this_nd.turb_e);
        }
        
        return;
    }

    pub fn get_nd_turb_edot_dfd1(& self, turb_edot : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>) {
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            turb_edot[i1].set_val(this_nd.turb_edot);
        }
        
        return;
    }

    pub fn eval_n_dfd1(&self, n_vec : &mut [DiffDoub1], d_nds : &mut [DiffDoub1], spt : & [f64]) {
        if self.this_type == 4 || self.this_type == 400 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]-spt[2]);
            d_nds[0].set_val(-1.0);
            d_nds[1].set_val(-1.0);
            d_nds[2].set_val(-1.0);
            
            n_vec[1].set_val(spt[0]);
            d_nds[3].set_val(1.0);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(spt[1]);
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(1.0);
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(spt[2]);
            d_nds[9].set_val(0.0);
            d_nds[10].set_val(0.0);
            d_nds[11].set_val(1.0);
        } else if self.this_type == 6 || self.this_type == 600 {
            n_vec[0].set_val(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
            d_nds[0].set_val(-0.5*(1.0-spt[2]));
            d_nds[1].set_val(-0.5*(1.0-spt[2]));
            d_nds[2].set_val(-0.5*(1.0-spt[0]-spt[1]));
            
            n_vec[1].set_val(0.5*spt[0]*(1.0-spt[2]));
            d_nds[3].set_val(0.5*(1.0-spt[2]));
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(-0.5*spt[0]);
            
            n_vec[2].set_val(0.5*spt[1]*(1.0-spt[2]));
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(0.5*(1.0-spt[2]));
            d_nds[8].set_val(-0.5*spt[1]);
            
            n_vec[3].set_val(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
            d_nds[9].set_val(-0.5*(1.0+spt[2]));
            d_nds[10].set_val(-0.5*(1.0+spt[2]));
            d_nds[11].set_val(0.5*(1.0-spt[0]-spt[1]));
            
            n_vec[4].set_val(0.5*spt[0]*(1.0+spt[2]));
            d_nds[12].set_val(0.5*(1.0+spt[2]));
            d_nds[13].set_val(0.0);
            d_nds[14].set_val(0.5*spt[0]);
            
            n_vec[5].set_val(0.5*spt[1]*(1.0+spt[2]));
            d_nds[15].set_val(0.0);
            d_nds[16].set_val(0.5*(1.0+spt[2]));
            d_nds[17].set_val(0.5*spt[1]);
        } else if self.this_type == 8 || self.this_type == 81 || self.this_type == 800 {
            n_vec[0].set_val(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
            d_nds[0].set_val(-0.125*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[1].set_val(-0.125*(1.0-spt[0])*(1.0-spt[2]));
            d_nds[2].set_val(-0.125*(1.0-spt[0])*(1.0-spt[1]));
            
            n_vec[1].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[3].set_val(0.125*(1.0-spt[1])*(1.0-spt[2]));
            d_nds[4].set_val(-0.125*(1.0+spt[0])*(1.0-spt[2]));
            d_nds[5].set_val(-0.125*(1.0+spt[0])*(1.0-spt[1]));
            
            n_vec[2].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[6].set_val(0.125*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[7].set_val(0.125*(1.0+spt[0])*(1.0-spt[2]));
            d_nds[8].set_val(-0.125*(1.0+spt[0])*(1.0+spt[1]));
            
            n_vec[3].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[9].set_val(-0.125*(1.0+spt[1])*(1.0-spt[2]));
            d_nds[10].set_val(0.125*(1.0-spt[0])*(1.0-spt[2]));
            d_nds[11].set_val(-0.125*(1.0-spt[0])*(1.0+spt[1]));
            
            n_vec[4].set_val(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[12].set_val(-0.125*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[13].set_val(-0.125*(1.0-spt[0])*(1.0+spt[2]));
            d_nds[14].set_val(0.125*(1.0-spt[0])*(1.0-spt[1]));
            
            n_vec[5].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[15].set_val(0.125*(1.0-spt[1])*(1.0+spt[2]));
            d_nds[16].set_val(-0.125*(1.0+spt[0])*(1.0+spt[2]));
            d_nds[17].set_val(0.125*(1.0+spt[0])*(1.0-spt[1]));
            
            n_vec[6].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[18].set_val(0.125*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[19].set_val(0.125*(1.0+spt[0])*(1.0+spt[2]));
            d_nds[20].set_val(0.125*(1.0+spt[0])*(1.0+spt[1]));
            
            n_vec[7].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[21].set_val(-0.125*(1.0+spt[1])*(1.0+spt[2]));
            d_nds[22].set_val(0.125*(1.0-spt[0])*(1.0+spt[2]));
            d_nds[23].set_val(0.125*(1.0-spt[0])*(1.0+spt[1]));
            
            if self.this_type == 81 {
                n_vec[8].set_val(1.0-spt[0]*spt[0]);
                d_nds[24].set_val(-2.0*spt[0]);
                d_nds[25].set_val(0.0);
                d_nds[26].set_val(0.0);
                
                n_vec[9].set_val(1.0-spt[1]*spt[1]);
                d_nds[27].set_val(0.0);
                d_nds[28].set_val(-2.0*spt[1]);
                d_nds[29].set_val(0.0);
                
                n_vec[10].set_val(1.0-spt[2]*spt[2]);
                d_nds[30].set_val(0.0);
                d_nds[31].set_val(0.0);
                d_nds[32].set_val(-2.0*spt[2]);
            }
        }
        else if self.this_type == 10 || self.this_type == 1000 {
            let p1 : f64 =  1.0 - spt[0] - spt[1] - spt[2];
            let p2 : f64 =  p1 - 0.5;
            n_vec[0].set_val(2.0 * p1 * p2);
            d_nds[0].set_val(2.0 * (-p2 - p1));
            d_nds[1].set_val(2.0 * (-p2 - p1));
            d_nds[2].set_val(2.0 * (-p2 - p1));
            
            n_vec[1].set_val(-2.0 * spt[0] * (0.5 - spt[0]));
            d_nds[3].set_val(-2.0 * (0.5 - spt[0] - spt[0]));
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(-2.0 * spt[1] * (0.5 - spt[1]));
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(-2.0 * (0.5 - spt[1] - spt[1]));
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(-2.0 * spt[2] * (0.5 - spt[2]));
            d_nds[9].set_val(0.0);
            d_nds[10].set_val(0.0);
            d_nds[11].set_val(-2.0 * (0.5 - spt[2] - spt[2]));
            
            n_vec[4].set_val(4.0 * spt[0] * p1);
            d_nds[12].set_val(4.0 * (p1 - spt[0]));
            d_nds[13].set_val(-4.0 * spt[0]);
            d_nds[14].set_val(-4.0 * spt[0]);
            
            n_vec[5].set_val(4.0 * spt[0] * spt[1]);
            d_nds[15].set_val(4.0 * spt[1]);
            d_nds[16].set_val(4.0 * spt[0]);
            d_nds[17].set_val(0.0);
            
            n_vec[6].set_val(4.0 * spt[1] * p1);
            d_nds[18].set_val(-4.0 * spt[1]);
            d_nds[19].set_val(4.0 * (p1 - spt[1]));
            d_nds[20].set_val(-4.0 * spt[1]);
            
            n_vec[7].set_val(4.0 * spt[2] * p1);
            d_nds[21].set_val(-4.0 * spt[2]);
            d_nds[22].set_val(-4.0 * spt[2]);
            d_nds[23].set_val(4.0 * (p1 - spt[2]));
            
            n_vec[8].set_val(4.0 * spt[0] * spt[2]);
            d_nds[24].set_val(4.0 * spt[2]);
            d_nds[25].set_val(0.0);
            d_nds[26].set_val(4.0 * spt[0]);
            
            n_vec[9].set_val(4.0 * spt[1] * spt[2]);
            d_nds[27].set_val(0.0);
            d_nds[28].set_val(4.0 * spt[2]);
            d_nds[29].set_val(4.0 * spt[1]);
        }
        else if self.this_type == 3 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]);
            d_nds[0].set_val(-1.0);
            d_nds[1].set_val(-1.0);
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(spt[0]);
            d_nds[3].set_val(1.0);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(spt[1]);
            d_nds[6].set_val(0.0);
            d_nds[7].set_val(1.0);
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(spt[0]*(1.0-spt[0]-spt[1]));
            d_nds[9].set_val(1.0-spt[0]-spt[1] - spt[0]);
            d_nds[10].set_val(-spt[0]);
            d_nds[11].set_val(0.0);
            
            n_vec[4].set_val(spt[0]*spt[1]);
            d_nds[12].set_val(spt[1]);
            d_nds[13].set_val(spt[0]);
            d_nds[14].set_val(0.0);
            
            n_vec[5].set_val(spt[1]*(1.0-spt[0]-spt[1]));
            d_nds[15].set_val(-spt[1]);
            d_nds[16].set_val(1.0-spt[0]-spt[1]-spt[1]);
            d_nds[17].set_val(0.0);
            
        } else if self.this_type == 41 {
            n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));
            d_nds[0].set_val(-0.25*(1.0-spt[1]));
            d_nds[1].set_val(-0.25*(1.0-spt[0]));
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
            d_nds[3].set_val(0.25*(1.0-spt[1]));
            d_nds[4].set_val(-0.25*(1.0+spt[0]));
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
            d_nds[6].set_val(0.25*(1.0+spt[1]));
            d_nds[7].set_val(0.25*(1.0+spt[0]));
            d_nds[8].set_val(0.0);
            
            n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
            d_nds[9].set_val(-0.25*(1.0+spt[1]));
            d_nds[10].set_val(0.25*(1.0-spt[0]));
            d_nds[11].set_val(0.0);
            
            n_vec[4].set_val(1.0-spt[0]*spt[0]);
            d_nds[12].set_val(-2.0*spt[0]);
            d_nds[13].set_val(0.0);
            d_nds[14].set_val(0.0);
            
            n_vec[5].set_val(1.0-spt[1]*spt[1]);
            d_nds[15].set_val(0.0);
            d_nds[16].set_val(-2.0*spt[1]);
            d_nds[17].set_val(0.0);
            
            n_vec[6].set_val((1.0-spt[0]*spt[0])*(1.0-spt[1]));
            d_nds[18].set_val(-2.0*spt[0]*(1.0-spt[1]));
            d_nds[19].set_val(-1.0+spt[0]*spt[0]);
            d_nds[20].set_val(0.0);
            
            n_vec[7].set_val((1.0-spt[1]*spt[1])*(1.0+spt[0]));
            d_nds[21].set_val(1.0-spt[1]*spt[1]);
            d_nds[22].set_val(-2.0*spt[1]*(1.0+spt[0]));
            d_nds[23].set_val(0.0);
            
            n_vec[8].set_val((1.0-spt[0]*spt[0])*(1.0+spt[1]));
            d_nds[24].set_val(-2.0*spt[0]*(1.0+spt[1]));
            d_nds[25].set_val(1.0-spt[0]*spt[0]);
            d_nds[26].set_val(0.0);
            
            n_vec[9].set_val((1.0-spt[1]*spt[1])*(1.0-spt[0]));
            d_nds[27].set_val(-1.0+spt[1]*spt[1]);
            d_nds[28].set_val(-2.0*spt[1]*(1.0-spt[0]));
            d_nds[29].set_val(0.0);
            
        } else if self.this_type == 2 {
            n_vec[0].set_val(0.5*(1.0 - spt[0]));
            d_nds[0].set_val(-0.5);
            d_nds[1].set_val(0.0);
            d_nds[2].set_val(0.0);
            
            n_vec[1].set_val(0.5*(1.0 + spt[0]));
            d_nds[3].set_val(0.5);
            d_nds[4].set_val(0.0);
            d_nds[5].set_val(0.0);
            
            n_vec[2].set_val(1.0 - spt[0]*spt[0]);
            d_nds[0].set_val(-2.0*spt[0]);
            d_nds[1].set_val(0.0);
            d_nds[2].set_val(0.0);
        }
        return;
    }

    pub fn get_ip_data_dfd1(&self, n_vec : &mut [DiffDoub1], d_ndx : &mut [DiffDoub1], det_j : &mut DiffDoub1, loc_nds : &mut Vec<DiffDoub1>, spt : & [f64]) {
        let i1 : usize;
        let i2 : usize;
        let mut n_cent = [DiffDoub1::new(); 11];
        let mut d_nds = [DiffDoub1::new(); 33];
        let mut d_nds_cent = [DiffDoub1::new(); 33];
        let mut j_mat = [DiffDoub1::new(); 9];
        let mut tmp_nds = [DiffDoub1::new(); 30];
        let mut j_cent = [DiffDoub1::new(); 9];
        let mut det_cent = DiffDoub1::new();
        let mut j_inv = [DiffDoub1::new(); 9];
        let mut j_inv_cent = [DiffDoub1::new(); 9];
        let mut x_vec = [DiffDoub1::new(); 3];
        let mut b_vec = [DiffDoub1::new(); 3];
        let mut z_dir = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        
        self.eval_n_dfd1(n_vec, &mut d_nds, spt);
        vec_to_ar_dfd1(&mut tmp_nds, loc_nds,  0,  30);
        mat_mul_ar_dfd1(&mut j_mat, &mut tmp_nds, &mut d_nds, 3, self.num_nds, 3);
        
        if self.this_type == 41 || self.this_type == 3 {
            z_dir.set_val_dfd1(& j_mat[0]);
            z_dir.mult(& j_mat[4]);
            tmp.set_val_dfd1(& j_mat[3]);
            tmp.mult(& j_mat[1]);
            z_dir.sub(& tmp);
            if z_dir.val > 0.0 {
                j_mat[8].set_val(1.0);
            } else {
                j_mat[8].set_val(-1.0);
            }
        } else if self.this_type == 2 {
            j_mat[4].set_val(1.0);
            if j_mat[0].val > 0.0 {
                j_mat[8].set_val(1.0);
            } else {
                j_mat[8].set_val(-1.0);
            }
        }
        
        get_det_inv_ar_dfd1(det_j, &mut j_inv, &mut j_mat, 3, 0, &mut x_vec, &mut b_vec);
        
        // mat_mul_dfd1(d_nds,d_nds,j_inv,self.n_dim,3,3);
        mat_mul_ar_dfd1(d_ndx, &mut d_nds, &mut j_inv, self.num_nds, 3, 3);
        
        if self.n_dim > self.num_nds {
            self.eval_n_dfd1(&mut n_cent, &mut  d_nds_cent, & self.s_cent);
            mat_mul_ar_dfd1(&mut j_cent, &mut  tmp_nds, &mut  d_nds_cent,  3,  self.num_nds,  3);
            
            if self.this_type == 41 || self.this_type == 3 {
                z_dir.set_val_dfd1(& j_cent[0]);
                z_dir.mult(& j_cent[4]);
                tmp.set_val_dfd1(& j_cent[3]);
                tmp.mult(& j_cent[1]);
                z_dir.sub(& tmp);
                if z_dir.val > 0.0 {
                    j_cent[8].set_val(1.0);
                }
                else {
                    j_cent[8].set_val(-1.0);
                }
            }
            else if self.this_type == 2 {
                j_cent[4].set_val(1.0);
                if j_cent[0].val > 0.0 {
                    j_cent[8].set_val(1.0);
                }
                else {
                    j_cent[8].set_val(-1.0);
                }
            }
            
            get_det_inv_ar_dfd1(&mut det_cent, &mut  j_inv_cent, &mut  j_cent,  3,  0, &mut  x_vec, &mut  b_vec);
            
            i1 = 3 * self.num_nds;
            i2 = self.n_dim - self.num_nds;
            mat_mul_ar_dfd1(&mut d_ndx[i1..], &mut d_nds[i1..], &mut  j_inv_cent,  i2,  3,  3);
        }
        return;
    }

    pub fn get_inst_ori_dfd1(& self, inst_ori_mat : &mut Vec<DiffDoub1>, loc_ori : &mut Vec<DiffDoub1>, glob_disp : &mut Vec<DiffDoub1>, stat : usize) {
        // stat = 1: nonlinear geometry, Element ori from nodal theta, no derivatives, 1st order diff for self.nodes
        // stat = 2: nonlinear geometry, 2nd order diff for DiffDoub1 version, 1st order for DiffDoub1 version
        let mut rot = [DiffDoub1::new(); 3];
        let mut nnds = DiffDoub1::new();
        let mut one = DiffDoub1::new();
        let mut tmp_ori = [DiffDoub1::new(); 9];
        let mut tmp_inst = [DiffDoub1::new(); 9];
        let mut tmp = DiffDoub1::new();
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut st_index : usize;
        let i_ori_size : usize =  (self.num_nds+1)*144;
        for i1 in 0..i_ori_size {
            inst_ori_mat[i1].set_val(0.0);
        }
        let is_diff : bool =  loc_ori[0].diff_type();
        
        nnds.set_val(0.0);
        one.set_val(1.0);
        rot[0].set_val(0.0);
        rot[1].set_val(0.0);
        rot[2].set_val(0.0);
        i2 = 3 * self.n_dim;
        for _i1 in 0..self.num_nds {
            rot[0].add(& glob_disp[i2]);
            rot[1].add(& glob_disp[i2 + self.n_dim]);
            rot[2].add(& glob_disp[i2 + 2 * self.n_dim]);
            nnds.add(& one);
            i2 += 1usize;
        }
        rot[0].dvd(& nnds);
        rot[1].dvd(& nnds);
        rot[2].dvd(& nnds);
        
        vec_to_ar_dfd1(&mut tmp_ori, loc_ori,  0,  9);
        
        if stat == 1 {
            d_orid_thet_dfd1(&mut tmp_inst, &mut  tmp_ori, &mut  rot,  0,  0);
            ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  0,  9);
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd1(& glob_disp[i2]);
                rot[1].set_val_dfd1(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd1(& glob_disp[i2+2*self.n_dim]);
                st_index = 144 * i1;
                d_orid_thet_dfd1(&mut tmp_inst, &mut  tmp_ori, &mut  rot,  0,  0);
                ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                for i2 in 1..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd1(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd1(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
        } else if is_diff {
            for i2 in 0..4 {
                st_index = 36*i2;
                d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                i3 = st_index;
                i4 = 9*i2;
                for _i5 in 0..9 {
                    tmp.set_val_dfd1(& inst_ori_mat[i3]);
                    inst_ori_mat[i4].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4 += 1usize;
                }
            }
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd1(& glob_disp[i2]);
                rot[1].set_val_dfd1(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd1(& glob_disp[i2+2*self.n_dim]);
                for i2 in 0..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd1(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd1(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
        } else {
            for i2 in 0..4 {
                st_index = 36*i2;
                d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                i3 = st_index;
                i4 = 9*i2;
                for _i5 in 0..9 {
                    tmp.set_val_dfd1(& inst_ori_mat[i3]);
                    inst_ori_mat[i4].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4 += 1usize;
                }
                for i6 in i2..4 {
                    st_index = 36*i2 + 9*i6;
                    d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, i6);
                    ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 36*i6 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd1(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd1(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                }
            }
            for i1 in 1..=self.num_nds {
                i2 = 3*self.n_dim + i1 - 1;
                rot[0].set_val_dfd1(& glob_disp[i2]);
                rot[1].set_val_dfd1(& glob_disp[i2+self.n_dim]);
                rot[2].set_val_dfd1(& glob_disp[i2+2*self.n_dim]);
                for i2 in 0..4 {
                    st_index = 144*i1 + 36*i2;
                    d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, 0);
                    ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                    i3 = st_index;
                    i4 = 144*i1 + 9*i2;
                    for _i5 in 0..9 {
                        tmp.set_val_dfd1(& inst_ori_mat[i3]);
                        inst_ori_mat[i4].set_val_dfd1(& tmp);
                        i3 += 1usize;
                        i4 += 1usize;
                    }
                    for i6 in i2..4 {
                        st_index = 144*i1 + 36*i2 + 9*i6;
                        d_orid_thet_dfd1(&mut tmp_inst, &mut tmp_ori, &mut rot, i2, i6);
                        ar_to_vec_dfd1(&mut tmp_inst, inst_ori_mat,  st_index,  st_index + 9);
                        i3 = st_index;
                        i4 = 144*i1 + 36*i6 + 9*i2;
                        for _i5 in 0..9 {
                            tmp.set_val_dfd1(& inst_ori_mat[i3]);
                            inst_ori_mat[i4].set_val_dfd1(& tmp);
                            i3 += 1usize;
                            i4 += 1usize;
                        }
                    }
                }
            }
        }
        
        return;
    }

    pub fn get_inst_disp_dfd1(& self, inst_disp : &mut [DiffDoub1], glob_disp : &mut Vec<DiffDoub1>, inst_ori_mat : &mut Vec<DiffDoub1>, loc_ori : &mut Vec<DiffDoub1>, x_glob : &mut Vec<DiffDoub1>, n_lgeom : bool, dv1 : usize, dv2 : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let mut nd_oind : usize;
        let mut dof : usize;
        let mut dof_oind : usize;
        let nd2 : usize;
        let mut dof2 : usize;
        let dof2_oind : usize;
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut nn_inv = DiffDoub1::new();
        let mut nn_inv2 = DiffDoub1::new();
        
        i2 = 6*self.n_dim;
        for i1 in 0..i2 {
            inst_disp[i1].set_val(0.0);
        }
        
        if !n_lgeom {
            if dv1 == MAX_INT && dv2 == MAX_INT {
                i7 = 3 * self.n_dim;
                for i1 in 0..3 {
                    i4 = i1 * self.n_dim;
                    for i2 in 0..self.num_nds {
                        i5 = i1 * 3;
                        i6 = i2;
                        for _i3 in 0..3 {
                            //i4 = i1 * self.n_dim + i2;
                            //i5 = i1 * 3 + i3;
                            //i6 = i3 * self.n_dim + i2;
                            tmp.set_val_dfd1(& loc_ori[i5]);
                            tmp.mult(& glob_disp[i6]);
                            inst_disp[i4].add(& tmp);
                            tmp.set_val_dfd1(& loc_ori[i5]);
                            tmp.mult(& glob_disp[i6 + i7]);
                            inst_disp[i4 + i7].add(& tmp);
                            i5 += 1usize;
                            i6  +=  self.n_dim;
                        }
                        i4 += 1usize;
                    }
                }
                i2 = 2 * self.num_nds * self.dof_per_nd;
                for _i1 in 0..self.num_int_dof {
                    nd = self.dof_table[i2];
                    dof = self.dof_table[i2 + 1];
                    i3 = dof * self.n_dim + nd;
                    inst_disp[i3].set_val_dfd1(& glob_disp[i3]);
                    i2  +=  2;
                }
            }
            else if (dv1 + dv2) >= MAX_INT {
                let dnz = dv1 + dv2 - MAX_INT;
                //dv1 = dv1 + dv2 - MAX_INT;
                nd = self.dof_table[2 * dnz];
                dof = self.dof_table[2 * dnz + 1];
                if dof < 3 {
                    if nd < self.num_nds {
                        i2 = nd;
                        i3 = dof;
                        for _i1 in 0..3 {
                            //i2 = i1 * self.n_dim + nd;
                            //i3 = i1 * 3 + dof;
                            inst_disp[i2].set_val_dfd1(& loc_ori[i3]);
                            i2  +=  self.n_dim;
                            i3  +=  3;
                        }
                    }
                    else {
                        i1 = dof * self.n_dim + nd;
                        inst_disp[i1].set_val(1.0);
                    }
                }
                else {
                    i2 = 3 * self.n_dim + nd;
                    i3 = dof - 3;
                    for _i1 in 0..3 {
                        //i2 = (i1 + 3) * self.n_dim + nd;
                        //i3 = i1 * 3 + (dof - 3);
                        inst_disp[i2].set_val_dfd1(& loc_ori[i3]);
                        i2  +=  self.n_dim;
                        i3  +=  3;
                    }
                }
            }
        }
        else {
            if dv1 == MAX_INT && dv2 == MAX_INT {
                for i1 in 0..3 {
                    for i2 in 0..self.num_nds {
                        i4 = i1 * self.n_dim + i2;
                        i5 = i2;
                        i6 = i2;
                        i7 = i1 * 3;
                        for _i3 in 0..3 {
                            tmp.set_val_dfd1(& glob_disp[i5]);
                            tmp.add(& x_glob[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp2.set_val_dfd1(& loc_ori[i7]);
                            tmp2.mult(& x_glob[i6]);
                            tmp.sub(& tmp2);
                            inst_disp[i4].add(& tmp);
                            i5  +=  self.n_dim;
                            i6  +=  self.num_nds;
                            i7 += 1usize;
                        }
                    }
                }
                
                i2 = self.num_nds * self.dof_per_nd;
                i3 = i2 + self.num_int_dof;
                for i1 in i2..i3 {
                    i4 = self.dof_table[2 * i1];
                    i5 = self.dof_table[2 * i1 + 1];
                    i6 = i5 * self.n_dim + i4;
                    inst_disp[i6].set_val_dfd1(& glob_disp[i6]);
                }
                
                for i1 in 0..self.num_nds {
                    i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                    i4 = 4 * self.n_dim + i1;
                    i5 = 5 * self.n_dim + i1;
                    nd_oind = 144 * (i1 + 1);
                    for i2 in 0..3 {
                        i6 = nd_oind + 3 + i2;
                        i7 = 6 + i2;
                        tmp.set_val_dfd1(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i3].add(& tmp);
                        i6 = nd_oind + 6 + i2;
                        i7 = i2;
                        tmp.set_val_dfd1(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i4].add(& tmp);
                        i6 = nd_oind + i2;
                        i7 = 3 + i2;
                        tmp.set_val_dfd1(& inst_ori_mat[i6]);
                        tmp.mult(& inst_ori_mat[i7]);
                        inst_disp[i5].add(& tmp);
                    }
                }
            }
            else if (dv1 + dv2) >= MAX_INT {
                let dnz = dv1 + dv2 - MAX_INT;
                //dv1 = dv1 + dv2 - MAX_INT;
                nd = self.dof_table[2 * dnz];
                dof = self.dof_table[2 * dnz + 1];
                if dof < 3 {
                    if nd < self.num_nds {
                        i1 = nd;// index in inst_disp
                        i2 = dof;//index in inst_ori_mat
                        inst_disp[i1].set_val_dfd1(& inst_ori_mat[i2]);
                        i1 = self.n_dim + nd;
                        i2 = 3 + dof;
                        inst_disp[i1].set_val_dfd1(& inst_ori_mat[i2]);
                        i1 = 2 * self.n_dim + nd;
                        i2 = 6 + dof;
                        inst_disp[i1].set_val_dfd1(& inst_ori_mat[i2]);
                    }
                    else {
                        i1 = dof * self.n_dim + nd;
                        inst_disp[i1].set_val(1.0);
                    }
                }
                else {// dof is rotation
                    nn_inv.set_val(1.0 / (self.num_nds as f64));
                    dof_oind = 36 * (dof - 2);
                    for i1 in 0..3 {
                        for i2 in 0..self.num_nds {
                            i4 = i1 * self.n_dim + i2;
                            i5 = i2;
                            i6 = i2;
                            i7 = dof_oind + i1 * 3;
                            for _i3 in 0..3 {
                                tmp.set_val_dfd1(& glob_disp[i5]);
                                tmp.add(& x_glob[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i5  +=  self.n_dim;
                                i6  +=  self.num_nds;
                                i7 += 1usize;
                            }
                        }
                    }
                    
                    for i1 in 0..self.num_nds {
                        i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                        i4 = 4 * self.n_dim + i1;
                        i5 = 5 * self.n_dim + i1;
                        nd_oind = 144 * (i1 + 1);
                        for i2 in 0..3 {
                            i6 = nd_oind + 3 + i2;
                            i7 = dof_oind + 6 + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i3].add(& tmp);
                            i6 = nd_oind + 6 + i2;
                            i7 = dof_oind + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i4].add(& tmp);
                            i6 = nd_oind + i2;
                            i7 = dof_oind + 3 + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv);
                            inst_disp[i5].add(& tmp);
                            if i1 == nd {
                                i6 = nd_oind + dof_oind + 3 + i2;
                                i7 = 6 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + 6 + i2;
                                i7 = i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + i2;
                                i7 = 3 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i5].add(& tmp);
                            }
                        }
                    }
                }
            }
            else {
                nd = self.dof_table[2 * dv1];
                dof = self.dof_table[2 * dv1 + 1];
                nd2 = self.dof_table[2 * dv2];
                dof2 = self.dof_table[2 * dv2 + 1];
                nn_inv.set_val(1.0 / (self.num_nds as f64));
                nn_inv2.set_val_dfd1(& nn_inv);
                nn_inv2.sqr();
                if dof > 2 && dof2 > 2 {
                    for i1 in 0..3 {
                        dof_oind = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
                        for i2 in 0..self.num_nds {
                            i4 = self.n_dim * i1 + i2;
                            i5 = i2;
                            i6 = i2;
                            i7 = dof_oind;
                            for _i3 in 0..3 {
                                tmp.set_val_dfd1(& glob_disp[i5]);
                                tmp.add(& x_glob[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv2);
                                inst_disp[i4].add(& tmp);
                                i5  +=  self.n_dim;
                                i6  +=  self.num_nds;
                                i7 += 1usize;
                            }
                        }
                    }
                    
                    dof_oind = 36 * (dof - 2);
                    dof2_oind = 9 * (dof2 - 2);
                    for i1 in 0..self.num_nds {
                        i3 = 3 * self.n_dim + i1;// indexes of thetax, y and z for Node i1
                        i4 = 4 * self.n_dim + i1;
                        i5 = 5 * self.n_dim + i1;
                        nd_oind = 144 * (i1 + 1);
                        for i2 in 0..3 {
                            i6 = nd_oind + 3 + i2;
                            i7 = dof_oind + dof2_oind + 6 + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i3].add(& tmp);
                            i6 = nd_oind + 6 + i2;
                            i7 = dof_oind + dof2_oind + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i4].add(& tmp);
                            i6 = nd_oind + i2;
                            i7 = dof_oind + dof2_oind + 3 + i2;
                            tmp.set_val_dfd1(& inst_ori_mat[i6]);
                            tmp.mult(& inst_ori_mat[i7]);
                            tmp.mult(& nn_inv2);
                            inst_disp[i5].add(& tmp);
                            if i1 == nd {
                                i6 = nd_oind + dof_oind + 3 + i2;
                                i7 = dof2_oind + 6 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + 6 + i2;
                                i7 = dof2_oind + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + i2;
                                i7 = dof2_oind + 3 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i5].add(& tmp);
                            }
                            if i1 == nd2 {
                                i6 = nd_oind + dof2_oind + 3 + i2;
                                i7 = dof_oind + 6 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof2_oind + 6 + i2;
                                i7 = dof_oind + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof2_oind + i2;
                                i7 = dof_oind + 3 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                tmp.mult(& nn_inv);
                                inst_disp[i5].add(& tmp);
                            }
                            if i1 == nd && i1 == nd2 {
                                i6 = nd_oind + dof_oind + dof2_oind + 3 + i2;
                                i7 = 6 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i3].add(& tmp);
                                i6 = nd_oind + dof_oind + dof2_oind + 6 + i2;
                                i7 = i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i4].add(& tmp);
                                i6 = nd_oind + dof_oind + dof2_oind + i2;
                                i7 = 3 + i2;
                                tmp.set_val_dfd1(& inst_ori_mat[i6]);
                                tmp.mult(& inst_ori_mat[i7]);
                                inst_disp[i5].add(& tmp);
                            }
                        }
                    }
                }
                else if dof < 3 && dof2 < 3 {
                    return;
                }
                else {
                    if dof > 2 {
                        i1 = dof;
                        dof = dof2;
                        dof2 = i1;
                        //i1 = nd;
                        nd = nd2;
                        //nd2 = i1;
                    }
                    i1 = 36 * (dof2 - 2) + dof;
                    tmp.set_val_dfd1(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = nd;
                    inst_disp[i2].add(& tmp);
                    i1 = 36 * (dof2 - 2) + 3 + dof;
                    tmp.set_val_dfd1(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = self.n_dim + nd;
                    inst_disp[i2].add(& tmp);
                    i1 = 36 * (dof2 - 2) + 6 + dof;
                    tmp.set_val_dfd1(& inst_ori_mat[i1]);
                    tmp.mult(& nn_inv);
                    i2 = 2 * self.n_dim + nd;
                    inst_disp[i2].add(& tmp);
                }
            }
        }

        return;
    }

    pub fn get_stress_prereq_dfd1(& self, pre : &mut DiffDoub1StressPrereq, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut offset = DiffDoub1::new();
        self.get_nd_crds_dfd1(&mut pre.glob_nds, nd_ar, dv_ar);
        self.get_loc_ori_dfd1(&mut pre.loc_ori, sec_ar, dv_ar);
        self.get_nd_disp_dfd1(&mut pre.glob_disp, nd_ar);
        self.get_nd_vel_dfd1(&mut pre.glob_vel, nd_ar);
        self.get_nd_acc_dfd1(&mut pre.glob_acc, nd_ar);
        self.get_nd_temp_dfd1(&mut pre.glob_temp, nd_ar);
        self.get_nd_tdot_dfd1(&mut pre.glob_tdot, nd_ar);
        if self.dof_per_nd == 6 {
            self.correct_orient_dfd1(&mut pre.loc_ori, &mut  pre.glob_nds);
            if self.this_type != 2 {
                self.get_layer_thk_z_dfd1(&mut pre.layer_thk, &mut  pre.layer_z, &mut  offset, sec_ar, dv_ar);
                self.get_layer_angle_dfd1(&mut pre.layer_ang, sec_ar, dv_ar);
                self.get_layer_q_dfd1(&mut pre.layer_q, sec_ar, mat_ar, dv_ar);
                self.get_layer_d_dfd1(&mut pre.layer_d, sec_ar, mat_ar, dv_ar);
                self.get_layer_th_exp_dfd1(&mut pre.layer_te, sec_ar, mat_ar, dv_ar);
                self.get_layer_einit_dfd1(&mut pre.layer_e0, sec_ar, dv_ar);
                self.get_layer_den_dfd1(&mut pre.layer_den, sec_ar, mat_ar, dv_ar);
                self.get_layer_cond_dfd1(&mut pre.layer_tc, sec_ar, mat_ar, dv_ar);
                self.get_layer_spec_heat_dfd1(&mut pre.layer_sh, sec_ar, mat_ar, dv_ar);
                self.get_abd_dfd1(&mut pre.cmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_q, &mut  pre.layer_ang, sec_ar);
                self.get_shell_damp_dfd1(&mut pre.dmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_d, &mut  pre.layer_ang, sec_ar);
                self.get_shell_exp_load_dfd1(&mut pre.therm_exp, &mut  pre.einit, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_q, &mut  pre.layer_te, &mut  pre.layer_e0, &mut  pre.layer_ang, sec_ar);
                self.get_shell_mass_dfd1(&mut pre.mmat, &mut  pre.layer_thk, &mut  pre.layer_z, &mut  pre.layer_den, sec_ar);
                self.get_shell_cond_dfd1(&mut pre.tcmat, &mut  pre.layer_thk, &mut  pre.layer_ang, &mut  pre.layer_tc, sec_ar);
                self.get_shell_spec_heat_dfd1(&mut pre.spec_heat, &mut  pre.layer_thk, &mut  pre.layer_sh, &mut  pre.layer_den, sec_ar);
            }
            else {
                self.get_beam_stiff_dfd1(&mut pre.cmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_damp_dfd1(&mut pre.dmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_exp_load_dfd1(&mut pre.therm_exp, &mut  pre.einit, sec_ar, mat_ar, dv_ar);
                self.get_beam_mass_dfd1(&mut pre.mmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_cond_dfd1(&mut pre.tcmat, sec_ar, mat_ar, dv_ar);
                self.get_beam_spec_heat_dfd1(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
            }
        }
        else if self.this_type == 21 {
            self.get_frc_fld_const_dfd1(&mut pre.frc_fld_coef, &mut  pre.frc_fld_exp, sec_ar, dv_ar);
            self.get_thrm_fld_const_dfd1(&mut pre.thrm_fld_coef, &mut  pre.ref_temp, sec_ar, dv_ar);
        }
        else if self.this_type == 1 {
            self.get_mass_per_el_dfd1(&mut pre.mass_per_el, sec_ar, dv_ar);
            self.get_specific_heat_dfd1(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
        }
        else {
            self.get_solid_stiff_dfd1(&mut pre.cmat, sec_ar, mat_ar, dv_ar);
            self.get_solid_damp_dfd1(&mut pre.dmat, sec_ar, mat_ar, dv_ar);
            self.get_thermal_exp_dfd1(&mut pre.therm_exp, &mut  pre.einit, sec_ar, mat_ar, dv_ar);
            self.get_density_dfd1(&mut pre.mmat[0],  0, sec_ar, mat_ar, dv_ar);
            self.get_conductivity_dfd1(&mut pre.tcmat, sec_ar, mat_ar, dv_ar);
            self.get_specific_heat_dfd1(&mut pre.spec_heat, sec_ar, mat_ar, dv_ar);
        }
        mat_mul_dfd1(&mut pre.loc_nds, &mut  pre.loc_ori, &mut  pre.glob_nds,  3,  3,  self.num_nds);
        
        
        return;
    }

    pub fn get_fluid_prereq_dfd1(&mut self, pre : &mut DiffDoub1FlPrereq, sec_ar : &mut Vec<Section>, fl_ar : &mut Vec<Fluid>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>, el_cent : & Vec<CentData>) {
        let i2 : usize;
        let mut sec_prop : f64;
        
        pre.num_nds = self.num_nds;
        self.get_nd_crds_dfd1(&mut pre.glob_nds, nd_ar, dv_ar);
        self.get_nd_disp_dfd1(&mut pre.glob_disp, nd_ar);
        i2 = 3 * self.num_nds;
        for i1 in 0..i2 {
            pre.glob_nds[i1].add(& pre.glob_disp[i1]);
        }
        self.get_nd_vel_dfd1(&mut pre.glob_vel, nd_ar);
        self.get_nd_fl_den_dfd1(&mut pre.fl_den, nd_ar);
        self.get_nd_fl_vel_dfd1(&mut pre.fl_vel, nd_ar);
        self.get_nd_temp_dfd1(&mut pre.fl_temp, nd_ar);
        self.get_nd_turb_e_dfd1(&mut pre.fl_turb_e, nd_ar);
        self.get_nd_fl_den_dot_dfd1(&mut pre.fl_den_dot, nd_ar);
        self.get_nd_fl_vdot_dfd1(&mut pre.fl_vel_dot, nd_ar);
        self.get_nd_tdot_dfd1(&mut pre.fl_tdot, nd_ar);
        self.get_nd_turb_edot_dfd1(&mut pre.fl_turb_edot, nd_ar);

        let ec : &CentData = &el_cent[self.label];
        for i in 0..3 {
            pre.f_per_mass[i].set_val_dfd1(&ec.f_dfd1[i]);
        }
        self.get_f_per_mass_dfd1(&mut pre.f_per_mass, dv_ar);
        pre.hg_per_vol.set_val_dfd1(&ec.hg_dfd1);
        self.get_gen_prop_dfd1(&mut pre.hg_per_vol, &mut CppStr::from("bodyHeatGen"), dv_ar);
        
        let this_sec = &sec_ar[self.sect_ptr];
        let this_fl = &fl_ar[this_sec.fl_ptr];

        pre.compressible = this_fl.compressible;
        
        sec_prop = this_fl.grad_turb_coef;
        pre.grad_turb_coef.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.grad_turb_coef, &mut CppStr::from("gradVTurbCoef"), dv_ar);
        
        sec_prop = this_fl.diss_turb_coef;
        pre.diss_turb_coef.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.diss_turb_coef, &mut CppStr::from("dissTurbCoef"), dv_ar);
        
        sec_prop = this_fl.viscosity;
        pre.ref_visc.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.ref_visc, &mut CppStr::from("viscosity"), dv_ar);
        
        sec_prop = this_fl.temp_vis_coef;
        pre.temp_vis_coef.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.temp_vis_coef, &mut CppStr::from("tempVisCoef"), dv_ar);
        
        sec_prop = this_fl.turb_vis_coef;
        pre.turb_vis_coef.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.turb_vis_coef, &mut CppStr::from("turbVisCoef"), dv_ar);
        
        sec_prop = this_fl.ref_enth;
        pre.ref_enth.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.ref_enth, &mut CppStr::from("refEnth"), dv_ar);
        
        sec_prop = this_fl.ref_den;
        pre.ref_den.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.ref_den, &mut CppStr::from("refDen"), dv_ar);
        
        sec_prop = this_fl.ref_temp;
        pre.ref_temp.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.ref_temp, &mut CppStr::from("refTemp"), dv_ar);
        
        sec_prop = this_fl.therm_cond;
        pre.therm_cond.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.therm_cond, &mut CppStr::from("thermCond"), dv_ar);

        sec_prop = this_fl.expansion;
        pre.expansion.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.expansion, &mut CppStr::from("thermExp"), dv_ar);
        
        sec_prop = this_fl.spec_heat;
        pre.spec_heat.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.spec_heat, &mut CppStr::from("specHeat"), dv_ar);
        
        sec_prop = this_fl.ideal_gas;
        pre.i_gconst.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.i_gconst, &mut CppStr::from("iGConst"), dv_ar);

        sec_prop = this_fl.bulk_modulus;
        pre.bulk_mod.set_val(sec_prop);
        self.get_gen_prop_dfd1(&mut pre.bulk_mod, &mut CppStr::from("bulkModulus"), dv_ar);

        if pre.compressible {
            pre.den_vis_coef.set_val_dfd1(&pre.ref_visc);
            pre.den_vis_coef.dvd(&pre.ref_den);
            pre.den_vis_coef.dvd(&pre.ref_temp);
            
            pre.den_pres_coef.set_val_dfd1(&pre.bulk_mod);
            pre.den_pres_coef.dvd(&pre.ref_den);
            
            pre.temp_pres_coef.set_val(3.0);
            pre.temp_pres_coef.mult(&pre.bulk_mod);
            pre.temp_pres_coef.mult(&pre.expansion);
        }
        else {
            pre.den_vis_coef.set_val(0.0);
            pre.den_pres_coef.set_val(0.0);
            pre.temp_pres_coef.set_val(0.0);
        }
        
        
        return;
    }

    pub fn get_volume_dfd1(& self, vol : &mut DiffDoub1, pre : &mut DiffDoub1StressPrereq, layer : usize, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut thk = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut s_tmp : [f64; 3] = [0f64; 3];
        let mut this_dv : &DesignVariable;
        
        if self.this_type == 2 {
            thk.set_val(sec_ar[self.sect_ptr].area);//rem: update to factor in dvars;
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    thk.add(& dv_val);
                }
            }
        } else if self.this_type == 3 || self.this_type == 41 {
            thk.set_val_dfd1(& pre.layer_thk[layer]);
        } else {
            thk.set_val(1.0);
        }
        
        vol.set_val(0.0);
        for i1 in 0..self.num_ip {
            i2 = 3 * i1;
            vec_to_ar(&mut s_tmp, &self.int_pts,  i2,  i2 + 3);
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, &mut  s_tmp);
            tmp.set_val(self.ip_wt[i1]);
            tmp.mult(& det_j);
            tmp.mult(& thk);
            vol.add(& tmp);
        }
        
        return;
    }

    pub fn get_section_def_dfd1(&mut self, sec_def : &mut [DiffDoub1], glob_disp : &mut Vec<DiffDoub1>, inst_ori_mat : &mut Vec<DiffDoub1>, loc_ori : &mut Vec<DiffDoub1>, x_glob : &mut Vec<DiffDoub1>, d_ndx : &mut [DiffDoub1], n_vec : &mut [DiffDoub1], n_lgeom : bool, dv1 : usize, dv2 : usize) {
        let i1 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut nd : usize;
        let dof : usize;
        let nd2 : usize;
        let dof2 : usize;
        let mut inst_disp = [DiffDoub1::new(); 60];
        let mut ux = [DiffDoub1::new(); 9];
        let mut rx = [DiffDoub1::new(); 9];
        let mut rot = [DiffDoub1::new(); 3];
        let mut tmp = DiffDoub1::new();
        
        if self.dof_per_nd != 6 {
            for i1 in 0..9 {
                sec_def[i1].set_val(0.0);
            }
            return;
        }
        
        self.get_inst_disp_dfd1(&mut inst_disp, glob_disp, inst_ori_mat, loc_ori, x_glob,  n_lgeom,  dv1,  dv2);
        
        if dv1 == MAX_INT && dv2 == MAX_INT {
            mat_mul_ar_dfd1(&mut ux, &mut inst_disp, d_ndx, 3, self.n_dim, 3);
            i1 = 3*self.n_dim;
            mat_mul_ar_dfd1(&mut rx, &mut inst_disp[i1..], d_ndx, 3, self.n_dim, 3);
            mat_mul_ar_dfd1(&mut rot, &mut inst_disp[i1..], n_vec, 3, self.n_dim, 1);
        } else if (dv1 + dv2) >= MAX_INT {
            let dnz = dv1 + dv2 - MAX_INT;
            //dv1 = dv1 + dv2 - MAX_INT;
            nd = self.dof_table[2*dnz];
            dof = self.dof_table[2*dnz+1];
            if dof < 3 {
                i3 = 0;
                i4 = nd;
                for i1 in 0..3 {
                    i5 = 3*nd;
                    for _i2 in 0..3 {
                        ux[i3].set_val_dfd1(& inst_disp[i4]);
                        ux[i3].mult(& d_ndx[i5]);
                        rx[i3].set_val(0.0);
                        i3 += 1usize;
                        i5 += 1usize;
                    }
                    rot[i1].set_val(0.0);
                    i4 +=  self.n_dim;
                }
            } else {
                i4 = 0;
                for i1 in 0..3 {
                    for i2 in 0..3 {
                        i5 = self.n_dim*i1;
                        i6 = i2;
                        i7 = self.n_dim*(i1+3);
                        ux[i4].set_val(0.0);
                        rx[i4].set_val(0.0);
                        for _i3 in 0..self.num_nds {
                            tmp.set_val_dfd1(& inst_disp[i5]);
                            tmp.mult(& d_ndx[i6]);
                            ux[i4].add(& tmp);
                            tmp.set_val_dfd1(& inst_disp[i7]);
                            tmp.mult(& d_ndx[i6]);
                            rx[i4].add(& tmp);
                            i5 += 1usize;
                            i6 +=  3;
                            i7 += 1usize;
                        }
                        i4 += 1usize;
                    }
                }
                for i1 in 0..3 {
                    rot[i1].set_val(0.0);
                    i3 = self.n_dim*(i1+3);
                    for i2 in 0..self.num_nds {
                        tmp.set_val_dfd1(& inst_disp[i3]);
                        tmp.mult(& n_vec[i2]);
                        rot[i1].add(& tmp);
                        i3 += 1usize;
                    }
                }
            }
        } else {
            nd = self.dof_table[2*dv1];
            dof = self.dof_table[2*dv1+1];
            nd2 = self.dof_table[2*dv2];
            dof2 = self.dof_table[2*dv2+1];
            if dof < 3 && dof2 < 3 {
                for i1 in 0..9 {
                    sec_def[i1].set_val(0.0);
                }
                return;
            } else if dof > 2 && dof2 > 2 {
                i4 = 0;
                for i1 in 0..3 {
                    for i2 in 0..3 {
                        i5 = self.n_dim*i1;
                        i6 = i2;
                        i7 = self.n_dim*(i1+3);
                        ux[i4].set_val(0.0);
                        for _i3 in 0..self.num_nds {
                            tmp.set_val_dfd1(& inst_disp[i5]);
                            tmp.mult(& d_ndx[i6]);
                            ux[i4].add(& tmp);
                            tmp.set_val_dfd1(& inst_disp[i7]);
                            tmp.mult(& d_ndx[i6]);
                            rx[i4].add(& tmp);
                            i5 += 1usize;
                            i6 +=  3;
                            i7 += 1usize;
                        }
                        i4 += 1usize;
                    }
                }
                for i1 in 0..3 {
                    rot[i1].set_val(0.0);
                    i3 = self.n_dim*(i1+3);
                    for i2 in 0..self.num_nds {
                        tmp.set_val_dfd1(& inst_disp[i3]);
                        tmp.mult(& n_vec[i2]);
                        rot[i1].add(& tmp);
                        i3 += 1usize;
                    }
                }
            } else {
                if dof > dof2 {
                    //i1 = dof;
                    //dof = dof2;
                    //dof2 = i1;
                    //i1 = nd;
                    nd = nd2;
                    //nd2 = i1;
                }
                i3 = 0;
                i4 = nd;
                for i1 in 0..3 {
                    i5 = 3*nd;
                    for _i2 in 0..3 {
                        ux[i3].set_val_dfd1(& inst_disp[i4]);
                        ux[i3].mult(& d_ndx[i5]);
                        rx[i3].set_val(0.0);
                        i3 += 1usize;
                        i5 += 1usize;
                    }
                    rot[i1].set_val(0.0);
                    i4 +=  self.n_dim;
                }
            }
        }
        
        if self.this_type == 2 {
            sec_def[0].set_val_dfd1(& ux[0]);
            sec_def[1].set_val_dfd1(& ux[3]);
            sec_def[1].sub(& rot[2]);
            sec_def[2].set_val_dfd1(& ux[6]);
            sec_def[2].add(& rot[1]);
            sec_def[3].set_val_dfd1(& rx[0]);
            sec_def[4].set_val_dfd1(& rx[3]);
            sec_def[5].set_val_dfd1(& rx[6]);
        } else {
            sec_def[0].set_val_dfd1(& ux[0]);
            sec_def[1].set_val_dfd1(& ux[4]);
            sec_def[2].set_val_dfd1(& ux[1]);
            sec_def[2].add(& ux[3]);
            sec_def[3].set_val_dfd1(& rx[3]);
            sec_def[4].set_val_dfd1(& rx[1]);
            sec_def[4].neg();
            sec_def[5].set_val_dfd1(& rx[4]);
            sec_def[5].sub(& rx[0]);
            sec_def[6].set_val_dfd1(& ux[6]);
            sec_def[6].add(& rot[1]);
            sec_def[7].set_val_dfd1(& ux[7]);
            sec_def[7].sub(& rot[0]);
            sec_def[8].set_val(2.0);
            sec_def[8].mult(& rot[2]);
            sec_def[8].sub(& ux[3]);
            sec_def[8].add(& ux[1]);
        }
        
        return;
    }

    pub fn get_solid_strain_dfd1(&mut self, strain : &mut [DiffDoub1], ux : &mut [DiffDoub1], d_ndx : &mut [DiffDoub1], loc_ori : &mut Vec<DiffDoub1>, dv1 : usize, dv2 : usize, n_lgeom : bool) {
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let nd : usize;
        let dof : usize;
        let nd2 : usize;
        let dof2 : usize;
        let mut ux_l = [DiffDoub1::new(); 9];
        let mut strn_mat = [DiffDoub1::new(); 9];
        let mut tmp_ori = [DiffDoub1::new(); 9];
        let mut tmp = DiffDoub1::new();
        
        for i1 in 0..9 {
            strn_mat[i1].set_val(0.0);
        }
        if dv1 == MAX_INT && dv2 == MAX_INT {
            vec_to_ar_dfd1(&mut tmp_ori, loc_ori,  0,  9);
            mat_mul_ar_dfd1(&mut ux_l, &mut tmp_ori, ux, 3, 3, 3);
            for i1 in 0..3 {
                i4 = 4*i1;
                i5 = 4*i1;
                for i2 in i1..3 {
                    strn_mat[i4].add(& ux_l[i4]);
                    strn_mat[i4].add(& ux_l[i5]);
                    if n_lgeom {
                        i6 = i1;
                        i7 = i2;
                        for _i3 in 0..3 {
                            tmp.set_val_dfd1(& ux[i6]);
                            tmp.mult(& ux[i7]);
                            strn_mat[i4].add(& tmp);
                            i6 +=  3;
                            i7 +=  3;
                        }
                    }
                    i4 += 1usize;
                    i5 +=  3;
                }
            }
        } else if (dv1 + dv2) >= MAX_INT {
            if dv1 < MAX_INT {
                nd = self.dof_table[2*dv1];
                dof = self.dof_table[2*dv1+1];
            } else {
                nd = self.dof_table[2*dv2];
                dof = self.dof_table[2*dv2+1];
            }
            for i1 in 0..3 {
                i4 = 4*i1;
                for i2 in i1..3 {
                    i5 = 3*i1 + dof;
                    i6 = 3*nd + i2;
                    tmp.set_val_dfd1(& loc_ori[i5]);
                    tmp.mult(& d_ndx[i6]);
                    strn_mat[i4].add(& tmp);
                    i5 = 3*i2 + dof;
                    i6 = 3*nd + i1;
                    tmp.set_val_dfd1(& loc_ori[i5]);
                    tmp.mult(& d_ndx[i6]);
                    strn_mat[i4].add(& tmp);
                    if n_lgeom {
                        i5 = 3*nd + i1;
                        i6 = 3*dof + i2;
                        tmp.set_val_dfd1(& d_ndx[i5]);
                        tmp.mult(& ux[i6]);
                        strn_mat[i4].add(& tmp);
                        i5 = 3*nd + i2;
                        i6 = 3*dof + i1;
                        tmp.set_val_dfd1(& d_ndx[i5]);
                        tmp.mult(& ux[i6]);
                        strn_mat[i4].add(& tmp);
                    }
                    i4 += 1usize;
                }
            }
        } else {
            if n_lgeom {
                nd = self.dof_table[2*dv1];
                dof = self.dof_table[2*dv1+1];
                nd2 = self.dof_table[2*dv2];
                dof2 = self.dof_table[2*dv2+1];
                if dof == dof2 {
                    for i1 in 0..3 {
                        i4 = 4*i1;
                        for i2 in i1..3 {
                            i5 = 3*nd + i1;
                            i6 = 3*nd2 + i2;
                            tmp.set_val_dfd1(& d_ndx[i5]);
                            tmp.mult(& d_ndx[i6]);
                            strn_mat[i4].add(& tmp);
                            i5 = 3*nd2 + i1;
                            i6 = 3*nd + i2;
                            tmp.set_val_dfd1(& d_ndx[i5]);
                            tmp.mult(& d_ndx[i6]);
                            strn_mat[i4].add(& tmp);
                            i4 += 1usize;
                        }
                    }
                }
            }
        }
        
        tmp.set_val(0.5);
        strain[0].set_val_dfd1(& strn_mat[0]);
        strain[0].mult(& tmp);
        strain[1].set_val_dfd1(& strn_mat[4]);
        strain[1].mult(& tmp);
        strain[2].set_val_dfd1(& strn_mat[8]);
        strain[2].mult(& tmp);
        strain[3].set_val_dfd1(& strn_mat[1]);
        strain[4].set_val_dfd1(& strn_mat[2]);
        strain[5].set_val_dfd1(& strn_mat[5]);
        
        return;
    }

    pub fn get_stress_strain_dfd1(&mut self, stress : &mut [DiffDoub1], strain : &mut [DiffDoub1], t_strain : &mut [DiffDoub1], spt : &mut [f64], layer : usize, n_lgeom : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i2 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ux = [DiffDoub1::new(); 9];
        let mut sec_def = [DiffDoub1::new(); 9];
        let mut sect_strn = [DiffDoub1::new(); 3];
        let mut adj_stn = [DiffDoub1::new(); 6];
        let mut tmp = DiffDoub1::new();
        let mut ip_temp = DiffDoub1::new();
        let mut tmp_ar = [DiffDoub1::new(); 60];

        for i in 0..6 {
            stress[i].set_val(0f64);
            strain[i].set_val(0f64);
            t_strain[i].set_val(0f64);
        }
        
        self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        
        ip_temp.set_val(0.0);
        for i1 in 0..self.num_nds {
            tmp.set_val_dfd1(& pre.glob_temp[i1]);
            tmp.mult(& n_vec[i1]);
            ip_temp.add(& tmp);
        }
        
        if self.this_type == 41 || self.this_type == 3 {
            if n_lgeom {
                self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            
            self.get_section_def_dfd1(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  MAX_INT,  MAX_INT);
            
            sect_strn[0].set_val_dfd1(& sec_def[0]);
            tmp.set_val_dfd1(& pre.layer_z[layer]);
            tmp.mult(& sec_def[3]);
            sect_strn[0].add(& tmp);
            
            sect_strn[1].set_val_dfd1(& sec_def[1]);
            tmp.set_val_dfd1(& pre.layer_z[layer]);
            tmp.mult(& sec_def[4]);
            sect_strn[1].sub(& tmp);
            
            sect_strn[2].set_val_dfd1(& sec_def[2]);
            tmp.set_val_dfd1(& pre.layer_z[layer]);
            tmp.mult(& sec_def[5]);
            sect_strn[2].add(& tmp);
            
            tmp.set_val_dfd1(& pre.layer_ang[layer]);
            tmp.neg();
            self.transform_strain_dfd1(strain, &mut  sect_strn, &mut  tmp);
            i2 = 3 * layer;
            for i1 in 0..3 {
                tmp.set_val_dfd1(& pre.layer_te[i2]);
                tmp.mult(& ip_temp);
                t_strain[i1].set_val_dfd1(&tmp);
                adj_stn[i1].set_val_dfd1(& strain[i1]);
                adj_stn[i1].sub(& tmp);
                adj_stn[i1].sub(& pre.layer_e0[i2]);
                i2 += 1usize;
            }
            mat_mul_ar_dfd1(stress, &mut  pre.layer_q[(9*layer)..], &mut  adj_stn,  3,  3,  1);
            
            tmp.set_val_dfd1(& strain[2]);
            strain[3].set_val_dfd1(& tmp);
            strain[2].set_val(0.0);

            tmp.set_val_dfd1(& t_strain[2]);
            t_strain[3].set_val_dfd1(&tmp);
            t_strain[2].set_val(0.0);
            
            tmp.set_val_dfd1(& stress[2]);
            stress[3].set_val_dfd1(& tmp);
            stress[2].set_val(0.0);
            
        } else if self.this_type != 2 {
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.glob_disp,  0,  60);
            mat_mul_ar_dfd1(&mut ux, &mut  tmp_ar, &mut  d_ndx,  3,  self.n_dim,  3);
            self.get_solid_strain_dfd1(strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  MAX_INT,  MAX_INT,  n_lgeom);
            for i1 in 0..6 {
                tmp.set_val_dfd1(& pre.therm_exp[i1]);
                tmp.mult(& ip_temp);
                t_strain[i1].set_val_dfd1(&tmp);
                adj_stn[i1].set_val_dfd1(& strain[i1]);
                adj_stn[i1].sub(& tmp);
                adj_stn[i1].sub(& pre.einit[i1]);
            }
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.cmat,  0,  36);
            mat_mul_ar_dfd1(stress, &mut tmp_ar, &mut adj_stn,  6,  6,  1);
        }
        return;
    }

    pub fn d_stress_straind_u_dfd1(&mut self, dsd_u : &mut Vec<DiffDoub1>, ded_u : &mut Vec<DiffDoub1>, dsd_t : &mut Vec<DiffDoub1>, spt : &mut [f64], layer : usize, n_lgeom : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ux = [DiffDoub1::new(); 9];
        let mut sec_def = [DiffDoub1::new(); 9];
        let mut sect_strn = [DiffDoub1::new(); 3];
        let mut d_strain = [DiffDoub1::new(); 6];
        let mut d_stress = [DiffDoub1::new(); 6];
        let mut cte = [DiffDoub1::new(); 6];
        let mut cten = [DiffDoub1::new(); 60];
        let mut tmp = DiffDoub1::new();
        let mut tmp_ar = [DiffDoub1::new(); 60];
        let mut tmp_ar2 = [DiffDoub1::new(); 6];
        
        let tot_dof : usize =  self.num_nds * self.dof_per_nd + self.num_int_dof;
        i2 = 6 * tot_dof;
        for i1 in 0..i2 {
            dsd_u[i1].set_val(0.0);
            ded_u[i1].set_val(0.0);
        }
        i2 = 6 * self.num_nds;
        for i1 in 0..i2 {
            dsd_t[i1].set_val(0.0);
        }
        
        if self.this_type == 41 || self.this_type == 3 {
            if n_lgeom {
                self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
            }
            
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
            for i1 in 0..tot_dof {
                self.get_section_def_dfd1(&mut sec_def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  i1,  MAX_INT);
                
                sect_strn[0].set_val_dfd1(& sec_def[0]);
                tmp.set_val_dfd1(& pre.layer_z[layer]);
                tmp.mult(& sec_def[3]);
                sect_strn[0].add(& tmp);
                
                sect_strn[1].set_val_dfd1(& sec_def[1]);
                tmp.set_val_dfd1(& pre.layer_z[layer]);
                tmp.mult(& sec_def[4]);
                sect_strn[1].sub(& tmp);
                
                sect_strn[2].set_val_dfd1(& sec_def[2]);
                tmp.set_val_dfd1(& pre.layer_z[layer]);
                tmp.mult(& sec_def[5]);
                sect_strn[2].add(& tmp);
                
                tmp.set_val_dfd1(& pre.layer_ang[layer]);
                tmp.neg();
                self.transform_strain_dfd1(&mut d_strain, &mut  sect_strn, &mut  tmp);
                mat_mul_ar_dfd1(&mut d_stress, &mut pre.layer_q[(9*layer)..], &mut  d_strain,  3,  3,  1);
                
                ded_u[i1].set_val_dfd1(& d_strain[0]);
                ded_u[i1 + tot_dof].set_val_dfd1(& d_strain[1]);
                ded_u[i1 + 3 * tot_dof].set_val_dfd1(& d_strain[2]);
                
                dsd_u[i1].set_val_dfd1(& d_stress[0]);
                dsd_u[tot_dof + i1].set_val_dfd1(& d_stress[1]);
                dsd_u[3 * tot_dof + i1].set_val_dfd1(& d_stress[2]);
            }
            mat_mul_ar_dfd1(&mut cte, &mut pre.layer_q[(9*layer)..], &mut pre.layer_te[(3*layer)..],  3,  3,  1);
            cte[0].neg();
            cte[1].neg();
            cte[2].neg();
            mat_mul_ar_dfd1(&mut cten, &mut  cte, &mut  n_vec,  3,  1,  self.num_nds);
            for i1 in 0..self.num_nds {
                dsd_t[i1].set_val_dfd1(& cten[i1]);
                dsd_t[i1 + self.num_nds].set_val_dfd1(& cten[i1 + self.num_nds]);
                dsd_t[i1 + 3 * self.num_nds].set_val_dfd1(& cten[i1 + 2 * self.num_nds]);
            }
        }
        else if self.this_type != 2 {
            self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.glob_disp,  0,  60);
            mat_mul_ar_dfd1(&mut ux, &mut  tmp_ar, &mut  d_ndx,  3,  self.n_dim,  3);
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.cmat,  0,  36);
            for i1 in 0..tot_dof {
                self.get_solid_strain_dfd1(&mut d_strain, &mut  ux, &mut  d_ndx, &mut  pre.loc_ori,  i1,  MAX_INT,  n_lgeom);
                mat_mul_ar_dfd1(&mut d_stress, &mut  tmp_ar, &mut  d_strain,  6,  6,  1);
                i3 = i1;
                for i2 in 0..6 {
                    ded_u[i3].set_val_dfd1(& d_strain[i2]);
                    dsd_u[i3].set_val_dfd1(& d_stress[i2]);
                    i3  +=  tot_dof;
                }
            }
            vec_to_ar_dfd1(&mut tmp_ar2, &mut  pre.therm_exp,  0,  6);
            mat_mul_ar_dfd1(&mut cte, &mut  tmp_ar, &mut  tmp_ar2,  6,  6,  1);
            for i1 in 0..6 {
                cte[i1].neg();
            }
            mat_mul_ar_dfd1(&mut tmp_ar, &mut  cte, &mut  n_vec,  6,  1,  self.num_nds);
            ar_to_vec_dfd1(&mut tmp_ar, dsd_t,  0,  60);
        }
        
        return;
    }

    pub fn get_def_frc_mom_dfd1(&mut self, def : &mut [DiffDoub1], frc_mom : &mut [DiffDoub1], spt : &mut [f64], n_lgeom : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut pt_temp = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp_ar = [DiffDoub1::new(); 36];
        
        if self.dof_per_nd != 6 {
            for i1 in 0..9 {
                def[i1].set_val(0.0);
            }
            return;
        }
        
        if n_lgeom {
            self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
        }
        
        self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        self.get_section_def_dfd1(def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  MAX_INT,  MAX_INT);
        
        pt_temp.set_val(0.0);
        for i1 in 0..self.num_nds {
            tmp.set_val_dfd1(& pre.glob_temp[i1]);
            tmp.mult(& n_vec[i1]);
            pt_temp.add(& tmp);
        }
        
        vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.cmat,  0,  36);
        mat_mul_ar_dfd1(frc_mom, &mut  tmp_ar, def,  self.def_dim,  self.def_dim,  1);
        
        for i1 in 0..6 {
            frc_mom[i1].sub(& pre.einit[i1]);
            tmp.set_val_dfd1(& pt_temp);
            tmp.mult(& pre.therm_exp[i1]);
            frc_mom[i1].sub(& tmp);
        }
        
        return;
    }

    pub fn d_def_frc_momd_u_dfd1(&mut self, d_defd_u : &mut Vec<DiffDoub1>, d_frc_momd_u : &mut Vec<DiffDoub1>, d_frc_momd_t : &mut Vec<DiffDoub1>, spt : &mut [f64], n_lgeom : bool, pre : &mut DiffDoub1StressPrereq) {
        let mut i2 : usize;
        let tot_dof : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut def = [DiffDoub1::new(); 9];
        let mut tmp_ar = [DiffDoub1::new(); 60];
        let mut tmp_ar2 = [DiffDoub1::new(); 6];
        
        if n_lgeom {
            self.get_inst_ori_dfd1(&mut pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_disp,  2);
        }
        
        self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        tot_dof = self.num_nds * self.dof_per_nd + self.num_int_dof;
        for i1 in 0..tot_dof {
            self.get_section_def_dfd1(&mut def, &mut  pre.glob_disp, &mut  pre.inst_ori, &mut  pre.loc_ori, &mut  pre.glob_nds, &mut  d_ndx, &mut  n_vec,  n_lgeom,  i1,  MAX_INT);
            i2 = i1;
            for i3 in 0..self.def_dim {
                d_defd_u[i2].set_val_dfd1(& def[i3]);
                i2  +=  tot_dof;
            }
        }
        
        mat_mul_dfd1(d_frc_momd_u, &mut  pre.cmat, d_defd_u,  self.def_dim,  self.def_dim,  tot_dof);
        
        vec_to_ar_dfd1(&mut tmp_ar2, &mut  pre.therm_exp,  0,  6);
        mat_mul_ar_dfd1(&mut tmp_ar, &mut  tmp_ar2, &mut  n_vec,  6,  1,  self.num_nds);
        ar_to_vec_dfd1(&mut tmp_ar, d_frc_momd_t,  0,  60);
        
        return;
    }

    pub fn get_flux_tgrad_dfd1(&mut self, flux : &mut [DiffDoub1], t_grad : &mut [DiffDoub1], spt : &mut [f64], layer : usize, pre : &mut DiffDoub1StressPrereq) {
        let i1 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut tmp_ar = [DiffDoub1::new(); 10];
        
        self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.glob_temp,  0,  10);
        mat_mul_ar_dfd1(t_grad, &mut  tmp_ar, &mut  d_ndx,  1,  self.num_nds,  3);
        if self.this_type == 41 || self.this_type == 3 {
            i1 = 9 * layer;
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.layer_tc,  i1,  i1 + 9);
            mat_mul_ar_dfd1(flux, &mut  tmp_ar, t_grad,  3,  3,  1);
        }
        else {
            vec_to_ar_dfd1(&mut tmp_ar, &mut  pre.tcmat,  0,  9);
            mat_mul_ar_dfd1(flux, &mut  tmp_ar, t_grad,  3,  3,  1);
        }
        flux[0].neg();
        flux[1].neg();
        flux[2].neg();
        
        return;
    }

    pub fn d_flux_tgradd_t_dfd1(&mut self, d_fd_t : &mut Vec<DiffDoub1>, d_tg : &mut Vec<DiffDoub1>, spt : &mut [f64], layer : usize, pre : &mut DiffDoub1StressPrereq) {
        let i1 : usize;
        let i2 : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut tmp_ar1 = [DiffDoub1::new(); 33];
        let mut tmp_ar2 = [DiffDoub1::new(); 9];
        let mut tmp_ar3 = [DiffDoub1::new(); 30];
        
        self.get_ip_data_dfd1(&mut n_vec, &mut  d_ndx, &mut  det_j, &mut  pre.loc_nds, spt);
        transpose_ar_dfd1(&mut tmp_ar1, &mut  d_ndx,  self.num_nds,  3);
        ar_to_vec_dfd1(&mut tmp_ar1, d_tg,  0,  33);
        if self.this_type == 41 || self.this_type == 3 {
            i1 = 9 * layer;
            vec_to_ar_dfd1(&mut tmp_ar2, &mut  pre.layer_tc,  i1,  i1 + 9);
            mat_mul_ar_dfd1(&mut tmp_ar3, &mut  tmp_ar2, &mut  tmp_ar1,  3,  3,  self.num_nds);
            ar_to_vec_dfd1(&mut tmp_ar3, d_fd_t,  0,  30);
        }
        else {
            mat_mul_dfd1(d_fd_t, &mut  pre.tcmat, d_tg,  3,  3,  self.num_nds);
        }
        i2 = 3 * self.num_nds;
        for i1 in 0..i2 {
            d_fd_t[i1].neg();
        }
        
        return;
    }

    pub fn get_cent_v_grad_dfd1(&mut self, v_grad : &mut [DiffDoub1], t_grad : &mut [DiffDoub1], pre : &mut DiffDoub1FlPrereq, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable> ) {
        // get fluid prerequisites before calling
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut spt = [0f64; 3];
       
        for i in 0..3 {
            spt[i] = self.s_cent[i];
        }

        // self.get_nd_crds_dfd1(&mut pre.glob_nds, nd_ar, dv_ar);
        // self.get_nd_disp_dfd1(&mut pre.glob_disp, nd_ar);
        // let i2 = 3 * self.num_nds;
        // for i1 in 0..i2 {
        //     pre.glob_nds[i1].add(& pre.glob_disp[i1]);
        // }
        // self.get_nd_fl_vel_dfd1(&mut pre.fl_vel, nd_ar);
        // self.get_nd_temp_dfd1(&mut pre.fl_temp, nd_ar);

        self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &mut spt);
        mat_mul_ar_dfd1(v_grad, &pre.fl_vel, &d_ndx, 3, self.num_nds, 3);
        mat_mul_ar_dfd1(t_grad, &pre.fl_temp, &d_ndx, 1, self.num_nds, 3);
        
    }

    pub fn put_vec_to_glob_mat_dfd1(&mut self, q_mat : &mut SparseMat, el_qvec : &mut Vec<DiffDoub1>, for_therm : bool, mat_row : usize, nd_ar : &mut Vec<Node>) {
        let nd_dof : usize =  self.num_nds*self.dof_per_nd;
        let tot_dof : usize =  nd_dof + self.num_int_dof;
        let mut nd : usize;
        let mut dof : usize;
        let mut glob_ind : usize;
        
        if for_therm {
            for i1 in 0..self.num_nds {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                q_mat.add_entry(mat_row,   glob_ind,   el_qvec[i1].val);
            }
        }
        else {
            for i1 in 0..tot_dof {
                if i1 < nd_dof {
                    nd = self.nodes[self.dof_table[2 * i1]];
                    dof = self.dof_table[2 * i1 + 1];
                    glob_ind = nd_ar[nd].dof_index[dof];
                    q_mat.add_entry(mat_row,  glob_ind,  el_qvec[i1].val);
                }
                else {
                    dof = i1 - nd_dof;
                    glob_ind = self.int_dof_index + dof;
                    q_mat.add_entry(mat_row,   glob_ind,   el_qvec[i1].val);
                }
            }
        }
        
        return;
    }

    //end dup
 
//end skip 
 
 
 
 
 
    pub fn get_el_vec(&mut self, el_vec : &mut Vec<f64>, glob_vec : &mut Vec<f64>, for_therm : bool, intnl : bool, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut nd : usize;
        let mut dof : usize;
        let mut glob_ind : usize;
        let nd_dof : usize;
        
        if for_therm {
            for i1 in 0..self.num_nds {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                el_vec[i1] = glob_vec[glob_ind];
            }
        }
        else {
            nd_dof = self.num_nds * self.dof_per_nd;
            i2 = 0;
            for i1 in 0..nd_dof {
                nd = self.nodes[self.dof_table[i2]];
                dof = self.dof_table[i2 + 1];
                glob_ind = nd_ar[nd].dof_index[dof];
                el_vec[i1] = glob_vec[glob_ind];
                i2  +=  2;
            }
            if intnl {
                for i1 in 0..self.num_int_dof {
                    el_vec[i1 + nd_dof] = glob_vec[i1 + self.int_dof_index];
                }
            }
        }
        
        return;
    }

    pub fn add_to_glob_vec(&mut self, el_vec : &mut Vec<f64>, glob_vec : &mut Vec<f64>, for_therm : bool, intnl : bool, nd_ar : &mut Vec<Node>) {
        let mut i2 : usize;
        let mut nd : usize;
        let mut dof : usize;
        let mut glob_ind : usize;
        let nd_dof : usize;
        
        if for_therm {
            for i1 in 0..self.num_nds {
                nd = self.nodes[i1];
                glob_ind = nd_ar[nd].sorted_rank;
                glob_vec[glob_ind]  +=  el_vec[i1];
            }
        }
        else {
            nd_dof = self.num_nds * self.dof_per_nd;
            i2 = 0;
            for i1 in 0..nd_dof {
                nd = self.nodes[self.dof_table[i2]];
                dof = self.dof_table[i2 + 1];
                glob_ind = nd_ar[nd].dof_index[dof];
                glob_vec[glob_ind]  +=  el_vec[i1];
                i2  +=  2;
            }
            if intnl {
                for i1 in 0..self.num_int_dof {
                    glob_vec[i1 + self.int_dof_index]  +=  el_vec[i1 + nd_dof];
                }
            }
        }
        
        return;
    }

}


