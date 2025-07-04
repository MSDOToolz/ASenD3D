use crate::element::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::list_ent::*;
use crate::section::*;
use crate::face::*;

//dup1

impl DiffDoub0StressPrereq {
    pub fn allocate_layers_dfd0(&mut self, num_layers : usize) {
        if num_layers != 0 {
            self.layer_z = vec![DiffDoub0::new(); num_layers];
            self.layer_thk = vec![DiffDoub0::new(); num_layers];
            self.layer_ang = vec![DiffDoub0::new(); num_layers];
            self.layer_q = vec![DiffDoub0::new(); 9 * num_layers];
            self.layer_d = vec![DiffDoub0::new(); 9 * num_layers];
            self.layer_te = vec![DiffDoub0::new(); 3 * num_layers];
            self.layer_e0 = vec![DiffDoub0::new(); 3 * num_layers];
            self.layer_den = vec![DiffDoub0::new(); num_layers];
            self.layer_tc = vec![DiffDoub0::new(); 9 * num_layers];
            self.layer_sh = vec![DiffDoub0::new(); num_layers];
            self.layer_de = vec![DiffDoub0::new(); 3 * num_layers];
            self.layer_diff = vec![DiffDoub0::new(); 9 * num_layers];
            self.layer_max_con = vec![DiffDoub0::new(); num_layers];
        }
        self.current_lay_len = num_layers;
        return;
    }

}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

impl DiffDoub1StressPrereq {
    pub fn allocate_layers_dfd1(&mut self, num_layers : usize) {
        if num_layers != 0 {
            self.layer_z = vec![DiffDoub1::new(); num_layers];
            self.layer_thk = vec![DiffDoub1::new(); num_layers];
            self.layer_ang = vec![DiffDoub1::new(); num_layers];
            self.layer_q = vec![DiffDoub1::new(); 9 * num_layers];
            self.layer_d = vec![DiffDoub1::new(); 9 * num_layers];
            self.layer_te = vec![DiffDoub1::new(); 3 * num_layers];
            self.layer_e0 = vec![DiffDoub1::new(); 3 * num_layers];
            self.layer_den = vec![DiffDoub1::new(); num_layers];
            self.layer_tc = vec![DiffDoub1::new(); 9 * num_layers];
            self.layer_sh = vec![DiffDoub1::new(); num_layers];
            self.layer_de = vec![DiffDoub1::new(); 3 * num_layers];
            self.layer_diff = vec![DiffDoub1::new(); 9 * num_layers];
            self.layer_max_con = vec![DiffDoub1::new(); num_layers];
        }
        self.current_lay_len = num_layers;
        return;
    }

}

//end dup
 
//end skip 
 
 
 
 
 
 
 
impl Element {
    pub fn initialize_type(&mut self, new_type : usize) {
        self.this_type = new_type;
        let mut i1 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        
        self.s_cent[0] = 0.0;
        self.s_cent[1] = 0.0;
        self.s_cent[2] = 0.0;
        if self.this_type == 4 || self.this_type == 400 {
            self.num_nds = 4;
            self.dof_per_nd = 3;
            self.num_int_dof = 0;
            self.n_dim = 4;
            self.def_dim = 6;
            self.num_ip = 1;
            self.num_faces = 4;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; i1 * 2];//int[i1][2];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = 0.25;
            self.int_pts[1] = 0.25;
            self.int_pts[2] = 0.25;
            self.ip_wt[0] = R_1O6;
            self.nd_spts = vec![0f64; 12];
            self.nd_spts[0] = 0.0;
            self.nd_spts[1] = 0.0;
            self.nd_spts[2] = 0.0;
            self.nd_spts[3] = 1.0;
            self.nd_spts[4] = 0.0;
            self.nd_spts[5] = 0.0;
            self.nd_spts[6] = 0.0;
            self.nd_spts[7] = 1.0;
            self.nd_spts[8] = 0.0;
            self.nd_spts[9] = 0.0;
            self.nd_spts[10] = 0.0;
            self.nd_spts[11] = 1.0;
        } else if self.this_type == 6 || self.this_type == 600 {
            self.num_nds = 6;
            self.dof_per_nd = 3;
            self.num_int_dof = 0;
            self.n_dim = 6;
            self.def_dim = 6;
            self.num_ip = 2;
            self.num_faces = 5;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; i1 * 2];// int[i1][2];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = R_1O3;
            self.int_pts[1] = R_1O3;
            self.int_pts[2] = -R_1ORT3;
            self.int_pts[3] = R_1O3;
            self.int_pts[4] = R_1O3;
            self.int_pts[5] = R_1ORT3;
            self.ip_wt[0] = 0.5;
            self.ip_wt[1] = 0.5;
            self.s_cent[0] = R_1O3;
            self.s_cent[1] = R_1O3;
            self.s_cent[2] = 0.0;
            self.nd_spts = vec![0f64; 18];
            self.nd_spts[0] = 0.0;
            self.nd_spts[1] = 0.0;
            self.nd_spts[2] = -1.0;
            self.nd_spts[3] = 1.0;
            self.nd_spts[4] = 0.0;
            self.nd_spts[5] = -1.0;
            self.nd_spts[6] = 0.0;
            self.nd_spts[7] = 1.0;
            self.nd_spts[8] = -1.0;
            self.nd_spts[9] = 0.0;
            self.nd_spts[10] = 0.0;
            self.nd_spts[11] = 1.0;
            self.nd_spts[12] = 1.0;
            self.nd_spts[13] = 0.0;
            self.nd_spts[14] = 1.0;
            self.nd_spts[15] = 0.0;
            self.nd_spts[16] = 1.0;
            self.nd_spts[17] = 1.0;
        } else if self.this_type == 8 || self.this_type == 81 || self.this_type == 800 {
            self.num_nds = 8;
            self.dof_per_nd = 3;
            if self.this_type == 8 || self.this_type == 800 {
                self.num_int_dof = 0;
                self.n_dim = 8;
            } else {
                self.num_int_dof = 9;
                self.n_dim = 11;
            }
            self.def_dim = 6;
            self.num_ip = 8;
            self.num_faces = 6;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; i1 * 2];//int[i1][2];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            let s_val : [f64; 2] = [-R_1ORT3,R_1ORT3];
            i4 = 0;
            i5 = 0;
            for i1 in 0..2 {
                for i2 in 0..2 {
                    for i3 in 0..2 {
                        self.int_pts[i5] = s_val[i3];
                        i5 += 1usize;
                        self.int_pts[i5] = s_val[i2];
                        i5 += 1usize;
                        self.int_pts[i5] = s_val[i1];
                        i5 += 1usize;
                        self.ip_wt[i4] = 1.0;
                        i4 += 1usize;
                    }
                }
            }
            if self.this_type == 81 {
                i3 = 24;
                for i1 in 8..11 {
                    for i2 in 0..3 {
                        self.dof_table[2*i3] = i1;
                        self.dof_table[2*i3+1] = i2;
                        i3 += 1usize;
                    }
                }
            }
            self.nd_spts = vec![0f64; 24];
            self.nd_spts[0] = -1.0;
            self.nd_spts[1] = -1.0;
            self.nd_spts[2] = -1.0;
            self.nd_spts[3] = 1.0;
            self.nd_spts[4] = -1.0;
            self.nd_spts[5] = -1.0;
            self.nd_spts[6] = 1.0;
            self.nd_spts[7] = 1.0;
            self.nd_spts[8] = -1.0;
            self.nd_spts[9] = -1.0;
            self.nd_spts[10] = 1.0;
            self.nd_spts[11] = -1.0;
            self.nd_spts[12] = -1.0;
            self.nd_spts[13] = -1.0;
            self.nd_spts[14] = 1.0;
            self.nd_spts[15] = 1.0;
            self.nd_spts[16] = -1.0;
            self.nd_spts[17] = 1.0;
            self.nd_spts[18] = 1.0;
            self.nd_spts[19] = 1.0;
            self.nd_spts[20] = 1.0;
            self.nd_spts[21] = -1.0;
            self.nd_spts[22] = 1.0;
            self.nd_spts[23] = 1.0;
        }
        else if self.this_type == 10 || self.this_type == 1000 {
            self.num_nds = 10;
            self.dof_per_nd = 3;
            self.num_int_dof = 0;
            self.n_dim = 10;
            self.def_dim = 6;
            self.num_ip = 4;
            self.num_faces = 4;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds * self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; 2 * i1];
            self.int_pts = vec![0f64; self.num_ip * 3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = R_TET1;
            self.int_pts[1] = R_TET1;
            self.int_pts[2] = R_TET1;
            self.int_pts[3] = R_TET2;
            self.int_pts[4] = R_TET1;
            self.int_pts[5] = R_TET1;
            self.int_pts[6] = R_TET1;
            self.int_pts[7] = R_TET2;
            self.int_pts[8] = R_TET1;
            self.int_pts[9] = R_TET1;
            self.int_pts[10] = R_TET1;
            self.int_pts[11] = R_TET2;
            self.ip_wt[0] = R_1O24;
            self.ip_wt[1] = R_1O24;
            self.ip_wt[2] = R_1O24;
            self.ip_wt[3] = R_1O24;
            self.s_cent[0] = 0.25;
            self.s_cent[1] = 0.25;
            self.s_cent[2] = 0.25;
            self.nd_spts = vec![0f64; 30];
            self.nd_spts[0] = 0.0;
            self.nd_spts[1] = 0.0;
            self.nd_spts[2] = 0.0;
            self.nd_spts[3] = 1.0;
            self.nd_spts[4] = 0.0;
            self.nd_spts[5] = 0.0;
            self.nd_spts[6] = 0.0;
            self.nd_spts[7] = 1.0;
            self.nd_spts[8] = 0.0;
            self.nd_spts[9] = 0.0;
            self.nd_spts[10] = 0.0;
            self.nd_spts[11] = 1.0;
            self.nd_spts[12] = 0.5;
            self.nd_spts[13] = 0.0;
            self.nd_spts[14] = 0.0;
            self.nd_spts[15] = 0.5;
            self.nd_spts[16] = 0.5;
            self.nd_spts[17] = 0.0;
            self.nd_spts[18] = 0.0;
            self.nd_spts[19] = 0.5;
            self.nd_spts[20] = 0.0;
            self.nd_spts[21] = 0.0;
            self.nd_spts[22] = 0.0;
            self.nd_spts[23] = 0.5;
            self.nd_spts[24] = 0.5;
            self.nd_spts[25] = 0.0;
            self.nd_spts[26] = 0.5;
            self.nd_spts[27] = 0.0;
            self.nd_spts[28] = 0.5;
            self.nd_spts[29] = 0.5;
        }
        else if self.this_type == 3 {
            self.num_nds = 3;
            self.dof_per_nd = 6;
            self.num_int_dof = 3;
            self.n_dim = 6;
            self.def_dim = 9;
            self.num_ip = 3;
            self.num_faces = 2;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; 2*i1];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = R_1O6;
            self.int_pts[1] = R_1O6;
            self.int_pts[2] = 0.0;
            self.int_pts[3] = R_2O3;
            self.int_pts[4] = R_1O6;
            self.int_pts[5] = 0.0;
            self.int_pts[6] = R_1O6;
            self.int_pts[7] = R_2O3;
            self.int_pts[8] = 0.0;
            self.ip_wt[0] = R_1O6;
            self.ip_wt[1] = R_1O6;
            self.ip_wt[2] = R_1O6;
            self.dof_table[36] = 3;
            self.dof_table[37] = 2;
            self.dof_table[38] = 4;
            self.dof_table[39] = 2;
            self.dof_table[40] = 5;
            self.dof_table[41] = 2;
            self.s_cent[0] = R_1O3;
            self.s_cent[1] = R_1O3;
            self.s_cent[2] = R_1O3;
        } else if self.this_type == 41 {
            self.num_nds = 4;
            self.dof_per_nd = 6;
            self.num_int_dof = 8;
            self.n_dim = 10;
            self.def_dim = 9;
            self.num_ip = 4;
            self.num_faces = 2;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; 2*i1];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = -R_1ORT3;
            self.int_pts[1] = -R_1ORT3;
            self.int_pts[2] = 0.0;
            self.int_pts[3] = R_1ORT3;
            self.int_pts[4] = -R_1ORT3;
            self.int_pts[5] = 0.0;
            self.int_pts[6] = -R_1ORT3;
            self.int_pts[7] = R_1ORT3;
            self.int_pts[8] = 0.0;
            self.int_pts[9] = R_1ORT3;
            self.int_pts[10] = R_1ORT3;
            self.int_pts[11] = 0.0;
            self.ip_wt[0] = 1.0;
            self.ip_wt[1] = 1.0;
            self.ip_wt[2] = 1.0;
            self.ip_wt[3] = 1.0;
            self.dof_table[48] = 4;
            self.dof_table[49] = 0;
            self.dof_table[50] = 5;
            self.dof_table[51] = 0;
            self.dof_table[52] = 4;
            self.dof_table[53] = 1;
            self.dof_table[54] = 5;
            self.dof_table[55] = 1;
            self.dof_table[56] = 6;
            self.dof_table[57] = 2;
            self.dof_table[58] = 7;
            self.dof_table[59] = 2;
            self.dof_table[60] = 8;
            self.dof_table[61] = 2;
            self.dof_table[62] = 9;
            self.dof_table[63] = 2;
        } else if self.this_type == 2 {
            self.num_nds = 2;
            self.dof_per_nd = 6;
            self.num_int_dof = 2;
            self.n_dim = 3;
            self.def_dim = 6;
            self.num_ip = 2;
            self.num_faces = 0;
            self.nodes = vec![0usize; self.num_nds];
            i1 = self.num_nds*self.dof_per_nd + self.num_int_dof;
            self.dof_table = vec![0usize; i1*2];
            self.int_pts = vec![0f64; self.num_ip*3];
            self.ip_wt = vec![0f64; self.num_ip];
            self.int_pts[0] = -R_1ORT3;
            self.int_pts[1] = 0.0;
            self.int_pts[2] = 0.0;
            self.int_pts[3] = R_1ORT3;
            self.int_pts[4] = 0.0;
            self.int_pts[5] = 0.0;
            self.ip_wt[0] = 1.0;
            self.ip_wt[1] = 1.0;
            self.dof_table[24] = 2;
            self.dof_table[25] = 1;
            self.dof_table[26] = 2;
            self.dof_table[27] = 2;
        }
        else if self.this_type == 21 {//force field
            self.num_nds = 2;
            self.dof_per_nd = 3;
            self.num_int_dof = 0;
            self.n_dim = 2;
            self.def_dim = 0;
            self.num_ip = 0;
            self.num_faces = 0;
            self.nodes = vec![0usize; self.num_nds];
            self.dof_table = vec![0usize; 12];
            self.int_pts = vec![0f64; 1];
            self.ip_wt = vec![0f64; 1];
        }
        else if self.this_type == 1 {// mass
            self.num_nds = 1;
            self.dof_per_nd = 3;
            self.num_int_dof = 0;
            self.n_dim = 1;
            self.def_dim = 0;
            self.num_ip = 0;
            self.num_faces = 0;
            self.nodes = vec![0usize; self.num_nds];
            self.dof_table = vec![0usize; 6];
            self.int_pts = vec![0f64; 1];
            self.int_pts = vec![0f64; 1];
        }
        
        i3 = 0;
        for i1 in 0..self.num_nds {
            for i2 in 0..self.dof_per_nd {
                self.dof_table[i3] = i1;
                self.dof_table[i3+1] = i2;
                i3 +=  2;
            }
        }
        
        if self.num_int_dof != 0 {
            self.internal_disp = vec![0f64; self.num_int_dof];
            self.int_prev_disp = vec![0f64; self.num_int_dof];
            self.internald_ldu = vec![0f64; self.num_int_dof];
            self.internal_adj = vec![0f64; self.num_int_dof];
            self.internal_ru = vec![DiffDoub1::new(); self.num_int_dof];
            i1 = (self.num_nds*self.dof_per_nd + self.num_int_dof)*self.num_int_dof;
            self.internal_mat = vec![0f64; i1];
        }
        
        self.int_dof_index = 0;
        
        self.sect_ptr = MAX_INT;
        
        return;
    }

    pub fn set_nodes(&mut self, new_nds : &mut [usize]) {
        for i1 in 0..self.num_nds {
            self.nodes[i1] = new_nds[i1];
        }
        return;
    }

    pub fn initialize_faces(&mut self, glob_fc_lst : &mut Vec<Face>, fi : &mut usize) {
        //fi = the current number of self.faces that have been written into glob_fc_lst
        let mut new_fc : &mut Face;
        if self.this_type == 4 || self.this_type == 400 {
            glob_fc_lst[*fi].num_nds = 3;
            new_fc  = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 1, self.nodes[1]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 3, self.nodes[3]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 1, self.nodes[1]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 3, self.nodes[3]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 3, self.nodes[3]);
            new_fc.set_node(2, 2, self.nodes[2]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        } else if self.this_type == 6 || self.this_type == 600 {
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 1, self.nodes[1]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 3, self.nodes[3]);
            new_fc.set_node(1, 4, self.nodes[4]);
            new_fc.set_node(2, 5, self.nodes[5]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 4, self.nodes[4]);
            new_fc.set_node(3, 3, self.nodes[3]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 1, self.nodes[1]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 5, self.nodes[5]);
            new_fc.set_node(3, 4, self.nodes[4]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 3, self.nodes[3]);
            new_fc.set_node(2, 5, self.nodes[5]);
            new_fc.set_node(3, 2, self.nodes[2]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        } else if self.this_type == 8 || self.this_type == 81 || self.this_type == 800 {
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 3, self.nodes[3]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 1, self.nodes[1]);
            new_fc.set_node(3, 0, self.nodes[0]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 4, self.nodes[4]);
            new_fc.set_node(1, 5, self.nodes[5]);
            new_fc.set_node(2, 6, self.nodes[6]);
            new_fc.set_node(3, 7, self.nodes[7]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 5, self.nodes[5]);
            new_fc.set_node(3, 4, self.nodes[4]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 1, self.nodes[1]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 6, self.nodes[6]);
            new_fc.set_node(3, 5, self.nodes[5]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 2, self.nodes[2]);
            new_fc.set_node(1, 3, self.nodes[3]);
            new_fc.set_node(2, 7, self.nodes[7]);
            new_fc.set_node(3, 6, self.nodes[6]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 3, self.nodes[3]);
            new_fc.set_node(1, 0, self.nodes[0]);
            new_fc.set_node(2, 4, self.nodes[4]);
            new_fc.set_node(3, 7, self.nodes[7]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        }
        else if self.this_type == 10 || self.this_type == 1000 {
            glob_fc_lst[*fi].num_nds = 6;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0,  0,  self.nodes[0]);
            new_fc.set_node(1,  2,  self.nodes[2]);
            new_fc.set_node(2,  1,  self.nodes[1]);
            new_fc.set_node(3,  6,  self.nodes[6]);
            new_fc.set_node(4,  5,  self.nodes[5]);
            new_fc.set_node(5,  4,  self.nodes[4]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 6;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0,  0,  self.nodes[0]);
            new_fc.set_node(1,  1,  self.nodes[1]);
            new_fc.set_node(2,  3,  self.nodes[3]);
            new_fc.set_node(3,  4,  self.nodes[4]);
            new_fc.set_node(4,  8,  self.nodes[8]);
            new_fc.set_node(5,  7,  self.nodes[7]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 6;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0,  1,  self.nodes[1]);
            new_fc.set_node(1,  2,  self.nodes[2]);
            new_fc.set_node(2,  3,  self.nodes[3]);
            new_fc.set_node(3,  5,  self.nodes[5]);
            new_fc.set_node(4,  9,  self.nodes[9]);
            new_fc.set_node(5,  8,  self.nodes[8]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 6;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0,  0,  self.nodes[0]);
            new_fc.set_node(1,  3,  self.nodes[3]);
            new_fc.set_node(2,  2,  self.nodes[2]);
            new_fc.set_node(3,  7,  self.nodes[7]);
            new_fc.set_node(4,  9,  self.nodes[9]);
            new_fc.set_node(5,  6,  self.nodes[6]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        }
        else if self.this_type == 3 {
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 2, self.nodes[2]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 3;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 2, self.nodes[2]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 0, self.nodes[0]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        } else if self.this_type == 41 {
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 0, self.nodes[0]);
            new_fc.set_node(1, 1, self.nodes[1]);
            new_fc.set_node(2, 2, self.nodes[2]);
            new_fc.set_node(3, 3, self.nodes[3]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
            glob_fc_lst[*fi].num_nds = 4;
            new_fc = &mut glob_fc_lst[*fi];
            new_fc.set_node(0, 3, self.nodes[3]);
            new_fc.set_node(1, 2, self.nodes[2]);
            new_fc.set_node(2, 1, self.nodes[1]);
            new_fc.set_node(3, 0, self.nodes[0]);
            new_fc.host_el = self.label;
            self.faces.push_back(*fi);
            *fi += 1usize;
        }
        
        return;
    }

    pub fn set_int_disp(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..self.num_int_dof {
            self.internal_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn set_int_prev_disp(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..self.num_int_dof {
            self.int_prev_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn advance_int_disp(&mut self) {
        for i1 in 0..self.num_int_dof {
            self.int_prev_disp[i1] = self.internal_disp[i1];
        }
        return;
    }

    pub fn backstep_int_disp(&mut self) {
        for i1 in 0..self.num_int_dof {
            self.internal_disp[i1] = self.int_prev_disp[i1];
        }
        return;
    }

    pub fn set_intd_ld_u(&mut self, globd_ld_u : &mut Vec<f64>) {
        let mut i2 : usize =  self.int_dof_index;
        for i1 in 0..self.num_int_dof {
            self.internald_ldu[i1] = globd_ld_u[i2];
            i2 += 1usize;
        }
        return;
    }

    pub fn get_num_layers(&mut self, sec_lst : &mut Vec<Section>) -> usize {
        return  sec_lst[self.sect_ptr].layers.len();
    }

    pub fn add_design_variable(&mut self, d_index : usize, coef : f64) {
        let mut dv = IDCapsule::new();
        dv.int_dat = d_index;
        dv.doub_dat = coef;
        self.design_vars.push_back(dv);
        return;
    }

    pub fn add_comp_dvar(&mut self, d_index : usize) {
        for dv in self.comp_dvars.iter_mut() {
            if *dv == d_index {
                return;
            }
        }
        self.comp_dvars.push_back(d_index);
        return;
    }

}


