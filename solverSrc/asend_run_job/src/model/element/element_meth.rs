use crate::model::element::*;
use crate::constants::*;
use crate::model::section::*;
use crate::model::face::*;

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
        let i1 : usize;

        self.nodes = vec![0usize; self.num_nds()];
        
        if self.num_int_dof() != 0 {
            self.internal_disp = vec![0f64; self.num_int_dof()];
            self.int_prev_disp = vec![0f64; self.num_int_dof()];
            self.internald_ldu = vec![0f64; self.num_int_dof()];
            self.internal_adj = vec![0f64; self.num_int_dof()];
            self.internal_ru = vec![DiffDoub1::new(); self.num_int_dof()];
            i1 = (self.num_nds()*self.dof_per_nd() + self.num_int_dof())*self.num_int_dof();
            self.internal_mat = vec![0f64; i1];
        }
        
        self.int_dof_index = 0;
        
        self.sect_ptr = MAX_INT;
        
        return;
    }

    pub fn num_nds(&self) -> usize {
        match self.this_type {
            21 => 2,
            41 => 4,
            81 => 8,
            _ => self.this_type,
        }
    }

    pub fn dof_per_nd(&self) -> usize {
        match self.this_type {
            2 => 6,
            3 => 6,
            41 => 6,
            _ => 3,
        }
    }

    pub fn n_dim(&self) -> usize {
        match self.this_type {
            2 => 3,
            21 => 2,
            3 => 6,
            41 => 10,
            81 => 11,
            _ => self.this_type,
        }
    }

    pub fn num_int_dof(&self) -> usize {
        match self.this_type {
            2 => 2,
            3 => 3,
            41 => 8,
            81 => 9,
            _ => 0,
        }
    }

    pub fn num_faces(&self) -> usize {
        match self.this_type {
            3 => 2,
            41 => 2,
            4 => 4,
            6 => 5,
            8 => 6,
            81 => 6,
            10 => 4,
            _ => 0,
        }
    }

    pub fn num_ip(&self) -> usize {
        match self.this_type {
            2 => 2,
            3 => 3,
            41 => 4,
            4 => 1,
            6 => 2,
            8 => 8,
            81 => 8,
            10 => 4,
            _ => 0
        }
    }

    pub fn def_dim(&self) -> usize {
        match self.this_type {
            3 => 9,
            41 => 9,
            _ => 6,
        }
    }

    pub fn dof_table(&self, ind : usize) -> usize {
        let dofi : usize;
        let col = match ind.checked_rem(2) {None => 0, Some(x) => x};
        match self.this_type {
            2 => {if ind < 24 {
                       dofi = ind/2;
                       return match col {
                           0 => dofi/6,
                           _ => match dofi.checked_rem(6) {None => 0, Some(x) => x},
                       };
                   }
                   else {
                       return match ind {
                           24 => 2,
                           25 => 1,
                           26 => 2,
                           _ => 2,
                       };
                   }},
            3 => {if ind < 36 {
                      dofi = ind/2;
                      return match col {
                          0 => dofi/6,
                          _ => match dofi.checked_rem(6) {None => 0, Some(x) => x},
                      };
                  }
                  else {
                      return match ind {
                          36 => 3,
                          37 => 2,
                          38 => 4,
                          39 => 2,
                          40 => 5,
                          _ => 2,
                      }
                  }},
            41 => {if ind < 48 {
                       dofi = ind/2;
                       return match col {
                           0 => dofi/6,
                           _ => match dofi.checked_rem(6) {None => 0, Some(x) => x},
                       };
                   }
                   else {
                       return match ind {
                           48 => 4,
                           49 => 0,
                           50 => 5,
                           51 => 0,
                           52 => 4,
                           53 => 1,
                           54 => 5,
                           55 => 1,
                           56 => 6,
                           57 => 2,
                           58 => 7,
                           59 => 2,
                           60 => 8,
                           61 => 2,
                           62 => 9,
                           _ => 2,
                       }
                   }},
            _ => {dofi = ind/2;
                  return match col {
                      0 => dofi/3,
                      _ => match dofi.checked_rem(3) {None => 0, Some(x) => x},
                  };},
        }
    }

    pub fn ip_crd(&self, crd : &mut [f64], ind : usize) {
        let xrow : usize;
        let yrow : usize;
        let zrow : usize;
        if self.this_type == 8 || self.this_type == 81 {
            xrow = match ind.checked_rem(2) {
                None => 0,
                Some(x) => x,
            };
            crd[0] = match xrow {
                0 => -R_1ORT3,
                _ => R_1ORT3,
            };
            yrow = match (ind/2).checked_rem(2) {
                None => 0,
                Some(x) => x,
            };
            crd[1] = match yrow {
                0 => -R_1ORT3,
                _ => R_1ORT3,
            };
            zrow = ind/4;
            crd[2] = match zrow {
                0 => -R_1ORT3,
                _ => R_1ORT3,
            };
            return;
        }
        match self.this_type {
            2 => {match ind {
                      0 => {crd[0] = -R_1ORT3;
                            crd[1] = 0.0;
                            crd[2] = 0.0;},
                      _ => {crd[0] = R_1ORT3;
                            crd[1] = 0.0;
                            crd[2] = 0.0;},
                  }},
            3 => {match ind {
                      0 => {crd[0] = R_1O6;
                            crd[1] = R_1O6;
                            crd[2] = 0.0;},
                      1 => {crd[0] = R_2O3;
                            crd[1] = R_1O6;
                            crd[2] = 0.0;}
                      _ => {crd[0] = R_1O6;
                            crd[1] = R_2O3;
                            crd[2] = 0.0;},
            }},
            41 => {xrow = match ind.checked_rem(2) {
                       None => 0,
                       Some(x) => x,
                   };
                   crd[0] = match xrow {
                       0 => -R_1ORT3,
                       _ => R_1ORT3,
                   };
                   yrow = ind/2;
                   crd[1] = match yrow {
                       0 => -R_1ORT3,
                       _ => R_1ORT3,
                   };
                   crd[2] = 0.0;},
            4 => {crd[0] = 0.25;
                  crd[1] = 0.25;
                  crd[2] = 0.25;},
            6 => {match ind {
                      0 => {crd[0] = R_1O3;
                            crd[1] = R_1O3;
                            crd[2] = -R_1ORT3;},
                      _ => {crd[0] = R_1O3;
                            crd[1] = R_1O3;
                            crd[2] = R_1ORT3;},
                   }},
            10 => {match ind {
                       0 => {crd[0] = R_TET1;
                             crd[1] = R_TET1;
                             crd[2] = R_TET1;},
                       1 => {crd[0] = R_TET2;
                             crd[1] = R_TET1;
                             crd[2] = R_TET1;},
                       2 => {crd[0] = R_TET1;
                             crd[1] = R_TET2;
                             crd[2] = R_TET1;},
                       _ => {crd[0] = R_TET1;
                             crd[1] = R_TET1;
                             crd[2] = R_TET2;},
                   }},
            _ => (),
        }
    }

    pub fn cent_s_crd(&self, crd : &mut [f64]) {
        match self.this_type {
            3 => {crd[0] = R_1O3;
                  crd[1] = R_1O3;
                  crd[2] = 0.0;},
            4 => {crd[0] = 0.25;
                  crd[1] = 0.25;
                  crd[2] = 0.25;},
            6 => {crd[0] = R_1O3;
                  crd[1] = R_1O3;
                  crd[2] = 0.0;},
            10 => {crd[0] = 0.25;
                   crd[1] = 0.25;
                   crd[2] = 0.25;},
            _ => {crd[0] = 0.0;
                  crd[1] = 0.0;
                  crd[2] = 0.0;},
        }
    }

    pub fn ip_wt(&self, ind : usize) -> f64 {
        match self.this_type {
            2 => 1.0,
            3 => R_1O6,
            41 => 1.0,
            4 => R_1O6,
            6 => 0.5,
            8 => 1.0,
            81 => 1.0,
            10 => R_1O24,
            _ => 0.0,
        }
    }

    pub fn set_nodes(&mut self, new_nds : &mut [usize]) {
        for i1 in 0..self.num_nds() {
            self.nodes[i1] = new_nds[i1];
        }
        return;
    }

    pub fn initialize_faces(&mut self, glob_fc_lst : &mut Vec<Face>, fi : &mut usize) {
        //fi = the current number of self.faces that have been written into glob_fc_lst
        let mut new_fc : &mut Face;
        if self.this_type == 4 {
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
        } else if self.this_type == 6 {
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
        } else if self.this_type == 8 || self.this_type == 81 {
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
        else if self.this_type == 10 {
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
        for i1 in 0..self.num_int_dof() {
            self.internal_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn set_int_prev_disp(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..self.num_int_dof() {
            self.int_prev_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn advance_int_disp(&mut self) {
        for i1 in 0..self.num_int_dof() {
            self.int_prev_disp[i1] = self.internal_disp[i1];
        }
        return;
    }

    pub fn backstep_int_disp(&mut self) {
        for i1 in 0..self.num_int_dof() {
            self.internal_disp[i1] = self.int_prev_disp[i1];
        }
        return;
    }

    pub fn set_intd_ld_u(&mut self, globd_ld_u : &mut Vec<f64>) {
        let mut i2 : usize =  self.int_dof_index;
        for i1 in 0..self.num_int_dof() {
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


