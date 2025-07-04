use crate::node::*;
use crate::list_ent::*;
use crate::diff_doub::*;
use crate::design_var::*;
use crate::cpp_str::CppStr;


impl Node {
    pub fn set_crd(&mut self, new_crd : &mut [f64]) {
        self.coord[0] = new_crd[0];
        self.coord[1] = new_crd[1];
        self.coord[2] = new_crd[2];
        return;
    }

    pub fn set_displacement(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..6 {
            self.displacement[i1] = new_disp[i1];
        }
        return;
    }

    pub fn set_velocity(&mut self, new_vel : &mut [f64]) {
        for i1 in 0..6 {
            self.velocity[i1] = new_vel[i1];
        }
        return;
    }

    pub fn set_acceleration(&mut self, new_acc : &mut [f64]) {
        for i1 in 0..6 {
            self.acceleration[i1] = new_acc[i1];
        }
        return;
    }

    pub fn add_to_displacement(&mut self, del_disp : &mut [f64]) {
        for i1 in 0..6 {
            self.displacement[i1]  +=  del_disp[i1];
        }
        return;
    }

    pub fn set_initial_disp(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..6 {
            self.initial_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn set_initial_vel(&mut self, new_vel : &mut [f64]) {
        for i1 in 0..6 {
            self.initial_vel[i1] = new_vel[i1];
        }
        return;
    }

    pub fn set_initial_acc(&mut self, new_acc : &mut [f64]) {
        for i1 in 0..6 {
            self.initial_acc[i1] = new_acc[i1];
        }
        return;
    }

    pub fn set_initial_flow(&mut self, new_flow : & [f64]) {
        self.initial_fl_den = new_flow[0];
        for i1 in 0..3 {
            self.initial_fl_vel[i1] = new_flow[i1+1];
        }
        self.initial_temp = new_flow[4];
        self.initial_turb_e = new_flow[5];
    }

    pub fn set_initial_flowdot(&mut self, new_fdot : &[f64]) {
        self.initial_fl_den_dot = new_fdot[0];
        for i1 in 0..3 {
            self.initial_fl_vel_dot[i1] = new_fdot[i1+1];
        }
        self.initial_tdot = new_fdot[4];
        self.initial_turb_edot = new_fdot[5];
    }

    pub fn set_prev_disp(&mut self, new_disp : &mut [f64]) {
        for i1 in 0..6 {
            self.prev_disp[i1] = new_disp[i1];
        }
        return;
    }

    pub fn set_prev_vel(&mut self, new_vel : &mut [f64]) {
        for i1 in 0..6 {
            self.prev_vel[i1] = new_vel[i1];
        }
        return;
    }

    pub fn set_prev_acc(&mut self, new_acc : &mut [f64]) {
        for i1 in 0..6 {
            self.prev_acc[i1] = new_acc[i1];
        }
        return;
    }

    pub fn initialize_disp(&mut self) {
        for i1 in 0..6 {
            self.prev_disp[i1] = self.initial_disp[i1];
            self.prev_vel[i1] = self.initial_vel[i1];
            self.prev_acc[i1] = self.initial_acc[i1];
            self.displacement[i1] = self.initial_disp[i1];
        }
        return;
    }

    pub fn initialize_temp(&mut self) {
        self.prev_temp = self.initial_temp;
        self.prev_tdot = self.initial_tdot;
        self.temperature = self.initial_temp;
        return;
    }

    pub fn initialize_fl_den(&mut self) {
        self.prev_fl_den = self.initial_fl_den;
        self.prev_fl_den_dot = self.initial_fl_den_dot;
        return;
    }

    pub fn update_vel_acc(&mut self, nm_beta : f64, nm_gamma : f64, del_t : f64) {
        let c1 : f64;
        let c2 : f64;
        c1 = 1.0 / (del_t * del_t * (nm_beta - nm_gamma));
        c2 = del_t * del_t * (0.5 + nm_beta - nm_gamma);
        for i1 in 0..6 {
            self.acceleration[i1] = c1 * (self.prev_disp[i1] - self.displacement[i1] + del_t * self.prev_vel[i1] + c2 * self.prev_acc[i1]);
            self.velocity[i1] = self.prev_vel[i1] + del_t * ((1.0 - nm_gamma) * self.prev_acc[i1] + nm_gamma * self.acceleration[i1]);
        }
        
        return;
    }

    pub fn update_tdot(&mut self, nm_gamma : f64, del_t : f64) {
        let c1 : f64;
        let c2 : f64;
        c1 = 1.0 / nm_gamma;
        c2 = 1.0 / del_t;
        self.temp_change_rate = c1 * (c2 * (self.temperature - self.prev_temp) - (1.0 - nm_gamma) * self.prev_tdot);
        return;
    }

    pub fn update_fl_den_dot(&mut self, nm_gamma : f64, del_t : f64) {
        let c1 : f64;
        let c2 : f64;
        c1 = 1.0 / nm_gamma;
        c2 = 1.0 / del_t;
        self.fl_den_dot = c1 * (c2 * (self.fl_den - self.prev_fl_den) - (1.0 - nm_gamma) * self.prev_fl_den_dot);
        return;
    }

    pub fn advance_disp(&mut self) {
        for i1 in 0..6 {
            self.prev_disp[i1] = self.displacement[i1];
            self.prev_vel[i1] = self.velocity[i1];
            self.prev_acc[i1] = self.acceleration[i1];
        }
        return;
    }

    pub fn advance_temp(&mut self) {
        self.prev_temp = self.temperature;
        self.prev_tdot = self.temp_change_rate;
        return;
    }

    pub fn advance_fl_den(&mut self) {
        self.prev_fl_den = self.fl_den;
        self.prev_fl_den_dot = self.fl_den_dot;
        return;
    }

    pub fn advance_flow(&mut self) {
        self.prev_fl_den = self.fl_den;
        self.prev_fl_den_dot = self.fl_den_dot;
        for i in 0..3 {
            self.prev_fl_vel[i] = self.fl_vel[i];
            self.prev_fl_vel_dot[i] = self.fl_vel_dot[i];
        }
        self.prev_temp = self.temperature;
        self.prev_tdot = self.temp_change_rate;
        self.prev_turb_e = self.turb_e;
        self.prev_turb_edot = self.turb_edot;
    }

    pub fn backstep_disp(&mut self) {
        for i1 in 0..6 {
            self.displacement[i1] = self.prev_disp[i1];
            self.velocity[i1] = self.prev_vel[i1];
            self.acceleration[i1] = self.prev_acc[i1];
        }
        return;
    }

    pub fn backstep_temp(&mut self) {
        self.temperature = self.prev_temp;
        self.temp_change_rate = self.prev_tdot;
        return;
    }

    pub fn backstep_fl_den(&mut self) {
        self.fl_den = self.prev_fl_den;
        self.fl_den_dot = self.prev_fl_den_dot;
        return;
    }

    pub fn backstep_flow(&mut self) {
        self.fl_den = self.prev_fl_den;
        self.fl_den_dot = self.prev_fl_den_dot;
        for i in 0..3 {
            self.fl_vel[i] = self.prev_fl_vel[i];
            self.fl_vel_dot[i] = self.prev_fl_vel_dot[i];
        }
        self.temperature = self.prev_temp;
        self.temp_change_rate = self.prev_tdot;
        self.turb_e = self.prev_turb_e;
        self.turb_edot = self.prev_turb_edot;
    }

    pub fn add_design_variable(&mut self, d_index : usize, coef : f64) {
        let mut dv = IDCapsule::new();
        dv.int_dat = d_index;
        dv.doub_dat = coef;
        self.d_var_lst.push_back(dv);
        return;
    }

    pub fn add_element(&mut self, el_index : usize, nd_in_el : usize) {
        for e in self.el_lst.iter() {
            if e.i1 == el_index {
                return;
            }
        }
        self.el_lst.push_back(DualInt {i1 : el_index, i2 : nd_in_el});
    }

    pub fn add_conn_nd(&mut self, nd_index : usize) {
        for nd in self.conn_nds.iter() {
            if *nd == nd_index {
                return;
            }
        }
        self.conn_nds.push_back(nd_index);
    }

    //dup1

    pub fn calc_crd_dfd0(&mut self, dv_ar : & Vec<DesignVariable>) {
        self.coord_dfd0[0].set_val(self.coord[0]);
        self.coord_dfd0[1].set_val(self.coord[1]);
        self.coord_dfd0[2].set_val(self.coord[2]);
        let mut d_index : usize;
        let mut d_val = DiffDoub0::new();
        let mut comp : usize;
        let mut cat : CppStr;
        let mut coef = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        for dv in self.d_var_lst.iter() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd0(&mut d_val);
            cat = this_dv.category.clone();
            comp = this_dv.component - 1;
            if cat.s == "nodeCoord" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                self.coord_dfd0[comp].add(& coef);
            }
        }
    }

    pub fn get_crd_dfd0(&self, crd_out : &mut [DiffDoub0], dv_ar : & Vec<DesignVariable>) {
        crd_out[0].set_val_dfd0(&self.coord_dfd0[0]);
        crd_out[1].set_val_dfd0(&self.coord_dfd0[1]);
        crd_out[2].set_val_dfd0(&self.coord_dfd0[2]);
    }

    pub fn get_def_crd_dfd0(&self, crd_out : &mut [DiffDoub0]) {
        for i in 0..3 {
            crd_out[i].set_val(self.displacement[i]);
            crd_out[i].add(&self.coord_dfd0[i]);
        }
    }

    pub fn get_disp_dfd0(&mut self, disp : &mut [DiffDoub0]) {
        for i1 in 0..6 {
            disp[i1].set_val(self.displacement[i1]);
        }
        return;
    }

    pub fn get_elastic_dvload_dfd0(&mut self, ld : &mut [DiffDoub0], dv_ar : & Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut cat : CppStr;
        let mut comp : usize;
        let mut this_dv : &DesignVariable;
        
        for i1 in 0..6 {
            ld[i1].set_val(0.0);
        }
        for dv in self.d_var_lst.iter_mut() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd0(&mut d_val);
            cat = this_dv.category.clone();
            comp = this_dv.component - 1;
            if cat.s == "elasticLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld[comp].add(& coef);
            }
        }
        return;
    }

    pub fn get_thermal_dvload_dfd0(&mut self, ld : &mut DiffDoub0, dv_ar : & Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        ld.set_val(0.0);
        for dv in self.d_var_lst.iter_mut() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd0(&mut d_val);
            cat = this_dv.category.clone();
            if cat.s == "thermalLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld.add(& coef);
            }
        }
        return;
    }

    pub fn get_diff_dvload_dfd0(&self, ld : &mut DiffDoub0, dv_ar : &Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        ld.set_val(0.0);
        for dv in self.d_var_lst.iter() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd0(&mut d_val);
            cat = this_dv.category.clone();
            if cat.s == "diffusionLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld.add(& coef);
            }
        }
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn calc_crd_dfd1(&mut self, dv_ar : & Vec<DesignVariable>) {
        self.coord_dfd1[0].set_val(self.coord[0]);
        self.coord_dfd1[1].set_val(self.coord[1]);
        self.coord_dfd1[2].set_val(self.coord[2]);
        let mut d_index : usize;
        let mut d_val = DiffDoub1::new();
        let mut comp : usize;
        let mut cat : CppStr;
        let mut coef = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        for dv in self.d_var_lst.iter() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd1(&mut d_val);
            cat = this_dv.category.clone();
            comp = this_dv.component - 1;
            if cat.s == "nodeCoord" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                self.coord_dfd1[comp].add(& coef);
            }
        }
    }

    pub fn get_crd_dfd1(&self, crd_out : &mut [DiffDoub1], dv_ar : & Vec<DesignVariable>) {
        crd_out[0].set_val_dfd1(&self.coord_dfd1[0]);
        crd_out[1].set_val_dfd1(&self.coord_dfd1[1]);
        crd_out[2].set_val_dfd1(&self.coord_dfd1[2]);
    }

    pub fn get_def_crd_dfd1(&self, crd_out : &mut [DiffDoub1]) {
        for i in 0..3 {
            crd_out[i].set_val(self.displacement[i]);
            crd_out[i].add(&self.coord_dfd1[i]);
        }
    }

    pub fn get_disp_dfd1(&mut self, disp : &mut [DiffDoub1]) {
        for i1 in 0..6 {
            disp[i1].set_val(self.displacement[i1]);
        }
        return;
    }

    pub fn get_elastic_dvload_dfd1(&mut self, ld : &mut [DiffDoub1], dv_ar : & Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut cat : CppStr;
        let mut comp : usize;
        let mut this_dv : &DesignVariable;
        
        for i1 in 0..6 {
            ld[i1].set_val(0.0);
        }
        for dv in self.d_var_lst.iter_mut() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd1(&mut d_val);
            cat = this_dv.category.clone();
            comp = this_dv.component - 1;
            if cat.s == "elasticLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld[comp].add(& coef);
            }
        }
        return;
    }

    pub fn get_thermal_dvload_dfd1(&mut self, ld : &mut DiffDoub1, dv_ar : & Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        ld.set_val(0.0);
        for dv in self.d_var_lst.iter_mut() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd1(&mut d_val);
            cat = this_dv.category.clone();
            if cat.s == "thermalLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld.add(& coef);
            }
        }
        return;
    }

    pub fn get_diff_dvload_dfd1(&self, ld : &mut DiffDoub1, dv_ar : &Vec<DesignVariable>) {
        let mut d_index : usize;
        let mut d_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        ld.set_val(0.0);
        for dv in self.d_var_lst.iter() {
            d_index = dv.int_dat;
            this_dv = &dv_ar[d_index];
            this_dv.get_value_dfd1(&mut d_val);
            cat = this_dv.category.clone();
            if cat.s == "diffusionLoad" {
                coef.set_val(dv.doub_dat);
                coef.mult(& d_val);
                ld.add(& coef);
            }
        }
    }

    //end dup
 
//end skip 
 
 
 
 
 
 
 
  
}


