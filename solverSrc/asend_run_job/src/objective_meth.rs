use crate::objective::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::list_ent::*;
use crate::node::*;
use crate::element::*;
use crate::design_var::*;
use crate::nd_el_set::*;
use crate::section::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;
use crate::fmath::*;
use crate::scratch::*;

use std::collections::linked_list::IterMut;

impl ObjectiveTerm {
    pub fn set_active_time(&mut self, new_at : &mut [f64]) {
        self.active_time[0] = new_at[0];
        self.active_time[1] = new_at[1];
        return;
    }

    pub fn allocate_obj(&mut self, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>) {
        let mut i1 : usize;
        if (self.q_len == 0) {
            if (self.el_set_ptr < max_int) {
                self.q_len = el_sets[self.el_set_ptr].labels.len();
            }
            else if (self.nd_set_ptr < max_int) {
                self.q_len = nd_sets[self.nd_set_ptr].labels.len();
            }
            if (self.q_len > 0) {
                self.q_vec = vec![0f64; self.q_len];
                if (self.optr.s == "powerNorm") {
                    self.tgt_vec = vec![0f64; self.q_len];
                } else if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                    self.el_vol_vec = vec![0f64; self.q_len];
                    self.tgt_vec = vec![0f64; 1];
                }
            }
        }
        
        for i1 in 0..self.q_len {
            self.q_vec[i1] = 0.0;
            if (self.optr.s == "powerNorm") {
                self.tgt_vec[i1] = 0.0;
            }
            else {
                self.el_vol_vec[i1] = 0.0;
            }
        }
        if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
            self.tgt_vec[0] = 0.0;
        }
        
        return;
    }

    pub fn allocate_obj_grad(&mut self) {
        let mut i1 : usize;
        if (self.q_len > 0 && self.d_qd_u.dim == 0) {
            self.d_qd_u.set_dim(self.q_len);
            self.d_qd_v.set_dim(self.q_len);
            self.d_qd_a.set_dim(self.q_len);
            self.d_qd_t.set_dim(self.q_len);
            self.d_qd_tdot.set_dim(self.q_len);
            self.d_qd_d.set_dim(self.q_len);
            self.d_vd_d.set_dim(self.q_len);
            if (self.optr.s == "powerNorm") {
                self.err_norm_vec = vec![0f64; self.q_len];
            }
        }
        self.d_qd_u.zero_all();
        self.d_qd_v.zero_all();
        self.d_qd_a.zero_all();
        self.d_qd_t.zero_all();
        self.d_qd_tdot.zero_all();
        self.d_qd_d.zero_all();
        self.d_vd_d.zero_all();
        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] = 0.0;
        }
        return;
    }

    pub fn get_power_norm(&mut self) -> f64 {
        let mut i1 : usize;
        let mut q_err : f64;
        let mut p_sum : f64 =  0.0;
        for i1 in 0..self.q_len {
            q_err = self.q_vec[i1] - self.tgt_vec[i1];
            p_sum  +=  powf(q_err, self.expnt);
        }
        return  self.coef * p_sum;
    }

    pub fn d_power_normd_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        self.d_qd_u.vector_multiply(d_ld_u, &mut  self.err_norm_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut  self.err_norm_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut  self.err_norm_vec,  true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut  self.err_norm_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut  self.err_norm_vec,  true);
        
        return;
    }

    pub fn d_power_normd_d(&mut self, d_ld_d : &mut Vec<f64>) {
        self.d_qd_d.vector_multiply(d_ld_d, &mut self.err_norm_vec, true);
        
        return;
    }

    pub fn get_vol_integral(&mut self) -> f64 {
        let mut i1 : usize;
        let mut v_int : f64 =  0.0;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
        }
        v_int  -=  self.tgt_vec[0];
        return  self.coef * powf(v_int, self.expnt);
    }

    pub fn d_vol_integrald_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        let mut i1 : usize;
        let mut v_int : f64 =  0.0;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
        }
        v_int  -=  self.tgt_vec[0];
        v_int = self.coef*self.expnt*powf(v_int, self.expnt - 1.0);
        
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1] *=  v_int;
        }
        
        self.d_qd_u.vector_multiply(d_ld_u, &mut  self.el_vol_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut  self.el_vol_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut  self.el_vol_vec,  true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut  self.el_vol_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut  self.el_vol_vec,  true);
        
        v_int = 1.0 / v_int;
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  v_int;
        }
        
        return;
    }

    pub fn d_vol_integrald_d(&mut self, d_ld_d : &mut Vec<f64>) {
        let mut i1 : usize;
        let mut v_int : f64 =  0.0;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
        }
        v_int  -=  self.tgt_vec[0];
        v_int = self.coef * self.expnt * powf(v_int, self.expnt - 1.0);
        
        for i1 in 0..self.q_len {
            self.q_vec[i1]  *=  v_int;
            self.el_vol_vec[i1]  *=  v_int;
        }
        
        self.d_qd_d.vector_multiply(d_ld_d, &mut  self.el_vol_vec,  true);
        self.d_vd_d.vector_multiply(d_ld_d, &mut  self.q_vec,  true);
        
        v_int = 1.0 / v_int;
        for i1 in 0..self.q_len {
            self.q_vec[i1]  *=  v_int;
            self.el_vol_vec[i1]  *=  v_int;
        }
        
        return;
    }

    pub fn get_vol_average(&mut self) -> f64 {
        let mut i1 : usize;
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let mut vol_avg : f64;
        let mut va_err : f64;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
            tot_vol  +=  self.el_vol_vec[i1];
        }
        vol_avg = v_int/tot_vol;
        va_err = vol_avg - self.tgt_vec[0];
        return  self.coef * powf(va_err, self.expnt);
    }

    pub fn d_vol_averaged_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        let mut i1 : usize;
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let mut vol_avg : f64;
        let mut va_err : f64;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
            tot_vol  +=  self.el_vol_vec[i1];
        }
        vol_avg = v_int / tot_vol;
        va_err = vol_avg - self.tgt_vec[0];
        va_err = self.coef * self.expnt * powf(va_err, self.expnt - 1.0)/tot_vol;
        
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  va_err;
        }
        
        self.d_qd_u.vector_multiply(d_ld_u, &mut  self.el_vol_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut  self.el_vol_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut  self.el_vol_vec,  true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut  self.el_vol_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut  self.el_vol_vec,  true);
        
        va_err = 1.0 / va_err;
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  va_err;
        }
        
        return;
    }

    pub fn d_vol_averaged_d(&mut self, d_ld_d : &mut Vec<f64>) {
        let mut i1 : usize;
        let mut col : usize;
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let mut vol_avg : f64;
        let mut va_err : f64;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
            tot_vol  +=  self.el_vol_vec[i1];
        }
        vol_avg = v_int / tot_vol;
        va_err = vol_avg - self.tgt_vec[0];
        va_err = self.coef * self.expnt * powf(va_err, self.expnt - 1.0) / tot_vol;
        
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  va_err;
            self.q_vec[i1]  *=  va_err;
        }
        
        self.d_qd_d.vector_multiply(d_ld_d, &mut  self.el_vol_vec,  true);
        self.d_vd_d.vector_multiply(d_ld_d, &mut  self.q_vec,  true);
        
        va_err = 1.0 / va_err;
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  va_err;
            self.q_vec[i1]  *=  va_err;
        }
        va_err = 1.0 / va_err;
        
        va_err = va_err * vol_avg;
        
        for i1 in 0..self.q_len {
            //let mut this_row = &self.d_vd_d.matrix[i1];
            for me in self.d_vd_d.matrix[i1].row_vec.iter_mut() {
                col = me.col;
                d_ld_d[col]  -=  va_err * me.value;
            }
        }
        
        return;
    }

    pub fn get_obj_val(&mut self, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>, st_pre : &mut DiffDoub0StressPrereq) {
        if (time < self.active_time[0] || time > self.active_time[1]) {
            return;
        }
        
        let mut i1 : usize;
        let mut fi : usize;
        let mut num_lay : usize;
        let mut tgt_len : usize;
        let mut tgt_val : f64;
        self.allocate_obj(nd_sets, el_sets);
        let mut q_ind : usize;
        let mut nd_data : [f64; 6] = [0f64; 6];
        let mut strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut se_den : f64;
        let mut def = [DiffDoub0::new(); 9];
        let mut frc_mom = [DiffDoub0::new(); 9];
        let mut flux = [DiffDoub0::new(); 3];
        let mut t_grad = [DiffDoub0::new(); 3];
        let mut e_vol = DiffDoub0::new();
        let mut e_den = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        
        let mut cat_list = CppStr::from("displacement velocity acceleration temperature tdot");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.nd_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Node Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            for ndi in nd_sets[self.nd_set_ptr].labels.iter_mut() {
                if (self.category.s == "displacement") {
                    self.q_vec[q_ind] = nd_ar[*ndi].displacement[self.component - 1];
                } else if (self.category.s == "velocity") {
                    self.q_vec[q_ind] = nd_ar[*ndi].velocity[self.component - 1];
                } else if (self.category.s == "acceleration") {
                    self.q_vec[q_ind] = nd_ar[*ndi].acceleration[self.component - 1];
                } else if (self.category.s == "temperature") {
                    self.q_vec[q_ind] = nd_ar[*ndi].temperature;
                } else if (self.category.s == "tdot") {
                    self.q_vec[q_ind] = nd_ar[*ndi].temp_change_rate;
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                tgt_len = self.tgt_vals.len();
                if (tgt_len == 0) {
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = 0.0;
                    }
                } else if (tgt_len == 1) {
                    tgt_val = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = tgt_val;
                    }
                } else {
                    q_ind = 0;
                    for tv in self.tgt_vals.iter_mut() {
                        self.tgt_vec[q_ind] = *tv;
                        q_ind += 1usize;
                    }
                }
                self.value += self.get_power_norm();
                return;
            } else if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                for i1 in 0..self.q_len {
                    self.el_vol_vec[i1] = 1.0;
                }
                tgt_len = self.tgt_vals.len();
                if (tgt_len == 0) {
                    self.tgt_vec[0] = 0.0;
                } else {
                    self.tgt_vec[0] = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                }
                if (self.optr.s == "volumeIntegral") {
                    self.value += self.get_vol_integral();
                    return;
                } else {
                    self.value += self.get_vol_average();
                    return;
                }
            }
        }
        cat_list = CppStr::from("stress strain strainEnergyDen");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.",err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar,  dv_ar);
                this_el.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                if (self.category.s == "stress") {
                    self.q_vec[q_ind] = stress[self.component - 1].val;
                } else if (self.category.s == "strain") {
                    self.q_vec[q_ind] = strain[self.component - 1].val;
                } else {
                    se_den = 0.0;
                    for i1 in 0..6 {
                        se_den  +=  stress[i1].val * strain[i1].val;
                    }
                    se_den  *=  0.5;
                    self.q_vec[q_ind] = se_den;
                }
                if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                    this_el.get_volume_dfd0(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                    self.el_vol_vec[q_ind] = e_vol.val;
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                tgt_len = self.tgt_vals.len();
                if (tgt_len == 0) {
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = 0.0;
                    }
                }
                else if (tgt_len == 1) {
                    tgt_val = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = tgt_val;
                    }
                }
                else {
                    q_ind = 0;
                    for tv in self.tgt_vals.iter_mut() {
                        self.tgt_vec[q_ind] = *tv;
                        q_ind += 1usize;
                    }
                }
                self.value += self.get_power_norm();
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.tgt_vals.len() == 0) {
                    self.tgt_vec[0] = 0.0;
                } else {
                    self.tgt_vec[0] = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                }
                if (self.optr.s == "volumeIntegral") {
                    self.value += self.get_vol_integral();
                    return;
                } else {
                    self.value += self.get_vol_average();
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("sectionDef sectionFrcMom");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre,  sec_ar, mat_ar, nd_ar, dv_ar);
                //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                this_el.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut s_cent,  n_lgeom, st_pre);
                if (self.category.s == "sectionDef") {
                    self.q_vec[q_ind] = def[self.component - 1].val;
                }
                else if (self.category.s == "sectionFrcMom") {
                    self.q_vec[q_ind] = frc_mom[self.component - 1].val;
                }
                if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                    num_lay = sec_ar[this_el.sect_ptr].layers.len();
                    if (num_lay > 1) {
                        e_vol.set_val(0.0);
                        for i1 in 0..num_lay {
                            this_el.get_volume_dfd0(&mut tmp, st_pre,  i1,  sec_ar, dv_ar);
                            e_vol.add(& tmp);
                        }
                    }
                    else {
                        this_el.get_volume_dfd0(&mut e_vol, st_pre,  0, sec_ar, dv_ar);
                    }
                    self.el_vol_vec[q_ind] = e_vol.val;
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                tgt_len = self.tgt_vals.len();
                if (tgt_len == 0) {
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = 0.0;
                    }
                }
                else if (tgt_len == 1) {
                    tgt_val = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = tgt_val;
                    }
                }
                else {
                    q_ind = 0;
                    for tv in self.tgt_vals.iter_mut() {
                        self.tgt_vec[q_ind] = *tv;
                        q_ind += 1usize;
                    }
                }
                self.value  += self.get_power_norm();
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.tgt_vals.len() == 0) {
                    self.tgt_vec[0] = 0.0;
                }
                else {
                    self.tgt_vec[0] = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                }
                if (self.optr.s == "volumeIntegral") {
                    self.value  += self.get_vol_integral();
                    return;
                }
                else {
                    self.value  += self.get_vol_average();
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("flux tempGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre,  sec_ar, mat_ar, nd_ar, dv_ar);
                this_el.get_flux_tgrad_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                if (self.category.s == "flux") {
                    self.q_vec[q_ind] = flux[self.component - 1].val;
                }
                else if (self.category.s == "tempGradient") {
                    self.q_vec[q_ind] = t_grad[self.component - 1].val;
                }
                if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                    this_el.get_volume_dfd0(&mut e_vol, st_pre, self.layer, sec_ar, dv_ar);
                    self.el_vol_vec[q_ind] = e_vol.val;
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                tgt_len = self.tgt_vals.len();
                if (tgt_len == 0) {
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = 0.0;
                    }
                }
                else if (tgt_len == 1) {
                    tgt_val = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                    for i1 in 0..self.q_len {
                        self.tgt_vec[i1] = tgt_val;
                    }
                }
                else {
                    q_ind = 0;
                    for tv in self.tgt_vals.iter_mut() {
                        self.tgt_vec[q_ind] = *tv;
                        q_ind += 1usize;
                    }
                }
                self.value  += self.get_power_norm();
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.tgt_vals.len() == 0) {
                    self.tgt_vec[0] = 0.0;
                }
                else {
                    self.tgt_vec[0] = match self.tgt_vals.front() {
                        None => 0.0f64,
                        Some(x) => *x,
                    };
                }
                if (self.optr.s == "volumeIntegral") {
                    self.value += self.get_vol_integral();
                    return;
                }
                else {
                    self.value += self.get_vol_average();
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("mass volume");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &Element;
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &el_ar[*eli];
                this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                this_el.get_volume_dfd0(&mut e_vol, st_pre, self.layer, sec_ar, dv_ar);
                self.el_vol_vec[q_ind] = e_vol.val;
                if (self.category.s == "volume") {
                    self.q_vec[q_ind] = 1.0;
                } else {
                    this_el.get_density_dfd0(&mut e_den,  self.layer,  sec_ar,  mat_ar, dv_ar);
                    self.q_vec[q_ind] = e_den.val;
                }
                q_ind += 1usize;
            }
            if (self.tgt_vals.len() == 0) {
                self.tgt_vec[0] = 0.0;
            } else {
                self.tgt_vec[0] = match self.tgt_vals.front() {
                    None => 0.0f64,
                    Some(x) => *x,
                };
            }
            self.value += self.get_vol_integral();
            return;
        }
        
        return;
    }

    pub fn getd_ld_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>, time : f64, n_lgeom : bool,
        nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>,
        st_pre : &mut DiffDoub0StressPrereq, scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        if (time < self.active_time[0] || time > self.active_time[1]) {
            return;
        }
        
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut fi : usize;
        self.allocate_obj_grad();
        let mut q_ind : usize;
        let mut dof_ind : usize;
        let mut curr_rank : usize;
        let mut strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut dsd_u = vec![DiffDoub0::new(); 288];
        let mut ded_u = vec![DiffDoub0::new(); 288];
        let mut dsd_t = vec![DiffDoub0::new(); 90];
        let mut dse_dend_u = [DiffDoub0::new(); 33];
        let mut dse_dend_t = [DiffDoub0::new(); 10];
        let mut def = [DiffDoub0::new(); 9];
        let mut frc_mom = [DiffDoub0::new(); 9];
        let mut flux = [DiffDoub0::new(); 3];
        let mut t_grad = [DiffDoub0::new(); 3];
        let mut d_fd_t = vec![DiffDoub0::new(); 30];
        let mut d_tgd_t = vec![DiffDoub0::new(); 30];
        let mut el_num_nds : usize;
        let mut el_dof_per_nd : usize;
        let mut el_num_int_dof : usize;
        let mut el_tot_dof : usize;
        let mut e_vol = DiffDoub0::new();
        let mut e_den = DiffDoub0::new();

        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        
        let mut cat_list = CppStr::from("displacement velocity acceleration temperature tdot");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.nd_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Node Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            for ndi in nd_sets[self.nd_set_ptr].labels.iter_mut() {
                dof_ind = nd_ar[*ndi].dof_index[self.component - 1];
                curr_rank = nd_ar[*ndi].sorted_rank;
                if (self.category.s == "displacement") {
                    self.d_qd_u.add_entry(q_ind, dof_ind, 1.0);
                }
                else if (self.category.s == "velocity") {
                    self.d_qd_v.add_entry(q_ind, dof_ind,  1.0);
                }
                else if (self.category.s == "acceleration") {
                    self.d_qd_a.add_entry(q_ind, dof_ind,  1.0);
                }
                else if (self.category.s == "temperature") {
                    self.d_qd_t.add_entry(q_ind, curr_rank, 1.0);
                }
                else if (self.category.s == "tdot") {
                    self.d_qd_tdot.add_entry(q_ind,   curr_rank,   1.0);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                return;
            }
            else if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("stress strain strainEnergyDen");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                this_el.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                this_el.d_stress_straind_u_dfd0(&mut dsd_u, &mut  ded_u, &mut  dsd_t, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                el_num_nds = this_el.num_nds;
                el_dof_per_nd = this_el.dof_per_nd;
                el_num_int_dof = this_el.num_int_dof;
                el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
                i1 = el_tot_dof * (self.component - 1);
                i2 = el_num_nds * (self.component - 1);
                if (self.category.s == "stress") {
                    sub_vec_dfd0(&mut scr_v1, &mut  dsd_u,  i1,  i1 + el_tot_dof);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                    sub_vec_dfd0(&mut scr_v1, &mut  dsd_t,  i2,  i2 + el_num_nds);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                }
                else if (self.category.s == "strain") {
                    sub_vec_dfd0(&mut scr_v1, &mut  ded_u,  i1,  i1 + el_tot_dof);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                }
                else {
                    for i2 in 0..el_tot_dof {
                        dse_dend_u[i2].set_val(0.0);
                        i3 = i2;
                        for i1 in 0..6 {
                            dse_dend_u[i2].val  +=  stress[i1].val * ded_u[i3].val + dsd_u[i3].val * strain[i1].val;
                            i3  +=  el_tot_dof;
                        }
                        dse_dend_u[i2].val  *=  0.5;
                    }
                    for i2 in 0..el_num_nds {
                        dse_dend_t[i2].set_val(0.0);
                        i3 = i2;
                        for i1 in 0..6 {
                            dse_dend_t[i2].val  +=  dsd_t[i3].val * strain[i1].val;
                            i3  +=  el_num_nds;
                        }
                        dse_dend_t[i2].val  *=  0.5;
                    }
                    ar_to_vec_dfd0(&mut dse_dend_u, &mut  scr_v1,  0,  33);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                    ar_to_vec_dfd0(&mut dse_dend_t, &mut  scr_v1,  0,  10);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("sectionDef sectionFrcMom");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                //this_el->d_stress_straind_u_dfd0(dsd_u, ded_u, dsd_t, spt, self.layer, n_lgeom, st_pre);
                this_el.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut s_cent,  n_lgeom, st_pre);
                this_el.d_def_frc_momd_u_dfd0(&mut ded_u, &mut  dsd_u, &mut  dsd_t, &mut s_cent,  n_lgeom, st_pre);
                el_num_nds = this_el.num_nds;
                el_dof_per_nd = this_el.dof_per_nd;
                el_num_int_dof = this_el.num_int_dof;
                el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
                i1 = el_tot_dof * (self.component - 1);
                i2 = el_num_nds * (self.component - 1);
                if (self.category.s == "sectionFrcMom") {
                    sub_vec_dfd0(&mut scr_v1, &mut  dsd_u,  i1,  i1 + el_tot_dof);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                    sub_vec_dfd0(&mut scr_v1, &mut  dsd_t,  i2,  i2 + el_num_nds);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                }
                else if (self.category.s == "sectionDef") {
                    sub_vec_dfd0(&mut scr_v1, &mut  ded_u,  i1,  i1 + el_tot_dof);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("flux tempGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                this_el.get_flux_tgrad_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                this_el.d_flux_tgradd_t_dfd0(&mut d_fd_t, &mut  d_tgd_t, &mut s_cent,  self.layer, st_pre);
                el_num_nds = this_el.num_nds;
                i1 = el_num_nds * (self.component - 1);
                if (self.category.s == "flux") {
                    sub_vec_dfd0(&mut scr_v1, &mut  d_fd_t,  i1,  i1 + el_num_nds);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                }
                else if (self.category.s == "tempGradient") {
                    sub_vec_dfd0(&mut scr_v1, &mut  d_tgd_t,  i1,  i1 + el_num_nds);
                    this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        return;
    }

    pub fn getd_ld_d(&mut self, d_ld_d : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : &mut Vec<DesignVariable>, st_pre : &mut DiffDoub1StressPrereq) {
        if (time < self.active_time[0] || time > self.active_time[1]) {
            return;
        }
        
        let mut i1 : usize;
        let mut num_lay : usize;
        let mut fi : usize;
        let mut q_ind : usize;
        let mut dv_val = DiffDoub0::new();
        let mut strain = [DiffDoub1::new(); 6];
        let mut stress = [DiffDoub1::new(); 6];
        let mut se_den : f64;
        let mut def = [DiffDoub1::new(); 9];
        let mut frc_mom = [DiffDoub1::new(); 9];
        let mut flux = [DiffDoub1::new(); 3];
        let mut t_grad = [DiffDoub1::new(); 3];
        let mut e_vol = DiffDoub1::new();
        let mut e_den = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        
        let mut cat_list = CppStr::from("stress strain strainEnergyDen");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut this_dv : &mut DesignVariable;
            let mut tmp_dvi = vec![0usize; dv_ar.len()];
            let mut dv_len = 0usize;
            let mut dvi = 0usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                dv_len = 0;
                for cdv in this_el.comp_dvars.iter() {
                    tmp_dvi[dv_len] = *cdv;
                    dv_len += 1usize;
                }
                //for dvi in this_el.comp_dvars.iter_mut() {
                for i in 0..dv_len {
                    dvi = tmp_dvi[i];
                    //this_dv = &mut dv_ar[dvi];
                    dv_ar[dvi].get_value_dfd0(&mut dv_val);
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val, 1.0);
                    this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_stress_strain_dfd1(&mut stress, &mut  strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                    if (self.category.s == "stress") {
                        self.d_qd_d.add_entry(q_ind,  dvi,   stress[self.component - 1].dval);
                    }
                    else if (self.category.s == "strain") {
                        self.d_qd_d.add_entry(q_ind,  dvi,   strain[self.component - 1].dval);
                    }
                    else {
                        se_den = 0.0;
                        for i1 in 0..6 {
                            se_den  +=  stress[i1].val * strain[i1].dval + stress[i1].dval * strain[i1].val;
                        }
                        se_den  *=  0.5;
                        self.d_qd_d.add_entry(q_ind, dvi, se_den);
                    }
                    if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                        this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                        self.d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
                    }
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val, 0.0);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_d(d_ld_d);
                    return;
                }
                else {
                    self.d_vol_averaged_d(d_ld_d);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("sectionDef sectionFrcMom");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut this_dv : &mut DesignVariable;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len = 0usize;
            let mut dvi = 0usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                dv_len = 0usize;
                for cdv in this_el.comp_dvars.iter() {
                    tmp_dv[dv_len] = *cdv;
                    dv_len += 1usize;
                }
                //for dvi in this_el.comp_dvars.iter_mut() {
                for i in 0..dv_len {
                    dvi = tmp_dv[i];
                    //this_dv = &mut dv_ar[*dvi];
                    dv_ar[dvi].get_value_dfd0(&mut dv_val);
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val, 1.0);
                    this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                    this_el.get_def_frc_mom_dfd1(&mut def, &mut  frc_mom, &mut s_cent, n_lgeom, st_pre);
                    if (self.category.s == "sectionFrcMom") {
                        self.d_qd_d.add_entry(q_ind, dvi,   frc_mom[self.component - 1].dval);
                    }
                    else if (self.category.s == "sectionDef") {
                        self.d_qd_d.add_entry(q_ind, dvi,   def[self.component - 1].dval);
                    }
                    if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                        num_lay = sec_ar[this_el.sect_ptr].layers.len();
                        if (num_lay > 1) {
                            e_vol.set_val(0.0);
                            for i1 in 0..num_lay {
                                this_el.get_volume_dfd1(&mut tmp, st_pre,  i1, sec_ar, dv_ar);
                                e_vol.add(& tmp);
                            }
                        }
                        else {
                            this_el.get_volume_dfd1(&mut e_vol, st_pre,  0, sec_ar, dv_ar);
                        }
                        self.d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
                    }
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val, 0.0);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_d(d_ld_d);
                    return;
                }
                else {
                    self.d_vol_averaged_d(d_ld_d);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("flux tempGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut this_dv : &mut DesignVariable;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len = 0usize;
            let mut dvi = 0usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &mut el_ar[*eli];
                for i in 0..3 {
                    s_cent[i] = this_el.s_cent[i];
                }
                dv_len = 0usize;
                for cdv in this_el.comp_dvars.iter() {
                    tmp_dv[dv_len] = *cdv;
                    dv_len += 1usize;
                }
                //for dvi in this_el.comp_dvars.iter_mut() {
                for i in 0..dv_len {
                    dvi = tmp_dv[i];
                    //this_dv = &mut dv_ar[*dvi];
                    dv_ar[dvi].get_value_dfd0(&mut dv_val);
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val,    1.0);
                    this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_flux_tgrad_dfd1(&mut flux, &mut  t_grad, &mut s_cent, self.layer, st_pre);
                    if (self.category.s == "flux") {
                        self.d_qd_d.add_entry(q_ind, dvi,   flux[self.component - 1].dval);
                    }
                    else if (self.category.s == "tempGradient") {
                        self.d_qd_d.add_entry(q_ind, dvi,   t_grad[self.component - 1].dval);
                    }
                    if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                        this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                        self.d_vd_d.add_entry(q_ind, dvi,  e_vol.dval);
                    }
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val,    0.0);
                }
                q_ind += 1usize;
            }
            if (self.optr.s == "powerNorm") {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = self.coef * self.expnt * powf((self.q_vec[i1] - self.tgt_vec[i1]), (self.expnt - 1.0));
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            if (self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage") {
                if (self.optr.s == "volumeIntegral") {
                    self.d_vol_integrald_d(d_ld_d);
                    return;
                }
                else {
                    self.d_vol_averaged_d(d_ld_d);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("mass volume");
        fi = cat_list.find(&self.category.s.as_str());
        if (fi < max_int) {
            if (self.el_set_ptr == max_int) {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &Element;
            let mut this_dv : &mut DesignVariable;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len = 0usize;
            let mut dvi = 0usize;
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                this_el = &el_ar[*eli];
                dv_len = 0usize;
                for cdv in this_el.comp_dvars.iter() {
                    tmp_dv[dv_len] = *cdv;
                    dv_len += 1usize;
                }
                //for dvi in this_el.comp_dvars.iter_mut() {
                for i in 0..dv_len {
                    dvi = tmp_dv[i];
                    //this_dv = &mut dv_ar[*dvi];
                    dv_ar[dvi].get_value_dfd0(&mut dv_val);
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val, 1.0);
                    this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                    self.d_vd_d.add_entry(q_ind, dvi,   e_vol.dval);
                    if(self.category.s == "mass") {
                        this_el.get_density_dfd1(&mut e_den,  self.layer, sec_ar, mat_ar, dv_ar);
                        self.d_qd_d.add_entry(q_ind, dvi,   e_den.dval);
                    }
                    dv_ar[dvi].diff_val.set_val_2(dv_val.val,  0.0);
                }
                q_ind += 1usize;
            }
            self.d_vol_integrald_d(d_ld_d);
            return;
        }
        
        return;
    }

}

impl Objective {
    pub fn clear_values(&mut self) {
        for tm in self.terms.iter_mut() {
            tm.value = 0.0;
        }
        return;
    }

    pub fn calculate_terms(&mut self, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>, st_pre : &mut DiffDoub0StressPrereq) {
        for tm in self.terms.iter_mut() {
            tm.get_obj_val(time,  n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
        }
        return;
    }

    pub fn calculated_ld_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>,
        st_pre : &mut DiffDoub0StressPrereq, scr : &mut IterMut<'_,FltScr>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        for tm in self.terms.iter_mut() {
            tm.getd_ld_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot,  time,  n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre, scr, scr_dfd);
        }
        return;
    }

    pub fn calculated_ld_d(&mut self, d_ld_d : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : &mut Vec<DesignVariable>, st_pre : &mut DiffDoub1StressPrereq) {
        for tm in self.terms.iter_mut() {
            tm.getd_ld_d(d_ld_d,  time,  n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
        }
        return;
    }

}


