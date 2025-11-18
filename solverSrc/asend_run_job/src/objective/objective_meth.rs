use crate::objective::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::node::*;
use crate::element::*;
use crate::design_var::*;
use crate::nd_el_set::*;
use crate::section::*;
use crate::list_ent::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;
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
        if self.q_len == 0 {
            if self.el_set_ptr < MAX_INT {
                self.q_len = el_sets[self.el_set_ptr].labels.len();
            }
            else if self.nd_set_ptr < MAX_INT {
                self.q_len = nd_sets[self.nd_set_ptr].labels.len();
            }
            if self.q_len > 0 {
                self.q_vec = vec![0f64; self.q_len];
                if self.optr.s == "powerNorm" || self.optr.s == "ks" {
                    self.tgt_vec = vec![0f64; self.q_len];
                } else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                    self.el_vol_vec = vec![0f64; self.q_len];
                    self.tgt_vec = vec![0f64; 1];
                }
            }
        }
        
        for i1 in 0..self.q_len {
            self.q_vec[i1] = 0.0;
            if self.optr.s == "powerNorm" {
                self.tgt_vec[i1] = 0.0;
            }
            else {
                self.el_vol_vec[i1] = 0.0;
            }
        }
        if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
            self.tgt_vec[0] = 0.0;
        }
        
        return;
    }

    pub fn allocate_obj_grad(&mut self) {
        if self.q_len > 0 && self.d_qd_u.dim == 0 {
            self.d_qd_u.set_dim(self.q_len);
            self.d_qd_v.set_dim(self.q_len);
            self.d_qd_a.set_dim(self.q_len);
            self.d_qd_c.set_dim(self.q_len);
            self.d_qd_cdot.set_dim(self.q_len);
            self.d_qd_t.set_dim(self.q_len);
            self.d_qd_tdot.set_dim(self.q_len);
            self.d_qd_d.set_dim(self.q_len);
            self.d_vd_d.set_dim(self.q_len);
            if self.optr.s == "powerNorm" || self.optr.s == "ks" {
                self.err_norm_vec = vec![0f64; self.q_len];
                self.d_tgtd_d.set_dim(self.q_len);
            }
        }
        self.d_qd_u.zero_all();
        self.d_qd_v.zero_all();
        self.d_qd_a.zero_all();

        self.d_qd_c.zero_all();
        self.d_qd_cdot.zero_all();
        self.d_qd_t.zero_all();
        self.d_qd_tdot.zero_all();
        self.d_qd_d.zero_all();
        self.d_vd_d.zero_all();
        if self.optr.s == "powerNorm" || self.optr.s == "ks" {
            for i1 in 0..self.q_len {
                self.err_norm_vec[i1] = 0.0;
            }
            self.d_tgtd_d.zero_all();
        }
        
    }

    pub fn get_tgt_vec(&mut self, el_ar : &Vec<Element>, el_sets : &Vec<Set>, sec_ar : &Vec<Section>, mat_ar : &Vec<Material>, dv_ar : &Vec<DesignVariable>) {
        let mut i1 : usize;
        let frnt : f64;
        let mut sec_i : usize;
        let mut mat_i : usize;
        let mut this_el : &Element;
        let mut tmp = DiffDoub0::new();
        let mut prop : &Vec<f64>;
        let mut p_vec = vec![DiffDoub0::new(); 3];
        let t_len = self.tgt_vals.len();

        if t_len == self.q_len {
            i1 = 0;
            for t in self.tgt_vals.iter() {
                self.tgt_vec[i1] = *t;
                i1 += 1;
            }
        }
        else if t_len == 1 {
            frnt = match self.tgt_vals.front() {
                None => 0f64,
                Some(x) => *x,
            };
            for i in 0..self.q_len {
                self.tgt_vec[i] = frnt;
            }
        }
        else if self.tgt_tag.s.len() > 0 {
            if self.el_set_ptr == MAX_INT {
                panic!("Error: objectives with property name specified as a target must have a valid element set identified");
            }
            i1 = 0;
            for ei in el_sets[self.el_set_ptr].labels.iter() {
                this_el = &el_ar[*ei];
                sec_i = this_el.sect_ptr;
                if this_el.this_type == 41 || this_el.this_type == 3 {
                    mat_i = sec_ar[sec_i].get_layer_mat_ptr(self.layer);
                }
                else {
                    mat_i = sec_ar[sec_i].mat_ptr;
                }
                prop = match mat_ar[mat_i].custom.get(&self.tgt_tag.s) {
                    None => panic!("Error: property '{}' identified in objective was not found in a corresponding material's custom property list.", self.tgt_tag.s),
                    Some(x) => x,
                };
                tmp.set_val(prop[0]);
                this_el.get_gen_prop_dfd0(&mut tmp, &mut self.tgt_tag, dv_ar);
                self.tgt_vec[i1] = tmp.val;
                i1 += 1;
            }
        }
        else {
            let prop_cats = String::from("stress strain strainEnergy mises");
            if prop_cats.contains(&self.category.s) {
                let mut def_prop = CppStr::new();
                let p_i : usize; 
                if self.category.s == "stress" {
                    if self.component > 4 {
                        def_prop.s = String::from("tensileStrength");
                        p_i = self.component - 1;
                    }
                    else {
                        def_prop.s = String::from("shearStrength");
                        p_i = self.component - 3;
                    }
                }
                else if self.category.s == "strain" {
                    if self.component > 4 {
                        def_prop.s = String::from("tensileStrain");
                        p_i = self.component - 1;
                    }
                    else {
                        def_prop.s = String::from("shearStrain");
                        p_i = self.component - 3;
                    }
                }
                else if self.category.s == "strainEnergy" {
                    def_prop.s = String::from("maxStrainEnergy");
                    p_i = 0;
                }
                else {
                    def_prop.s = String::from("misesStrength");
                    p_i = 0;
                }
                i1 = 0;
                for ei in el_sets[self.el_set_ptr].labels.iter() {
                    this_el = &el_ar[*ei];
                    sec_i = this_el.sect_ptr;
                    if this_el.this_type == 41 || this_el.this_type == 3 {
                        mat_i = sec_ar[sec_i].get_layer_mat_ptr(self.layer);
                    }
                    else {
                        mat_i = sec_ar[sec_i].mat_ptr;
                    }
                    if mat_ar[mat_i].custom.contains_key(&def_prop.s) {
                        prop = match mat_ar[mat_i].custom.get(&def_prop.s) {
                            None => panic!(""),
                            Some(x) => x,
                        };
                        for i in 0..3 {
                            p_vec[i].set_val(prop[i]);
                        }
                        this_el.get_gen_pvec_dfd0(&mut p_vec, &mut def_prop, dv_ar);
                        self.tgt_vec[i1] = p_vec[p_i].val;
                    }
                    else {
                        self.tgt_vec[i1] = match self.optr.s.as_str() {
                            "ks" => 1.0,
                            &_ => 0.0,
                        };
                    }
                    i1 += 1;
                }
            }
            else { 
                let t_val = match self.optr.s.as_str() {
                    "ks" => 1.0f64,
                    &_ => 0.0f64,
                };
                for i in 0..self.q_len {
                    self.tgt_vec[i] = t_val;
                }
            }
        }
    }

    pub fn get_d_tgtd_d(&mut self, el_ar : &Vec<Element>, el_sets : &Vec<Set>, sec_ar : &Vec<Section>, mat_ar : &Vec<Material>, dv_ar : &mut Vec<DesignVariable>) {
        let mut i1 : usize;
        let frnt : f64;
        let mut sec_i : usize;
        let mut mat_i : usize;
        let mut this_el : &Element;
        let mut dv_val : f64;
        let mut tmp = DiffDoub1::new();
        let mut prop : &Vec<f64>;
        let mut p_vec = vec![DiffDoub1::new(); 3];

        self.d_tgtd_d.zero_all();

        if self.tgt_tag.s.len() > 0 {
            if self.el_set_ptr == MAX_INT {
                panic!("Error: objectives with property name specified as a target must have a valid element set identified");
            }
            i1 = 0;
            for ei in el_sets[self.el_set_ptr].labels.iter() {
                this_el = &el_ar[*ei];
                sec_i = this_el.sect_ptr;
                if this_el.this_type == 41 || this_el.this_type == 3 {
                    mat_i = sec_ar[sec_i].get_layer_mat_ptr(self.layer);
                }
                else {
                    mat_i = sec_ar[sec_i].mat_ptr;
                }
                prop = match mat_ar[mat_i].custom.get(&self.tgt_tag.s) {
                    None => panic!("Error: property '{}' identified in objective was not found in a corresponding material's custom property list.", self.tgt_tag.s),
                    Some(x) => x,
                };
                for dv in this_el.comp_dvars.iter() {
                    dv_val = dv_ar[*dv].value.val;
                    dv_ar[*dv].diff_val.set_val_2(dv_val, 1.0);
                    tmp.set_val(prop[0]);
                    this_el.get_gen_prop_dfd1(&mut tmp, &mut self.tgt_tag, dv_ar);
                    self.d_tgtd_d.add_entry(i1, *dv, tmp.dval);
                    dv_ar[*dv].diff_val.set_val_2(dv_val, 0.0);
                }
                i1 += 1;
            }
        }
        else if self.tgt_vals.len() == 0 {
            let prop_cats = String::from("stress strain strainEnergy mises");
            if prop_cats.contains(&self.category.s) {
                let mut def_prop = CppStr::new();
                let p_i : usize; 
                if self.category.s == "stress" {
                    if self.component > 4 {
                        def_prop.s = String::from("tensileStrength");
                        p_i = self.component - 1;
                    }
                    else {
                        def_prop.s = String::from("shearStrength");
                        p_i = self.component - 3;
                    }
                }
                else if self.category.s == "strain" {
                    if self.component > 4 {
                        def_prop.s = String::from("tensileStrain");
                        p_i = self.component - 1;
                    }
                    else {
                        def_prop.s = String::from("shearStrain");
                        p_i = self.component - 3;
                    }
                }
                else if self.category.s == "strainEnergy" {
                    def_prop.s = String::from("maxStrainEnergy");
                    p_i = 0;
                }
                else {
                    def_prop.s = String::from("misesStrength");
                    p_i = 0;
                }
                i1 = 0;
                for ei in el_sets[self.el_set_ptr].labels.iter() {
                    this_el = &el_ar[*ei];
                    sec_i = this_el.sect_ptr;
                    if this_el.this_type == 41 || this_el.this_type == 3 {
                        mat_i = sec_ar[sec_i].get_layer_mat_ptr(self.layer);
                    }
                    else {
                        mat_i = sec_ar[sec_i].mat_ptr;
                    }
                    if mat_ar[mat_i].custom.contains_key(&def_prop.s) {
                        prop = match mat_ar[mat_i].custom.get(&def_prop.s) {
                            None => panic!(""),
                            Some(x) => x,
                        };
                        for dv in this_el.comp_dvars.iter() {
                            dv_val = dv_ar[*dv].value.val;
                            dv_ar[*dv].diff_val.set_val_2(dv_val, 1.0);
                            for i in 0..3 {
                                p_vec[i].set_val(prop[i]);
                            }
                            this_el.get_gen_pvec_dfd1(&mut p_vec, &mut def_prop, dv_ar);
                            self.d_tgtd_d.add_entry(i1, *dv, p_vec[p_i].dval);
                            dv_ar[*dv].diff_val.set_val_2(dv_val, 0.0);
                        }
                    }
                    i1 += 1;
                }

            }
        }

    }

    pub fn get_power_norm(&mut self) -> f64 {
        let mut q_err : f64;
        let mut p_sum : f64 =  0.0;
        for i1 in 0..self.q_len {
            q_err = self.q_vec[i1] - self.tgt_vec[i1];
            p_sum  +=  powf(q_err, self.expnt);
        }
        return self.coef * powf(p_sum, 1f64/self.expnt);
    }

    pub fn d_power_normd_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        let mut q_err : f64;
        let mut p_sum = 0.0f64;
        for i1 in 0..self.q_len {
            q_err = self.q_vec[i1] - self.tgt_vec[i1];
            p_sum += powf(q_err, self.expnt);
        }
        p_sum = powf(p_sum, 1.0f64/self.expnt - 1.0f64);
        p_sum *= self.coef;

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] *= p_sum;
        }

        self.d_qd_u.vector_multiply(d_ld_u, &mut self.err_norm_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut self.err_norm_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut self.err_norm_vec,  true);
        self.d_qd_c.vector_multiply(d_ld_c, &mut self.err_norm_vec, true);
        self.d_qd_cdot.vector_multiply(d_ld_cdot, &mut self.err_norm_vec, true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut self.err_norm_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut self.err_norm_vec,  true);

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] /= p_sum;
        }
        
        return;
    }

    pub fn d_power_normd_d(&mut self, d_ld_d : &mut Vec<f64>) {
        let mut q_err : f64;
        let mut p_sum = 0.0f64;
        for i1 in 0..self.q_len {
            q_err = self.q_vec[i1] - self.tgt_vec[i1];
            p_sum += powf(q_err, self.expnt);
        }
        p_sum = powf(p_sum, 1.0f64/self.expnt - 1.0f64);
        p_sum *= self.coef;

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] *= p_sum;
        }

        self.d_qd_d.vector_multiply(d_ld_d, &mut self.err_norm_vec, true);

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] /= p_sum;
        }
        return;
    }

    pub fn get_ks(&mut self) -> f64 {
        let mut p_sum = 0.0f64;
        let mut exp : f64;
        for i1 in 0..self.q_len {
            exp = (self.expnt*self.q_vec[i1]).exp();
            p_sum += exp;
        }
        (self.coef/self.expnt) * p_sum.ln()
    }

    pub fn d_ksd_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        let mut p_sum = 0.0f64;
        for i1 in 0..self.q_len {
            p_sum += self.err_norm_vec[i1];
        }
        p_sum = self.coef/p_sum;

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] *= p_sum;
        }

        self.d_qd_u.vector_multiply(d_ld_u, &mut self.err_norm_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut self.err_norm_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut self.err_norm_vec,  true);
        self.d_qd_c.vector_multiply(d_ld_c, &mut self.err_norm_vec, true);
        self.d_qd_cdot.vector_multiply(d_ld_cdot, &mut self.err_norm_vec, true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut self.err_norm_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut self.err_norm_vec,  true);

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] /= p_sum;
        }
    }

    pub fn d_ksd_d(&mut self, d_ld_d : &mut Vec<f64>) {
        let mut p_sum = 0.0f64;
        for i1 in 0..self.q_len {
            p_sum += self.err_norm_vec[i1];
        }
        p_sum = self.coef/p_sum;

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] *= p_sum;
        }

        self.d_qd_d.vector_multiply(d_ld_d, &mut self.err_norm_vec, true);

        for i1 in 0..self.q_len {
            self.err_norm_vec[i1] /= p_sum;
        }        
    }

    pub fn get_vol_integral(&mut self) -> f64 {
        let mut v_int : f64 =  0.0;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
        }
        v_int  -=  self.tgt_vec[0];
        return  self.coef * powf(v_int, self.expnt);
    }

    pub fn d_vol_integrald_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
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
        self.d_qd_c.vector_multiply(d_ld_c, &mut self.el_vol_vec, true);
        self.d_qd_cdot.vector_multiply(d_ld_cdot, &mut self.el_vol_vec, true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut  self.el_vol_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut  self.el_vol_vec,  true);
        
        v_int = 1.0 / v_int;
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  v_int;
        }
        
        return;
    }

    pub fn d_vol_integrald_d(&mut self, d_ld_d : &mut Vec<f64>) {
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
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let vol_avg : f64;
        let va_err : f64;
        for i1 in 0..self.q_len {
            v_int  +=  self.q_vec[i1] * self.el_vol_vec[i1];
            tot_vol  +=  self.el_vol_vec[i1];
        }
        vol_avg = v_int/tot_vol;
        va_err = vol_avg - self.tgt_vec[0];
        return  self.coef * powf(va_err, self.expnt);
    }

    pub fn d_vol_averaged_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>) {
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let vol_avg : f64;
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
        
        self.d_qd_u.vector_multiply(d_ld_u, &mut self.el_vol_vec,  true);
        self.d_qd_v.vector_multiply(d_ld_v, &mut self.el_vol_vec,  true);
        self.d_qd_a.vector_multiply(d_ld_a, &mut self.el_vol_vec,  true);
        self.d_qd_c.vector_multiply(d_ld_c, &mut self.el_vol_vec, true);
        self.d_qd_cdot.vector_multiply(d_ld_cdot, &mut self.el_vol_vec, true);
        self.d_qd_t.vector_multiply(d_ld_t, &mut self.el_vol_vec,  true);
        self.d_qd_tdot.vector_multiply(d_ld_tdot, &mut self.el_vol_vec,  true);
        
        va_err = 1.0 / va_err;
        for i1 in 0..self.q_len {
            self.el_vol_vec[i1]  *=  va_err;
        }
        
    }

    pub fn d_vol_averaged_d(&mut self, d_ld_d : &mut Vec<f64>) {
        let mut col : usize;
        let mut v_int : f64 =  0.0;
        let mut tot_vol : f64 =  0.0;
        let vol_avg : f64;
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
        if time < self.active_time[0] || time > self.active_time[1] {
            return;
        }
        
        let mut fi : usize;
        let mut num_lay : usize;
        let tgt_len : usize;
        let tgt_val : f64;
        self.allocate_obj(nd_sets, el_sets);
        self.get_tgt_vec(el_ar, el_sets, sec_ar, mat_ar, dv_ar);
        let mut q_ind : usize;
        let mut strain = [DiffDoub0::new(); 6];
        let mut t_strain = [DiffDoub0::new(); 6];
        let mut d_strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut se_den : f64;
        let mut def = [DiffDoub0::new(); 9];
        let mut frc_mom = [DiffDoub0::new(); 9];
        let mut flux = [DiffDoub0::new(); 3];
        let mut t_grad = [DiffDoub0::new(); 3];
        let mut e_vol = DiffDoub0::new();
        let mut e_den = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        
        let mut cat_list = CppStr::from("displacement velocity acceleration concentration cdot temperature tdot");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.nd_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Node Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            for ndi in nd_sets[self.nd_set_ptr].labels.iter_mut() {
                self.q_vec[q_ind] = match self.category.s.as_str() {
                    "displacement" => nd_ar[*ndi].displacement[self.component - 1],
                    "velocity" => nd_ar[*ndi].velocity[self.component - 1],
                    "acceleration" => nd_ar[*ndi].acceleration[self.component - 1],
                    "concentration" => nd_ar[*ndi].fl_den,
                    "cdot" => nd_ar[*ndi].fl_den_dot,
                    "temperature" => nd_ar[*ndi].temperature,
                    "tdot" => nd_ar[*ndi].temp_change_rate,
                    &_ => 0.0,
                };
                
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                self.value += self.get_power_norm();
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    self.q_vec[i1] /= self.tgt_vec[i1];
                }
                self.value += self.get_ks();
            } 
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                for i1 in 0..self.q_len {
                    self.el_vol_vec[i1] = 1.0;
                }
                if self.optr.s == "volumeIntegral" {
                    self.value += self.get_vol_integral();
                    return;
                } else {
                    self.value += self.get_vol_average();
                    return;
                }
            }
        }

        cat_list = CppStr::from("stress strain strainEnergyDen mises tsaiWu");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.",err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar,  dv_ar);
                    this_el.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut t_strain, &mut d_strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                    if self.category.s == "stress" {
                        self.q_vec[q_ind] = stress[self.component - 1].val;
                    } 
                    else if self.category.s == "strain" {
                        self.q_vec[q_ind] = strain[self.component - 1].val;
                    } 
                    else if self.category.s == "strainEnergy" {
                        se_den = 0.0;
                        for i1 in 0..6 {
                            se_den  +=  stress[i1].val * strain[i1].val;
                        }
                        se_den  *=  0.5;
                        self.q_vec[q_ind] = se_den;
                    }
                    else if self.category.s == "mises" {
                        this_el.get_mises_dfd0(&mut tmp, &stress);
                        self.q_vec[q_ind] = tmp.val; 
                    }
                    else if self.category.s == "tsaiWu" {
                        this_el.get_tsai_wu_dfd0(&mut tmp, &stress, self.layer, sec_ar, mat_ar, dv_ar);
                        self.q_vec[q_ind] = tmp.val;
                    }

                    if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                        this_el.get_volume_dfd0(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                        self.el_vol_vec[q_ind] = e_vol.val;
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                self.value += self.get_power_norm();
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    self.q_vec[i1] /= self.tgt_vec[i1];
                }
                self.value += self.get_ks();
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
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
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre,  sec_ar, mat_ar, nd_ar, dv_ar);
                    //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                    this_el.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut s_cent,  n_lgeom, st_pre);
                    if self.category.s == "sectionDef" {
                        self.q_vec[q_ind] = def[self.component - 1].val;
                    }
                    else if self.category.s == "sectionFrcMom" {
                        self.q_vec[q_ind] = frc_mom[self.component - 1].val;
                    }
                    if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                        num_lay = sec_ar[this_el.sect_ptr].layers.len();
                        if num_lay > 1 {
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
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                self.value  += self.get_power_norm();
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    self.q_vec[i1] /= self.tgt_vec[i1];
                }
                self.value += self.get_ks();
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.value  += self.get_vol_integral();
                    return;
                }
                else {
                    self.value  += self.get_vol_average();
                    return;
                }
            }
        }

        cat_list = CppStr::from("massFlux concGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre,  sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_mass_flux_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                    if self.category.s == "massFlux" {
                        self.q_vec[q_ind] = flux[self.component - 1].val;
                    }
                    else if self.category.s == "concGradient" {
                        self.q_vec[q_ind] = t_grad[self.component - 1].val;
                    }
                    if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                        this_el.get_volume_dfd0(&mut e_vol, st_pre, self.layer, sec_ar, dv_ar);
                        self.el_vol_vec[q_ind] = e_vol.val;
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                self.value  += self.get_power_norm();
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    self.q_vec[i1] /= self.tgt_vec[i1];
                }
                self.value += self.get_ks();
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.value += self.get_vol_integral();
                    return;
                }
                else {
                    self.value += self.get_vol_average();
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("flux tempGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre,  sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_flux_tgrad_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                    if self.category.s == "flux" {
                        self.q_vec[q_ind] = flux[self.component - 1].val;
                    }
                    else if self.category.s == "tempGradient" {
                        self.q_vec[q_ind] = t_grad[self.component - 1].val;
                    }
                    if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                        this_el.get_volume_dfd0(&mut e_vol, st_pre, self.layer, sec_ar, dv_ar);
                        self.el_vol_vec[q_ind] = e_vol.val;
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                self.value  += self.get_power_norm();
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    self.q_vec[i1] /= self.tgt_vec[i1];
                }
                self.value += self.get_ks();
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
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
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &Element;
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &el_ar[*eli];
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_volume_dfd0(&mut e_vol, st_pre, self.layer, sec_ar, dv_ar);
                    self.el_vol_vec[q_ind] = e_vol.val;
                    if self.category.s == "volume" {
                        self.q_vec[q_ind] = 1.0;
                    } else {
                        this_el.get_density_dfd0(&mut e_den,  self.layer,  sec_ar,  mat_ar, dv_ar);
                        self.q_vec[q_ind] = e_den.val;
                    }
                }
                q_ind += 1usize;
            }
            self.value += self.get_vol_integral();
            return;
        }
        
        return;
    }

    pub fn getd_ld_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>, time : f64, n_lgeom : bool,
        nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>,
        st_pre : &mut DiffDoub0StressPrereq, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        if time < self.active_time[0] || time > self.active_time[1] {
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
        let mut t_strain = [DiffDoub0::new(); 6];
        let mut d_strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut dsd_u = vec![DiffDoub0::new(); 288];
        let mut ded_u = vec![DiffDoub0::new(); 288];
        let mut dsd_t = vec![DiffDoub0::new(); 90];
        let mut dsd_c = vec![DiffDoub0::new(); 90];
        let mut dse_dend_u = [DiffDoub0::new(); 33];
        let mut dse_dend_c = [DiffDoub0::new(); 10];
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
        let mut dq_mat : &mut SparseMat;
        let mut exp : f64;

        let mut scr_v1 = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors"),
            Some(x) => &mut x.dat,
        };
        
        let mut cat_list = CppStr::from("displacement velocity acceleration concentration cdot temperature tdot");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.nd_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Node Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            for ndi in nd_sets[self.nd_set_ptr].labels.iter_mut() {
                dof_ind = nd_ar[*ndi].dof_index[self.component - 1];
                curr_rank = nd_ar[*ndi].sorted_rank;
                match self.category.s.as_str() {
                    "displacement" => self.d_qd_u.add_entry(q_ind, dof_ind, 1.0),
                    "velocity" => self.d_qd_v.add_entry(q_ind, dof_ind,  1.0),
                    "acceleration" => self.d_qd_a.add_entry(q_ind, dof_ind,  1.0),
                    "concentration" => self.d_qd_c.add_entry(q_ind, curr_rank, 1.0),
                    "cdot" => self.d_qd_cdot.add_entry(q_ind, curr_rank, 1.0),
                    "temperature" => self.d_qd_t.add_entry(q_ind, curr_rank, 1.0),
                    "tdot" => self.d_qd_tdot.add_entry(q_ind,   curr_rank,   1.0),
                    &_ => panic!("Error: unrecognized objective category '{}' in getd_ld_u()", self.category.s),
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "ks" {
                dq_mat = match self.category.s.as_str() {
                    "displacement" => &mut self.d_qd_u,
                    "velocity" => &mut self.d_qd_v,
                    "acceleration" => &mut self.d_qd_a,
                    "concentration" => &mut self.d_qd_c,
                    "cdot" => &mut self.d_qd_cdot,
                    "temperature" => &mut self.d_qd_t,
                    "tdot" => &mut self.d_qd_tdot,
                    &_ => panic!("Error: unrecognized objective category '{}' in getd_ld_u()", self.category.s),
                };
                for i1 in 0..self.q_len {
                    for me in dq_mat.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_ksd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("stress strain strainEnergy mises tsaiWu");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut t_strain, &mut d_strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                    this_el.d_stress_straind_u_dfd0(&mut dsd_u, &mut  ded_u, &mut  dsd_t, &mut dsd_c, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                    el_num_nds = this_el.num_nds;
                    el_dof_per_nd = this_el.dof_per_nd;
                    el_num_int_dof = this_el.num_int_dof;
                    el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
                    i1 = el_tot_dof * (self.component - 1);
                    i2 = el_num_nds * (self.component - 1);
                    if self.category.s == "stress" {
                        sub_vec_dfd0(&mut scr_v1, &mut  dsd_u,  i1,  i1 + el_tot_dof);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                        sub_vec_dfd0(&mut scr_v1, &mut  dsd_t,  i2,  i2 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                    else if self.category.s == "strain" {
                        sub_vec_dfd0(&mut scr_v1, &mut  ded_u,  i1,  i1 + el_tot_dof);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                    }
                    else if self.category.s == "strainEnergy" {
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
                    else if self.category.s == "mises" {
                        this_el.d_mises_du_dfd0(&mut dse_dend_u, &mut dse_dend_t, &mut dse_dend_c, &stress, &dsd_u, &dsd_t, &dsd_c);

                        ar_to_vec_dfd0(&mut dse_dend_u, &mut  scr_v1,  0,  33);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                        ar_to_vec_dfd0(&mut dse_dend_c, &mut  scr_v1,  0,  10);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_c, &mut  scr_v1,  true,  q_ind, nd_ar);
                        ar_to_vec_dfd0(&mut dse_dend_t, &mut  scr_v1,  0,  10);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                    else if self.category.s == "tsaiWu" {
                        this_el.d_tw_du_dfd0(&mut dse_dend_u, &mut dse_dend_t, &mut dse_dend_c, &stress, &dsd_u, &dsd_t, &dsd_c, self.layer, sec_ar, mat_ar, dv_ar);

                        ar_to_vec_dfd0(&mut dse_dend_u, &mut  scr_v1,  0,  33);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                        ar_to_vec_dfd0(&mut dse_dend_c, &mut  scr_v1,  0,  10);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_c, &mut  scr_v1,  true,  q_ind, nd_ar);
                        ar_to_vec_dfd0(&mut dse_dend_t, &mut  scr_v1,  0,  10);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    for me in self.d_qd_u.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    for me in self.d_qd_c.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    for me in self.d_qd_t.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_ksd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("sectionDef sectionFrcMom");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                    //this_el->d_stress_straind_u_dfd0(dsd_u, ded_u, dsd_t, spt, self.layer, n_lgeom, st_pre);
                    this_el.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut s_cent,  n_lgeom, st_pre);
                    this_el.d_def_frc_momd_u_dfd0(&mut ded_u, &mut  dsd_u, &mut  dsd_t, &mut dsd_c, &mut s_cent,  n_lgeom, st_pre);
                    el_num_nds = this_el.num_nds;
                    el_dof_per_nd = this_el.dof_per_nd;
                    el_num_int_dof = this_el.num_int_dof;
                    el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
                    i1 = el_tot_dof * (self.component - 1);
                    i2 = el_num_nds * (self.component - 1);
                    if self.category.s == "sectionFrcMom" {
                        sub_vec_dfd0(&mut scr_v1, &mut  dsd_u,  i1,  i1 + el_tot_dof);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                        sub_vec_dfd0(&mut scr_v1, &mut  dsd_t,  i2,  i2 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                    else if self.category.s == "sectionDef" {
                        sub_vec_dfd0(&mut scr_v1, &mut  ded_u,  i1,  i1 + el_tot_dof);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_u, &mut  scr_v1,  false,  q_ind, nd_ar);
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    for me in self.d_qd_u.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    for me in self.d_qd_c.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    for me in self.d_qd_t.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_ksd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }

        cat_list = CppStr::from("massFlux concGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_mass_flux_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                    this_el.d_mass_flux_dt_dfd0(&mut d_fd_t, &mut  d_tgd_t, &mut s_cent,  self.layer, st_pre);
                    el_num_nds = this_el.num_nds;
                    i1 = el_num_nds * (self.component - 1);
                    if self.category.s == "massFlux" {
                        sub_vec_dfd0(&mut scr_v1, &mut  d_fd_t,  i1,  i1 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_c, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                    else if self.category.s == "concGradient" {
                        sub_vec_dfd0(&mut scr_v1, &mut  d_tgd_t,  i1,  i1 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_c, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    for me in self.d_qd_c.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_ksd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        cat_list = CppStr::from("flux tempGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
                    this_el = &mut el_ar[*eli];
                    for i in 0..3 {
                        s_cent[i] = this_el.s_cent[i];
                    }
                    this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                    this_el.get_flux_tgrad_dfd0(&mut flux, &mut  t_grad, &mut s_cent,  self.layer, st_pre);
                    this_el.d_flux_tgradd_t_dfd0(&mut d_fd_t, &mut  d_tgd_t, &mut s_cent,  self.layer, st_pre);
                    el_num_nds = this_el.num_nds;
                    i1 = el_num_nds * (self.component - 1);
                    if self.category.s == "flux" {
                        sub_vec_dfd0(&mut scr_v1, &mut  d_fd_t,  i1,  i1 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                    else if self.category.s == "tempGradient" {
                        sub_vec_dfd0(&mut scr_v1, &mut  d_tgd_t,  i1,  i1 + el_num_nds);
                        this_el.put_vec_to_glob_mat_dfd0(&mut self.d_qd_t, &mut  scr_v1,  true,  q_ind, nd_ar);
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    for me in self.d_qd_t.matrix[i1].row_vec.iter_mut() {
                        me.value /= self.tgt_vec[i1];
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_ksd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
                else {
                    self.d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot);
                    return;
                }
            }
        }
        
        return;
    }

    pub fn getd_ld_d(&mut self, d_ld_d : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : &mut Vec<DesignVariable>, st_pre : &mut DiffDoub1StressPrereq) {
        if time < self.active_time[0] || time > self.active_time[1] {
            return;
        }
        
        let mut num_lay : usize;
        let mut fi : usize;
        let mut q_ind : usize;
        let mut dv_val = DiffDoub0::new();
        let mut strain = [DiffDoub1::new(); 6];
        let mut t_strain = [DiffDoub1::new(); 6];
        let mut d_strain = [DiffDoub1::new(); 6];
        let mut stress = [DiffDoub1::new(); 6];
        let mut se_den : f64;
        let mut def = [DiffDoub1::new(); 9];
        let mut frc_mom = [DiffDoub1::new(); 9];
        let mut flux = [DiffDoub1::new(); 3];
        let mut t_grad = [DiffDoub1::new(); 3];
        let mut e_vol = DiffDoub1::new();
        let mut e_den = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut mul_fact : f64;

        self.get_d_tgtd_d(el_ar, el_sets, sec_ar, mat_ar, dv_ar);
        
        let mut cat_list = CppStr::from("stress strain strainEnergy mises tsaiWu");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut tmp_dvi = vec![0usize; dv_ar.len()];
            let mut dv_len : usize;
            let mut dvi : usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter() {
                if el_ar[*eli].is_active {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                        this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                        this_el.get_stress_strain_dfd1(&mut stress, &mut  strain, &mut t_strain, &mut d_strain, &mut s_cent,  self.layer,  n_lgeom, st_pre);
                        if self.category.s == "stress" {
                            self.d_qd_d.add_entry(q_ind,  dvi,   stress[self.component - 1].dval);
                        }
                        else if self.category.s == "strain" {
                            self.d_qd_d.add_entry(q_ind,  dvi,   strain[self.component - 1].dval);
                        }
                        else if self.category.s == "strainEnergy" {
                            se_den = 0.0;
                            for i1 in 0..6 {
                                se_den  +=  stress[i1].val * strain[i1].dval + stress[i1].dval * strain[i1].val;
                            }
                            se_den  *=  0.5;
                            self.d_qd_d.add_entry(q_ind, dvi, se_den);
                        }
                        else if self.category.s == "mises" {
                            this_el.get_mises_dfd1(&mut tmp, &stress);
                            self.d_qd_d.add_entry(q_ind, dvi, tmp.dval);
                        }
                        else if self.category.s == "tsaiWu" {
                            this_el.get_tsai_wu_dfd1(&mut tmp, &stress, self.layer, sec_ar, mat_ar, dv_ar);
                            self.d_qd_d.add_entry(q_ind, dvi, tmp.dval);
                        }

                        if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                            this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                            self.d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
                        }

                        dv_ar[dvi].diff_val.set_val_2(dv_val.val, 0.0);
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    mul_fact = 1.0/self.tgt_vec[i1];
                    for me in self.d_qd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    mul_fact = -self.q_vec[i1]/(self.tgt_vec[i1]*self.tgt_vec[i1]);
                    for me in self.d_tgtd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_qd_d.add_matrix(&mut self.d_tgtd_d);
                self.d_ksd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
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
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len : usize;
            let mut dvi : usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                        this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                        //this_el->get_stress_strain_dfd0(stress, strain, spt, self.layer, n_lgeom, st_pre);
                        this_el.get_def_frc_mom_dfd1(&mut def, &mut  frc_mom, &mut s_cent, n_lgeom, st_pre);
                        if self.category.s == "sectionFrcMom" {
                            self.d_qd_d.add_entry(q_ind, dvi,   frc_mom[self.component - 1].dval);
                        }
                        else if self.category.s == "sectionDef" {
                            self.d_qd_d.add_entry(q_ind, dvi,   def[self.component - 1].dval);
                        }
                        if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                            num_lay = sec_ar[this_el.sect_ptr].layers.len();
                            if num_lay > 1 {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    mul_fact = 1.0/self.tgt_vec[i1];
                    for me in self.d_qd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    mul_fact = -self.q_vec[i1]/(self.tgt_vec[i1]*self.tgt_vec[i1]);
                    for me in self.d_tgtd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_qd_d.add_matrix(&mut self.d_tgtd_d);
                self.d_ksd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
                    self.d_vol_integrald_d(d_ld_d);
                    return;
                }
                else {
                    self.d_vol_averaged_d(d_ld_d);
                    return;
                }
            }
        }

        cat_list = CppStr::from("massFlux concGradient");
        fi = cat_list.find(&self.category.s.as_str());
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len : usize;
            let mut dvi : usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                        this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                        this_el.get_mass_flux_dfd1(&mut flux, &mut  t_grad, &mut s_cent, self.layer, st_pre);
                        if self.category.s == "massFlux" {
                            self.d_qd_d.add_entry(q_ind, dvi, flux[self.component - 1].dval);
                        }
                        else if self.category.s == "concGradient" {
                            self.d_qd_d.add_entry(q_ind, dvi, t_grad[self.component - 1].dval);
                        }
                        if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                            this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                            self.d_vd_d.add_entry(q_ind, dvi,  e_vol.dval);
                        }
                        dv_ar[dvi].diff_val.set_val_2(dv_val.val,    0.0);
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    mul_fact = 1.0/self.tgt_vec[i1];
                    for me in self.d_qd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    mul_fact = -self.q_vec[i1]/(self.tgt_vec[i1]*self.tgt_vec[i1]);
                    for me in self.d_tgtd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_qd_d.add_matrix(&mut self.d_tgtd_d);
                self.d_ksd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
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
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &mut Element;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len : usize;
            let mut dvi : usize;
            let mut s_cent : [f64; 3] = [0.0f64; 3];
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                        this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                        this_el.get_flux_tgrad_dfd1(&mut flux, &mut  t_grad, &mut s_cent, self.layer, st_pre);
                        if self.category.s == "flux" {
                            self.d_qd_d.add_entry(q_ind, dvi,   flux[self.component - 1].dval);
                        }
                        else if self.category.s == "tempGradient" {
                            self.d_qd_d.add_entry(q_ind, dvi,   t_grad[self.component - 1].dval);
                        }
                        if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                            this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                            self.d_vd_d.add_entry(q_ind, dvi,  e_vol.dval);
                        }
                        dv_ar[dvi].diff_val.set_val_2(dv_val.val,    0.0);
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                    }
                }
                q_ind += 1usize;
            }
            if self.optr.s == "powerNorm" {
                for i1 in 0..self.q_len {
                    self.err_norm_vec[i1] = powf(self.q_vec[i1] - self.tgt_vec[i1], self.expnt - 1.0);
                }
                self.d_power_normd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "ks" {
                for i1 in 0..self.q_len {
                    mul_fact = 1.0/self.tgt_vec[i1];
                    for me in self.d_qd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    mul_fact = -self.q_vec[i1]/(self.tgt_vec[i1]*self.tgt_vec[i1]);
                    for me in self.d_tgtd_d.matrix[i1].row_vec.iter_mut() {
                        me.value *= mul_fact;
                    }
                    self.err_norm_vec[i1] = (self.expnt*self.q_vec[i1]).exp();
                }
                self.d_qd_d.add_matrix(&mut self.d_tgtd_d);
                self.d_ksd_d(d_ld_d);
                return;
            }
            else if self.optr.s == "volumeIntegral" || self.optr.s == "volumeAverage" {
                if self.optr.s == "volumeIntegral" {
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
        if fi < MAX_INT {
            if self.el_set_ptr == MAX_INT {
                let mut err_str = format!("Error: Objective terms of category '{}' must have a valid Element Set specified.\n", self.category.s);
                err_str = format!("{} Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.", err_str);
                panic!("{}", err_str);
            }
            q_ind = 0;
            let mut this_el : &Element;
            let mut tmp_dv = vec![0usize; dv_ar.len()];
            let mut dv_len : usize;
            let mut dvi : usize;
            for eli in el_sets[self.el_set_ptr].labels.iter_mut() {
                if el_ar[*eli].is_active {
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
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                        this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
                        this_el.get_volume_dfd1(&mut e_vol, st_pre,  self.layer, sec_ar, dv_ar);
                        self.d_vd_d.add_entry(q_ind, dvi,   e_vol.dval);
                        if self.category.s == "mass" {
                            this_el.get_density_dfd1(&mut e_den,  self.layer, sec_ar, mat_ar, dv_ar);
                            self.d_qd_d.add_entry(q_ind, dvi,   e_den.dval);
                        }
                        dv_ar[dvi].diff_val.set_val_2(dv_val.val,  0.0);
                        for nd in nd_ar.iter_mut() {
                            nd.calc_crd_dfd1(dv_ar);
                        }
                    }
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

    pub fn calculated_ld_u(&mut self, d_ld_u : &mut Vec<f64>, d_ld_v : &mut Vec<f64>, d_ld_a : &mut Vec<f64>, d_ld_c : &mut Vec<f64>, d_ld_cdot : &mut Vec<f64>, d_ld_t : &mut Vec<f64>, d_ld_tdot : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, nd_sets : &mut Vec<Set>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>,
        st_pre : &mut DiffDoub0StressPrereq, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        for tm in self.terms.iter_mut() {
            tm.getd_ld_u(d_ld_u, d_ld_v, d_ld_a, d_ld_c, d_ld_cdot, d_ld_t, d_ld_tdot, time, n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre, scr_dfd);
        }
        return;
    }

    pub fn calculated_ld_d(&mut self, d_ld_d : &mut Vec<f64>, time : f64, n_lgeom : bool, nd_ar : &mut Vec<Node>, el_ar : &mut Vec<Element>, el_sets : &mut Vec<Set>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : &mut Vec<DesignVariable>, st_pre : &mut DiffDoub1StressPrereq) {
        for tm in self.terms.iter_mut() {
            tm.getd_ld_d(d_ld_d,  time,  n_lgeom, nd_ar, el_ar, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
        }
        return;
    }

}


