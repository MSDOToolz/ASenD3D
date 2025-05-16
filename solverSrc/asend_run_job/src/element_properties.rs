use crate::element::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::design_var::*;
use crate::node::*;
use crate::section::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;


impl Element {

    //dup1
    pub fn get_gen_prop_dfd0(&self, prop : &mut DiffDoub0, prop_key : &mut CppStr, dv_ar : & Vec<DesignVariable>) {
        let mut d_vind : usize;
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut this_dv : & DesignVariable;
        
        for dv in self.design_vars.iter() {
            d_vind = dv.int_dat;
            this_dv = & dv_ar[d_vind];
            if this_dv.category.s == prop_key.s {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                prop.add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_f_per_mass_dfd0(&self, f_vec : &mut Vec<DiffDoub0>, dv_ar : &Vec<DesignVariable>) {
        let mut d_vind : usize;
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut this_dv : & DesignVariable;
        let mut comp : usize;
        
        for dv in self.design_vars.iter() {
            d_vind = dv.int_dat;
            this_dv = & dv_ar[d_vind];
            if this_dv.category.s == "bodyForce" {
                this_dv.get_value_dfd0(&mut dv_val);
                comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                //prop.add(& dv_val);
                f_vec[comp].add(&dv_val);
            }
        }
    }

    pub fn get_layer_thk_z_dfd0(&self, lay_thk : &mut Vec<DiffDoub0>, lay_z : &mut Vec<DiffDoub0>, z_offset : &mut DiffDoub0, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        //z_offset = 1: upper z surface is reference plane
        //z_offset = -1: lower z surface is reference plane
        let mut layi : usize =  0;
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut tot_thk = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut z_crd = DiffDoub0::new();
        let mut z_next = DiffDoub0::new();
        let mut z_mid = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        tot_thk.set_val(0.0);
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            lay_thk[layi].set_val(lay.thickness);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thickness" && this_dv.layer == layi {
                    this_dv.get_value_dfd0(&mut dv_val);
                    coef.set_val(dv.doub_dat);
                    dv_val.mult(& coef);
                    lay_thk[layi].add(& dv_val);
                }
            }
            tot_thk.add(& lay_thk[layi]);
            layi += 1usize;
        }
        
        z_offset.set_val(this_sec.z_offset);
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "zOffset" {
                this_dv.get_value_dfd0(&mut dv_val);
                coef.set_val(dv.doub_dat);
                dv_val.mult(& coef);
                z_offset.add(& dv_val);
            }
        }
        
        tmp.set_val(1.0);
        tmp.add(z_offset);
        z_crd.set_val(-0.5);
        z_crd.mult(& tot_thk);
        z_crd.mult(& tmp);
        
        layi = 0;
        for _lay in this_sec.layers.iter() {
            z_next.set_val_dfd0(& z_crd);
            z_next.add(& lay_thk[layi]);
            tmp.set_val_dfd0(& z_crd);
            tmp.add(& z_next);
            z_mid.set_val(0.5);
            z_mid.mult(& tmp);
            lay_z[layi].set_val_dfd0(& z_mid);
            layi += 1usize;
            z_crd.set_val_dfd0(& z_next);
        }
        
        return;
    }

    pub fn get_layer_q_dfd0(&self, lay_q : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub0::new(); 3];
        let mut poisson_dv = [DiffDoub0::new(); 3];
        let mut shear_mod_dv = [DiffDoub0::new(); 3];
        let mut smat = [DiffDoub0::new(); 9];
        let mut qmat = [DiffDoub0::new(); 9];
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut x_vec = [DiffDoub0::new(); 3];
        let mut b_vec = [DiffDoub0::new(); 3];
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        let mut tmp = DiffDoub0::new();
        
        i2 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "modulus" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd0(&mut dv_val);
                        dv_val.mult(& coef);
                        modulus_dv[dv_comp - 1].add(& dv_val);
                    }
                    else if dv_cat.s == "poissonRatio" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd0(&mut dv_val);
                        dv_val.mult(& coef);
                        poisson_dv[dv_comp - 1].add(& dv_val);
                    }
                    else if dv_cat.s == "shearModulus" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd0(&mut dv_val);
                        dv_val.mult(& coef);
                        shear_mod_dv[dv_comp - 1].add(& dv_val);
                    }
                }
            }
            
            for i1 in 0..9 {
                smat[i1].set_val(0.0);
            }
            smat[0].set_val(1.0);
            smat[0].dvd(& modulus_dv[0]);
            smat[1].set_val_dfd0(& poisson_dv[0]);
            smat[1].neg();
            smat[1].dvd(& modulus_dv[0]);
            tmp.set_val_dfd0(& smat[1]);
            smat[3].set_val_dfd0(& tmp);
            smat[4].set_val(1.0);
            smat[4].dvd(& modulus_dv[1]);
            smat[8].set_val(1.0);
            smat[8].dvd(& shear_mod_dv[0]);
            get_det_inv_ar_dfd0(&mut coef, &mut  qmat, &mut  smat,  3,  0, &mut  x_vec, &mut  b_vec);
            
            for i1 in 0..9 {
                lay_q[i2].set_val_dfd0(& qmat[i1]);
                i2 += 1usize;
            }
            
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_d_dfd0(&self, lay_d : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut damp_mat_dv = [DiffDoub0::new(); 36];
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        let mut tmp = DiffDoub0::new();
        
        i2 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..36 {
                damp_mat_dv[i1].set_val(this_mat.damping[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "dampingMat" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd0(&mut dv_val);
                        dv_val.mult(& coef);
                        damp_mat_dv[dv_comp].add(& dv_val);
                    }
                }
            }
            
            for i3 in 1..6 {
                i5 = 6 * i3;// lower tri term
                i6 = i3;// upper tri term
                for _i4 in 0..i3 {
                    tmp.set_val_dfd0(& damp_mat_dv[i6]);
                    damp_mat_dv[i5].set_val_dfd0(& tmp);
                    i5 += 1usize;
                    i6  +=  6;
                }
            }
            
            lay_d[i2].set_val_dfd0(& damp_mat_dv[0]);
            lay_d[i2 + 1].set_val_dfd0(& damp_mat_dv[1]);
            lay_d[i2 + 2].set_val_dfd0(& damp_mat_dv[3]);
            lay_d[i2 + 3].set_val_dfd0(& damp_mat_dv[6]);
            lay_d[i2 + 4].set_val_dfd0(& damp_mat_dv[7]);
            lay_d[i2 + 5].set_val_dfd0(& damp_mat_dv[9]);
            lay_d[i2 + 6].set_val_dfd0(& damp_mat_dv[18]);
            lay_d[i2 + 7].set_val_dfd0(& damp_mat_dv[19]);
            lay_d[i2 + 8].set_val_dfd0(& damp_mat_dv[21]);
            
            layi += 1usize;
            i2  +=  9;
        }
        return;
    }

    pub fn get_layer_angle_dfd0(&self, lay_ang : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut angle = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            angle.set_val(lay.angle);
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "angle" {
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd0(&mut dv_val);
                        dv_val.mult(& coef);
                        angle.add(& dv_val);
                    }
                }
            }
            lay_ang[layi].set_val_dfd0(& angle);
            
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_th_exp_dfd0(&self, lay_th_exp : &mut Vec<DiffDoub0>, lay_diff_exp : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut t_exp_dv = [DiffDoub0::new(); 6];
        let mut d_exp_dv = [DiffDoub0::new(); 6];
        let mut dv_comp : usize;
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        
        layi = 0;
        i2 = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..6 {
                t_exp_dv[i1].set_val(this_mat.expansion[i1]);
                d_exp_dv[i1].set_val(this_mat.diff_exp[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermalExp" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    t_exp_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "diffExp" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    d_exp_dv[dv_comp].add(& dv_val);
                }
            }
            lay_th_exp[i2].set_val_dfd0(& t_exp_dv[0]);
            lay_diff_exp[i2].set_val_dfd0(& d_exp_dv[0]);
            i2 += 1usize;
            lay_th_exp[i2].set_val_dfd0(& t_exp_dv[1]);
            lay_diff_exp[i2].set_val_dfd0(& d_exp_dv[1]);
            i2 += 1usize;
            lay_th_exp[i2].set_val_dfd0(& t_exp_dv[3]);
            lay_diff_exp[i2].set_val_dfd0(& d_exp_dv[3]);
            i2 += 1usize;
            layi += 1usize;
        }
        return;
    }

    pub fn get_layer_einit_dfd0(&self, lay_einit : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut e0_dv = [DiffDoub0::new(); 6];
        let mut dv_comp : usize;
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        i2 = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for _lay in this_sec.layers.iter() {
            for i1 in 0..6 {
                e0_dv[i1].set_val(0.0);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "initialStrain" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    e0_dv[dv_comp].add(& dv_val);
                }
            }
            lay_einit[i2].set_val_dfd0(& e0_dv[0]);
            i2 += 1usize;
            lay_einit[i2].set_val_dfd0(& e0_dv[1]);
            i2 += 1usize;
            lay_einit[i2].set_val_dfd0(& e0_dv[3]);
            i2 += 1usize;
            layi += 1usize;
        }
        return;
    }

    pub fn get_layer_den_dfd0(&self, layer_den : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut mat_den : f64;
        let mut den_dv = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            mat_den = mat_ar[lay.mat_ptr].density;
            den_dv.set_val(mat_den);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "density" && this_dv.layer == layi {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    den_dv.add(& dv_val);
                }
            }
            layer_den[layi].set_val_dfd0(& den_dv);
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_cond_dfd0(&self, lay_cond : &mut Vec<DiffDoub0>, lay_diff : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i1 : usize;
        let mut layi : usize;
        let mut cond_dv = [DiffDoub0::new(); 6];
        let mut diff_dv = [DiffDoub0::new(); 6];
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        
        i1 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i2 in 0..6 {
                cond_dv[i2].set_val(this_mat.conductivity[i2]);
                diff_dv[i2].set_val(this_mat.diffusivity[i2]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermalCond" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "diffusivity" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv[dv_comp].add(& dv_val);
                }
            }
            lay_cond[i1].set_val_dfd0(& cond_dv[0]);
            lay_cond[i1 + 1].set_val_dfd0(& cond_dv[3]);
            lay_cond[i1 + 2].set_val_dfd0(& cond_dv[4]);
            lay_cond[i1 + 3].set_val_dfd0(& cond_dv[3]);
            lay_cond[i1 + 4].set_val_dfd0(& cond_dv[1]);
            lay_cond[i1 + 5].set_val_dfd0(& cond_dv[5]);
            lay_cond[i1 + 6].set_val_dfd0(& cond_dv[4]);
            lay_cond[i1 + 7].set_val_dfd0(& cond_dv[5]);
            lay_cond[i1 + 8].set_val_dfd0(& cond_dv[2]);

            lay_diff[i1].set_val_dfd0(& diff_dv[0]);
            lay_diff[i1 + 1].set_val_dfd0(& diff_dv[3]);
            lay_diff[i1 + 2].set_val_dfd0(& diff_dv[4]);
            lay_diff[i1 + 3].set_val_dfd0(& diff_dv[3]);
            lay_diff[i1 + 4].set_val_dfd0(& diff_dv[1]);
            lay_diff[i1 + 5].set_val_dfd0(& diff_dv[5]);
            lay_diff[i1 + 6].set_val_dfd0(& diff_dv[4]);
            lay_diff[i1 + 7].set_val_dfd0(& diff_dv[5]);
            lay_diff[i1 + 8].set_val_dfd0(& diff_dv[2]);

            i1  +=  9;
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_spec_heat_dfd0(&self, lay_sh : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut mat_sh : f64;
        let mut sh_dv = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            mat_sh = mat_ar[lay.mat_ptr].spec_heat;
            sh_dv.set_val(mat_sh);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "specHeat" && this_dv.layer == layi {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    sh_dv.add(& dv_val);
                }
            }
            lay_sh[layi].set_val_dfd0(& sh_dv);
            layi += 1usize;
        }
        
        return;
    }

    pub fn transform_strain_dfd0(&self, stn_new : &mut [DiffDoub0], stn_orig : &mut [DiffDoub0], angle : &mut DiffDoub0) {
        let mut angle_rad = DiffDoub0::new();
        let mut a11 = DiffDoub0::new();
        let mut a12 = DiffDoub0::new();
        let mut a21 = DiffDoub0::new();
        let mut a22 = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut te = [DiffDoub0::new(); 9];
        
        
        angle_rad.set_val(R_PI_0180);
        angle_rad.mult(angle);
        a11.set_val_dfd0(& angle_rad);
        a11.cs();
        a21.set_val_dfd0(& angle_rad);
        a21.sn();
        a12.set_val_dfd0(& a21);
        a12.neg();
        a22.set_val_dfd0(& a11);
        
        te[0].set_val_dfd0(& a11);
        te[0].sqr();
        te[1].set_val_dfd0(& a12);
        te[1].sqr();
        te[2].set_val_dfd0(& a11);
        te[2].mult(& a12);
        te[3].set_val_dfd0(& a21);
        te[3].sqr();
        te[4].set_val_dfd0(& a22);
        te[4].sqr();
        te[5].set_val_dfd0(& a22);
        te[5].mult(& a21);
        te[6].set_val(2.0);
        te[6].mult(& a11);
        te[6].mult(& a21);
        te[7].set_val(2.0);
        te[7].mult(& a12);
        te[7].mult(& a22);
        te[8].set_val_dfd0(& a11);
        te[8].mult(& a22);
        tmp.set_val_dfd0(& a12);
        tmp.mult(& a21);
        te[8].add(& tmp);
        
        mat_mul_ar_dfd0(stn_new, &mut te, stn_orig,  3,  3,  1);
        
        return;
    }

    pub fn transform_q_dfd0(&self, q_new : &mut [DiffDoub0], q_orig : &mut [DiffDoub0], angle : &mut DiffDoub0) {
        let mut angle_rad = DiffDoub0::new();
        let mut a11 = DiffDoub0::new();
        let mut a12 = DiffDoub0::new();
        let mut a21 = DiffDoub0::new();
        let mut a22 = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut ts = [DiffDoub0::new(); 9];
        let mut te = [DiffDoub0::new(); 9];
        let mut te_inv = [DiffDoub0::new(); 9];
        let mut x_vec = [DiffDoub0::new(); 3];
        let mut b_vec = [DiffDoub0::new(); 3];
        
        
        angle_rad.set_val(R_PI_0180);
        angle_rad.mult(angle);
        a11.set_val_dfd0(& angle_rad);
        a11.cs();
        a21.set_val_dfd0(& angle_rad);
        a21.sn();
        a12.set_val_dfd0(& a21);
        a12.neg();
        a22.set_val_dfd0(& a11);
        
        ts[0].set_val_dfd0(& a11);
        ts[0].sqr();
        ts[1].set_val_dfd0(& a12);
        ts[1].sqr();
        ts[2].set_val(2.0);
        ts[2].mult(& a11);
        ts[2].mult(& a12);
        ts[3].set_val_dfd0(& a21);
        ts[3].sqr();
        ts[4].set_val_dfd0(& a22);
        ts[4].sqr();
        ts[5].set_val(2.0);
        ts[5].mult(& a22);
        ts[5].mult(& a21);
        ts[6].set_val_dfd0(& a11);
        ts[6].mult(& a21);
        ts[7].set_val_dfd0(& a12);
        ts[7].mult(& a22);
        ts[8].set_val_dfd0(& a11);
        ts[8].mult(& a22);
        tmp.set_val_dfd0(& a12);
        tmp.mult(& a21);
        ts[8].add(& tmp);
        
        for i1 in 0..9 {
            te[i1].set_val_dfd0(& ts[i1]);
        }
        tmp.set_val(0.5);
        te[2].mult(& tmp);
        te[5].mult(& tmp);
        tmp.set_val(2.0);
        te[6].mult(& tmp);
        te[7].mult(& tmp);
        
        get_det_inv_ar_dfd0(&mut coef, &mut  te_inv, &mut  te,  3,  0, &mut  x_vec, &mut  b_vec);
        
        mat_mul_ar_dfd0(&mut te, q_orig, &mut te_inv,  3,  3,  3);
        mat_mul_ar_dfd0(q_new, &mut  ts, &mut  te,  3,  3,  3);
        
        return;
    }

    pub fn get_solid_stiff_dfd0(&self, cmat : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub0::new(); 3];
        let mut poisson_dv = [DiffDoub0::new(); 3];
        let mut shear_mod_dv = [DiffDoub0::new(); 3];
        let mut smat = [DiffDoub0::new(); 36];
        let mut ctmp = [DiffDoub0::new(); 36];
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut x_vec = [DiffDoub0::new(); 6];
        let mut b_vec = [DiffDoub0::new(); 6];
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub0::new();
        
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_mat.stiffness[0] > 0.0 {
            for i1 in 0..36 {
                cmat[i1].set_val(this_mat.stiffness[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "stiffnessMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    cmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& cmat[i5]);
                    cmat[i4].set_val_dfd0(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                dv_cat = this_dv.category.clone();
                if dv_cat.s == "modulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    modulus_dv[dv_comp - 1].add(& dv_val);
                }
                else if dv_cat.s == "poissonRatio" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    poisson_dv[dv_comp - 1].add(& dv_val);
                }
                else if dv_cat.s == "shearModulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    shear_mod_dv[dv_comp - 1].add(& dv_val);
                }
            }
            for i1 in 0..36 {
                smat[i1].set_val(0.0);
            }
            
            smat[0].set_val(1.0);
            smat[0].dvd(& modulus_dv[0]);
            smat[1].set_val_dfd0(& poisson_dv[0]);
            smat[1].neg();
            smat[1].dvd(& modulus_dv[0]);
            smat[2].set_val_dfd0(& poisson_dv[1]);
            smat[2].neg();
            smat[2].dvd(& modulus_dv[0]);
            tmp.set_val_dfd0(& smat[1]);
            smat[6].set_val_dfd0(& tmp);
            smat[7].set_val(1.0);
            smat[7].dvd(& modulus_dv[1]);
            smat[8].set_val_dfd0(& poisson_dv[2]);
            smat[8].neg();
            smat[8].dvd(& modulus_dv[1]);
            tmp.set_val_dfd0(& smat[2]);
            smat[12].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& smat[8]);
            smat[13].set_val_dfd0(& tmp);
            smat[14].set_val(1.0);
            smat[14].dvd(& modulus_dv[2]);
            smat[21].set_val(1.0);
            smat[21].dvd(& shear_mod_dv[0]);
            smat[28].set_val(1.0);
            smat[28].dvd(& shear_mod_dv[1]);
            smat[35].set_val(1.0);
            smat[35].dvd(& shear_mod_dv[2]);
            
            get_det_inv_ar_dfd0(&mut coef, &mut  ctmp, &mut  smat,  6,  0, &mut  x_vec, &mut  b_vec);
            ar_to_vec_dfd0(&mut ctmp, cmat,  0,  36);
        }
        
        return;
    }

    pub fn get_abd_dfd0(&self, cmat : &mut Vec<DiffDoub0>, lay_thk : &mut Vec<DiffDoub0>, lay_z : &mut Vec<DiffDoub0>, lay_q : &mut Vec<DiffDoub0>, lay_ang : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let num_lay : usize;
        let mut z_max = DiffDoub0::new();
        let mut z_min = DiffDoub0::new();
        let mut thk = DiffDoub0::new();
        let mut qmat = [DiffDoub0::new(); 9];
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        
        for i1 in 0..81 {
            cmat[i1].set_val(0.0);
        }
        
        num_lay = sec_ar[self.sect_ptr].layers.len();
        i2 = 0;
        for i1 in 0..num_lay {
            thk.set_val(0.5);
            thk.mult(& lay_thk[i1]);
            z_min.set_val_dfd0(& lay_z[i1]);
            z_min.sub(& thk);
            z_max.set_val_dfd0(& lay_z[i1]);
            z_max.add(& thk);
            self.transform_q_dfd0(&mut qmat, &mut lay_q[i2..], &mut lay_ang[i1]);
            
            // a matrix portion
            tmp.set_val_dfd0(& z_max);
            tmp.sub(& z_min);
            for i3 in 0..3 {
                i5 = 9 * i3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    tmp2.set_val_dfd0(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            
            // b matrix portion
            tmp.set_val_dfd0(& z_max);
            tmp.sqr();
            tmp2.set_val_dfd0(& z_min);
            tmp2.sqr();
            tmp.sub(& tmp2);
            tmp2.set_val(0.5);
            tmp.mult(& tmp2);
            for i3 in 0..3 {
                i5 = 9 * i3 + 3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    //i6 = i3*3 + i4;
                    //i5 = i3*9 + (i4 + 3)
                    tmp2.set_val_dfd0(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            
            // d matrix portion				
            tmp.set_val_dfd0(& z_max);
            tmp.sqr();
            tmp.mult(& z_max);
            tmp2.set_val_dfd0(& z_min);
            tmp2.sqr();
            tmp2.mult(& z_min);
            tmp.sub(& tmp2);
            tmp2.set_val(R_1O3);
            tmp.mult(& tmp2);
            for i3 in 0..3 {
                i5 = 9 * (i3 + 3) + 3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    tmp2.set_val_dfd0(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            i2  +=  9;
        }
        
        for i1 in 1..6 {
            i3 = 9 * i1;
            i4 = i1;
            for _i2 in 0..i1 {
                //i3 = i1*9 + i2
                //i4 = i2*9 + i1
                tmp.set_val_dfd0(& cmat[i4]);
                cmat[i3].set_val_dfd0(& tmp);
                i3 += 1usize;
                i4  +=  9;
            }
        }
        
        tmp.set_val(1.0);
        tmp.mult(& cmat[20]);
        cmat[60].set_val_dfd0(& tmp);
        cmat[70].set_val_dfd0(& tmp);
        cmat[80].set_val_dfd0(& tmp);
        
        return;
    }

    pub fn get_beam_stiff_dfd0(&self, cmat : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub0::new(); 3];
        let mut shear_mod_dv = [DiffDoub0::new(); 3];
        let mut area_dv = DiffDoub0::new();
        let mut idv = [DiffDoub0::new(); 5];
        let mut jdv = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub0::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.stiffness[0] > 0.0 {
            for i1 in 0..36 {
                cmat[i1].set_val(this_sec.stiffness[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "stiffnessMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    cmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& cmat[i5]);
                    cmat[i4].set_val_dfd0(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            jdv.set_val(this_sec.polar_moment);
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "modulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    modulus_dv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "shearModulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    shear_mod_dv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if this_dv.category.s == "areaMoment" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            
            for i1 in 0..6 {
                i3 = 7 * i1;
                for _i2 in i1..6 {
                    cmat[i3].set_val(0.0);
                    i3 += 1usize;
                }
            }
            
            cmat[0].set_val_dfd0(& modulus_dv[0]);
            cmat[0].mult(& area_dv);
            cmat[4].set_val_dfd0(& modulus_dv[0]);
            cmat[4].mult(& idv[0]);
            cmat[5].set_val_dfd0(& modulus_dv[0]);
            cmat[5].neg();
            cmat[5].mult(& idv[1]);
            cmat[7].set_val_dfd0(& shear_mod_dv[0]);
            cmat[7].mult(& area_dv);
            cmat[9].set_val_dfd0(& shear_mod_dv[0]);
            cmat[9].neg();
            cmat[9].mult(& idv[0]);
            cmat[14].set_val_dfd0(& shear_mod_dv[1]);
            cmat[14].mult(& area_dv);
            cmat[15].set_val_dfd0(& shear_mod_dv[2]);
            cmat[15].mult(& idv[1]);
            cmat[21].set_val_dfd0(& shear_mod_dv[0]);
            cmat[21].mult(& jdv);
            cmat[28].set_val_dfd0(& modulus_dv[0]);
            cmat[28].mult(& idv[2]);
            cmat[29].set_val_dfd0(& modulus_dv[0]);
            cmat[29].neg();
            cmat[29].mult(& idv[4]);
            cmat[35].set_val_dfd0(& modulus_dv[0]);
            cmat[35].mult(& idv[3]);
            
            for i1 in 1..6 {
                i3 = 6 * i1;
                i4 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& cmat[i4]);
                    cmat[i3].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        return;
    }

    pub fn get_thermal_exp_dfd0(&self, th_exp : &mut Vec<DiffDoub0>, diff_exp : &mut Vec<DiffDoub0>, einit : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_comp : usize;
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec = &sec_ar[self.sect_ptr];
        let this_mat = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..6 {
            th_exp[i1].set_val(this_mat.expansion[i1]);
            diff_exp[i1].set_val(this_mat.diff_exp[i1]);
            einit[i1].set_val(0.0);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "thermalExp" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& tmp);
                th_exp[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "diffExp" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& tmp);
                diff_exp[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "initialStrain" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& tmp);
                einit[dv_comp].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_shell_exp_load_dfd0(&self, exp_ld : &mut Vec<DiffDoub0>, diff_ld : &mut Vec<DiffDoub0>, e0_ld : &mut Vec<DiffDoub0>, lay_thk : &mut Vec<DiffDoub0>, lay_z : &mut Vec<DiffDoub0>, 
        lay_q : &mut Vec<DiffDoub0>, lay_th_exp : &mut Vec<DiffDoub0>, lay_diff_exp : &mut Vec<DiffDoub0>, lay_einit : &mut Vec<DiffDoub0>, 
        lay_ang : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {

        let num_lay : usize;
        let mut qi : usize;
        let mut exi : usize;
        let mut sect_q = [DiffDoub0::new(); 9];
        let mut sect_te = [DiffDoub0::new(); 3];
        let mut sect_de = [DiffDoub0::new(); 3];
        let mut sect_e0 = [DiffDoub0::new(); 3];
        let mut qte_prod = [DiffDoub0::new(); 3];
        let mut qde_prod = [DiffDoub0::new(); 3];
        let mut qe0_prod = [DiffDoub0::new(); 3];
        let mut this_q = [DiffDoub0::new(); 9];
        let mut this_te = [DiffDoub0::new(); 3];
        let mut this_de = [DiffDoub0::new(); 3];
        let mut this_e0 = [DiffDoub0::new(); 3];
        let mut z_min = DiffDoub0::new();
        let mut z_max = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        
        for i1 in 0..6 {
            exp_ld[i1].set_val(0.0);
            e0_ld[i1].set_val(0.0);
        }
        
        num_lay = sec_ar[self.sect_ptr].layers.len();
        qi = 0;
        exi = 0;
        for layi in 0..num_lay {
            vec_to_ar_dfd0(&mut this_q, lay_q,  qi,  qi + 9);
            self.transform_q_dfd0(&mut sect_q, &mut  this_q, &mut lay_ang[layi]);
            
            vec_to_ar_dfd0(&mut this_te, lay_th_exp,  exi,  exi + 3);
            self.transform_strain_dfd0(&mut sect_te, &mut  this_te, &mut lay_ang[layi]);

            vec_to_ar_dfd0(&mut this_de, lay_diff_exp, exi, exi + 3);
            self.transform_strain_dfd0(&mut sect_de, &mut this_de, &mut lay_ang[layi]);
            
            vec_to_ar_dfd0(&mut this_e0, lay_einit,  exi,  exi + 3);
            self.transform_strain_dfd0(&mut sect_e0, &mut  this_e0, &mut  lay_ang[layi]);
            
            mat_mul_ar_dfd0(&mut qte_prod, &mut  sect_q, &mut  sect_te,  3,  3,  1);
            mat_mul_ar_dfd0(&mut qde_prod, &mut sect_q, &mut sect_de, 3, 3, 1);
            mat_mul_ar_dfd0(&mut qe0_prod, &mut  sect_q, &mut  sect_e0,  3,  3,  1);
            
            for i1 in 0..3 {
                tmp.set_val_dfd0(& qte_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                exp_ld[i1].add(& tmp);

                tmp.set_val_dfd0(& qde_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                diff_ld[i1].add(& tmp);

                tmp.set_val_dfd0(& qe0_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                e0_ld[i1].add(& tmp);
            }
            
            tmp.set_val(0.5);
            tmp.mult(& lay_thk[layi]);
            z_min.set_val_dfd0(& lay_z[layi]);
            z_min.sub(& tmp);
            z_min.sqr();
            z_max.set_val_dfd0(& lay_z[layi]);
            z_max.add(& tmp);
            z_max.sqr();
            tmp.set_val(0.5);
            tmp2.set_val_dfd0(& z_max);
            tmp2.sub(& z_min);
            tmp.mult(& tmp2);// tmp = 0.5*(z_max^2 - z_min^2)
            for i1 in 0..3 {
                tmp2.set_val_dfd0(& qte_prod[i1]);
                tmp2.mult(& tmp);
                exp_ld[i1 + 3].add(& tmp2);

                tmp2.set_val_dfd0(& qde_prod[i1]);
                tmp2.mult(& tmp);
                diff_ld[i1 + 3].add(& tmp2);

                tmp2.set_val_dfd0(& qe0_prod[i1]);
                tmp2.mult(& tmp);
                e0_ld[i1 + 3].add(& tmp2);
            }
            qi  +=  9;
            exi  +=  3;
        }
        
        return;
    }

    pub fn get_beam_exp_load_dfd0(&self, exp_ld : &mut Vec<DiffDoub0>, diff_ld : &mut Vec<DiffDoub0>, e0_ld : &mut Vec<DiffDoub0>,
         sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {

        let mut i2 : usize;
        let mut dv_val = DiffDoub0::new();
        let mut cat : CppStr;
        let mut cat_list : CppStr;
        let mut dv_comp : usize;
        let mut mod_dv = [DiffDoub0::new(); 3];
        let mut shr_mod_dv = [DiffDoub0::new(); 3];
        let mut te_coef_dv = [DiffDoub0::new(); 6];
        let mut de_coef_dv = [DiffDoub0::new(); 6];
        let mut e0_dv = [DiffDoub0::new(); 6];
        let mut area_dv = DiffDoub0::new();
        let mut idv = [DiffDoub0::new(); 5];
        let mut qmat = [DiffDoub0::new(); 9];
        let mut qte = [DiffDoub0::new(); 3];
        let mut qde = [DiffDoub0::new(); 3];
        let mut qe0 = [DiffDoub0::new(); 3];
        let mut dedgu = [DiffDoub0::new(); 18];
        let mut tmp_exp = [DiffDoub0::new(); 6];
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub0::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];

        if this_sec.exp_load_coef[0] > 0.0 || this_sec.exp_load_coef[0] > 0.0 {
            for i1 in 0..6 {
                exp_ld[i1].set_val(this_sec.exp_load_coef[i1]);
                diff_ld[i1].set_val(this_sec.diff_load_coef[i1]);
                e0_ld[i1].set_val(0.0);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "thermalExp" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    exp_ld[dv_comp].add(& dv_val);
                }
                else if cat.s == "diffExp" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_ld[dv_comp].add(& dv_val);
                }
                else if cat.s == "initialStrain" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    e0_ld[dv_comp].add(& dv_val);
                }
            }
        }
        else {
            for i1 in 0..3 {
                mod_dv[i1].set_val(this_mat.modulus[i1]);
                shr_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
                te_coef_dv[i1].set_val(this_mat.expansion[i1]);
                te_coef_dv[i1 + 3].set_val(this_mat.expansion[i1]);
                de_coef_dv[i1].set_val(this_mat.diff_exp[i1]);
                de_coef_dv[i1 + 3].set_val(this_mat.diff_exp[i1]);
                e0_dv[i1].set_val(0.0);
                e0_dv[i1 + 3].set_val(0.0);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            cat_list = CppStr::from("modulus shearModulus thermalExp diffExp initialStrain area areaMoment");
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                i2 = cat_list.find(&cat.s.as_str());
                if i2 < MAX_INT {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    if cat.s == "modulus" {
                        mod_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "shearModulus" {
                        shr_mod_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "thermalExp" {
                        te_coef_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "diffExp" {
                        de_coef_dv[dv_comp].add( & dv_val);
                    }
                    else if cat.s == "initialStrain" {
                        e0_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "area" {
                        area_dv.add(& dv_val);
                    }
                    else if cat.s == "areaMoment" {
                        idv[dv_comp].add(& dv_val);
                    }
                }
            }
            for i1 in 1..8 {
                qmat[i1].set_val(0.0);
            }
            qmat[0].set_val_dfd0(& mod_dv[0]);
            qmat[4].set_val_dfd0(& shr_mod_dv[0]);
            qmat[8].set_val_dfd0(& shr_mod_dv[1]);
            
            tmp.set_val_dfd0(& te_coef_dv[3]);
            te_coef_dv[1].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& te_coef_dv[4]);
            te_coef_dv[2].set_val_dfd0(& tmp);
            mat_mul_ar_dfd0(&mut qte, &mut  qmat, &mut  te_coef_dv,  3,  3,  1);

            tmp.set_val_dfd0(& de_coef_dv[3]);
            de_coef_dv[1].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& de_coef_dv[4]);
            de_coef_dv[2].set_val_dfd0(& tmp);
            mat_mul_ar_dfd0(&mut qde, &mut  qmat, &mut  de_coef_dv,  3,  3,  1);
            
            tmp.set_val_dfd0(& e0_dv[3]);
            e0_dv[1].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& e0_dv[4]);
            e0_dv[2].set_val_dfd0(& tmp);
            mat_mul_ar_dfd0(&mut qe0, &mut  qmat, &mut  e0_dv,  3,  3,  1);
            
            for i1 in 0..18 {
                dedgu[i1].set_val(0.0);
            }
            dedgu[0].set_val_dfd0(& area_dv);
            dedgu[4].set_val_dfd0(& area_dv);
            dedgu[8].set_val_dfd0(& area_dv);
            dedgu[10].set_val_dfd0(& idv[0]);
            dedgu[10].neg();
            dedgu[11].set_val_dfd0(& idv[1]);
            dedgu[12].set_val_dfd0(& idv[0]);
            dedgu[15].set_val_dfd0(& idv[1]);
            dedgu[15].neg();
            
            mat_mul_ar_dfd0(&mut tmp_exp, &mut  dedgu, &mut  qte,  6,  3,  1);
            ar_to_vec_dfd0(&mut tmp_exp, exp_ld,  0,  6);

            mat_mul_ar_dfd0(&mut tmp_exp, &mut  dedgu, &mut  qde,  6,  3,  1);
            ar_to_vec_dfd0(&mut tmp_exp, diff_ld,  0,  6);
            
            mat_mul_ar_dfd0(&mut tmp_exp, &mut  dedgu, &mut  qe0,  6,  3,  1);
            ar_to_vec_dfd0(&mut tmp_exp,  e0_ld,  0,  6);
        }
        
        return;
    }

    pub fn get_density_dfd0(&self, den : &mut DiffDoub0, layer : usize, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut dv_lay : usize;
        let mut coef = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let this_mat : &Material;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        if self.this_type == 3 || self.this_type == 41 {
            this_mat = &mat_ar[this_sec.get_layer_mat_ptr(layer)];
            den.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                dv_lay = this_dv.layer;
                if cat.s == "density" && dv_lay == layer {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    den.add(& dv_val);
                }
            }
        } else {
            this_mat = &mat_ar[this_sec.mat_ptr];
            den.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "density" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    den.add(& dv_val);
                }
            }
        }
        
        return;
    }

    pub fn get_shell_mass_dfd0(&self, mmat : &mut Vec<DiffDoub0>, lay_thk : &mut Vec<DiffDoub0>, lay_z : &mut Vec<DiffDoub0>, lay_den : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {
        let i1 : usize;
        let mut tmp = DiffDoub0::new();
        let mut tmp2 = DiffDoub0::new();
        let mut z_min = DiffDoub0::new();
        let mut z_min2 = DiffDoub0::new();
        let mut z_max = DiffDoub0::new();
        let mut z_max2 = DiffDoub0::new();
        
        for i1 in 0..36 {
            mmat[i1].set_val(0.0);
        }
        
        i1 = sec_ar[self.sect_ptr].layers.len();
        for layi in 0..i1 {
            tmp.set_val_dfd0(& lay_den[layi]);
            tmp.mult(& lay_thk[layi]);
            mmat[0].add(& tmp);
            mmat[7].add(& tmp);
            mmat[14].add(& tmp);
            
            tmp.set_val(0.5);
            tmp.mult(& lay_thk[layi]);
            z_min.set_val_dfd0(& lay_z[layi]);
            z_min.sub(& tmp);
            z_min2.set_val_dfd0(& z_min);
            z_min2.sqr();
            z_max.set_val_dfd0(& lay_z[layi]);
            z_max.add(& tmp);
            z_max2.set_val_dfd0(& z_max);
            z_max2.sqr();
            tmp.set_val_dfd0(& z_max2);
            tmp.sub(& z_min2);// tmp == z_max^2 - z_min^2
            tmp2.set_val(0.5);
            tmp2.mult(& lay_den[layi]);
            tmp2.mult(& tmp);// tmp2 = 0.5*rho*(z_max^2 - z_min^2)
            mmat[4].add(& tmp2);
            mmat[24].add(& tmp2);
            mmat[9].sub(& tmp2);
            mmat[19].sub(& tmp2);
            
            z_max2.mult(& z_max);
            z_min2.mult(& z_min);
            tmp.set_val_dfd0(& z_max2);
            tmp.sub(& z_min2);// tmp == z_max^3 - z_min^3
            tmp2.set_val(R_1O3);
            tmp2.mult(& lay_den[layi]);
            tmp2.mult(& tmp);
            mmat[21].add(& tmp2);
            mmat[28].add(& tmp2);
            mmat[35].add(& tmp2);
        }
        
        return;
    }

    pub fn get_beam_mass_dfd0(&self, mmat : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut dv_comp : usize;
        
        
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut d_cat : CppStr;
        
        
        let mut den_dv = DiffDoub0::new();
        let mut area_dv = DiffDoub0::new();
        let mut idv = [DiffDoub0::new(); 5];
        let mut jdv = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.mass[0] > 0.0 {
            for i1 in 0..36 {
                mmat[i1].set_val(this_sec.mass[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "massMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    mmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i3 = 6 * i1;// lower tri term
                i4 = i1;// upper tri term
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& mmat[i4]);
                    mmat[i3].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        else {
            for i1 in 0..36 {
                mmat[i1].set_val(0.0);
            }
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            area_dv.set_val(this_sec.area);
            jdv.set_val(this_sec.polar_moment);
            den_dv.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                d_cat = this_dv.category.clone();
                if d_cat.s == "density" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    den_dv.add(& dv_val);
                }
                else if d_cat.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if d_cat.s == "areaMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    dv_comp = this_dv.component;
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if d_cat.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& area_dv);
            mmat[0].set_val_dfd0(& tmp);
            mmat[7].set_val_dfd0(& tmp);
            mmat[14].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& idv[0]);
            mmat[4].set_val_dfd0(& tmp);
            mmat[9].set_val_dfd0(& tmp);
            mmat[9].neg();
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& idv[1]);
            mmat[5].set_val_dfd0(& tmp);
            mmat[5].neg();
            mmat[15].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& jdv);
            mmat[21].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& idv[2]);
            mmat[28].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& idv[4]);
            mmat[29].set_val_dfd0(& tmp);
            tmp.set_val_dfd0(& den_dv);
            tmp.mult(& idv[3]);
            mmat[35].set_val_dfd0(& tmp);
            
            for i1 in 3..6 {
                i3 = 6 * i1;// lower tri term
                i4 = i1;// upper tri term
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& mmat[i4]);
                    mmat[i3].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        return;
    }

    pub fn get_solid_damp_dfd0(&self, dmat : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        
        let mut temp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..36 {
            dmat[i1].set_val(this_mat.damping[i1]);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "dampingMat" {
                this_dv.get_value_dfd0(&mut dv_val);
                dv_comp = this_dv.component;
                temp.set_val(dv.doub_dat);
                dv_val.mult(& temp);
                dmat[dv_comp].add(& dv_val);
            }
        }
        
        for i1 in 1..6 {
            i3 = 6 * i1;// lower tri term
            i4 = i1;// upper tri term
            for _i2 in 0..i1 {
                temp.set_val_dfd0(& dmat[i4]);
                dmat[i3].set_val_dfd0(& temp);
                i3 += 1usize;
                i4  +=  6;
            }
        }
        
        return;
    }

    pub fn get_shell_damp_dfd0(&self, dmat : &mut Vec<DiffDoub0>, lay_thk : &mut Vec<DiffDoub0>, lay_z : &mut Vec<DiffDoub0>, lay_d : &mut Vec<DiffDoub0>, lay_ang : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {
        self.get_abd_dfd0(dmat, lay_thk, lay_z, lay_d, lay_ang, sec_ar);
        dmat[60].set_val(0.0);
        dmat[70].set_val(0.0);
        dmat[80].set_val(0.0);
        return;
    }

    pub fn get_beam_damp_dfd0(&self, dmat : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_comp : usize;
        let mut dmat_dv = [DiffDoub0::new(); 36];
        let mut area_dv = DiffDoub0::new();
        let mut idv = [DiffDoub0::new(); 5];
        let mut jdv = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub0::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.damping[0] > 0.0 {
            for i1 in 0..36 {
                dmat[i1].set_val(this_sec.damping[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "dampingMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    dmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& dmat[i5]);
                    dmat[i4].set_val_dfd0(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..36 {
                dmat_dv[i1].set_val(this_mat.damping[i1]);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            jdv.set_val(this_sec.polar_moment);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "dampingMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    dmat_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if this_dv.category.s == "areaMoment" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            
            for i1 in 0..6 {
                i3 = 7 * i1;
                for _i2 in i1..6 {
                    dmat[i3].set_val(0.0);
                    i3 += 1usize;
                }
            }
            
            dmat[0].set_val_dfd0(& dmat_dv[0]);
            dmat[0].mult(& area_dv);
            dmat[4].set_val_dfd0(& dmat_dv[0]);
            dmat[4].mult(& idv[0]);
            dmat[5].set_val_dfd0(& dmat_dv[0]);
            dmat[5].neg();
            dmat[5].mult(& idv[1]);
            dmat[7].set_val_dfd0(& dmat_dv[21]);
            dmat[7].mult(& area_dv);
            dmat[9].set_val_dfd0(& dmat_dv[21]);
            dmat[9].neg();
            dmat[9].mult(& idv[0]);
            dmat[14].set_val_dfd0(& dmat_dv[28]);// 28
            dmat[14].mult(& area_dv);
            dmat[15].set_val_dfd0(& dmat_dv[35]);// 35
            dmat[15].mult(& idv[1]);
            dmat[21].set_val_dfd0(& dmat_dv[21]);// 21
            dmat[21].mult(& jdv);
            dmat[28].set_val_dfd0(& dmat_dv[0]);
            dmat[28].mult(& idv[2]);
            dmat[29].set_val_dfd0(& dmat_dv[0]);
            dmat[29].neg();
            dmat[29].mult(& idv[4]);
            dmat[35].set_val_dfd0(& dmat_dv[0]);
            dmat[35].mult(& idv[3]);
            
            for i1 in 1..6 {
                i3 = 6 * i1;
                i4 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd0(& dmat[i4]);
                    dmat[i3].set_val_dfd0(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        
        return;
    }

    pub fn get_conductivity_dfd0(&self, t_cond : &mut Vec<DiffDoub0>, diff : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cond_dv = [DiffDoub0::new(); 6];
        let mut diff_dv = [DiffDoub0::new(); 6];
        
        let mut temp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..6 {
            cond_dv[i1].set_val(this_mat.conductivity[i1]);
            diff_dv[i1].set_val(this_mat.diffusivity[i1]);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "thermalCond" {
                dv_comp = this_dv.component - 1;
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& temp);
                cond_dv[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "diffusivity" {
                dv_comp = this_dv.component - 1;
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& temp);
                diff_dv[dv_comp].add(& dv_val);
            }
        }
        
        t_cond[0] = cond_dv[0];
        t_cond[1] = cond_dv[3];
        t_cond[2] = cond_dv[4];
        t_cond[3] = cond_dv[3];
        t_cond[4] = cond_dv[1];
        t_cond[5] = cond_dv[5];
        t_cond[6] = cond_dv[4];
        t_cond[7] = cond_dv[5];
        t_cond[8] = cond_dv[2];

        diff[0] = diff_dv[0];
        diff[1] = diff_dv[3];
        diff[2] = diff_dv[4];
        diff[3] = diff_dv[3];
        diff[4] = diff_dv[1];
        diff[5] = diff_dv[5];
        diff[6] = diff_dv[4];
        diff[7] = diff_dv[5];
        diff[8] = diff_dv[2];
        
        return;
    }

    pub fn get_shell_cond_dfd0(&self, t_cond : &mut Vec<DiffDoub0>, diff : &mut Vec<DiffDoub0>, lay_thk : &mut Vec<DiffDoub0>, lay_ang : &mut Vec<DiffDoub0>, 
        lay_cond : &mut Vec<DiffDoub0>, lay_diff : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {

        let mut i2 : usize;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let mut layer_mat = [DiffDoub0::new(); 9];
        let mut al_mat = [DiffDoub0::new(); 9];
        let mut al_t = [DiffDoub0::new(); 9];
        let mut tmp = [DiffDoub0::new(); 9];
        let mut this_cond = [DiffDoub0::new(); 9];
        let mut this_diff = [DiffDoub0::new(); 9];
        let mut tmp2 = DiffDoub0::new();
        
        for i1 in 0..9 {
            t_cond[i1].set_val(0.0);
            diff[i1].set_val(0.0);
            al_mat[i1].set_val(0.0);
        }
        al_mat[8].set_val(1.0);
        
        for i1 in 0..num_lay {
            al_mat[0].set_val(R_PI_0180);
            al_mat[0].mult(& lay_ang[i1]);
            al_mat[0].cs();
            tmp2.set_val_dfd0(& al_mat[0]);
            al_mat[4].set_val_dfd0(& tmp2);
            al_mat[3].set_val(R_PI_0180);
            al_mat[3].mult(& lay_ang[i1]);
            al_mat[3].sn();
            tmp2.set_val_dfd0(& al_mat[3]);
            al_mat[1].set_val_dfd0(& tmp2);
            al_mat[1].neg();
            transpose_ar_dfd0(&mut al_t, &mut  al_mat,  3,  3);
            i2 = 9 * i1;
            
            vec_to_ar_dfd0(&mut this_cond, lay_cond,  i2,  i2 + 9);
            mat_mul_ar_dfd0(&mut tmp, &mut  this_cond, &mut  al_t,  3,  3,  3);
            mat_mul_ar_dfd0(&mut layer_mat, &mut  al_mat, &mut  tmp,  3,  3,  3);
            for i2 in 0..9 {
                layer_mat[i2].mult(& lay_thk[i1]);
                t_cond[i2].add(& layer_mat[i2]);
            }

            vec_to_ar_dfd0(&mut this_diff, lay_diff, i2,  i2 + 9);
            mat_mul_ar_dfd0(&mut tmp, &mut  this_diff, &mut  al_t,  3,  3,  3);
            mat_mul_ar_dfd0(&mut layer_mat, &mut  al_mat, &mut  tmp,  3,  3,  3);
            for i2 in 0..9 {
                layer_mat[i2].mult(& lay_thk[i1]);
                diff[i2].add(& layer_mat[i2]);
            }
        }
        
        return;
    }

    pub fn get_beam_cond_dfd0(&self, t_cond : &mut Vec<DiffDoub0>, diff : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut cond_dv = DiffDoub0::new();
        let mut diff_dv = DiffDoub0::new();
        let mut area_dv = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];

        if this_sec.conductivity > 0.0 || this_sec.diffusivity > 0.0 {
            cond_dv.set_val(this_sec.conductivity);
            diff_dv.set_val(this_sec.diffusivity);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermCond" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv.add(& dv_val);
                }
                else if this_dv.category.s == "diffusivity" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv.add(& dv_val);
                }
            }
        }
        else {
            cond_dv.set_val(this_mat.conductivity[0]);
            diff_dv.set_val(this_mat.diffusivity[0]);
            area_dv.set_val(this_sec.area);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "thermCond" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv.add(& dv_val);
                }
                else if cat.s == "diffusivity" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv.add(& dv_val);
                }
                else if cat.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    area_dv.add(& dv_val);
                }
            }
            cond_dv.mult(& area_dv);
            diff_dv.mult(& area_dv);
        }
        
        for i1 in 1..8 {
            t_cond[i1].set_val(0.0);
            diff[i1].set_val(0.0);
        }
        
        t_cond[0].set_val_dfd0(& cond_dv);
        t_cond[4].set_val_dfd0(& cond_dv);
        t_cond[8].set_val_dfd0(& cond_dv);

        diff[0].set_val_dfd0(& diff_dv);
        diff[4].set_val_dfd0(& diff_dv);
        diff[8].set_val_dfd0(& diff_dv);
        
        return;
    }

    pub fn get_specific_heat_dfd0(&self, spec_heat : &mut DiffDoub0, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut temp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        spec_heat.set_val(this_mat.spec_heat);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "specHeat" {
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& temp);
                spec_heat.add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_shell_spec_heat_dfd0(&self, spec_heat : &mut DiffDoub0, lay_thk : &mut Vec<DiffDoub0>, lay_sh : &mut Vec<DiffDoub0>, lay_den : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>) {
        let num_lay : usize;
        let mut tmp = DiffDoub0::new();
        
        spec_heat.set_val(0.0);
        num_lay = sec_ar[self.sect_ptr].layers.len();
        for i1 in 0..num_lay {
            tmp.set_val_dfd0(& lay_den[i1]);
            tmp.mult(& lay_sh[i1]);
            tmp.mult(& lay_thk[i1]);
            spec_heat.add(& tmp);
        }
        
        return;
    }

    pub fn get_beam_spec_heat_dfd0(&self, spec_heat : &mut DiffDoub0, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut density_dv = DiffDoub0::new();
        let mut area_dv = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut dv_val = DiffDoub0::new();
        
        let mut this_dv : &DesignVariable;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        let sec_sh : f64 =  this_sec.spec_heat;
        let mat_sh : f64 =  this_mat.spec_heat;
        let mat_den : f64 =  this_mat.density;
        
        if sec_sh > 0.0 {
            spec_heat.set_val(sec_sh);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "specHeat" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    spec_heat.add(& dv_val);
                }
            }
        }
        else {
            spec_heat.set_val(mat_sh);
            density_dv.set_val(mat_den);
            area_dv.set_val(this_sec.area);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "specHeatDV" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    spec_heat.add(& dv_val);
                }
                else if cat.s == "density" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    density_dv.add(& dv_val);
                }
                else if cat.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd0(&mut dv_val);
                    dv_val.mult(& tmp);
                    area_dv.add(& dv_val);
                }
            }
            spec_heat.mult(& density_dv);
            spec_heat.mult(& area_dv);
        }
        
        return;
    }

    pub fn get_nd_crds_dfd0(&self, x_glob : &mut Vec<DiffDoub0>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut nd_crd = [DiffDoub0::new(); 3];
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            this_nd.get_crd_dfd0(&mut nd_crd, dv_ar);
            x_glob[i1].set_val_dfd0(& nd_crd[0]);
            x_glob[i1+self.num_nds].set_val_dfd0(& nd_crd[1]);
            x_glob[i1+2*self.num_nds].set_val_dfd0(& nd_crd[2]);
        }
        
        return;
    }

    pub fn get_loc_ori_dfd0(&self, loc_ori : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i1 : usize;
        let mut dv_cat : CppStr;
        let mut rot = [DiffDoub0::new(); 3];
        let mut dv_val = DiffDoub0::new();
        let mut coef = DiffDoub0::new();
        let mut ori_copy = [DiffDoub0::new(); 9];
        let mut tmp_ori = [DiffDoub0::new(); 9];
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for i1 in 0..9 {
            ori_copy[i1].set_val(this_sec.orientation[i1]);
        }
        
        rot[0].set_val(0.0);
        rot[1].set_val(0.0);
        rot[2].set_val(0.0);
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            dv_cat = this_dv.category.clone();
            if dv_cat.s == "orientation" {
                i1 = this_dv.component - 1;
                coef.set_val(dv.doub_dat);
                this_dv.get_value_dfd0(&mut dv_val);
                dv_val.mult(& coef);
                rot[i1].add(& dv_val);
            }
        }
        
        rotate_orient_dfd0(&mut tmp_ori, &mut  ori_copy, &mut  rot);
        ar_to_vec_dfd0(&mut tmp_ori, loc_ori,  0,  9);
        
        return;
    }

    pub fn correct_orient_dfd0(&self, loc_ori : &mut Vec<DiffDoub0>, x_glob : &mut Vec<DiffDoub0>) {
        let mut v1 = [DiffDoub0::new(); 3];
        let mut v2 = [DiffDoub0::new(); 3];
        let mut v3 = [DiffDoub0::new(); 3];
        let mut rot = [DiffDoub0::new(); 3];
        let mut ori_copy = [DiffDoub0::new(); 9];
        let mut tmp_ori = [DiffDoub0::new(); 9];
        
        let mut dp = DiffDoub0::new();
        let mut magv3 = DiffDoub0::new();
        let mut mag_cp = DiffDoub0::new();
        let mut theta = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        
        for i1 in 0..9 {
            ori_copy[i1].set_val_dfd0(& loc_ori[i1]);
        }
        
        if self.this_type == 3 || self.this_type == 41 {
            if self.this_type == 3 {
                v1[0].set_val_dfd0(& x_glob[1]);
                v1[0].sub(& x_glob[0]);
                v1[1].set_val_dfd0(& x_glob[4]);
                v1[1].sub(& x_glob[3]);
                v1[2].set_val_dfd0(& x_glob[7]);
                v1[2].sub(& x_glob[6]);
                
                v2[0].set_val_dfd0(& x_glob[2]);
                v2[0].sub(& x_glob[0]);
                v2[1].set_val_dfd0(& x_glob[5]);
                v2[1].sub(& x_glob[3]);
                v2[2].set_val_dfd0(& x_glob[8]);
                v2[2].sub(& x_glob[6]);
            } else {
                v1[0].set_val_dfd0(& x_glob[2]);
                v1[0].sub(& x_glob[0]);
                v1[1].set_val_dfd0(& x_glob[6]);
                v1[1].sub(& x_glob[4]);
                v1[2].set_val_dfd0(& x_glob[10]);
                v1[2].sub(& x_glob[8]);
                
                v2[0].set_val_dfd0(& x_glob[3]);
                v2[0].sub(& x_glob[1]);
                v2[1].set_val_dfd0(& x_glob[7]);
                v2[1].sub(& x_glob[5]);
                v2[2].set_val_dfd0(& x_glob[11]);
                v2[2].sub(& x_glob[9]);
            }
            cross_prod_dfd0(&mut v3, &mut  v1, &mut  v2);
            
            dp.set_val_dfd0(& v3[0]);
            dp.mult(& loc_ori[6]);
            tmp.set_val_dfd0(& v3[1]);
            tmp.mult(& loc_ori[7]);
            dp.add(& tmp);
            tmp.set_val_dfd0(& v3[2]);
            tmp.mult(& loc_ori[8]);
            dp.add(& tmp);
            
            if dp.val < 0.0 {
                v3[0].neg();
                v3[1].neg();
                v3[2].neg();
            }
            
            magv3.set_val_dfd0(& v3[0]);
            magv3.sqr();
            tmp.set_val_dfd0(& v3[1]);
            tmp.sqr();
            magv3.add(& tmp);
            tmp.set_val_dfd0(& v3[2]);
            tmp.sqr();
            magv3.add(& tmp);
            magv3.sqt();
            
            cross_prod_dfd0(&mut rot, &mut loc_ori[6..], &mut v3);
            mag_cp.set_val_dfd0(& rot[0]);
            mag_cp.sqr();
            tmp.set_val_dfd0(& rot[1]);
            tmp.sqr();
            mag_cp.add(& tmp);
            tmp.set_val_dfd0(& rot[2]);
            tmp.sqr();
            mag_cp.add(& tmp);
            mag_cp.sqt();
            if mag_cp.val < 1e-12 {
                return;
            }
            
            theta.set_val_dfd0(& mag_cp);
            theta.dvd(& magv3);
            theta.asn();
            
            tmp.set_val_dfd0(& theta);
            tmp.dvd(& mag_cp);
            
            rot[0].mult(& tmp);
            rot[1].mult(& tmp);
            rot[2].mult(& tmp);
            
            rotate_orient_dfd0(&mut tmp_ori, &mut  ori_copy, &mut  rot);
            ar_to_vec_dfd0(&mut tmp_ori, loc_ori,  0,  9);
        }
        
        return;
    }

    pub fn get_frc_fld_const_dfd0(&self, coef : &mut Vec<DiffDoub0>, exp : &mut Vec<DiffDoub0>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        coef[0].set_val(this_sec.pot_coef);
        coef[1].set_val(this_sec.damp_coef);
        exp[0].set_val(this_sec.pot_exp);
        exp[1].set_val(this_sec.damp_exp);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "potFldCoef" {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[0].add(& dv_val);
            }
            else if cat.s == "dampFldCoef" {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[1].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_thrm_fld_const_dfd0(&self, coef : &mut Vec<DiffDoub0>, ref_t : &mut DiffDoub0, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        coef[0].set_val(this_sec.cond_coef);
        coef[1].set_val(this_sec.rad_coef);
        ref_t.set_val(this_sec.ref_temp);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "condCoef" {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[0].add(& dv_val);
            }
            else if cat.s == "radCoef" {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[1].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_mass_per_el_dfd0(&self, mass_per_el : &mut DiffDoub0, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec = &sec_ar[self.sect_ptr];
        mass_per_el.set_val(this_sec.mass_per_el);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "massPerEl" {
                this_dv.get_value_dfd0(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                mass_per_el.add(& dv_val);
            }
        }
        
        return;
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1
    pub fn get_gen_prop_dfd1(&self, prop : &mut DiffDoub1, prop_key : &mut CppStr, dv_ar : & Vec<DesignVariable>) {
        let mut d_vind : usize;
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut this_dv : & DesignVariable;
        
        for dv in self.design_vars.iter() {
            d_vind = dv.int_dat;
            this_dv = & dv_ar[d_vind];
            if this_dv.category.s == prop_key.s {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                prop.add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_f_per_mass_dfd1(&self, f_vec : &mut Vec<DiffDoub1>, dv_ar : &Vec<DesignVariable>) {
        let mut d_vind : usize;
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut this_dv : & DesignVariable;
        let mut comp : usize;
        
        for dv in self.design_vars.iter() {
            d_vind = dv.int_dat;
            this_dv = & dv_ar[d_vind];
            if this_dv.category.s == "bodyForce" {
                this_dv.get_value_dfd1(&mut dv_val);
                comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                //prop.add(& dv_val);
                f_vec[comp].add(&dv_val);
            }
        }
    }

    pub fn get_layer_thk_z_dfd1(&self, lay_thk : &mut Vec<DiffDoub1>, lay_z : &mut Vec<DiffDoub1>, z_offset : &mut DiffDoub1, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        //z_offset = 1: upper z surface is reference plane
        //z_offset = -1: lower z surface is reference plane
        let mut layi : usize =  0;
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut tot_thk = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut z_crd = DiffDoub1::new();
        let mut z_next = DiffDoub1::new();
        let mut z_mid = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        tot_thk.set_val(0.0);
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            lay_thk[layi].set_val(lay.thickness);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thickness" && this_dv.layer == layi {
                    this_dv.get_value_dfd1(&mut dv_val);
                    coef.set_val(dv.doub_dat);
                    dv_val.mult(& coef);
                    lay_thk[layi].add(& dv_val);
                }
            }
            tot_thk.add(& lay_thk[layi]);
            layi += 1usize;
        }
        
        z_offset.set_val(this_sec.z_offset);
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "zOffset" {
                this_dv.get_value_dfd1(&mut dv_val);
                coef.set_val(dv.doub_dat);
                dv_val.mult(& coef);
                z_offset.add(& dv_val);
            }
        }
        
        tmp.set_val(1.0);
        tmp.add(z_offset);
        z_crd.set_val(-0.5);
        z_crd.mult(& tot_thk);
        z_crd.mult(& tmp);
        
        layi = 0;
        for _lay in this_sec.layers.iter() {
            z_next.set_val_dfd1(& z_crd);
            z_next.add(& lay_thk[layi]);
            tmp.set_val_dfd1(& z_crd);
            tmp.add(& z_next);
            z_mid.set_val(0.5);
            z_mid.mult(& tmp);
            lay_z[layi].set_val_dfd1(& z_mid);
            layi += 1usize;
            z_crd.set_val_dfd1(& z_next);
        }
        
        return;
    }

    pub fn get_layer_q_dfd1(&self, lay_q : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub1::new(); 3];
        let mut poisson_dv = [DiffDoub1::new(); 3];
        let mut shear_mod_dv = [DiffDoub1::new(); 3];
        let mut smat = [DiffDoub1::new(); 9];
        let mut qmat = [DiffDoub1::new(); 9];
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut x_vec = [DiffDoub1::new(); 3];
        let mut b_vec = [DiffDoub1::new(); 3];
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        let mut tmp = DiffDoub1::new();
        
        i2 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "modulus" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd1(&mut dv_val);
                        dv_val.mult(& coef);
                        modulus_dv[dv_comp - 1].add(& dv_val);
                    }
                    else if dv_cat.s == "poissonRatio" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd1(&mut dv_val);
                        dv_val.mult(& coef);
                        poisson_dv[dv_comp - 1].add(& dv_val);
                    }
                    else if dv_cat.s == "shearModulus" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd1(&mut dv_val);
                        dv_val.mult(& coef);
                        shear_mod_dv[dv_comp - 1].add(& dv_val);
                    }
                }
            }
            
            for i1 in 0..9 {
                smat[i1].set_val(0.0);
            }
            smat[0].set_val(1.0);
            smat[0].dvd(& modulus_dv[0]);
            smat[1].set_val_dfd1(& poisson_dv[0]);
            smat[1].neg();
            smat[1].dvd(& modulus_dv[0]);
            tmp.set_val_dfd1(& smat[1]);
            smat[3].set_val_dfd1(& tmp);
            smat[4].set_val(1.0);
            smat[4].dvd(& modulus_dv[1]);
            smat[8].set_val(1.0);
            smat[8].dvd(& shear_mod_dv[0]);
            get_det_inv_ar_dfd1(&mut coef, &mut  qmat, &mut  smat,  3,  0, &mut  x_vec, &mut  b_vec);
            
            for i1 in 0..9 {
                lay_q[i2].set_val_dfd1(& qmat[i1]);
                i2 += 1usize;
            }
            
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_d_dfd1(&self, lay_d : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut damp_mat_dv = [DiffDoub1::new(); 36];
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        let mut tmp = DiffDoub1::new();
        
        i2 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..36 {
                damp_mat_dv[i1].set_val(this_mat.damping[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "dampingMat" {
                        dv_comp = this_dv.component;
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd1(&mut dv_val);
                        dv_val.mult(& coef);
                        damp_mat_dv[dv_comp].add(& dv_val);
                    }
                }
            }
            
            for i3 in 1..6 {
                i5 = 6 * i3;// lower tri term
                i6 = i3;// upper tri term
                for _i4 in 0..i3 {
                    tmp.set_val_dfd1(& damp_mat_dv[i6]);
                    damp_mat_dv[i5].set_val_dfd1(& tmp);
                    i5 += 1usize;
                    i6  +=  6;
                }
            }
            
            lay_d[i2].set_val_dfd1(& damp_mat_dv[0]);
            lay_d[i2 + 1].set_val_dfd1(& damp_mat_dv[1]);
            lay_d[i2 + 2].set_val_dfd1(& damp_mat_dv[3]);
            lay_d[i2 + 3].set_val_dfd1(& damp_mat_dv[6]);
            lay_d[i2 + 4].set_val_dfd1(& damp_mat_dv[7]);
            lay_d[i2 + 5].set_val_dfd1(& damp_mat_dv[9]);
            lay_d[i2 + 6].set_val_dfd1(& damp_mat_dv[18]);
            lay_d[i2 + 7].set_val_dfd1(& damp_mat_dv[19]);
            lay_d[i2 + 8].set_val_dfd1(& damp_mat_dv[21]);
            
            layi += 1usize;
            i2  +=  9;
        }
        return;
    }

    pub fn get_layer_angle_dfd1(&self, lay_ang : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut dv_ind : usize;
        let mut angle = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            angle.set_val(lay.angle);
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.layer == layi {
                    dv_cat = this_dv.category.clone();
                    if dv_cat.s == "angle" {
                        coef.set_val(dv.doub_dat);
                        this_dv.get_value_dfd1(&mut dv_val);
                        dv_val.mult(& coef);
                        angle.add(& dv_val);
                    }
                }
            }
            lay_ang[layi].set_val_dfd1(& angle);
            
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_th_exp_dfd1(&self, lay_th_exp : &mut Vec<DiffDoub1>, lay_diff_exp : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut t_exp_dv = [DiffDoub1::new(); 6];
        let mut d_exp_dv = [DiffDoub1::new(); 6];
        let mut dv_comp : usize;
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        
        layi = 0;
        i2 = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i1 in 0..6 {
                t_exp_dv[i1].set_val(this_mat.expansion[i1]);
                d_exp_dv[i1].set_val(this_mat.diff_exp[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermalExp" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    t_exp_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "diffExp" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    d_exp_dv[dv_comp].add(& dv_val);
                }
            }
            lay_th_exp[i2].set_val_dfd1(& t_exp_dv[0]);
            lay_diff_exp[i2].set_val_dfd1(& d_exp_dv[0]);
            i2 += 1usize;
            lay_th_exp[i2].set_val_dfd1(& t_exp_dv[1]);
            lay_diff_exp[i2].set_val_dfd1(& d_exp_dv[1]);
            i2 += 1usize;
            lay_th_exp[i2].set_val_dfd1(& t_exp_dv[3]);
            lay_diff_exp[i2].set_val_dfd1(& d_exp_dv[3]);
            i2 += 1usize;
            layi += 1usize;
        }
        return;
    }

    pub fn get_layer_einit_dfd1(&self, lay_einit : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i2 : usize;
        let mut layi : usize;
        let mut e0_dv = [DiffDoub1::new(); 6];
        let mut dv_comp : usize;
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        i2 = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for _lay in this_sec.layers.iter() {
            for i1 in 0..6 {
                e0_dv[i1].set_val(0.0);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "initialStrain" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    e0_dv[dv_comp].add(& dv_val);
                }
            }
            lay_einit[i2].set_val_dfd1(& e0_dv[0]);
            i2 += 1usize;
            lay_einit[i2].set_val_dfd1(& e0_dv[1]);
            i2 += 1usize;
            lay_einit[i2].set_val_dfd1(& e0_dv[3]);
            i2 += 1usize;
            layi += 1usize;
        }
        return;
    }

    pub fn get_layer_den_dfd1(&self, layer_den : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut mat_den : f64;
        let mut den_dv = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            mat_den = mat_ar[lay.mat_ptr].density;
            den_dv.set_val(mat_den);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "density" && this_dv.layer == layi {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    den_dv.add(& dv_val);
                }
            }
            layer_den[layi].set_val_dfd1(& den_dv);
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_cond_dfd1(&self, lay_cond : &mut Vec<DiffDoub1>, lay_diff : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i1 : usize;
        let mut layi : usize;
        let mut cond_dv = [DiffDoub1::new(); 6];
        let mut diff_dv = [DiffDoub1::new(); 6];
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        let mut this_mat : &Material;
        
        i1 = 0;
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            this_mat = &mat_ar[lay.mat_ptr];
            for i2 in 0..6 {
                cond_dv[i2].set_val(this_mat.conductivity[i2]);
                diff_dv[i2].set_val(this_mat.diffusivity[i2]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermalCond" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "diffusivity" && this_dv.layer == layi {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv[dv_comp].add(& dv_val);
                }
            }
            lay_cond[i1].set_val_dfd1(& cond_dv[0]);
            lay_cond[i1 + 1].set_val_dfd1(& cond_dv[3]);
            lay_cond[i1 + 2].set_val_dfd1(& cond_dv[4]);
            lay_cond[i1 + 3].set_val_dfd1(& cond_dv[3]);
            lay_cond[i1 + 4].set_val_dfd1(& cond_dv[1]);
            lay_cond[i1 + 5].set_val_dfd1(& cond_dv[5]);
            lay_cond[i1 + 6].set_val_dfd1(& cond_dv[4]);
            lay_cond[i1 + 7].set_val_dfd1(& cond_dv[5]);
            lay_cond[i1 + 8].set_val_dfd1(& cond_dv[2]);

            lay_diff[i1].set_val_dfd1(& diff_dv[0]);
            lay_diff[i1 + 1].set_val_dfd1(& diff_dv[3]);
            lay_diff[i1 + 2].set_val_dfd1(& diff_dv[4]);
            lay_diff[i1 + 3].set_val_dfd1(& diff_dv[3]);
            lay_diff[i1 + 4].set_val_dfd1(& diff_dv[1]);
            lay_diff[i1 + 5].set_val_dfd1(& diff_dv[5]);
            lay_diff[i1 + 6].set_val_dfd1(& diff_dv[4]);
            lay_diff[i1 + 7].set_val_dfd1(& diff_dv[5]);
            lay_diff[i1 + 8].set_val_dfd1(& diff_dv[2]);

            i1  +=  9;
            layi += 1usize;
        }
        
        return;
    }

    pub fn get_layer_spec_heat_dfd1(&self, lay_sh : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut layi : usize;
        let mut mat_sh : f64;
        let mut sh_dv = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        layi = 0;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for lay in this_sec.layers.iter() {
            mat_sh = mat_ar[lay.mat_ptr].spec_heat;
            sh_dv.set_val(mat_sh);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "specHeat" && this_dv.layer == layi {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    sh_dv.add(& dv_val);
                }
            }
            lay_sh[layi].set_val_dfd1(& sh_dv);
            layi += 1usize;
        }
        
        return;
    }

    pub fn transform_strain_dfd1(&self, stn_new : &mut [DiffDoub1], stn_orig : &mut [DiffDoub1], angle : &mut DiffDoub1) {
        let mut angle_rad = DiffDoub1::new();
        let mut a11 = DiffDoub1::new();
        let mut a12 = DiffDoub1::new();
        let mut a21 = DiffDoub1::new();
        let mut a22 = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut te = [DiffDoub1::new(); 9];
        
        
        angle_rad.set_val(R_PI_0180);
        angle_rad.mult(angle);
        a11.set_val_dfd1(& angle_rad);
        a11.cs();
        a21.set_val_dfd1(& angle_rad);
        a21.sn();
        a12.set_val_dfd1(& a21);
        a12.neg();
        a22.set_val_dfd1(& a11);
        
        te[0].set_val_dfd1(& a11);
        te[0].sqr();
        te[1].set_val_dfd1(& a12);
        te[1].sqr();
        te[2].set_val_dfd1(& a11);
        te[2].mult(& a12);
        te[3].set_val_dfd1(& a21);
        te[3].sqr();
        te[4].set_val_dfd1(& a22);
        te[4].sqr();
        te[5].set_val_dfd1(& a22);
        te[5].mult(& a21);
        te[6].set_val(2.0);
        te[6].mult(& a11);
        te[6].mult(& a21);
        te[7].set_val(2.0);
        te[7].mult(& a12);
        te[7].mult(& a22);
        te[8].set_val_dfd1(& a11);
        te[8].mult(& a22);
        tmp.set_val_dfd1(& a12);
        tmp.mult(& a21);
        te[8].add(& tmp);
        
        mat_mul_ar_dfd1(stn_new, &mut te, stn_orig,  3,  3,  1);
        
        return;
    }

    pub fn transform_q_dfd1(&self, q_new : &mut [DiffDoub1], q_orig : &mut [DiffDoub1], angle : &mut DiffDoub1) {
        let mut angle_rad = DiffDoub1::new();
        let mut a11 = DiffDoub1::new();
        let mut a12 = DiffDoub1::new();
        let mut a21 = DiffDoub1::new();
        let mut a22 = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut ts = [DiffDoub1::new(); 9];
        let mut te = [DiffDoub1::new(); 9];
        let mut te_inv = [DiffDoub1::new(); 9];
        let mut x_vec = [DiffDoub1::new(); 3];
        let mut b_vec = [DiffDoub1::new(); 3];
        
        
        angle_rad.set_val(R_PI_0180);
        angle_rad.mult(angle);
        a11.set_val_dfd1(& angle_rad);
        a11.cs();
        a21.set_val_dfd1(& angle_rad);
        a21.sn();
        a12.set_val_dfd1(& a21);
        a12.neg();
        a22.set_val_dfd1(& a11);
        
        ts[0].set_val_dfd1(& a11);
        ts[0].sqr();
        ts[1].set_val_dfd1(& a12);
        ts[1].sqr();
        ts[2].set_val(2.0);
        ts[2].mult(& a11);
        ts[2].mult(& a12);
        ts[3].set_val_dfd1(& a21);
        ts[3].sqr();
        ts[4].set_val_dfd1(& a22);
        ts[4].sqr();
        ts[5].set_val(2.0);
        ts[5].mult(& a22);
        ts[5].mult(& a21);
        ts[6].set_val_dfd1(& a11);
        ts[6].mult(& a21);
        ts[7].set_val_dfd1(& a12);
        ts[7].mult(& a22);
        ts[8].set_val_dfd1(& a11);
        ts[8].mult(& a22);
        tmp.set_val_dfd1(& a12);
        tmp.mult(& a21);
        ts[8].add(& tmp);
        
        for i1 in 0..9 {
            te[i1].set_val_dfd1(& ts[i1]);
        }
        tmp.set_val(0.5);
        te[2].mult(& tmp);
        te[5].mult(& tmp);
        tmp.set_val(2.0);
        te[6].mult(& tmp);
        te[7].mult(& tmp);
        
        get_det_inv_ar_dfd1(&mut coef, &mut  te_inv, &mut  te,  3,  0, &mut  x_vec, &mut  b_vec);
        
        mat_mul_ar_dfd1(&mut te, q_orig, &mut te_inv,  3,  3,  3);
        mat_mul_ar_dfd1(q_new, &mut  ts, &mut  te,  3,  3,  3);
        
        return;
    }

    pub fn get_solid_stiff_dfd1(&self, cmat : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub1::new(); 3];
        let mut poisson_dv = [DiffDoub1::new(); 3];
        let mut shear_mod_dv = [DiffDoub1::new(); 3];
        let mut smat = [DiffDoub1::new(); 36];
        let mut ctmp = [DiffDoub1::new(); 36];
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut x_vec = [DiffDoub1::new(); 6];
        let mut b_vec = [DiffDoub1::new(); 6];
        let mut dv_cat : CppStr;
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub1::new();
        
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_mat.stiffness[0] > 0.0 {
            for i1 in 0..36 {
                cmat[i1].set_val(this_mat.stiffness[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "stiffnessMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    cmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& cmat[i5]);
                    cmat[i4].set_val_dfd1(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                dv_cat = this_dv.category.clone();
                if dv_cat.s == "modulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    modulus_dv[dv_comp - 1].add(& dv_val);
                }
                else if dv_cat.s == "poissonRatio" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    poisson_dv[dv_comp - 1].add(& dv_val);
                }
                else if dv_cat.s == "shearModulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    shear_mod_dv[dv_comp - 1].add(& dv_val);
                }
            }
            for i1 in 0..36 {
                smat[i1].set_val(0.0);
            }
            
            smat[0].set_val(1.0);
            smat[0].dvd(& modulus_dv[0]);
            smat[1].set_val_dfd1(& poisson_dv[0]);
            smat[1].neg();
            smat[1].dvd(& modulus_dv[0]);
            smat[2].set_val_dfd1(& poisson_dv[1]);
            smat[2].neg();
            smat[2].dvd(& modulus_dv[0]);
            tmp.set_val_dfd1(& smat[1]);
            smat[6].set_val_dfd1(& tmp);
            smat[7].set_val(1.0);
            smat[7].dvd(& modulus_dv[1]);
            smat[8].set_val_dfd1(& poisson_dv[2]);
            smat[8].neg();
            smat[8].dvd(& modulus_dv[1]);
            tmp.set_val_dfd1(& smat[2]);
            smat[12].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& smat[8]);
            smat[13].set_val_dfd1(& tmp);
            smat[14].set_val(1.0);
            smat[14].dvd(& modulus_dv[2]);
            smat[21].set_val(1.0);
            smat[21].dvd(& shear_mod_dv[0]);
            smat[28].set_val(1.0);
            smat[28].dvd(& shear_mod_dv[1]);
            smat[35].set_val(1.0);
            smat[35].dvd(& shear_mod_dv[2]);
            
            get_det_inv_ar_dfd1(&mut coef, &mut  ctmp, &mut  smat,  6,  0, &mut  x_vec, &mut  b_vec);
            ar_to_vec_dfd1(&mut ctmp, cmat,  0,  36);
        }
        
        return;
    }

    pub fn get_abd_dfd1(&self, cmat : &mut Vec<DiffDoub1>, lay_thk : &mut Vec<DiffDoub1>, lay_z : &mut Vec<DiffDoub1>, lay_q : &mut Vec<DiffDoub1>, lay_ang : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let num_lay : usize;
        let mut z_max = DiffDoub1::new();
        let mut z_min = DiffDoub1::new();
        let mut thk = DiffDoub1::new();
        let mut qmat = [DiffDoub1::new(); 9];
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        
        for i1 in 0..81 {
            cmat[i1].set_val(0.0);
        }
        
        num_lay = sec_ar[self.sect_ptr].layers.len();
        i2 = 0;
        for i1 in 0..num_lay {
            thk.set_val(0.5);
            thk.mult(& lay_thk[i1]);
            z_min.set_val_dfd1(& lay_z[i1]);
            z_min.sub(& thk);
            z_max.set_val_dfd1(& lay_z[i1]);
            z_max.add(& thk);
            self.transform_q_dfd1(&mut qmat, &mut lay_q[i2..], &mut lay_ang[i1]);
            
            // a matrix portion
            tmp.set_val_dfd1(& z_max);
            tmp.sub(& z_min);
            for i3 in 0..3 {
                i5 = 9 * i3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    tmp2.set_val_dfd1(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            
            // b matrix portion
            tmp.set_val_dfd1(& z_max);
            tmp.sqr();
            tmp2.set_val_dfd1(& z_min);
            tmp2.sqr();
            tmp.sub(& tmp2);
            tmp2.set_val(0.5);
            tmp.mult(& tmp2);
            for i3 in 0..3 {
                i5 = 9 * i3 + 3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    //i6 = i3*3 + i4;
                    //i5 = i3*9 + (i4 + 3)
                    tmp2.set_val_dfd1(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            
            // d matrix portion				
            tmp.set_val_dfd1(& z_max);
            tmp.sqr();
            tmp.mult(& z_max);
            tmp2.set_val_dfd1(& z_min);
            tmp2.sqr();
            tmp2.mult(& z_min);
            tmp.sub(& tmp2);
            tmp2.set_val(R_1O3);
            tmp.mult(& tmp2);
            for i3 in 0..3 {
                i5 = 9 * (i3 + 3) + 3;
                i6 = 3 * i3;
                for _i4 in 0..3 {
                    tmp2.set_val_dfd1(& tmp);
                    tmp2.mult(& qmat[i6]);
                    cmat[i5].add(& tmp2);
                    i5 += 1usize;
                    i6 += 1usize;
                }
            }
            i2  +=  9;
        }
        
        for i1 in 1..6 {
            i3 = 9 * i1;
            i4 = i1;
            for _i2 in 0..i1 {
                //i3 = i1*9 + i2
                //i4 = i2*9 + i1
                tmp.set_val_dfd1(& cmat[i4]);
                cmat[i3].set_val_dfd1(& tmp);
                i3 += 1usize;
                i4  +=  9;
            }
        }
        
        tmp.set_val(1.0);
        tmp.mult(& cmat[20]);
        cmat[60].set_val_dfd1(& tmp);
        cmat[70].set_val_dfd1(& tmp);
        cmat[80].set_val_dfd1(& tmp);
        
        return;
    }

    pub fn get_beam_stiff_dfd1(&self, cmat : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_ind : usize;
        let mut dv_comp : usize;
        let mut modulus_dv = [DiffDoub1::new(); 3];
        let mut shear_mod_dv = [DiffDoub1::new(); 3];
        let mut area_dv = DiffDoub1::new();
        let mut idv = [DiffDoub1::new(); 5];
        let mut jdv = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub1::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.stiffness[0] > 0.0 {
            for i1 in 0..36 {
                cmat[i1].set_val(this_sec.stiffness[i1]);
            }
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "stiffnessMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    cmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& cmat[i5]);
                    cmat[i4].set_val_dfd1(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..3 {
                modulus_dv[i1].set_val(this_mat.modulus[i1]);
                shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            jdv.set_val(this_sec.polar_moment);
            for dv in self.design_vars.iter() {
                dv_ind = dv.int_dat;
                this_dv = &dv_ar[dv_ind];
                if this_dv.category.s == "modulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    modulus_dv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "shearModulus" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    shear_mod_dv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if this_dv.category.s == "areaMoment" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            
            for i1 in 0..6 {
                i3 = 7 * i1;
                for _i2 in i1..6 {
                    cmat[i3].set_val(0.0);
                    i3 += 1usize;
                }
            }
            
            cmat[0].set_val_dfd1(& modulus_dv[0]);
            cmat[0].mult(& area_dv);
            cmat[4].set_val_dfd1(& modulus_dv[0]);
            cmat[4].mult(& idv[0]);
            cmat[5].set_val_dfd1(& modulus_dv[0]);
            cmat[5].neg();
            cmat[5].mult(& idv[1]);
            cmat[7].set_val_dfd1(& shear_mod_dv[0]);
            cmat[7].mult(& area_dv);
            cmat[9].set_val_dfd1(& shear_mod_dv[0]);
            cmat[9].neg();
            cmat[9].mult(& idv[0]);
            cmat[14].set_val_dfd1(& shear_mod_dv[1]);
            cmat[14].mult(& area_dv);
            cmat[15].set_val_dfd1(& shear_mod_dv[2]);
            cmat[15].mult(& idv[1]);
            cmat[21].set_val_dfd1(& shear_mod_dv[0]);
            cmat[21].mult(& jdv);
            cmat[28].set_val_dfd1(& modulus_dv[0]);
            cmat[28].mult(& idv[2]);
            cmat[29].set_val_dfd1(& modulus_dv[0]);
            cmat[29].neg();
            cmat[29].mult(& idv[4]);
            cmat[35].set_val_dfd1(& modulus_dv[0]);
            cmat[35].mult(& idv[3]);
            
            for i1 in 1..6 {
                i3 = 6 * i1;
                i4 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& cmat[i4]);
                    cmat[i3].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        return;
    }

    pub fn get_thermal_exp_dfd1(&self, th_exp : &mut Vec<DiffDoub1>, diff_exp : &mut Vec<DiffDoub1>, einit : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_comp : usize;
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec = &sec_ar[self.sect_ptr];
        let this_mat = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..6 {
            th_exp[i1].set_val(this_mat.expansion[i1]);
            diff_exp[i1].set_val(this_mat.diff_exp[i1]);
            einit[i1].set_val(0.0);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "thermalExp" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& tmp);
                th_exp[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "diffExp" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& tmp);
                diff_exp[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "initialStrain" {
                dv_comp = this_dv.component - 1;
                tmp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& tmp);
                einit[dv_comp].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_shell_exp_load_dfd1(&self, exp_ld : &mut Vec<DiffDoub1>, diff_ld : &mut Vec<DiffDoub1>, e0_ld : &mut Vec<DiffDoub1>, lay_thk : &mut Vec<DiffDoub1>, lay_z : &mut Vec<DiffDoub1>, 
        lay_q : &mut Vec<DiffDoub1>, lay_th_exp : &mut Vec<DiffDoub1>, lay_diff_exp : &mut Vec<DiffDoub1>, lay_einit : &mut Vec<DiffDoub1>, 
        lay_ang : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {

        let num_lay : usize;
        let mut qi : usize;
        let mut exi : usize;
        let mut sect_q = [DiffDoub1::new(); 9];
        let mut sect_te = [DiffDoub1::new(); 3];
        let mut sect_de = [DiffDoub1::new(); 3];
        let mut sect_e0 = [DiffDoub1::new(); 3];
        let mut qte_prod = [DiffDoub1::new(); 3];
        let mut qde_prod = [DiffDoub1::new(); 3];
        let mut qe0_prod = [DiffDoub1::new(); 3];
        let mut this_q = [DiffDoub1::new(); 9];
        let mut this_te = [DiffDoub1::new(); 3];
        let mut this_de = [DiffDoub1::new(); 3];
        let mut this_e0 = [DiffDoub1::new(); 3];
        let mut z_min = DiffDoub1::new();
        let mut z_max = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        
        for i1 in 0..6 {
            exp_ld[i1].set_val(0.0);
            e0_ld[i1].set_val(0.0);
        }
        
        num_lay = sec_ar[self.sect_ptr].layers.len();
        qi = 0;
        exi = 0;
        for layi in 0..num_lay {
            vec_to_ar_dfd1(&mut this_q, lay_q,  qi,  qi + 9);
            self.transform_q_dfd1(&mut sect_q, &mut  this_q, &mut lay_ang[layi]);
            
            vec_to_ar_dfd1(&mut this_te, lay_th_exp,  exi,  exi + 3);
            self.transform_strain_dfd1(&mut sect_te, &mut  this_te, &mut lay_ang[layi]);

            vec_to_ar_dfd1(&mut this_de, lay_diff_exp, exi, exi + 3);
            self.transform_strain_dfd1(&mut sect_de, &mut this_de, &mut lay_ang[layi]);
            
            vec_to_ar_dfd1(&mut this_e0, lay_einit,  exi,  exi + 3);
            self.transform_strain_dfd1(&mut sect_e0, &mut  this_e0, &mut  lay_ang[layi]);
            
            mat_mul_ar_dfd1(&mut qte_prod, &mut  sect_q, &mut  sect_te,  3,  3,  1);
            mat_mul_ar_dfd1(&mut qde_prod, &mut sect_q, &mut sect_de, 3, 3, 1);
            mat_mul_ar_dfd1(&mut qe0_prod, &mut  sect_q, &mut  sect_e0,  3,  3,  1);
            
            for i1 in 0..3 {
                tmp.set_val_dfd1(& qte_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                exp_ld[i1].add(& tmp);

                tmp.set_val_dfd1(& qde_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                diff_ld[i1].add(& tmp);

                tmp.set_val_dfd1(& qe0_prod[i1]);
                tmp.mult(& lay_thk[layi]);
                e0_ld[i1].add(& tmp);
            }
            
            tmp.set_val(0.5);
            tmp.mult(& lay_thk[layi]);
            z_min.set_val_dfd1(& lay_z[layi]);
            z_min.sub(& tmp);
            z_min.sqr();
            z_max.set_val_dfd1(& lay_z[layi]);
            z_max.add(& tmp);
            z_max.sqr();
            tmp.set_val(0.5);
            tmp2.set_val_dfd1(& z_max);
            tmp2.sub(& z_min);
            tmp.mult(& tmp2);// tmp = 0.5*(z_max^2 - z_min^2)
            for i1 in 0..3 {
                tmp2.set_val_dfd1(& qte_prod[i1]);
                tmp2.mult(& tmp);
                exp_ld[i1 + 3].add(& tmp2);

                tmp2.set_val_dfd1(& qde_prod[i1]);
                tmp2.mult(& tmp);
                diff_ld[i1 + 3].add(& tmp2);

                tmp2.set_val_dfd1(& qe0_prod[i1]);
                tmp2.mult(& tmp);
                e0_ld[i1 + 3].add(& tmp2);
            }
            qi  +=  9;
            exi  +=  3;
        }
        
        return;
    }

    pub fn get_beam_exp_load_dfd1(&self, exp_ld : &mut Vec<DiffDoub1>, diff_ld : &mut Vec<DiffDoub1>, e0_ld : &mut Vec<DiffDoub1>,
         sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {

        let mut i2 : usize;
        let mut dv_val = DiffDoub1::new();
        let mut cat : CppStr;
        let mut cat_list : CppStr;
        let mut dv_comp : usize;
        let mut mod_dv = [DiffDoub1::new(); 3];
        let mut shr_mod_dv = [DiffDoub1::new(); 3];
        let mut te_coef_dv = [DiffDoub1::new(); 6];
        let mut de_coef_dv = [DiffDoub1::new(); 6];
        let mut e0_dv = [DiffDoub1::new(); 6];
        let mut area_dv = DiffDoub1::new();
        let mut idv = [DiffDoub1::new(); 5];
        let mut qmat = [DiffDoub1::new(); 9];
        let mut qte = [DiffDoub1::new(); 3];
        let mut qde = [DiffDoub1::new(); 3];
        let mut qe0 = [DiffDoub1::new(); 3];
        let mut dedgu = [DiffDoub1::new(); 18];
        let mut tmp_exp = [DiffDoub1::new(); 6];
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub1::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];

        if this_sec.exp_load_coef[0] > 0.0 || this_sec.exp_load_coef[0] > 0.0 {
            for i1 in 0..6 {
                exp_ld[i1].set_val(this_sec.exp_load_coef[i1]);
                diff_ld[i1].set_val(this_sec.diff_load_coef[i1]);
                e0_ld[i1].set_val(0.0);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "thermalExp" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    exp_ld[dv_comp].add(& dv_val);
                }
                else if cat.s == "diffExp" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_ld[dv_comp].add(& dv_val);
                }
                else if cat.s == "initialStrain" {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    e0_ld[dv_comp].add(& dv_val);
                }
            }
        }
        else {
            for i1 in 0..3 {
                mod_dv[i1].set_val(this_mat.modulus[i1]);
                shr_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
                te_coef_dv[i1].set_val(this_mat.expansion[i1]);
                te_coef_dv[i1 + 3].set_val(this_mat.expansion[i1]);
                de_coef_dv[i1].set_val(this_mat.diff_exp[i1]);
                de_coef_dv[i1 + 3].set_val(this_mat.diff_exp[i1]);
                e0_dv[i1].set_val(0.0);
                e0_dv[i1 + 3].set_val(0.0);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            cat_list = CppStr::from("modulus shearModulus thermalExp diffExp initialStrain area areaMoment");
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                i2 = cat_list.find(&cat.s.as_str());
                if i2 < MAX_INT {
                    dv_comp = this_dv.component - 1;
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    if cat.s == "modulus" {
                        mod_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "shearModulus" {
                        shr_mod_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "thermalExp" {
                        te_coef_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "diffExp" {
                        de_coef_dv[dv_comp].add( & dv_val);
                    }
                    else if cat.s == "initialStrain" {
                        e0_dv[dv_comp].add(& dv_val);
                    }
                    else if cat.s == "area" {
                        area_dv.add(& dv_val);
                    }
                    else if cat.s == "areaMoment" {
                        idv[dv_comp].add(& dv_val);
                    }
                }
            }
            for i1 in 1..8 {
                qmat[i1].set_val(0.0);
            }
            qmat[0].set_val_dfd1(& mod_dv[0]);
            qmat[4].set_val_dfd1(& shr_mod_dv[0]);
            qmat[8].set_val_dfd1(& shr_mod_dv[1]);
            
            tmp.set_val_dfd1(& te_coef_dv[3]);
            te_coef_dv[1].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& te_coef_dv[4]);
            te_coef_dv[2].set_val_dfd1(& tmp);
            mat_mul_ar_dfd1(&mut qte, &mut  qmat, &mut  te_coef_dv,  3,  3,  1);

            tmp.set_val_dfd1(& de_coef_dv[3]);
            de_coef_dv[1].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& de_coef_dv[4]);
            de_coef_dv[2].set_val_dfd1(& tmp);
            mat_mul_ar_dfd1(&mut qde, &mut  qmat, &mut  de_coef_dv,  3,  3,  1);
            
            tmp.set_val_dfd1(& e0_dv[3]);
            e0_dv[1].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& e0_dv[4]);
            e0_dv[2].set_val_dfd1(& tmp);
            mat_mul_ar_dfd1(&mut qe0, &mut  qmat, &mut  e0_dv,  3,  3,  1);
            
            for i1 in 0..18 {
                dedgu[i1].set_val(0.0);
            }
            dedgu[0].set_val_dfd1(& area_dv);
            dedgu[4].set_val_dfd1(& area_dv);
            dedgu[8].set_val_dfd1(& area_dv);
            dedgu[10].set_val_dfd1(& idv[0]);
            dedgu[10].neg();
            dedgu[11].set_val_dfd1(& idv[1]);
            dedgu[12].set_val_dfd1(& idv[0]);
            dedgu[15].set_val_dfd1(& idv[1]);
            dedgu[15].neg();
            
            mat_mul_ar_dfd1(&mut tmp_exp, &mut  dedgu, &mut  qte,  6,  3,  1);
            ar_to_vec_dfd1(&mut tmp_exp, exp_ld,  0,  6);

            mat_mul_ar_dfd1(&mut tmp_exp, &mut  dedgu, &mut  qde,  6,  3,  1);
            ar_to_vec_dfd1(&mut tmp_exp, diff_ld,  0,  6);
            
            mat_mul_ar_dfd1(&mut tmp_exp, &mut  dedgu, &mut  qe0,  6,  3,  1);
            ar_to_vec_dfd1(&mut tmp_exp,  e0_ld,  0,  6);
        }
        
        return;
    }

    pub fn get_density_dfd1(&self, den : &mut DiffDoub1, layer : usize, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut dv_lay : usize;
        let mut coef = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let this_mat : &Material;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        if self.this_type == 3 || self.this_type == 41 {
            this_mat = &mat_ar[this_sec.get_layer_mat_ptr(layer)];
            den.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                dv_lay = this_dv.layer;
                if cat.s == "density" && dv_lay == layer {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    den.add(& dv_val);
                }
            }
        } else {
            this_mat = &mat_ar[this_sec.mat_ptr];
            den.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "density" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    den.add(& dv_val);
                }
            }
        }
        
        return;
    }

    pub fn get_shell_mass_dfd1(&self, mmat : &mut Vec<DiffDoub1>, lay_thk : &mut Vec<DiffDoub1>, lay_z : &mut Vec<DiffDoub1>, lay_den : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {
        let i1 : usize;
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut z_min = DiffDoub1::new();
        let mut z_min2 = DiffDoub1::new();
        let mut z_max = DiffDoub1::new();
        let mut z_max2 = DiffDoub1::new();
        
        for i1 in 0..36 {
            mmat[i1].set_val(0.0);
        }
        
        i1 = sec_ar[self.sect_ptr].layers.len();
        for layi in 0..i1 {
            tmp.set_val_dfd1(& lay_den[layi]);
            tmp.mult(& lay_thk[layi]);
            mmat[0].add(& tmp);
            mmat[7].add(& tmp);
            mmat[14].add(& tmp);
            
            tmp.set_val(0.5);
            tmp.mult(& lay_thk[layi]);
            z_min.set_val_dfd1(& lay_z[layi]);
            z_min.sub(& tmp);
            z_min2.set_val_dfd1(& z_min);
            z_min2.sqr();
            z_max.set_val_dfd1(& lay_z[layi]);
            z_max.add(& tmp);
            z_max2.set_val_dfd1(& z_max);
            z_max2.sqr();
            tmp.set_val_dfd1(& z_max2);
            tmp.sub(& z_min2);// tmp == z_max^2 - z_min^2
            tmp2.set_val(0.5);
            tmp2.mult(& lay_den[layi]);
            tmp2.mult(& tmp);// tmp2 = 0.5*rho*(z_max^2 - z_min^2)
            mmat[4].add(& tmp2);
            mmat[24].add(& tmp2);
            mmat[9].sub(& tmp2);
            mmat[19].sub(& tmp2);
            
            z_max2.mult(& z_max);
            z_min2.mult(& z_min);
            tmp.set_val_dfd1(& z_max2);
            tmp.sub(& z_min2);// tmp == z_max^3 - z_min^3
            tmp2.set_val(R_1O3);
            tmp2.mult(& lay_den[layi]);
            tmp2.mult(& tmp);
            mmat[21].add(& tmp2);
            mmat[28].add(& tmp2);
            mmat[35].add(& tmp2);
        }
        
        return;
    }

    pub fn get_beam_mass_dfd1(&self, mmat : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut dv_comp : usize;
        
        
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut d_cat : CppStr;
        
        
        let mut den_dv = DiffDoub1::new();
        let mut area_dv = DiffDoub1::new();
        let mut idv = [DiffDoub1::new(); 5];
        let mut jdv = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.mass[0] > 0.0 {
            for i1 in 0..36 {
                mmat[i1].set_val(this_sec.mass[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "massMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    mmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i3 = 6 * i1;// lower tri term
                i4 = i1;// upper tri term
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& mmat[i4]);
                    mmat[i3].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        else {
            for i1 in 0..36 {
                mmat[i1].set_val(0.0);
            }
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            area_dv.set_val(this_sec.area);
            jdv.set_val(this_sec.polar_moment);
            den_dv.set_val(this_mat.density);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                d_cat = this_dv.category.clone();
                if d_cat.s == "density" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    den_dv.add(& dv_val);
                }
                else if d_cat.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if d_cat.s == "areaMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    dv_comp = this_dv.component;
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if d_cat.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& area_dv);
            mmat[0].set_val_dfd1(& tmp);
            mmat[7].set_val_dfd1(& tmp);
            mmat[14].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& idv[0]);
            mmat[4].set_val_dfd1(& tmp);
            mmat[9].set_val_dfd1(& tmp);
            mmat[9].neg();
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& idv[1]);
            mmat[5].set_val_dfd1(& tmp);
            mmat[5].neg();
            mmat[15].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& jdv);
            mmat[21].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& idv[2]);
            mmat[28].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& idv[4]);
            mmat[29].set_val_dfd1(& tmp);
            tmp.set_val_dfd1(& den_dv);
            tmp.mult(& idv[3]);
            mmat[35].set_val_dfd1(& tmp);
            
            for i1 in 3..6 {
                i3 = 6 * i1;// lower tri term
                i4 = i1;// upper tri term
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& mmat[i4]);
                    mmat[i3].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        return;
    }

    pub fn get_solid_damp_dfd1(&self, dmat : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        
        let mut temp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..36 {
            dmat[i1].set_val(this_mat.damping[i1]);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "dampingMat" {
                this_dv.get_value_dfd1(&mut dv_val);
                dv_comp = this_dv.component;
                temp.set_val(dv.doub_dat);
                dv_val.mult(& temp);
                dmat[dv_comp].add(& dv_val);
            }
        }
        
        for i1 in 1..6 {
            i3 = 6 * i1;// lower tri term
            i4 = i1;// upper tri term
            for _i2 in 0..i1 {
                temp.set_val_dfd1(& dmat[i4]);
                dmat[i3].set_val_dfd1(& temp);
                i3 += 1usize;
                i4  +=  6;
            }
        }
        
        return;
    }

    pub fn get_shell_damp_dfd1(&self, dmat : &mut Vec<DiffDoub1>, lay_thk : &mut Vec<DiffDoub1>, lay_z : &mut Vec<DiffDoub1>, lay_d : &mut Vec<DiffDoub1>, lay_ang : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {
        self.get_abd_dfd1(dmat, lay_thk, lay_z, lay_d, lay_ang, sec_ar);
        dmat[60].set_val(0.0);
        dmat[70].set_val(0.0);
        dmat[80].set_val(0.0);
        return;
    }

    pub fn get_beam_damp_dfd1(&self, dmat : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut dv_comp : usize;
        let mut dmat_dv = [DiffDoub1::new(); 36];
        let mut area_dv = DiffDoub1::new();
        let mut idv = [DiffDoub1::new(); 5];
        let mut jdv = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        let mut tmp = DiffDoub1::new();
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        if this_sec.damping[0] > 0.0 {
            for i1 in 0..36 {
                dmat[i1].set_val(this_sec.damping[i1]);
            }
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "dampingMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    dmat[dv_comp].add(& dv_val);
                }
            }
            for i1 in 1..6 {
                i4 = 6 * i1;
                i5 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& dmat[i5]);
                    dmat[i4].set_val_dfd1(& tmp);
                    i4 += 1usize;
                    i5  +=  6;
                }
            }
        }
        else {
            for i1 in 0..36 {
                dmat_dv[i1].set_val(this_mat.damping[i1]);
            }
            area_dv.set_val(this_sec.area);
            for i1 in 0..5 {
                idv[i1].set_val(this_sec.area_moment[i1]);
            }
            jdv.set_val(this_sec.polar_moment);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "dampingMat" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    dmat_dv[dv_comp].add(& dv_val);
                }
                else if this_dv.category.s == "area" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    area_dv.add(& dv_val);
                }
                else if this_dv.category.s == "areaMoment" {
                    dv_comp = this_dv.component;
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    idv[dv_comp - 1].add(& dv_val);
                }
                else if this_dv.category.s == "polarMoment" {
                    coef.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& coef);
                    jdv.add(& dv_val);
                }
            }
            
            for i1 in 0..6 {
                i3 = 7 * i1;
                for _i2 in i1..6 {
                    dmat[i3].set_val(0.0);
                    i3 += 1usize;
                }
            }
            
            dmat[0].set_val_dfd1(& dmat_dv[0]);
            dmat[0].mult(& area_dv);
            dmat[4].set_val_dfd1(& dmat_dv[0]);
            dmat[4].mult(& idv[0]);
            dmat[5].set_val_dfd1(& dmat_dv[0]);
            dmat[5].neg();
            dmat[5].mult(& idv[1]);
            dmat[7].set_val_dfd1(& dmat_dv[21]);
            dmat[7].mult(& area_dv);
            dmat[9].set_val_dfd1(& dmat_dv[21]);
            dmat[9].neg();
            dmat[9].mult(& idv[0]);
            dmat[14].set_val_dfd1(& dmat_dv[28]);// 28
            dmat[14].mult(& area_dv);
            dmat[15].set_val_dfd1(& dmat_dv[35]);// 35
            dmat[15].mult(& idv[1]);
            dmat[21].set_val_dfd1(& dmat_dv[21]);// 21
            dmat[21].mult(& jdv);
            dmat[28].set_val_dfd1(& dmat_dv[0]);
            dmat[28].mult(& idv[2]);
            dmat[29].set_val_dfd1(& dmat_dv[0]);
            dmat[29].neg();
            dmat[29].mult(& idv[4]);
            dmat[35].set_val_dfd1(& dmat_dv[0]);
            dmat[35].mult(& idv[3]);
            
            for i1 in 1..6 {
                i3 = 6 * i1;
                i4 = i1;
                for _i2 in 0..i1 {
                    tmp.set_val_dfd1(& dmat[i4]);
                    dmat[i3].set_val_dfd1(& tmp);
                    i3 += 1usize;
                    i4  +=  6;
                }
            }
        }
        
        
        return;
    }

    pub fn get_conductivity_dfd1(&self, t_cond : &mut Vec<DiffDoub1>, diff : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cond_dv = [DiffDoub1::new(); 6];
        let mut diff_dv = [DiffDoub1::new(); 6];
        
        let mut temp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut dv_comp : usize;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        for i1 in 0..6 {
            cond_dv[i1].set_val(this_mat.conductivity[i1]);
            diff_dv[i1].set_val(this_mat.diffusivity[i1]);
        }
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "thermalCond" {
                dv_comp = this_dv.component - 1;
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& temp);
                cond_dv[dv_comp].add(& dv_val);
            }
            else if this_dv.category.s == "diffusivity" {
                dv_comp = this_dv.component - 1;
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& temp);
                diff_dv[dv_comp].add(& dv_val);
            }
        }
        
        t_cond[0] = cond_dv[0];
        t_cond[1] = cond_dv[3];
        t_cond[2] = cond_dv[4];
        t_cond[3] = cond_dv[3];
        t_cond[4] = cond_dv[1];
        t_cond[5] = cond_dv[5];
        t_cond[6] = cond_dv[4];
        t_cond[7] = cond_dv[5];
        t_cond[8] = cond_dv[2];

        diff[0] = diff_dv[0];
        diff[1] = diff_dv[3];
        diff[2] = diff_dv[4];
        diff[3] = diff_dv[3];
        diff[4] = diff_dv[1];
        diff[5] = diff_dv[5];
        diff[6] = diff_dv[4];
        diff[7] = diff_dv[5];
        diff[8] = diff_dv[2];
        
        return;
    }

    pub fn get_shell_cond_dfd1(&self, t_cond : &mut Vec<DiffDoub1>, diff : &mut Vec<DiffDoub1>, lay_thk : &mut Vec<DiffDoub1>, lay_ang : &mut Vec<DiffDoub1>, 
        lay_cond : &mut Vec<DiffDoub1>, lay_diff : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {

        let mut i2 : usize;
        let num_lay : usize =  sec_ar[self.sect_ptr].layers.len();
        let mut layer_mat = [DiffDoub1::new(); 9];
        let mut al_mat = [DiffDoub1::new(); 9];
        let mut al_t = [DiffDoub1::new(); 9];
        let mut tmp = [DiffDoub1::new(); 9];
        let mut this_cond = [DiffDoub1::new(); 9];
        let mut this_diff = [DiffDoub1::new(); 9];
        let mut tmp2 = DiffDoub1::new();
        
        for i1 in 0..9 {
            t_cond[i1].set_val(0.0);
            diff[i1].set_val(0.0);
            al_mat[i1].set_val(0.0);
        }
        al_mat[8].set_val(1.0);
        
        for i1 in 0..num_lay {
            al_mat[0].set_val(R_PI_0180);
            al_mat[0].mult(& lay_ang[i1]);
            al_mat[0].cs();
            tmp2.set_val_dfd1(& al_mat[0]);
            al_mat[4].set_val_dfd1(& tmp2);
            al_mat[3].set_val(R_PI_0180);
            al_mat[3].mult(& lay_ang[i1]);
            al_mat[3].sn();
            tmp2.set_val_dfd1(& al_mat[3]);
            al_mat[1].set_val_dfd1(& tmp2);
            al_mat[1].neg();
            transpose_ar_dfd1(&mut al_t, &mut  al_mat,  3,  3);
            i2 = 9 * i1;
            
            vec_to_ar_dfd1(&mut this_cond, lay_cond,  i2,  i2 + 9);
            mat_mul_ar_dfd1(&mut tmp, &mut  this_cond, &mut  al_t,  3,  3,  3);
            mat_mul_ar_dfd1(&mut layer_mat, &mut  al_mat, &mut  tmp,  3,  3,  3);
            for i2 in 0..9 {
                layer_mat[i2].mult(& lay_thk[i1]);
                t_cond[i2].add(& layer_mat[i2]);
            }

            vec_to_ar_dfd1(&mut this_diff, lay_diff, i2,  i2 + 9);
            mat_mul_ar_dfd1(&mut tmp, &mut  this_diff, &mut  al_t,  3,  3,  3);
            mat_mul_ar_dfd1(&mut layer_mat, &mut  al_mat, &mut  tmp,  3,  3,  3);
            for i2 in 0..9 {
                layer_mat[i2].mult(& lay_thk[i1]);
                diff[i2].add(& layer_mat[i2]);
            }
        }
        
        return;
    }

    pub fn get_beam_cond_dfd1(&self, t_cond : &mut Vec<DiffDoub1>, diff : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut cond_dv = DiffDoub1::new();
        let mut diff_dv = DiffDoub1::new();
        let mut area_dv = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];

        if this_sec.conductivity > 0.0 || this_sec.diffusivity > 0.0 {
            cond_dv.set_val(this_sec.conductivity);
            diff_dv.set_val(this_sec.diffusivity);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "thermCond" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv.add(& dv_val);
                }
                else if this_dv.category.s == "diffusivity" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv.add(& dv_val);
                }
            }
        }
        else {
            cond_dv.set_val(this_mat.conductivity[0]);
            diff_dv.set_val(this_mat.diffusivity[0]);
            area_dv.set_val(this_sec.area);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "thermCond" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    cond_dv.add(& dv_val);
                }
                else if cat.s == "diffusivity" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    diff_dv.add(& dv_val);
                }
                else if cat.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    area_dv.add(& dv_val);
                }
            }
            cond_dv.mult(& area_dv);
            diff_dv.mult(& area_dv);
        }
        
        for i1 in 1..8 {
            t_cond[i1].set_val(0.0);
            diff[i1].set_val(0.0);
        }
        
        t_cond[0].set_val_dfd1(& cond_dv);
        t_cond[4].set_val_dfd1(& cond_dv);
        t_cond[8].set_val_dfd1(& cond_dv);

        diff[0].set_val_dfd1(& diff_dv);
        diff[4].set_val_dfd1(& diff_dv);
        diff[8].set_val_dfd1(& diff_dv);
        
        return;
    }

    pub fn get_specific_heat_dfd1(&self, spec_heat : &mut DiffDoub1, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut temp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        spec_heat.set_val(this_mat.spec_heat);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            if this_dv.category.s == "specHeat" {
                temp.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& temp);
                spec_heat.add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_shell_spec_heat_dfd1(&self, spec_heat : &mut DiffDoub1, lay_thk : &mut Vec<DiffDoub1>, lay_sh : &mut Vec<DiffDoub1>, lay_den : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>) {
        let num_lay : usize;
        let mut tmp = DiffDoub1::new();
        
        spec_heat.set_val(0.0);
        num_lay = sec_ar[self.sect_ptr].layers.len();
        for i1 in 0..num_lay {
            tmp.set_val_dfd1(& lay_den[i1]);
            tmp.mult(& lay_sh[i1]);
            tmp.mult(& lay_thk[i1]);
            spec_heat.add(& tmp);
        }
        
        return;
    }

    pub fn get_beam_spec_heat_dfd1(&self, spec_heat : &mut DiffDoub1, sec_ar : &mut Vec<Section>, mat_ar : &mut Vec<Material>, dv_ar : & Vec<DesignVariable>) {
        let mut cat : CppStr;
        let mut density_dv = DiffDoub1::new();
        let mut area_dv = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut dv_val = DiffDoub1::new();
        
        let mut this_dv : &DesignVariable;
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        let this_mat : &Material = &mat_ar[this_sec.mat_ptr];
        let sec_sh : f64 =  this_sec.spec_heat;
        let mat_sh : f64 =  this_mat.spec_heat;
        let mat_den : f64 =  this_mat.density;
        
        if sec_sh > 0.0 {
            spec_heat.set_val(sec_sh);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                if this_dv.category.s == "specHeat" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    spec_heat.add(& dv_val);
                }
            }
        }
        else {
            spec_heat.set_val(mat_sh);
            density_dv.set_val(mat_den);
            area_dv.set_val(this_sec.area);
            for dv in self.design_vars.iter() {
                this_dv = &dv_ar[dv.int_dat];
                cat = this_dv.category.clone();
                if cat.s == "specHeatDV" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    spec_heat.add(& dv_val);
                }
                else if cat.s == "density" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    density_dv.add(& dv_val);
                }
                else if cat.s == "area" {
                    tmp.set_val(dv.doub_dat);
                    this_dv.get_value_dfd1(&mut dv_val);
                    dv_val.mult(& tmp);
                    area_dv.add(& dv_val);
                }
            }
            spec_heat.mult(& density_dv);
            spec_heat.mult(& area_dv);
        }
        
        return;
    }

    pub fn get_nd_crds_dfd1(&self, x_glob : &mut Vec<DiffDoub1>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
        let mut nd_crd = [DiffDoub1::new(); 3];
        let mut this_nd : &Node;
        
        for i1 in 0..self.num_nds {
            this_nd = &nd_ar[self.nodes[i1]];
            this_nd.get_crd_dfd1(&mut nd_crd, dv_ar);
            x_glob[i1].set_val_dfd1(& nd_crd[0]);
            x_glob[i1+self.num_nds].set_val_dfd1(& nd_crd[1]);
            x_glob[i1+2*self.num_nds].set_val_dfd1(& nd_crd[2]);
        }
        
        return;
    }

    pub fn get_loc_ori_dfd1(&self, loc_ori : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut i1 : usize;
        let mut dv_cat : CppStr;
        let mut rot = [DiffDoub1::new(); 3];
        let mut dv_val = DiffDoub1::new();
        let mut coef = DiffDoub1::new();
        let mut ori_copy = [DiffDoub1::new(); 9];
        let mut tmp_ori = [DiffDoub1::new(); 9];
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        for i1 in 0..9 {
            ori_copy[i1].set_val(this_sec.orientation[i1]);
        }
        
        rot[0].set_val(0.0);
        rot[1].set_val(0.0);
        rot[2].set_val(0.0);
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            dv_cat = this_dv.category.clone();
            if dv_cat.s == "orientation" {
                i1 = this_dv.component - 1;
                coef.set_val(dv.doub_dat);
                this_dv.get_value_dfd1(&mut dv_val);
                dv_val.mult(& coef);
                rot[i1].add(& dv_val);
            }
        }
        
        rotate_orient_dfd1(&mut tmp_ori, &mut  ori_copy, &mut  rot);
        ar_to_vec_dfd1(&mut tmp_ori, loc_ori,  0,  9);
        
        return;
    }

    pub fn correct_orient_dfd1(&self, loc_ori : &mut Vec<DiffDoub1>, x_glob : &mut Vec<DiffDoub1>) {
        let mut v1 = [DiffDoub1::new(); 3];
        let mut v2 = [DiffDoub1::new(); 3];
        let mut v3 = [DiffDoub1::new(); 3];
        let mut rot = [DiffDoub1::new(); 3];
        let mut ori_copy = [DiffDoub1::new(); 9];
        let mut tmp_ori = [DiffDoub1::new(); 9];
        
        let mut dp = DiffDoub1::new();
        let mut magv3 = DiffDoub1::new();
        let mut mag_cp = DiffDoub1::new();
        let mut theta = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        
        for i1 in 0..9 {
            ori_copy[i1].set_val_dfd1(& loc_ori[i1]);
        }
        
        if self.this_type == 3 || self.this_type == 41 {
            if self.this_type == 3 {
                v1[0].set_val_dfd1(& x_glob[1]);
                v1[0].sub(& x_glob[0]);
                v1[1].set_val_dfd1(& x_glob[4]);
                v1[1].sub(& x_glob[3]);
                v1[2].set_val_dfd1(& x_glob[7]);
                v1[2].sub(& x_glob[6]);
                
                v2[0].set_val_dfd1(& x_glob[2]);
                v2[0].sub(& x_glob[0]);
                v2[1].set_val_dfd1(& x_glob[5]);
                v2[1].sub(& x_glob[3]);
                v2[2].set_val_dfd1(& x_glob[8]);
                v2[2].sub(& x_glob[6]);
            } else {
                v1[0].set_val_dfd1(& x_glob[2]);
                v1[0].sub(& x_glob[0]);
                v1[1].set_val_dfd1(& x_glob[6]);
                v1[1].sub(& x_glob[4]);
                v1[2].set_val_dfd1(& x_glob[10]);
                v1[2].sub(& x_glob[8]);
                
                v2[0].set_val_dfd1(& x_glob[3]);
                v2[0].sub(& x_glob[1]);
                v2[1].set_val_dfd1(& x_glob[7]);
                v2[1].sub(& x_glob[5]);
                v2[2].set_val_dfd1(& x_glob[11]);
                v2[2].sub(& x_glob[9]);
            }
            cross_prod_dfd1(&mut v3, &mut  v1, &mut  v2);
            
            dp.set_val_dfd1(& v3[0]);
            dp.mult(& loc_ori[6]);
            tmp.set_val_dfd1(& v3[1]);
            tmp.mult(& loc_ori[7]);
            dp.add(& tmp);
            tmp.set_val_dfd1(& v3[2]);
            tmp.mult(& loc_ori[8]);
            dp.add(& tmp);
            
            if dp.val < 0.0 {
                v3[0].neg();
                v3[1].neg();
                v3[2].neg();
            }
            
            magv3.set_val_dfd1(& v3[0]);
            magv3.sqr();
            tmp.set_val_dfd1(& v3[1]);
            tmp.sqr();
            magv3.add(& tmp);
            tmp.set_val_dfd1(& v3[2]);
            tmp.sqr();
            magv3.add(& tmp);
            magv3.sqt();
            
            cross_prod_dfd1(&mut rot, &mut loc_ori[6..], &mut v3);
            mag_cp.set_val_dfd1(& rot[0]);
            mag_cp.sqr();
            tmp.set_val_dfd1(& rot[1]);
            tmp.sqr();
            mag_cp.add(& tmp);
            tmp.set_val_dfd1(& rot[2]);
            tmp.sqr();
            mag_cp.add(& tmp);
            mag_cp.sqt();
            if mag_cp.val < 1e-12 {
                return;
            }
            
            theta.set_val_dfd1(& mag_cp);
            theta.dvd(& magv3);
            theta.asn();
            
            tmp.set_val_dfd1(& theta);
            tmp.dvd(& mag_cp);
            
            rot[0].mult(& tmp);
            rot[1].mult(& tmp);
            rot[2].mult(& tmp);
            
            rotate_orient_dfd1(&mut tmp_ori, &mut  ori_copy, &mut  rot);
            ar_to_vec_dfd1(&mut tmp_ori, loc_ori,  0,  9);
        }
        
        return;
    }

    pub fn get_frc_fld_const_dfd1(&self, coef : &mut Vec<DiffDoub1>, exp : &mut Vec<DiffDoub1>, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        coef[0].set_val(this_sec.pot_coef);
        coef[1].set_val(this_sec.damp_coef);
        exp[0].set_val(this_sec.pot_exp);
        exp[1].set_val(this_sec.damp_exp);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "potFldCoef" {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[0].add(& dv_val);
            }
            else if cat.s == "dampFldCoef" {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[1].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_thrm_fld_const_dfd1(&self, coef : &mut Vec<DiffDoub1>, ref_t : &mut DiffDoub1, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec : &Section = &sec_ar[self.sect_ptr];
        coef[0].set_val(this_sec.cond_coef);
        coef[1].set_val(this_sec.rad_coef);
        ref_t.set_val(this_sec.ref_temp);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "condCoef" {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[0].add(& dv_val);
            }
            else if cat.s == "radCoef" {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                coef[1].add(& dv_val);
            }
        }
        
        return;
    }

    pub fn get_mass_per_el_dfd1(&self, mass_per_el : &mut DiffDoub1, sec_ar : &mut Vec<Section>, dv_ar : & Vec<DesignVariable>) {
        let mut dv_val = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut cat : CppStr;
        let mut this_dv : &DesignVariable;
        
        let this_sec = &sec_ar[self.sect_ptr];
        mass_per_el.set_val(this_sec.mass_per_el);
        
        for dv in self.design_vars.iter() {
            this_dv = &dv_ar[dv.int_dat];
            cat = this_dv.category.clone();
            if cat.s == "massPerEl" {
                this_dv.get_value_dfd1(&mut dv_val);
                tmp.set_val(dv.doub_dat);
                dv_val.mult(& tmp);
                mass_per_el.add(& dv_val);
            }
        }
        
        return;
    }

    //end dup
 
//end skip 
 
 
 
 
 
 
 
 
 
}


