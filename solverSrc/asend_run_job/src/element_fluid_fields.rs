use crate::element::*;
use crate::diff_doub::*;
use crate::design_var::*;
use crate::node::*;
use crate::section::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;

impl Element {
    //dup1

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

    pub fn get_fluid_prereq_dfd0(&mut self, pre : &mut DiffDoub0FlPrereq, sec_ar : &mut Vec<Section>, fl_ar : &mut Vec<Fluid>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
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

        for i in 0..3 {
            pre.f_per_mass[i].set_val(self.body_force[i]);
        }
        self.get_f_per_mass_dfd0(&mut pre.f_per_mass, dv_ar);
        pre.hg_per_vol.set_val(self.body_heat_gen);
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

            pre.den_pres_coef.set_val(0.0);
            pre.temp_pres_coef.set_val(0.0);
            
            pre.turb_vis_coef.set_val(0.0);
            pre.temp_vis_coef.set_val(0.0);
        }
        else {
            pre.den_vis_coef.set_val(0.0);

            pre.den_pres_coef.set_val_dfd0(&pre.bulk_mod);
            pre.den_pres_coef.dvd(&pre.ref_den);

            pre.temp_pres_coef.set_val(3.0);
            pre.temp_pres_coef.mult(&pre.bulk_mod);
            pre.temp_pres_coef.mult(&pre.expansion);
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

        self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &mut spt);
        mat_mul_ar_dfd0(v_grad, &pre.fl_vel, &d_ndx, 3, self.num_nds, 3);
        mat_mul_ar_dfd0(t_grad, &pre.fl_temp, &d_ndx, 1, self.num_nds, 3);
        
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

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

    pub fn get_fluid_prereq_dfd1(&mut self, pre : &mut DiffDoub1FlPrereq, sec_ar : &mut Vec<Section>, fl_ar : &mut Vec<Fluid>, nd_ar : &mut Vec<Node>, dv_ar : & Vec<DesignVariable>) {
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

        for i in 0..3 {
            pre.f_per_mass[i].set_val(self.body_force[i]);
        }
        self.get_f_per_mass_dfd1(&mut pre.f_per_mass, dv_ar);
        pre.hg_per_vol.set_val(self.body_heat_gen);
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

            pre.den_pres_coef.set_val(0.0);
            pre.temp_pres_coef.set_val(0.0);
            
            pre.turb_vis_coef.set_val(0.0);
            pre.temp_vis_coef.set_val(0.0);
        }
        else {
            pre.den_vis_coef.set_val(0.0);

            pre.den_pres_coef.set_val_dfd1(&pre.bulk_mod);
            pre.den_pres_coef.dvd(&pre.ref_den);

            pre.temp_pres_coef.set_val(3.0);
            pre.temp_pres_coef.mult(&pre.bulk_mod);
            pre.temp_pres_coef.mult(&pre.expansion);
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

        self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &mut spt);
        mat_mul_ar_dfd1(v_grad, &pre.fl_vel, &d_ndx, 3, self.num_nds, 3);
        mat_mul_ar_dfd1(t_grad, &pre.fl_temp, &d_ndx, 1, self.num_nds, 3);
        
    }

    //end dup
 
//end skip 
 
 
 
}