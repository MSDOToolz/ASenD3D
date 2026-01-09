use crate::fluid_domain::cell::*;
use crate::fluid_domain::face::*;
use crate::fluid_domain::node::*;
use crate::fluid_domain::load::*;
use crate::matrix_functions::*;
use crate::constants::*;

impl Cell {

    pub fn get_eqn_prereq(&self, pre : &mut EqnPrereq, nd_ar : &Vec<Node>, sd_ar : &Vec<SubDomain>, fl_ar : &Vec<Fluid>, dv_ar : &Vec<DesignVariable>) {
        let mut dfcrd = [DiffDoub1::new(); 3];
        let mut nd : &Node;

        let mut k : usize;
        for i in 0..4 {
            nd = &nd_ar[self.nodes[i]];
            pre.fl_den[i].set_val_dfd1(&nd.fl_den);
            pre.fl_den_dot[i].set_val_dfd1(&nd.fl_den_dot);
            pre.temp[i].set_val_dfd1(&nd.temperature);
            pre.temp_dot[i].set_val_dfd1(&nd.temp_dot);
            pre.turb[i].set_val_dfd1(&nd.turb_e);
            pre.turb_dot[i].set_val_dfd1(&nd.turb_e_dot);

            nd.get_def_crd(&mut dfcrd);
            k = i;
            for j in 0..3 {
                pre.def_coord[k].set_val_dfd1(&dfcrd[j]);
                pre.fl_vel[k].set_val_dfd1(&nd.fl_vel[j]);
                pre.fl_vel_dot[k].set_val_dfd1(&nd.fl_vel_dot[j]);
                pre.vrel[k].set_val_dfd1(&nd.fl_vel[j]);
                pre.vrel[k].sub(&nd.velocity[j]);
                k += 4;
            }
        }

        let fld = &fl_ar[sd_ar[self.sub_dom_pt].fluid_ptr];

        pre.viscosity.set_val(fld.viscosity);
        pre.conductivity.set_val(fld.therm_cond);
        pre.expansion.set_val(fld.expansion);
        pre.spec_heat.set_val(fld.spec_heat);
        pre.ideal_gas.set_val(fld.ideal_gas);
        pre.bulk_mod.set_val(fld.bulk_modulus);
        pre.ref_temp.set_val(fld.ref_temp);
        pre.ref_pres.set_val(fld.ref_pres);
        pre.ref_den.set_val(fld.ref_den);
        pre.ref_enth.set_val(fld.ref_enth);
        pre.temp_vis_coef.set_val(fld.temp_vis_coef);
        pre.turb_vis_coef.set_val(fld.turb_vis_coef);
        pre.grad_turb_coef.set_val(fld.grad_turb_coef);
        pre.diss_turb_coef.set_val(fld.diss_turb_coef);
        pre.compressible = fld.compressible;

        let mut dv : &DesignVariable;
        let mut tmp = DiffDoub1::new();
        for dvi in self.dvars.iter() {
            dv = &dv_ar[dvi.int_dat];
            tmp.set_val(dvi.doub_dat);
            tmp.mult(&dv.diff_val);
            match dv.category.s.as_str() {
                "viscosity" => pre.viscosity.add(&tmp),
                "conductivity" => pre.conductivity.add(&tmp),
                "expansion" => pre.expansion.add(&tmp),
                "specHeat" => pre.spec_heat.add(&tmp),
                "idealGasConst" => pre.ideal_gas.add(&tmp),
                "bulkModulus" => pre.bulk_mod.add(&tmp),
                "refTemp" => pre.ref_temp.add(&tmp),
                "refPres" => pre.ref_pres.add(&tmp),
                "refDen" => pre.ref_den.add(&tmp),
                "refEnth" => pre.ref_enth.add(&tmp),
                "tempVisCoef" => pre.temp_vis_coef.add(&tmp),
                "turbVisCoef" => pre.turb_vis_coef.add(&tmp),
                "gradTurbCoef" => pre.grad_turb_coef.add(&tmp),
                "dissTurbCoef" => pre.diss_turb_coef.add(&tmp),
                &_ => (),
            }
        }


    }

    pub fn update_volume(&mut self, pre : &EqnPrereq) {
        let mut vmat = [DiffDoub1::new(); 9];
        
        vmat[0].set_val_dfd1(&pre.def_coord[1]);
        vmat[0].sub(&pre.def_coord[0]);

        vmat[1].set_val_dfd1(&pre.def_coord[5]);
        vmat[1].sub(&pre.def_coord[4]);

        vmat[2].set_val_dfd1(&pre.def_coord[9]);
        vmat[2].sub(&pre.def_coord[8]);

        vmat[3].set_val_dfd1(&pre.def_coord[2]);
        vmat[3].sub(&pre.def_coord[0]);

        vmat[4].set_val_dfd1(&pre.def_coord[6]);
        vmat[4].sub(&pre.def_coord[4]);

        vmat[5].set_val_dfd1(&pre.def_coord[10]);
        vmat[5].sub(&pre.def_coord[8]);

        vmat[6].set_val_dfd1(&pre.def_coord[3]);
        vmat[6].sub(&pre.def_coord[0]);

        vmat[7].set_val_dfd1(&pre.def_coord[7]);
        vmat[7].sub(&pre.def_coord[4]);

        vmat[8].set_val_dfd1(&pre.def_coord[11]);
        vmat[8].sub(&pre.def_coord[8]);

        q_rfactor_ar_dfd1(&mut vmat, 3, 0, 2, 0, 2, 0);
        self.volume.set_val(R_1O6);
        self.volume.mult(&vmat[0]);
        self.volume.mult(&vmat[4]);
        self.volume.mult(&vmat[8]);

    }

    pub fn update_cent_data(&self, pre : &EqnPrereq, cdat : &mut Vec<CellData>) {
        let cd = &mut cdat[self.label];

        cd.den.set_val(0.0);
        cd.temp.set_val(0.0);
        cd.turb.set_val(0.0);
        for j in 0..3 {
            cd.vel[j].set_val(0.0);
            cd.v_rel[j].set_val(0.0);
        }

        let mut k : usize;
        for i in 0..4 {
            cd.den.add(&pre.fl_den[i]);
            cd.temp.add(&pre.temp[i]);
            cd.turb.add(&pre.turb[i]);
            k = i;
            for j in 0..3 {
                cd.vel[j].add(&pre.fl_vel[k]);
                cd.v_rel[j].add(&pre.vrel[k]);
                k += 4;
            }
        }

        let mut tmp = DiffDoub1::new();
        tmp.set_val(0.25);
        
        cd.den.mult(&tmp);
        cd.temp.mult(&tmp);
        cd.turb.mult(&tmp);
        for j in 0..3 {
            cd.vel[j].mult(&tmp);
            cd.v_rel[j].mult(&tmp);
        }

        let mut dnds = [DiffDoub1::new(); 12];
        
        dnds[0].set_val(-1.0);
        dnds[1].set_val(-1.0);
        dnds[2].set_val(-1.0);
        dnds[3].set_val(1.0);
        dnds[7].set_val(1.0);
        dnds[11].set_val(1.0);

        let mut dxds = [DiffDoub1::new(); 9];
        let mut det = DiffDoub1::new();
        let mut dsdx = [DiffDoub1::new(); 9];
        let mut x = [DiffDoub1::new(); 3];
        let mut b = [DiffDoub1::new(); 3];

        mat_mul_ar_dfd1(&mut dxds, &pre.def_coord, &dnds, 3, 4, 3);
        get_det_inv_ar_dfd1(&mut det, &mut dsdx, &mut dxds, 3, 0, &mut x, &mut b);

        let mut dndx = [DiffDoub1::new(); 12];
        mat_mul_ar_dfd1(&mut dndx, &dnds, &dsdx, 4, 3, 3);
        
        mat_mul_ar_dfd1(&mut cdat[self.label].v_grad, &pre.fl_vel, &dndx, 3, 4, 3);
        mat_mul_ar_dfd1(&mut cdat[self.label].t_grad, &pre.temp, &dndx, 1, 4, 3);

    }

    pub fn unsteady(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq) {
        let mut den = DiffDoub1::new();
        let mut den_dot = DiffDoub1::new();
        let mut vel = [DiffDoub1::new(); 3];
        let mut vel_dot = [DiffDoub1::new(); 3];
        let mut temp = DiffDoub1::new();
        let mut temp_dot = DiffDoub1::new();
        let mut turb = DiffDoub1::new();
        let mut turb_dot = DiffDoub1::new();

        let mut enth = DiffDoub1::new();
        let mut enth_dot = DiffDoub1::new();

        let mut k : usize;
        for i in 0..4 {
            den.add(&pre.fl_den[i]);
            den_dot.add(&pre.fl_den_dot[i]);
            temp.add(&pre.temp[i]);
            temp_dot.add(&pre.temp_dot[i]);
            turb.add(&pre.turb[i]);
            turb_dot.add(&pre.turb_dot[i]);
            k = i;
            for j in 0..3 {
                vel[j].add(&pre.fl_vel[k]);
                vel_dot[j].add(&pre.fl_vel_dot[k]);
                k += 4;
            }
        }
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        tmp.set_val(0.25);
        den.mult(&tmp);
        den_dot.mult(&tmp);
        temp.mult(&tmp);
        temp_dot.mult(&tmp);
        turb.mult(&tmp);
        turb_dot.mult(&tmp);
        for i in 0..3 {
            vel[i].mult(&tmp);
            vel_dot[i].mult(&tmp);
        }
        
        let mut r_vec = [DiffDoub1::new(); 6];
        
        // mass

        r_vec[0].set_val_dfd1(&den_dot);
        r_vec[0].mult(&self.volume);

        // momentum

        for i in 0..3 {
            tmp.set_val(0.0);
            tmp2.set_val_dfd1(&den_dot);
            tmp2.mult(&vel[i]);
            tmp.add(&tmp2);
            tmp2.set_val_dfd1(&den);
            tmp2.mult(&vel_dot[i]);
            tmp.add(&tmp2);
            tmp.mult(&self.volume);
            r_vec[i+1].set_val_dfd1(&tmp);
        }

        // energy

        tmp.set_val_dfd1(&pre.spec_heat);
        tmp.mult(&temp);
        enth.set_val_dfd1(&pre.ref_enth);
        enth.add(&tmp);

        enth_dot.set_val_dfd1(&pre.spec_heat);
        enth_dot.mult(&temp_dot);

        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&vel[i]);
            tmp2.sqr();
            tmp.add(&tmp2);
        }
        tmp2.set_val(0.5);
        tmp.mult(&tmp2);

        tmp.add(&enth);
        tmp.add(&turb);
        tmp.mult(&den_dot);
        r_vec[4].set_val_dfd1(&tmp);

        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&vel[i]);
            tmp2.mult(&vel_dot[i]);
            tmp.add(&tmp2);
        }

        tmp.add(&enth_dot);
        tmp.add(&turb_dot);
        tmp.mult(&den);
        r_vec[4].add(&tmp);

        r_vec[4].mult(&self.volume);

        // turbulence

        tmp.set_val(0.0);
        tmp2.set_val_dfd1(&den_dot);
        tmp2.mult(&turb);
        tmp.add(&tmp2);
        tmp2.set_val_dfd1(&den);
        tmp2.mult(&turb_dot);
        tmp.add(&tmp2);
        tmp.mult(&self.volume);

        r_vec[5].set_val_dfd1(&tmp);

        // construct output

        k = 0;
        for _i in 0..4 {
            for j in 0..6 {
                r_out[k].add(&r_vec[j]);
                k += 1;
            }
        }

    }

    pub fn turbulence(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, cdat : &Vec<CellData>) {
        let mut den = DiffDoub1::new();
        let mut turb = DiffDoub1::new();
        let mut v_mag = DiffDoub1::new();

        let mut k : usize;
        for i in 0..4 {
            den.add(&pre.fl_den[i]);
            turb.add(&pre.turb[i]);
        }
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        tmp.set_val(0.25);
        den.mult(&tmp);
        turb.mult(&tmp);

        let vg = &cdat[self.label].v_grad;
        for i in 0..9 {
            tmp.set_val_dfd1(&vg[i]);
            tmp.sqr();
            v_mag.add(&tmp);
        }
        v_mag.sqt();

        tmp.set_val_dfd1(&pre.grad_turb_coef);
        tmp.mult(&v_mag);
        tmp.mult(&den);
        tmp.neg();
        tmp2.set_val_dfd1(&tmp);

        tmp.set_val_dfd1(&pre.diss_turb_coef);
        tmp.mult(&den);
        tmp.mult(&turb);
        tmp2.add(&tmp);

        tmp2.mult(&self.volume);

        k = 5;
        for _i in 0..4 {
            r_out[k].add(&tmp2);
            k += 6;
        }

    }

    pub fn convective(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, fcdat : &[FaceData], fc_ar : &Vec<Face>) {
        let mut k = 0usize;
        for i in 0..4 {
            fc_ar[self.faces[i]].convective(&mut r_out[k..k+6], pre, &fcdat[i]);
            k += 6;
        }
    }

    pub fn viscous(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, fcdat : &[FaceData], fc_ar : &Vec<Face>) {
        let mut k = 0usize;
        for i in 0..4 {
            fc_ar[self.faces[i]].viscous(&mut r_out[k..k+6], pre, &fcdat[i]);
            k += 6;
        }
    }

    pub fn pressure(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, fcdat : &[FaceData], fc_ar : &Vec<Face>) {
        let mut k = 0usize;
        for i in 0..4 {
            fc_ar[self.faces[i]].pressure(&mut r_out[k..k+6], pre, &fcdat[i]);
            k += 6;
        }
    }

    pub fn heat_flux(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, fcdat : &[FaceData], fc_ar : &Vec<Face>) {
        let mut k = 0usize;
        for i in 0..4 {
            fc_ar[self.faces[i]].heat_flux(&mut r_out[k..k+6], pre, &fcdat[i]);
            k += 6;
        }
    }

    pub fn grav_r(&self, r_vec : &mut [DiffDoub1], den : &DiffDoub1, vel : &[DiffDoub1], ld : &[f64]) {
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        for i in 0..3 {
            tmp.set_val(ld[i]);
            tmp.mult(&den);
            tmp.mult(&self.volume);
            tmp.neg();
            r_vec[i+1].add(&tmp);
        }
        for i in 0..3 {
            tmp.set_val(ld[i]);
            tmp.mult(&vel[i]);
            tmp2.add(&tmp);
        }
        tmp2.mult(&den);
        tmp2.mult(&self.volume);
        tmp2.neg();
        r_vec[4].add(&tmp2);

    }

    pub fn cent_acc(&self, ld : &mut [f64], cent : &[DiffDoub1], this_ld : &Load) {
        ld[0] = cent[0].val - this_ld.center[0];
        ld[1] = cent[1].val - this_ld.center[1];
        ld[2] = cent[2].val - this_ld.center[2];

        let ax = & this_ld.axis;
        let dp = ld[0]*ax[0] + ld[1]*ax[1] + ld[2]*ax[2];

        let ang2 = this_ld.angular_vel*this_ld.angular_vel;
        for i in 0..3 {
            ld[i] -= dp*ax[i];
            ld[i] *= ang2;
        }
    }

    pub fn load(&self, r_out : &mut Vec<DiffDoub1>, pre : &EqnPrereq, ld_ar : &Vec<Load>, time : f64) {
        if self.loads.is_empty() {
            return;
        }

        let mut r_vec = [DiffDoub1::new(); 6];
        let mut this_ld : &Load;
        let mut ld = [0f64; 3];
        let mut den = DiffDoub1::new();
        let mut vel = [DiffDoub1::new(); 3];
        let mut cent = [DiffDoub1::new(); 3];

        let mut k : usize;
        for i in 0..4 {
            den.add(&pre.fl_den[i]);
            k = i;
            for j in 0..3 {
                vel[j].add(&pre.fl_vel[k]);
                cent[j].add(&pre.def_coord[k]);
                k += 4;
            }
        }
        let mut tmp = DiffDoub1::new();
        tmp.set_val(0.25);
        den.mult(&tmp);
        for j in 0..3 {
            vel[j].mult(&tmp);
            cent[j].mult(&tmp);
        }

        for ldi in self.loads.iter() {
            this_ld = &ld_ar[*ldi];
            this_ld.get_current_ld(&mut ld, time);
            match this_ld.this_type.s.as_str() {
                "gravitational" => self.grav_r(&mut r_vec, &den, &vel, &ld),
                
                "centrifugal" => {self.cent_acc(&mut ld, &cent, this_ld);
                                  self.grav_r(&mut r_vec, &den, &vel, &ld);},
                
                "volHeatGen" => {tmp.set_val(-ld[0]);
                                 tmp.mult(&self.volume);
                                 r_vec[4].add(&tmp);},
                &_ => (),
            }
        }

        k = 0;
        for _i in 0..4 {
            for j in 0..6 {
                r_out[k].add(&r_vec[j]);
                k += 1;
            }
        }

    }

    pub fn put_to_glob_r(&self, glob_r : &mut Vec<DiffDoub1>, cell_r : &mut Vec<DiffDoub1>, nodes : &Vec<Node>, pre : &EqnPrereq, fcdat : &mut [FaceData], fc_ar : &Vec<Face>, cdat : &Vec<CellData>, ld_ar : &Vec<Load>, time : f64) {
        for i in 0..24 {
            cell_r[i].set_val(0.0);
        }

        for i in 0..4 {
            fc_ar[self.faces[i]].update_area_norm(&mut fcdat[i], nodes);
            fc_ar[self.faces[i]].get_flux_vars(&mut fcdat[i], pre, fc_ar, cdat);
        }

        self.unsteady(cell_r, pre);
        self.convective(cell_r, pre, fcdat, fc_ar);
        self.viscous(cell_r, pre, fcdat, fc_ar);
        self.pressure(cell_r, pre, fcdat, fc_ar);
        self.heat_flux(cell_r, pre, fcdat, fc_ar);
        self.turbulence(cell_r, pre, cdat);
        self.load(cell_r, pre, ld_ar, time);

        let mut this_nd : &Node;
        let mut k = 0usize;
        let mut gind : usize;
        for i in 0..4 {
            this_nd = &nodes[self.nodes[i]];
            gind = 6*this_nd.sorted_rank;
            for _j in 0..6 {
                glob_r[gind].add(&cell_r[k]);
                gind += 1;
                k += 1;
            }
        }
    }

    pub fn put_to_glob_mat(&self, glob_mat : &mut SparseMat, cell_r : &mut Vec<DiffDoub1>, col : usize, nodes : &Vec<Node>, pre : &EqnPrereq, fcdat : &mut [FaceData], fc_ar : &Vec<Face>, cdat : &Vec<CellData>, ld_ar : &Vec<Load>, time : f64) {
        for i in 0..24 {
            cell_r[i].set_val(0.0);
        }

        for i in 0..4 {
            fc_ar[self.faces[i]].update_area_norm(&mut fcdat[i], nodes);
            fc_ar[self.faces[i]].get_flux_vars(&mut fcdat[i], pre, fc_ar, cdat);
        }

        self.unsteady(cell_r, pre);
        self.convective(cell_r, pre, fcdat, fc_ar);
        self.viscous(cell_r, pre, fcdat, fc_ar);
        self.pressure(cell_r, pre, fcdat, fc_ar);
        self.heat_flux(cell_r, pre, fcdat, fc_ar);
        self.turbulence(cell_r, pre, cdat);
        self.load(cell_r, pre, ld_ar, time);

        let mut this_nd : &Node;
        let mut k = 0usize;
        let mut gind : usize;
        for i in 0..4 {
            this_nd = &nodes[self.nodes[i]];
            gind = 6*this_nd.sorted_rank;
            for _j in 0..6 {
                glob_mat.add_entry(gind, col, cell_r[k].dval);
                gind += 1;
                k += 1;
            }
        }
    }

    pub fn put_unsteady_to_mat(&self, glob_mat : &mut SparseMat, cell_r : &mut Vec<DiffDoub1>, col : usize, mul_fact : f64, nodes : &Vec<Node>, pre : &EqnPrereq) {
        for i in 0..24 {
            cell_r[i].set_val(0.0);
        }

        self.unsteady(cell_r, pre);

        let mut this_nd : &Node;
        let mut k = 0usize;
        let mut gind : usize;
        for i in 0..4 {
            this_nd = &nodes[self.nodes[i]];
            gind = 6*this_nd.sorted_rank;
            for _j in 0..6 {
                glob_mat.add_entry(gind, col, mul_fact*cell_r[k].dval);
                gind += 1;
                k += 1;
            }
        }
    }

}