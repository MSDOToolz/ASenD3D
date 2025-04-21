use crate::face::*;
use crate::fmath::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::element::{DiffDoub0FlPrereq, DiffDoub1FlPrereq};
use crate::node::*;
use crate::design_var::*;
use crate::matrix_functions::*;


impl Face {

    pub fn initialize_int_pts(&mut self) {
        if self.num_nds == 3 {
            self.int_pts = vec![R_1O6, R_1O6, R_2O3, R_1O6, R_1O6, R_2O3];
            self.num_ip = 3;
        }
        else if self.num_nds == 4 {
            self.int_pts = vec![-R_1ORT3, -R_1ORT3, R_1ORT3, -R_1ORT3, -R_1ORT3, R_1ORT3, R_1ORT3, R_1ORT3];
            self.num_ip = 4;
        }
    }

    pub fn set_node(&mut self, place : usize, loc_nd : usize, glob_nd : usize) {
        self.loc_nodes[place] = loc_nd;
        self.glob_nodes[place] = glob_nd;
        return;
    }

    pub fn sorted_nodes(&mut self, srt_nds : &mut [usize]) {
        let i3 : usize;
        let mut i4 : usize;
        let mut swap : usize;
        for i1 in 0..self.num_nds {
            srt_nds[i1] = self.glob_nodes[i1];
        }
        i3 = self.num_nds - 1;
        for _i1 in 0..i3 {
            for i2 in 0..i3 {
                i4 = i2 + 1;
                if srt_nds[i4] < srt_nds[i2] {
                    swap = srt_nds[i2];
                    srt_nds[i2] = srt_nds[i4];
                    srt_nds[i4] = swap;
                }
            }
        }
        return;
    }

    pub fn get_low_nd(&mut self) -> usize {
        let mut low_nd : usize =  self.glob_nodes[0];
        for i1 in 1..self.num_nds {
            if self.glob_nodes[i1] < low_nd {
                low_nd = self.glob_nodes[i1];
            }
        }
        return  low_nd;
    }

    pub fn get_centroid(&self, cent : &mut [f64], nd_ar : &Vec<Node>) {
        for j in 0..3 {
            cent[j] = 0f64;
        }

        let mut crd : &[f64];
        for i in 0..self.num_nds {
            crd = &nd_ar[self.glob_nodes[i]].coord;
            for j in 0..3 {
                cent[j] += crd[j];
            }
        }
        
        let nn_inv = 1.0f64/(self.num_nds as f64);
        for j in 0..3 {
            cent[j] *= nn_inv;
        }
    }

    pub fn spt_in_face(&self, pt : &[f64]) -> bool {
        // pt = s coordinate vector from centroid of element
        if self.num_nds == 4 || self.num_nds == 8 {
            if fabs(pt[0]) > 1.0 {
                return false;
            }
            if fabs(pt[1]) > 1.0 {
                return false;
            }
            return true;
        }
        if pt[0] < -R_1O3 {
            return false;
        }
        if pt[1] < -R_1O3 {
            return false;
        }
        if pt[0] + pt[1] > R_1O3 {
            return false;
        }
        return true;
    }

    pub fn get_proj_dist(&self, s_crd : &mut [f64], pt : &[f64], nd_ar : &Vec<Node>) -> f64 {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut cent = [0f64; 3];
        let mut dx_ds1 = [0f64; 3];
        let mut dx_ds2 = [0f64; 3];
        let mut dist = 0f64;

        self.get_centroid(&mut cent, nd_ar);

        if self.num_nds == 4 || self.num_nds == 8 {
            let crd1 = &nd_ar[self.glob_nodes[1]].coord;
            let crd2 = &nd_ar[self.glob_nodes[2]].coord;
            let crd3 = &nd_ar[self.glob_nodes[3]].coord;
            for i in 0..3 {
                dx_ds1[i] = 0.5f64*(crd1[i] + crd2[i]) - cent[i];
                dx_ds2[i] = 0.5f64*(crd2[i] + crd3[i]) - cent[i];
            }
        }
        else {
            let crd0 = &nd_ar[self.glob_nodes[0]].coord;
            let crd1 = &nd_ar[self.glob_nodes[1]].coord;
            let crd2 = &nd_ar[self.glob_nodes[2]].coord;
            for i in 0..3 {
                dx_ds1[i] = crd1[i] - crd0[i];
                dx_ds2[i] = crd2[i] - crd0[i];
            }
        }
        
        let mut mat = [dx_ds1[0], dx_ds2[0], dx_ds1[1], dx_ds2[1], dx_ds1[2], dx_ds2[2]];
        let mut bvec = [pt[0] - cent[0], pt[1] - cent[1], pt[2] - cent[2]];
        let mut soln = [0.0f64; 2];

        q_rfactor_ar(&mut mat,2, 0, 2, 0, 1, 0);
        solveq_rx_eqb_ar(&mut soln, &mut mat, &mut bvec, 2, 0, 2, 0, 1, 0);

        if !self.spt_in_face(&soln) {
            let mut mul_fac = 0.5f64;
            while mul_fac < 0.98 {
                while !self.spt_in_face(&soln) {
                    soln[0] *= mul_fac;
                    soln[1] *= mul_fac;
                }
                soln[0] /= mul_fac;
                soln[1] /= mul_fac;
                mul_fac = mul_fac.sqrt();
            }
        }
        
        if self.num_nds == 4 || self.num_nds ==8 {
            s_crd[0] = soln[0];
            s_crd[1] = soln[1];
        }
        else {
            s_crd[0] = soln[0] + R_1O3;
            s_crd[1] = soln[1] + R_1O3;
        }
        
        
        let mut proj : f64;
        let mut dp = 0.0f64;
        for i in 0..3 {
            proj = (cent[i] + soln[0]*dx_ds1[i] + soln[1]*dx_ds2[i]) - pt[i];
            dp += proj*proj;
        }

        dp.sqrt()
    }

    //dup1

    pub fn get_area_normal_dfd0(&self, area : &mut DiffDoub0, norm : &mut [DiffDoub0], nd_crd : & Vec<DiffDoub0>, el_nn : usize) {
        let mut v1 = [DiffDoub0::new(); 3];
        let mut v2 = [DiffDoub0::new(); 3];
        // let mut tmp_v = [DiffDoub0::new(); 3];
        let mut tmp = DiffDoub0::new();
        let mut shft : usize;
        
        if self.num_nds == 4 {
            // nd_ar[self.glob_nodes[2]].get_def_crd_dfd0(&mut v1);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd0(&mut tmp_v);
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd0(&nd_crd[2 + shft]);
                v1[i].sub(&nd_crd[0 + shft]);
                v2[i].set_val_dfd0(&nd_crd[3 + shft]);
                v2[i].sub(&nd_crd[1 + shft]);   
            }
            // v1[0].sub(& tmp_v[0]);
            // v1[1].sub(& tmp_v[1]);
            // v1[2].sub(& tmp_v[2]);
            // nd_ar[self.glob_nodes[3]].get_def_crd_dfd0(&mut v2);
            // nd_ar[self.glob_nodes[1]].get_def_crd_dfd0(&mut tmp_v);
            // v2[0].sub(& tmp_v[0]);
            // v2[1].sub(& tmp_v[1]);
            // v2[2].sub(& tmp_v[2]);
        }
        else {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd0(&nd_crd[1 + shft]);
                v1[i].sub(&nd_crd[0 + shft]);
                v2[i].set_val_dfd0(&nd_crd[2 + shft]);
                v2[i].sub(&nd_crd[0 + shft]);
            }
            // nd_ar[self.glob_nodes[1]].get_def_crd_dfd0(&mut v1);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd0(&mut tmp_v);
            // v1[0].sub(& tmp_v[0]);
            // v1[1].sub(& tmp_v[1]);
            // v1[2].sub(& tmp_v[2]);
            // nd_ar[self.glob_nodes[2]].get_def_crd_dfd0(&mut v2);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd0(&mut tmp_v);
            // v2[0].sub(& tmp_v[0]);
            // v2[1].sub(& tmp_v[1]);
            // v2[2].sub(& tmp_v[2]);
        }
        
        cross_prod_dfd0(norm, &mut  v1, &mut  v2);
        area.set_val_dfd0(& norm[0]);
        area.sqr();
        tmp.set_val_dfd0(& norm[1]);
        tmp.sqr();
        area.add(& tmp);
        tmp.set_val_dfd0(& norm[2]);
        tmp.sqr();
        area.add(& tmp);
        area.sqt();
        
        tmp.set_val(1.0);
        tmp.dvd(& area);
        norm[0].mult(& tmp);
        norm[1].mult(& tmp);
        norm[2].mult(& tmp);
        
        tmp.set_val(0.5);
        area.mult(& tmp);
        
        return;
    }

    pub fn get_basis_dfd0(&self, n_vec : &mut [DiffDoub0], spt : &[f64]) {
        if self.num_nds == 3 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]);
            n_vec[1].set_val(spt[0]);
            n_vec[2].set_val(spt[1]);
        }
        else if self.num_nds == 4 {
            n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));  
            n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
            n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
            n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
        }
    }

    pub fn get_rf_convective_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq) {
        let mut rho = DiffDoub0::new();
        let mut vel = [DiffDoub0::new(); 3];
        let mut fl_vel = [DiffDoub0::new(); 3];
        let mut v_rel = [DiffDoub0::new(); 3];
        let mut temp = DiffDoub0::new();
        let mut turb_e = DiffDoub0::new();
        let mut enth = DiffDoub0::new();

        let mut rho_v_n = DiffDoub0::new();
        let mut r_sum = [DiffDoub0::new(); 6];

        let mut area = DiffDoub0::new();
        let mut norm = [DiffDoub0::new(); 3];
        self.get_area_normal_dfd0(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub0::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);
        
        let mut n_vec = [DiffDoub0::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();

        for _ip in 0..self.num_ip {
            self.get_basis_dfd0(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd0(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd0(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
                td.set_val_dfd0(&pre.fl_turb_e[nd]);
                td.mult(&n_vec[j]);
                turb_e.add(&td);
            }
            for k in 0..3 {
                vel[k].set_val(0.0);
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd0(&pre.glob_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    vel[k].add(&td);
                    td.set_val_dfd0(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
                v_rel[k].set_val_dfd0(&fl_vel[k]);
                v_rel[k].sub(&vel[k]);
            }

            enth.set_val_dfd0(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);

            rho_v_n.set_val(0.0);
            for j in 0..3 {
                td.set_val_dfd0(&v_rel[j]);
                td.mult(&norm[j]);
                rho_v_n.add(&td);
            }
            rho_v_n.mult(&rho);
            rho_v_n.mult(&da);
            
            // mass
            r_sum[0].add(&rho_v_n);

            // momentum
            for j in 0..3 {
                td.set_val_dfd0(&rho_v_n);
                td.mult(&fl_vel[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&fl_vel[j]);
                td2.sqr();
                td.add(&td2);
            }
            td2.set_val(0.5);
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_v_n);
            r_sum[4].add(&td);

            // turbulence
            td.set_val_dfd0(&turb_e);
            td.mult(&rho_v_n);
            r_sum[5].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_viscous_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub0::new();
        let mut temp = DiffDoub0::new();
        let mut fl_vel = [DiffDoub0::new(); 3];
        let mut turb_e = DiffDoub0::new();
        let mut vis = DiffDoub0::new();

        let mut r_sum = [DiffDoub0::new(); 6];

        let mut area = DiffDoub0::new();
        let mut norm = [DiffDoub0::new(); 3];
        self.get_area_normal_dfd0(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub0::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);

        let mut n_vec = [DiffDoub0::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        let mut oth_hst : usize;

        let mut grad_v = [DiffDoub0::new(); 9];
        if self.twin_id < MAX_INT {
            td.set_val(0.5);
            for i in 0..9 {
                grad_v[i].set_val_dfd0(&ec_vec[self.host_el].v_grad_dfd0[i]);
                oth_hst = fc_vec[self.twin_id].host_el;
                grad_v[i].add(&ec_vec[oth_hst].v_grad_dfd0[i]);
                grad_v[i].mult(&td);
            }
        }
        else {
            for i in 0..9 {
                grad_v[i].set_val_dfd0(&ec_vec[self.host_el].v_grad_dfd0[i]);
            }
        }

        let mut gv_n_da = [DiffDoub0::new(); 3];
        mat_mul_ar_dfd0(&mut gv_n_da, &grad_v, &norm, 3, 3, 1);
        for i in 0..3 {
            gv_n_da[i].mult(&da);
        }

        for _ip in 0..self.num_ip {
            self.get_basis_dfd0(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd0(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd0(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
                td.set_val_dfd0(&pre.fl_turb_e[nd]);
                td.mult(&n_vec[j]);
                turb_e.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd0(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
            }

            // calculate viscosity
            vis.set_val_dfd0(&pre.ref_temp);
            vis.add(&temp);
            vis.mult(&rho);
            vis.mult(&pre.den_vis_coef);
            td.set_val_dfd0(&pre.turb_vis_coef);
            td.mult(&turb_e);
            vis.add(&td);
            td.set_val_dfd0(&pre.temp_vis_coef);
            td.mult(&temp);
            vis.add(&td);
            

            // momentum
            for j in 0..3 {
                td.set_val_dfd0(&vis);
                td.mult(&gv_n_da[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&fl_vel[j]);
                td2.mult(&gv_n_da[j]);
                td.add(&td2);
            }
            td.mult(&vis);
            r_sum[4].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }   
    }

    pub fn get_rf_pressure_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub0::new();
        let mut temp = DiffDoub0::new();
        let mut fl_vel = [DiffDoub0::new(); 3];
        let mut pres = DiffDoub0::new();

        let mut r_sum = [DiffDoub0::new(); 6];

        let mut area = DiffDoub0::new();
        let mut norm = [DiffDoub0::new(); 3];
        self.get_area_normal_dfd0(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub0::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);
        let mut n_da = [DiffDoub0::new(); 3];
        for i in 0..3 {
            n_da[i].set_val_dfd0(&norm[i]);
            n_da[i].mult(&da);
        }

        let mut n_vec = [DiffDoub0::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        let mut oth_hst : usize;

        let mut grad_t = [DiffDoub0::new(); 3];
        if self.twin_id < MAX_INT {
            td.set_val(0.5);
            for i in 0..3 {
                grad_t[i].set_val_dfd0(&ec_vec[self.host_el].t_grad_dfd0[i]);
                oth_hst = fc_vec[self.twin_id].host_el;
                grad_t[i].add(&ec_vec[oth_hst].t_grad_dfd0[i]);
                grad_t[i].mult(&td);
            }
        }
        else {
            for i in 0..9 {
                grad_t[i].set_val_dfd0(&ec_vec[self.host_el].t_grad_dfd0[i]);
            }
        }

        for _ip in 0..self.num_ip {
            self.get_basis_dfd0(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd0(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd0(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd0(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
            }

            // calculate pressure
            if pre.compressible {
                pres.set_val_dfd0(&temp);
                pres.add(&pre.ref_temp);
                pres.mult(&pre.i_gconst);
                pres.mult(&rho);
            }
            else {
                pres.set_val_dfd0(&pre.ref_pres);
                td.set_val_dfd0(&rho);
                td.sub(&pre.ref_den);
                td.mult(&pre.den_pres_coef);
                pres.add(&td);
                td.set_val_dfd0(&pre.temp_pres_coef);
                td.mult(&temp);
                pres.add(&td);
            }
            

            // momentum
            for j in 0..3 {
                td.set_val_dfd0(&pres);
                td.mult(&n_da[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&fl_vel[j]);
                td2.mult(&n_da[j]);
                td.add(&td2);
            }
            td.mult(&pres);
            r_sum[4].add(&td);
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&grad_t[j]);
                td2.mult(&n_da[j]);
                td.add(&td2);
            }
            td.mult(&pre.therm_cond);
            td.neg();
            r_sum[4].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }          
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn get_area_normal_dfd1(&self, area : &mut DiffDoub1, norm : &mut [DiffDoub1], nd_crd : & Vec<DiffDoub1>, el_nn : usize) {
        let mut v1 = [DiffDoub1::new(); 3];
        let mut v2 = [DiffDoub1::new(); 3];
        // let mut tmp_v = [DiffDoub1::new(); 3];
        let mut tmp = DiffDoub1::new();
        let mut shft : usize;
        
        if self.num_nds == 4 {
            // nd_ar[self.glob_nodes[2]].get_def_crd_dfd1(&mut v1);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd1(&mut tmp_v);
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd1(&nd_crd[2 + shft]);
                v1[i].sub(&nd_crd[0 + shft]);
                v2[i].set_val_dfd1(&nd_crd[3 + shft]);
                v2[i].sub(&nd_crd[1 + shft]);   
            }
            // v1[0].sub(& tmp_v[0]);
            // v1[1].sub(& tmp_v[1]);
            // v1[2].sub(& tmp_v[2]);
            // nd_ar[self.glob_nodes[3]].get_def_crd_dfd1(&mut v2);
            // nd_ar[self.glob_nodes[1]].get_def_crd_dfd1(&mut tmp_v);
            // v2[0].sub(& tmp_v[0]);
            // v2[1].sub(& tmp_v[1]);
            // v2[2].sub(& tmp_v[2]);
        }
        else {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd1(&nd_crd[1 + shft]);
                v1[i].sub(&nd_crd[0 + shft]);
                v2[i].set_val_dfd1(&nd_crd[2 + shft]);
                v2[i].sub(&nd_crd[0 + shft]);
            }
            // nd_ar[self.glob_nodes[1]].get_def_crd_dfd1(&mut v1);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd1(&mut tmp_v);
            // v1[0].sub(& tmp_v[0]);
            // v1[1].sub(& tmp_v[1]);
            // v1[2].sub(& tmp_v[2]);
            // nd_ar[self.glob_nodes[2]].get_def_crd_dfd1(&mut v2);
            // nd_ar[self.glob_nodes[0]].get_def_crd_dfd1(&mut tmp_v);
            // v2[0].sub(& tmp_v[0]);
            // v2[1].sub(& tmp_v[1]);
            // v2[2].sub(& tmp_v[2]);
        }
        
        cross_prod_dfd1(norm, &mut  v1, &mut  v2);
        area.set_val_dfd1(& norm[0]);
        area.sqr();
        tmp.set_val_dfd1(& norm[1]);
        tmp.sqr();
        area.add(& tmp);
        tmp.set_val_dfd1(& norm[2]);
        tmp.sqr();
        area.add(& tmp);
        area.sqt();
        
        tmp.set_val(1.0);
        tmp.dvd(& area);
        norm[0].mult(& tmp);
        norm[1].mult(& tmp);
        norm[2].mult(& tmp);
        
        tmp.set_val(0.5);
        area.mult(& tmp);
        
        return;
    }

    pub fn get_basis_dfd1(&self, n_vec : &mut [DiffDoub1], spt : &[f64]) {
        if self.num_nds == 3 {
            n_vec[0].set_val(1.0-spt[0]-spt[1]);
            n_vec[1].set_val(spt[0]);
            n_vec[2].set_val(spt[1]);
        }
        else if self.num_nds == 4 {
            n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));  
            n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
            n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
            n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
        }
    }

    pub fn get_rf_convective_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq) {
        let mut rho = DiffDoub1::new();
        let mut vel = [DiffDoub1::new(); 3];
        let mut fl_vel = [DiffDoub1::new(); 3];
        let mut v_rel = [DiffDoub1::new(); 3];
        let mut temp = DiffDoub1::new();
        let mut turb_e = DiffDoub1::new();
        let mut enth = DiffDoub1::new();

        let mut rho_v_n = DiffDoub1::new();
        let mut r_sum = [DiffDoub1::new(); 6];

        let mut area = DiffDoub1::new();
        let mut norm = [DiffDoub1::new(); 3];
        self.get_area_normal_dfd1(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub1::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);
        
        let mut n_vec = [DiffDoub1::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();

        for _ip in 0..self.num_ip {
            self.get_basis_dfd1(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd1(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd1(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
                td.set_val_dfd1(&pre.fl_turb_e[nd]);
                td.mult(&n_vec[j]);
                turb_e.add(&td);
            }
            for k in 0..3 {
                vel[k].set_val(0.0);
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd1(&pre.glob_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    vel[k].add(&td);
                    td.set_val_dfd1(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
                v_rel[k].set_val_dfd1(&fl_vel[k]);
                v_rel[k].sub(&vel[k]);
            }

            enth.set_val_dfd1(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);

            rho_v_n.set_val(0.0);
            for j in 0..3 {
                td.set_val_dfd1(&v_rel[j]);
                td.mult(&norm[j]);
                rho_v_n.add(&td);
            }
            rho_v_n.mult(&rho);
            rho_v_n.mult(&da);
            
            // mass
            r_sum[0].add(&rho_v_n);

            // momentum
            for j in 0..3 {
                td.set_val_dfd1(&rho_v_n);
                td.mult(&fl_vel[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&fl_vel[j]);
                td2.sqr();
                td.add(&td2);
            }
            td2.set_val(0.5);
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_v_n);
            r_sum[4].add(&td);

            // turbulence
            td.set_val_dfd1(&turb_e);
            td.mult(&rho_v_n);
            r_sum[5].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_viscous_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub1::new();
        let mut temp = DiffDoub1::new();
        let mut fl_vel = [DiffDoub1::new(); 3];
        let mut turb_e = DiffDoub1::new();
        let mut vis = DiffDoub1::new();

        let mut r_sum = [DiffDoub1::new(); 6];

        let mut area = DiffDoub1::new();
        let mut norm = [DiffDoub1::new(); 3];
        self.get_area_normal_dfd1(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub1::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);

        let mut n_vec = [DiffDoub1::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        let mut oth_hst : usize;

        let mut grad_v = [DiffDoub1::new(); 9];
        if self.twin_id < MAX_INT {
            td.set_val(0.5);
            for i in 0..9 {
                grad_v[i].set_val_dfd1(&ec_vec[self.host_el].v_grad_dfd1[i]);
                oth_hst = fc_vec[self.twin_id].host_el;
                grad_v[i].add(&ec_vec[oth_hst].v_grad_dfd1[i]);
                grad_v[i].mult(&td);
            }
        }
        else {
            for i in 0..9 {
                grad_v[i].set_val_dfd1(&ec_vec[self.host_el].v_grad_dfd1[i]);
            }
        }

        let mut gv_n_da = [DiffDoub1::new(); 3];
        mat_mul_ar_dfd1(&mut gv_n_da, &grad_v, &norm, 3, 3, 1);
        for i in 0..3 {
            gv_n_da[i].mult(&da);
        }

        for _ip in 0..self.num_ip {
            self.get_basis_dfd1(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd1(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd1(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
                td.set_val_dfd1(&pre.fl_turb_e[nd]);
                td.mult(&n_vec[j]);
                turb_e.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd1(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
            }

            // calculate viscosity
            vis.set_val_dfd1(&pre.ref_temp);
            vis.add(&temp);
            vis.mult(&rho);
            vis.mult(&pre.den_vis_coef);
            td.set_val_dfd1(&pre.turb_vis_coef);
            td.mult(&turb_e);
            vis.add(&td);
            td.set_val_dfd1(&pre.temp_vis_coef);
            td.mult(&temp);
            vis.add(&td);
            

            // momentum
            for j in 0..3 {
                td.set_val_dfd1(&vis);
                td.mult(&gv_n_da[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&fl_vel[j]);
                td2.mult(&gv_n_da[j]);
                td.add(&td2);
            }
            td.mult(&vis);
            r_sum[4].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }   
    }

    pub fn get_rf_pressure_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub1::new();
        let mut temp = DiffDoub1::new();
        let mut fl_vel = [DiffDoub1::new(); 3];
        let mut pres = DiffDoub1::new();

        let mut r_sum = [DiffDoub1::new(); 6];

        let mut area = DiffDoub1::new();
        let mut norm = [DiffDoub1::new(); 3];
        self.get_area_normal_dfd1(&mut area, &mut norm, &pre.glob_nds, pre.num_nds);
        let mut da = DiffDoub1::new();
        da.set_val(1.0f64/(self.num_ip as f64));
        da.mult(&area);
        let mut n_da = [DiffDoub1::new(); 3];
        for i in 0..3 {
            n_da[i].set_val_dfd1(&norm[i]);
            n_da[i].mult(&da);
        }

        let mut n_vec = [DiffDoub1::new(); 8];
        let mut si = 0usize;
        let mut nd : usize;
        let mut shft : usize;
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        let mut oth_hst : usize;

        let mut grad_t = [DiffDoub1::new(); 3];
        if self.twin_id < MAX_INT {
            td.set_val(0.5);
            for i in 0..3 {
                grad_t[i].set_val_dfd1(&ec_vec[self.host_el].t_grad_dfd1[i]);
                oth_hst = fc_vec[self.twin_id].host_el;
                grad_t[i].add(&ec_vec[oth_hst].t_grad_dfd1[i]);
                grad_t[i].mult(&td);
            }
        }
        else {
            for i in 0..9 {
                grad_t[i].set_val_dfd1(&ec_vec[self.host_el].t_grad_dfd1[i]);
            }
        }

        for _ip in 0..self.num_ip {
            self.get_basis_dfd1(&mut n_vec,&self.int_pts[si..si+2]);
            rho.set_val(0.0);
            temp.set_val(0.0);
            for j in 0..self.num_nds {
                nd = self.loc_nodes[j];
                td.set_val_dfd1(&pre.fl_den[nd]);
                td.mult(&n_vec[j]);
                rho.add(&td);
                td.set_val_dfd1(&pre.fl_temp[nd]);
                td.mult(&n_vec[j]);
                temp.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                shft = k*pre.num_nds;
                for j in 0..self.num_nds {
                    nd = self.loc_nodes[j];
                    td.set_val_dfd1(&pre.fl_vel[nd + shft]);
                    td.mult(&n_vec[j]);
                    fl_vel[k].add(&td);
                }
            }

            // calculate pressure
            if pre.compressible {
                pres.set_val_dfd1(&temp);
                pres.add(&pre.ref_temp);
                pres.mult(&pre.i_gconst);
                pres.mult(&rho);
            }
            else {
                pres.set_val_dfd1(&pre.ref_pres);
                td.set_val_dfd1(&rho);
                td.sub(&pre.ref_den);
                td.mult(&pre.den_pres_coef);
                pres.add(&td);
                td.set_val_dfd1(&pre.temp_pres_coef);
                td.mult(&temp);
                pres.add(&td);
            }
            

            // momentum
            for j in 0..3 {
                td.set_val_dfd1(&pres);
                td.mult(&n_da[j]);
                r_sum[j+1].add(&td);
            }

            // energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&fl_vel[j]);
                td2.mult(&n_da[j]);
                td.add(&td2);
            }
            td.mult(&pres);
            r_sum[4].add(&td);
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&grad_t[j]);
                td2.mult(&n_da[j]);
                td.add(&td2);
            }
            td.mult(&pre.therm_cond);
            td.neg();
            r_sum[4].add(&td);

            si += 2;
        }

        let mut ri = 0usize;
        for _i in 0..pre.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }          
    }

    //end dup
 
//end skip 
 
 
 
 
 
}

impl FacePtList {
    pub fn add_face(&mut self, new_i : usize) {
        self.fc_list.push_back(new_i);
        return;
    }

    pub fn add_if_absent(&mut self, new_i : usize, glob_faces : &mut Vec<Face>) -> bool {
        let new_num_nds : usize;
        let mut new_srtd : [usize; 8] = [0usize; 8];
        let mut this_num_nds : usize;
        let mut this_srtd : [usize; 8] = [0usize; 8];
        let mut all_match : bool;
        
        new_num_nds = glob_faces[new_i].num_nds;
        glob_faces[new_i].sorted_nodes(&mut new_srtd);
        for fi in self.fc_list.iter_mut() {
            this_num_nds = glob_faces[*fi].num_nds;
            if this_num_nds == new_num_nds {
                glob_faces[*fi].sorted_nodes(&mut this_srtd);
                all_match = true;
                for i1 in 0..this_num_nds {
                    if this_srtd[i1] != new_srtd[i1] {
                        all_match = false;
                    }
                }
                if all_match {
                    glob_faces[new_i].on_surf = false;
                    glob_faces[*fi].on_surf = false;
                    glob_faces[new_i].twin_id = *fi;
                    glob_faces[*fi].twin_id = new_i;
                    return  false;
                }
            }
        }
        
        self.fc_list.push_back(new_i);
        return  true;
    }

}


