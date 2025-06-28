use crate::element::*;
use crate::node::*;
use crate::diff_doub::*;
use crate::face::*;
use crate::list_ent::*;
use crate::scratch::*;

use std::collections::linked_list::IterMut;

impl Element {
    //dup1

    pub fn get_rf_unsteady_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq) {
        let mut rho = DiffDoub0::new();
        let mut rho_dot = DiffDoub0::new();
        let mut fl_vel = [DiffDoub0::new(); 3];
        let mut vel_dot = [DiffDoub0::new(); 3];
        let mut temp = DiffDoub0::new();
        let mut t_dot = DiffDoub0::new();
        let mut turb_e = DiffDoub0::new();
        let mut turb_e_dot = DiffDoub0::new();
        let mut enth = DiffDoub0::new();
        let mut enth_dot = DiffDoub0::new();

        let mut si = 0usize;
        let mut vi : usize;
        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut r_sum = [DiffDoub0::new(); 6];
        let mut ipw_j = DiffDoub0::new();
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        let mut td3 = DiffDoub0::new();

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            rho_dot.set_val(0.0);
            temp.set_val(0.0);
            t_dot.set_val(0.0);
            turb_e.set_val(0.0);
            turb_e_dot.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den_dot[j]);
                rho_dot.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_temp[j]);
                temp.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_tdot[j]);
                t_dot.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_edot[j]);
                turb_e_dot.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                vel_dot[k].set_val(0.0);
                vi = self.num_nds*k;
                for j in 0..self.num_nds {
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel[vi]);
                    fl_vel[k].add(&td);
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel_dot[vi]);
                    vel_dot[k].add(&td);
                    vi += 1;
                }
            }
            enth.set_val_dfd0(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);
            enth_dot.set_val_dfd0(&pre.spec_heat);
            enth_dot.mult(&t_dot);

            //mass
            td.set_val_dfd0(&rho_dot);
            td.mult(&ipw_j);
            r_sum[0].add(&td);

            //momentum
            for j in 0..3 {
                td.set_val_dfd0(&rho_dot);
                td.mult(&fl_vel[j]);
                td2.set_val_dfd0(&rho);
                td2.mult(&vel_dot[j]);
                td.add(&td2);
                td.mult(&ipw_j);
                r_sum[j+1].add(&td);
            }
            
            //energy
            td.set_val(0.5);
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd0(&fl_vel[j]);
                td3.sqr();
                td2.add(&td3);
            }
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_dot); // td = rho_dot(0.5v*v + h + tau)
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd0(&fl_vel[j]);
                td3.mult(&vel_dot[j]);
                td2.add(&td3);
            }
            td2.add(&enth_dot);
            td2.add(&turb_e_dot);
            td2.mult(&rho);  //td2 = rho(v*v_dot + h_dot + tau_dot)
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[4].add(&td);

            //turbulence
            td.set_val_dfd0(&rho_dot);
            td.mult(&turb_e);
            td2.set_val_dfd0(&rho);
            td2.mult(&turb_e_dot);
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[5].add(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_convective_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_convective_dfd0(r_vec,pre);
        }
    }

    pub fn get_rf_viscous_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub0::new();
        let mut turb_e = DiffDoub0::new();
        let mut grad_v = [DiffDoub0::new(); 9];
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();

        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ipw_j = DiffDoub0::new();
        let mut r_sum = DiffDoub0::new();
        let mut si = 0usize;
        let mut gvi : usize;
        let mut nvi : usize;
        let mut gni : usize;

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
            }
            gvi = 0;
            for j in 0..3 {
                for k in 0..3 {
                    grad_v[gvi].set_val(0.0);
                    nvi = j*self.num_nds;
                    gni = k;
                    for _m  in 0..self.num_nds {
                        //nvi = j*self.num_nds + m;
                        //gni = m*3 + k;
                        td.set_val_dfd0(&pre.fl_vel[nvi]);
                        td.mult(&d_ndx[gni]);
                        grad_v[gvi].add(&td);
                        nvi += 1;
                        gni += 3;
                    }
                    gvi += 1;
                }
            }
            td.set_val(0.0);
            for j in 0..9 {
                td2.set_val_dfd0(&grad_v[j]);
                td2.sqr();
                td.add(&td2);
            }
            td.sqt(); // td = mag(grad_v)
            
            td2.set_val_dfd0(&pre.grad_turb_coef);
            td2.mult(&rho);
            td2.mult(&td);
            td2.mult(&ipw_j);
            r_sum.sub(&td2);

            td2.set_val_dfd0(&pre.diss_turb_coef);
            td2.mult(&rho);
            td2.mult(&turb_e);
            td2.mult(&ipw_j);
            r_sum.add(&td2);

            si += 3;
        }

        let mut ri = 5usize;
        for j in 0..self.num_nds {
            r_vec[ri].add(&r_sum);
            ri += 6;
        }

        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_viscous_dfd0(r_vec, pre, fc_vec, ec_vec);
        }

    }

    pub fn get_rf_pressure_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &DiffDoub0FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_pressure_dfd0(r_vec,pre, fc_vec, ec_vec);
        }
    }

    pub fn get_rf_load_dfd0(&self, r_vec : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq) {
        let mut rho = DiffDoub0::new();
        let mut vel = [DiffDoub0::new(); 3];

        let mut n_vec = [DiffDoub0::new(); 11];
        let mut d_ndx = [DiffDoub0::new(); 33];
        let mut det_j = DiffDoub0::new();
        let mut ipw_j = DiffDoub0::new();
        let mut si = 0usize;
        let mut nvi : usize;
        let mut r_sum = [DiffDoub0::new(); 6];
        let mut td = DiffDoub0::new();
        let mut td2 = DiffDoub0::new();
        
        for ip in 0..self.num_ip {
            self.get_ip_data_dfd0(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd0(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
            }
            nvi = 0;
            for k in 0..3 {
                vel[k].set_val(0.0);
                for j in 0..self.num_nds {
                    //nvi = k*self.num_nds + j;
                    td.set_val_dfd0(&n_vec[j]);
                    td.mult(&pre.fl_vel[nvi]);
                    vel[k].add(&td);
                    nvi += 1;
                }
            }

            //momentum
            for j in 1..4 {
                td.set_val_dfd0(&pre.f_per_mass[j-1]);
                td.mult(&rho);
                td.mult(&ipw_j);
                r_sum[j].sub(&td);
            }

            //energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd0(&vel[j]);
                td2.mult(&pre.f_per_mass[j]);
                td.add(&td2);
            }
            td.mult(&rho);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            td.set_val_dfd0(&pre.hg_per_vol);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_dfd0(&self, r_glob : &mut Vec<DiffDoub0>, pre : &mut DiffDoub0FlPrereq, dynamic : bool, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub0Scr>) {
        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_dfd0"),
            Some(x) => &mut x.dat,
        };
        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd0(r_vec, pre);
        }
        self.get_rf_convective_dfd0( r_vec, pre, fc_vec);
        self.get_rf_viscous_dfd0(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_pressure_dfd0(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_load_dfd0(r_vec,pre);

        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                r_glob[rg+j].add(&r_vec[ri]);
                ri += 1;
            }
        }
    }

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn get_rf_unsteady_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq) {
        let mut rho = DiffDoub1::new();
        let mut rho_dot = DiffDoub1::new();
        let mut fl_vel = [DiffDoub1::new(); 3];
        let mut vel_dot = [DiffDoub1::new(); 3];
        let mut temp = DiffDoub1::new();
        let mut t_dot = DiffDoub1::new();
        let mut turb_e = DiffDoub1::new();
        let mut turb_e_dot = DiffDoub1::new();
        let mut enth = DiffDoub1::new();
        let mut enth_dot = DiffDoub1::new();

        let mut si = 0usize;
        let mut vi : usize;
        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut r_sum = [DiffDoub1::new(); 6];
        let mut ipw_j = DiffDoub1::new();
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        let mut td3 = DiffDoub1::new();

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            rho_dot.set_val(0.0);
            temp.set_val(0.0);
            t_dot.set_val(0.0);
            turb_e.set_val(0.0);
            turb_e_dot.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den_dot[j]);
                rho_dot.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_temp[j]);
                temp.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_tdot[j]);
                t_dot.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_edot[j]);
                turb_e_dot.add(&td);
            }
            for k in 0..3 {
                fl_vel[k].set_val(0.0);
                vel_dot[k].set_val(0.0);
                vi = self.num_nds*k;
                for j in 0..self.num_nds {
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel[vi]);
                    fl_vel[k].add(&td);
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel_dot[vi]);
                    vel_dot[k].add(&td);
                    vi += 1;
                }
            }
            enth.set_val_dfd1(&pre.spec_heat);
            enth.mult(&temp);
            enth.add(&pre.ref_enth);
            enth_dot.set_val_dfd1(&pre.spec_heat);
            enth_dot.mult(&t_dot);

            //mass
            td.set_val_dfd1(&rho_dot);
            td.mult(&ipw_j);
            r_sum[0].add(&td);

            //momentum
            for j in 0..3 {
                td.set_val_dfd1(&rho_dot);
                td.mult(&fl_vel[j]);
                td2.set_val_dfd1(&rho);
                td2.mult(&vel_dot[j]);
                td.add(&td2);
                td.mult(&ipw_j);
                r_sum[j+1].add(&td);
            }
            
            //energy
            td.set_val(0.5);
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd1(&fl_vel[j]);
                td3.sqr();
                td2.add(&td3);
            }
            td.mult(&td2);
            td.add(&enth);
            td.add(&turb_e);
            td.mult(&rho_dot); // td = rho_dot(0.5v*v + h + tau)
            td2.set_val(0.0);
            for j in 0..3 {
                td3.set_val_dfd1(&fl_vel[j]);
                td3.mult(&vel_dot[j]);
                td2.add(&td3);
            }
            td2.add(&enth_dot);
            td2.add(&turb_e_dot);
            td2.mult(&rho);  //td2 = rho(v*v_dot + h_dot + tau_dot)
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[4].add(&td);

            //turbulence
            td.set_val_dfd1(&rho_dot);
            td.mult(&turb_e);
            td2.set_val_dfd1(&rho);
            td2.mult(&turb_e_dot);
            td.add(&td2);
            td.mult(&ipw_j);
            r_sum[5].add(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_convective_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_convective_dfd1(r_vec,pre);
        }
    }

    pub fn get_rf_viscous_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        let mut rho = DiffDoub1::new();
        let mut turb_e = DiffDoub1::new();
        let mut grad_v = [DiffDoub1::new(); 9];
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();

        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ipw_j = DiffDoub1::new();
        let mut r_sum = DiffDoub1::new();
        let mut si = 0usize;
        let mut gvi : usize;
        let mut nvi : usize;
        let mut gni : usize;

        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            turb_e.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_turb_e[j]);
                turb_e.add(&td);
            }
            gvi = 0;
            for j in 0..3 {
                for k in 0..3 {
                    grad_v[gvi].set_val(0.0);
                    nvi = j*self.num_nds;
                    gni = k;
                    for _m  in 0..self.num_nds {
                        //nvi = j*self.num_nds + m;
                        //gni = m*3 + k;
                        td.set_val_dfd1(&pre.fl_vel[nvi]);
                        td.mult(&d_ndx[gni]);
                        grad_v[gvi].add(&td);
                        nvi += 1;
                        gni += 3;
                    }
                    gvi += 1;
                }
            }
            td.set_val(0.0);
            for j in 0..9 {
                td2.set_val_dfd1(&grad_v[j]);
                td2.sqr();
                td.add(&td2);
            }
            td.sqt(); // td = mag(grad_v)
            
            td2.set_val_dfd1(&pre.grad_turb_coef);
            td2.mult(&rho);
            td2.mult(&td);
            td2.mult(&ipw_j);
            r_sum.sub(&td2);

            td2.set_val_dfd1(&pre.diss_turb_coef);
            td2.mult(&rho);
            td2.mult(&turb_e);
            td2.mult(&ipw_j);
            r_sum.add(&td2);

            si += 3;
        }

        let mut ri = 5usize;
        for j in 0..self.num_nds {
            r_vec[ri].add(&r_sum);
            ri += 6;
        }

        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
        }

    }

    pub fn get_rf_pressure_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &DiffDoub1FlPrereq, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>) {
        for fi in self.faces.iter() {
            fc_vec[*fi].get_rf_pressure_dfd1(r_vec,pre, fc_vec, ec_vec);
        }
    }

    pub fn get_rf_load_dfd1(&self, r_vec : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq) {
        let mut rho = DiffDoub1::new();
        let mut vel = [DiffDoub1::new(); 3];

        let mut n_vec = [DiffDoub1::new(); 11];
        let mut d_ndx = [DiffDoub1::new(); 33];
        let mut det_j = DiffDoub1::new();
        let mut ipw_j = DiffDoub1::new();
        let mut si = 0usize;
        let mut nvi : usize;
        let mut r_sum = [DiffDoub1::new(); 6];
        let mut td = DiffDoub1::new();
        let mut td2 = DiffDoub1::new();
        
        for ip in 0..self.num_ip {
            self.get_ip_data_dfd1(&mut n_vec, &mut d_ndx, &mut det_j, &mut pre.glob_nds, &self.int_pts[si..si+3]);
            ipw_j.set_val(self.ip_wt[ip]);
            ipw_j.mult(&det_j);
            rho.set_val(0.0);
            for j in 0..self.num_nds {
                td.set_val_dfd1(&n_vec[j]);
                td.mult(&pre.fl_den[j]);
                rho.add(&td);
            }
            nvi = 0;
            for k in 0..3 {
                vel[k].set_val(0.0);
                for j in 0..self.num_nds {
                    //nvi = k*self.num_nds + j;
                    td.set_val_dfd1(&n_vec[j]);
                    td.mult(&pre.fl_vel[nvi]);
                    vel[k].add(&td);
                    nvi += 1;
                }
            }

            //momentum
            for j in 1..4 {
                td.set_val_dfd1(&pre.f_per_mass[j-1]);
                td.mult(&rho);
                td.mult(&ipw_j);
                r_sum[j].sub(&td);
            }

            //energy
            td.set_val(0.0);
            for j in 0..3 {
                td2.set_val_dfd1(&vel[j]);
                td2.mult(&pre.f_per_mass[j]);
                td.add(&td2);
            }
            td.mult(&rho);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            td.set_val_dfd1(&pre.hg_per_vol);
            td.mult(&ipw_j);
            r_sum[4].sub(&td);

            si += 3;
        }

        let mut ri = 0usize;
        for _i in 0..self.num_nds {
            for j in 0..6 {
                r_vec[ri].add(&r_sum[j]);
                ri += 1;
            }
        }

    }

    pub fn get_rf_dfd1(&self, r_glob : &mut Vec<DiffDoub1>, pre : &mut DiffDoub1FlPrereq, dynamic : bool, fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {
        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_dfd1"),
            Some(x) => &mut x.dat,
        };
        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd1(r_vec, pre);
        }
        self.get_rf_convective_dfd1( r_vec, pre, fc_vec);
        self.get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_pressure_dfd1(r_vec, pre, fc_vec, ec_vec);
        self.get_rf_load_dfd1(r_vec,pre);

        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                r_glob[rg+j].add(&r_vec[ri]);
                ri += 1;
            }
        }
    }

    //end dup
 
//end skip 
 
 
 

    pub fn get_rf_mat_col(&self, glob_mat : &mut SparseMat, pre : &mut DiffDoub1FlPrereq, dynamic : bool, glob_col : usize, prim_var : bool, mult_fact : f64,
        fc_vec : &Vec<Face>, ec_vec : &Vec<CentData>, nd_ar : &Vec<Node>, scr_dfd : &mut IterMut<'_,DiffDoub1Scr>) {

        let r_vec = match scr_dfd.next() {
            None => panic!("Error: ran out of scratch vectors, get_rf_mat()"),
            Some(x) => &mut x.dat,
        };

        for i in r_vec.iter_mut() {
            i.set_val(0.0);
        }
        if dynamic {
            self.get_rf_unsteady_dfd1(r_vec, pre);
        }
        if prim_var {
            self.get_rf_convective_dfd1(r_vec, pre, fc_vec);
            self.get_rf_viscous_dfd1(r_vec, pre, fc_vec, ec_vec);
            self.get_rf_pressure_dfd1(r_vec, pre, fc_vec, ec_vec);
            self.get_rf_load_dfd1(r_vec,pre);
        }
        
        let mut ri = 0usize;
        let mut rg : usize;
        for nd in self.nodes.iter() {
            rg = 6*nd_ar[*nd].sorted_rank;
            for j in 0..6 {
                glob_mat.add_entry(rg+j, glob_col, mult_fact*r_vec[ri].dval);
                ri += 1;
            }
        }
    }

}