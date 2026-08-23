use crate::model::particle_source::ParticleSource;
use crate::constants::MAX_INT;
use crate::list_ent::{QuadFloat, DualFloat};
use crate::model::element::Element;
use crate::model::node::Node;
use crate::nd_el_set::Set;
use crate::diff_doub::*;
use crate::matrix_functions::*;

impl ParticleSource {
    pub fn is_active(&self, time : f64) -> bool {
        self.active_time[0] <= time && self.active_time[1] >= time
    }

    pub fn get_coord(&self, crd : &mut [f64], time : f64) {
        let mut prev = QuadFloat::new();
        match self.coord.front() {
            None => (),
            Some(x) => {prev.f1 = x.f1;
                                    prev.f2 = x.f2;
                                    prev.f3 = x.f3;
                                    prev.f4 = x.f4;},
        }
        for ent in self.coord.iter() {
            if ent.f1 > time {
                crd[0] = prev.f2 + (ent.f2 - prev.f2)*(time - prev.f1)/(ent.f1 - prev.f1);
                crd[1] = prev.f3 + (ent.f3 - prev.f3)*(time - prev.f1)/(ent.f1 - prev.f1);
                crd[2] = prev.f4 + (ent.f4 - prev.f4)*(time - prev.f1)/(ent.f1 - prev.f1);
                return;
            }
            prev.f1 = ent.f1;
            prev.f2 = ent.f2;
            prev.f3 = ent.f3;
            prev.f4 = ent.f4;
        }
        crd[0] = prev.f2;
        crd[1] = prev.f3;
        crd[2] = prev.f4;
    }

    pub fn get_mean_vel(&self, vel : &mut [f64], time : f64) {
        let mut prev = QuadFloat::new();
        match self.mean_vel.front() {
            None => (),
            Some(x) => {prev.f1 = x.f1;
                                    prev.f2 = x.f2;
                                    prev.f3 = x.f3;
                                    prev.f4 = x.f4;},
        }
        for ent in self.mean_vel.iter() {
            if ent.f1 > time {
                vel[0] = prev.f2 + (ent.f2 - prev.f2)*(time - prev.f1)/(ent.f1 - prev.f1);
                vel[1] = prev.f3 + (ent.f3 - prev.f3)*(time - prev.f1)/(ent.f1 - prev.f1);
                vel[2] = prev.f4 + (ent.f4 - prev.f4)*(time - prev.f1)/(ent.f1 - prev.f1);
                return;
            }
            prev.f1 = ent.f1;
            prev.f2 = ent.f2;
            prev.f3 = ent.f3;
            prev.f4 = ent.f4;
        }
        vel[0] = prev.f2;
        vel[1] = prev.f3;
        vel[2] = prev.f4;
    }

    pub fn get_temperature(&self, time : f64) -> f64 {
        let mut prev = DualFloat::new();
        match self.temp.front() {
            None => (),
            Some(x) => {prev.f1 = x.f1;
                                    prev.f2 = x.f2;},
        }
        for ent in self.temp.iter() {
            if ent.f1 > time {
                return prev.f2 + (ent.f2 - prev.f2)*(time - prev.f1)/(ent.f1 - prev.f1);
            }
            prev.f1 = ent.f1;
            prev.f2 = ent.f2;
        }
        prev.f2
    }

    pub fn get_frequency(&self, time : f64) -> f64 {
        let mut prev = DualFloat::new();
        match self.frequency.front() {
            None => (),
            Some(x) => {prev.f1 = x.f1;
                                    prev.f2 = x.f2;},
        }
        for ent in self.frequency.iter() {
            if ent.f1 > time {
                return prev.f2 + (ent.f2 - prev.f2)*(time - prev.f1)/(ent.f1 - prev.f1);
            }
            prev.f1 = ent.f1;
            prev.f2 = ent.f2;
        }
        prev.f2
    }

    pub fn out_of_bounds(&self, crd : &[f64]) -> bool {
        if crd[0] < self.x_range[0] {
            return true;
        }
        if crd[0] > self.x_range[1] {
            return true;
        }
        if crd[1] < self.y_range[0] {
            return true;
        }
        if crd[1] > self.y_range[1] {
            return true;
        }
        if crd[2] < self.z_range[0] {
            return true;
        }
        if crd[2] > self.z_range[1] {
            return true;
        }
        false
    }

    pub fn get_dir_cos(&self, a_mat : &mut [DiffDoub1], n1_crd : &mut [DiffDoub1], nd_ar : &mut Vec<Node>) {
        let mut n2_crd = [DiffDoub1::new(); 3];
        let mut n3_crd = [DiffDoub1::new(); 3];
        let mut v1 = [DiffDoub1::new(); 3];
        let mut v2 = [DiffDoub1::new(); 3];
        let mut v3 = [DiffDoub1::new(); 3];
        let mut mag = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();

        if self.ref_nodes_i[0] == MAX_INT {
            for i in 0..9 {
                a_mat[i].set_val(0.0);
            }
            for i in 0..3 {
                a_mat[i*4].set_val(1.0);
                n1_crd[i].set_val(0.0);
            }
            return;
        }

        nd_ar[self.ref_nodes_i[0]].get_def_crd_vel(n1_crd);
        nd_ar[self.ref_nodes_i[1]].get_def_crd_vel(&mut n2_crd);
        nd_ar[self.ref_nodes_i[2]].get_def_crd_vel(&mut n3_crd);

        for i in 0..3 {
            v1[i].set_val_dfd1(&n2_crd[i]);
            v1[i].sub(&n1_crd[i]);
        }
        
        for i in 0..3 {
            tmp2.set_val_dfd1(&v1[i]);
            tmp2.sqr();
            tmp.add(&tmp2);
        }
        tmp.sqt();
        mag.set_val(1.0);
        mag.dvd(&tmp);
        for i in 0..3 {
            a_mat[i].set_val_dfd1(&mag);
            a_mat[i].mult(&v1[i]);
        }

        for i in 0..3 {
            v2[i].set_val_dfd1(&n3_crd[i]);
            v2[i].sub(&n1_crd[i]);
        }
        cross_prod_dfd1(&mut v3, &v1, &v2);
        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&v3[i]);
            tmp2.sqr();
            tmp.add(&tmp2);
        }
        tmp.sqt();
        mag.set_val(1.0);
        mag.dvd(&tmp);
        for i in 0..3 {
            a_mat[i+6].set_val_dfd1(&mag);
            a_mat[i+6].mult(&v3[i]);
        }

        cross_prod_dfd1(&mut v2, &a_mat[6..9], &a_mat[0..3]);
        for i in 0..3 {
            a_mat[i+3].set_val_dfd1(&v2[i]);
        }

    }

    pub fn get_global_crd(&self, glob_crd : &mut [DiffDoub1], a_mat : &[DiffDoub1], n1_crd : &[DiffDoub1], time : f64) {
        glob_crd[0].set_val(0.0);
        glob_crd[1].set_val(0.0);
        glob_crd[2].set_val(0.0);
        let mut loc_crd = [DiffDoub1::new(); 3];
        let mut crd = [0f64; 3];
        let mut vel = [0f64; 3];
        let mut tmp = DiffDoub1::new();
        self.get_coord(&mut crd, time);
        self.get_mean_vel(&mut vel, time);

        for i in 0..3 {
            loc_crd[i].set_val(crd[i]);
            if self.vel_in_local {
                loc_crd[i].dval = vel[i];
            }
        }

        let mut k = 0usize;
        for i in 0..3 {
            for j in 0..3 {
                tmp.set_val_dfd1(&a_mat[k]);
                tmp.mult(&loc_crd[i]);
                glob_crd[j].add(&tmp);
                k += 1;
            }
        }

        glob_crd[0].add(&n1_crd[0]);
        glob_crd[1].add(&n1_crd[1]);
        glob_crd[2].add(&n1_crd[2]);

        if !self.vel_in_local {
            for i in 0..3 {
                glob_crd[i].dval = vel[i];
            }
        }
    }

    pub fn get_local_crd(&self, loc_crd : &mut [f64], glob_crd : &mut [DiffDoub0], a_mat : &[DiffDoub1], n1_crd : &[DiffDoub1]) {
        glob_crd[0].val -= n1_crd[0].val;
        glob_crd[1].val -= n1_crd[1].val;
        glob_crd[2].val -= n1_crd[2].val;

        loc_crd[0] = 0.0;
        loc_crd[1] = 0.0;
        loc_crd[2] = 0.0;

        let mut k = 0usize;
        for i in 0..3 {
            for j in 0..3 {
                loc_crd[i] += a_mat[k].val*glob_crd[j].val;
                k += 1;
            }
        }

    }

    pub fn deact_ob_els(&self, el_ar : &mut Vec<Element>, nd_ar : &mut Vec<Node>, el_sets : &Vec<Set>) {
        let mut crd = [DiffDoub0::new(); 3];
        let mut loc_crd = [0f64; 3];
        let mut a_mat = [DiffDoub1::new(); 9];
        a_mat[0].set_val(1.0);
        a_mat[4].set_val(1.0);
        a_mat[8].set_val(1.0);
        let mut n1_crd = [DiffDoub1::new(); 3];

        if self.ref_nodes_i[0] < MAX_INT {
            self.get_dir_cos(&mut a_mat, &mut n1_crd, nd_ar);
        }
        for el in el_sets[self.elset_pt].labels.iter() {
            if el_ar[*el].this_type == 1 {
                nd_ar[el_ar[*el].nodes[0]].get_def_crd_dfd0(&mut crd);
                self.get_local_crd(&mut loc_crd, &mut crd, &a_mat, &n1_crd);
                if self.out_of_bounds(&loc_crd) {
                    el_ar[*el].is_active = false;
                }
            }
        }
    }

    pub fn release_if_clear(&mut self, time : f64, del_t : f64, el_ar : &mut Vec<Element>, nd_ar : &mut Vec<Node>, el_sets : &Vec<Set>) {
        if self.swapset_pt < MAX_INT {
            self.release_swap(time, del_t, el_ar, nd_ar, el_sets);
        }
        else {
            self.release_std(time, del_t, el_ar, nd_ar, el_sets);
        }
    }

    pub fn release_std(&mut self, time : f64, del_t : f64, el_ar : &mut Vec<Element>, nd_ar : &mut Vec<Node>, el_sets : &Vec<Set>) {
        let mut crd = [DiffDoub0::new(); 3];
        let mut ndi : usize;
        let mut nd : &mut Node;
        let mut a_mat = [DiffDoub1::new(); 9];
        a_mat[0].set_val(1.0);
        a_mat[4].set_val(1.0);
        a_mat[8].set_val(1.0);
        let mut n1_crd = [DiffDoub1::new(); 3];
        let mut glob_crd = [DiffDoub1::new(); 3];
        let mut loc_crd = [0f64; 3];

        self.since_release += del_t;
        if self.is_active(time) {
            let freq = self.get_frequency(time);
            if self.since_release >= freq {
                if self.ref_nodes_i[0] < MAX_INT {
                    self.get_dir_cos(&mut a_mat, &mut n1_crd, nd_ar);
                }
                self.get_global_crd(&mut glob_crd, &a_mat, &n1_crd, time);
                //--------------------------------------------------------
                // todo: come up with a way of generating random velocity
                //-------------------------------------------------------
                for el in el_sets[self.elset_pt].labels.iter() {
                    if el_ar[*el].this_type == 1 {
                        ndi = el_ar[*el].nodes[0];
                        nd = &mut nd_ar[ndi];
                        nd.get_def_crd_dfd0(&mut crd);
                        self.get_local_crd(&mut loc_crd, &mut crd, &a_mat, &n1_crd);

                        if self.out_of_bounds(&loc_crd) {
                            el_ar[*el].is_active = true;
                            for i in 0..3 {
                                nd.prev_disp[i] = glob_crd[i].val - nd.coord_dfd0[i].val;
                                nd.pp_disp[i] = nd.prev_disp[i] - del_t*glob_crd[i].dval;
                                nd.displacement[i] = nd.prev_disp[i] + del_t*glob_crd[i].dval;
                                nd.prev_vel[i] = glob_crd[i].dval;
                                nd.velocity[i] = glob_crd[i].dval;
                                nd.prev_acc[i] = 0f64;
                                nd.acceleration[i] = 0f64;
                            }
                            nd.prev_temp = self.get_temperature(time);
                            nd.pp_temp = nd.prev_temp;
                            nd.temperature = nd.prev_temp;
                            nd.prev_tdot = 0.0;
                            nd.temp_change_rate = 0.0;
                            self.since_release -= freq;
                            return;
                        }
                    }
                }
            }
        }
    }

    pub fn release_swap(&mut self, time : f64, del_t : f64, el_ar : &mut Vec<Element>, nd_ar : &mut Vec<Node>, el_sets : &Vec<Set>) {
        let mut near_el = MAX_INT;
        let mut near_dist = 1.0e+100;
        let mut near_nd : &Node;
        let mut near_crd = [DiffDoub0::new(); 3];
        let mut del : f64;
        let mut dist : f64;
        let mut a_mat = [DiffDoub1::new(); 9];
        a_mat[0].set_val(1.0);
        a_mat[4].set_val(1.0);
        a_mat[8].set_val(1.0);
        let mut n1_crd = [DiffDoub1::new(); 3];
        let mut glob_crd = [DiffDoub1::new(); 3];
        let mut loc_crd = [0f64; 3];
        let half_len = 0.5*self.spacing*((self.refine_lev - 1) as f64);
        let mut corner_crd = [0f64; 3];
        let mut ct : usize;
        let n1 : usize;
        let mut n2 : usize;
        let mut ins_el : usize;

        if self.is_active(time) {
            if self.ref_nodes_i[0] < MAX_INT {
                self.get_dir_cos(&mut a_mat, &mut n1_crd, nd_ar);
            }
            self.get_global_crd(&mut glob_crd, &a_mat, &n1_crd, time);

            for el in el_sets[self.swapset_pt].labels.iter() {
                near_nd = &nd_ar[el_ar[*el].nodes[0]];
                near_nd.get_def_crd_dfd0(&mut near_crd);
                
                dist = 0.0;
                for i in 0..3 {
                    del = near_crd[i].val - glob_crd[i].val;
                    dist += del*del;
                }
                dist = dist.sqrt();
                
                if dist < near_dist {
                    near_el = *el;
                    near_dist = dist;
                }
            }

            if near_dist < self.swap_dist {
                n1 = el_ar[near_el].nodes[0];
                nd_ar[n1].get_def_crd_dfd0(&mut near_crd);
                for i in 0..3 {
                    corner_crd[i] = near_crd[i].val - half_len;
                }

                for e in self.insert_els.iter_mut() {
                    *e = MAX_INT;
                }

                ct = 0;
                for el in el_sets[self.elset_pt].labels.iter() {
                    if ct < self.insert_els.len() {
                        nd_ar[el_ar[*el].nodes[0]].get_def_crd_dfd0(&mut near_crd);
                        self.get_local_crd(&mut loc_crd, &mut near_crd, &a_mat, &n1_crd);
                        if self.out_of_bounds(&loc_crd) {
                            self.insert_els[ct] = *el;
                            ct += 1;
                        }
                    }   
                }

                ct = 0;
                for i in 0..self.refine_lev {
                    loc_crd[0] = corner_crd[0] + (i as f64)*self.spacing;
                    for j in 0..self.refine_lev {
                        loc_crd[1] = corner_crd[1] + (j as f64)*self.spacing;
                        for k in 0..self.refine_lev {
                            loc_crd[2] = corner_crd[2] + (k as f64)*self.spacing;
                            ins_el = self.insert_els[ct];
                            if ins_el < MAX_INT {
                                el_ar[ins_el].is_active = true;
                                n2 = el_ar[ins_el].nodes[0];
                                for m in 0..3 {
                                    nd_ar[n2].displacement[m] = loc_crd[m] - nd_ar[n2].coord_dfd0[m].val;
                                    nd_ar[n2].velocity[m] = nd_ar[n1].velocity[m];
                                    nd_ar[n2].acceleration[m] = nd_ar[n1].acceleration[m];
                                    nd_ar[n2].prev_disp[m] = nd_ar[n2].displacement[m] - del_t*nd_ar[n2].velocity[m] + 0.5*del_t*del_t*nd_ar[n2].acceleration[m];
                                    nd_ar[n2].pp_disp[m] = nd_ar[n2].displacement[m] - 2.0*del_t*nd_ar[n2].velocity[m] + 2.0*del_t*del_t*nd_ar[n2].acceleration[m];
                                    nd_ar[n2].prev_vel[m] = nd_ar[n1].prev_vel[m];
                                    nd_ar[n2].prev_acc[m] = nd_ar[n1].prev_acc[m];
                                }
                                nd_ar[n2].temperature = nd_ar[n1].temperature;
                                nd_ar[n2].temp_change_rate = nd_ar[n1].temp_change_rate;
                                nd_ar[n2].fl_den = nd_ar[n1].fl_den;
                                nd_ar[n2].fl_den_dot = nd_ar[n1].fl_den_dot;
                                
                                nd_ar[n2].prev_temp = nd_ar[n1].prev_temp;
                                nd_ar[n2].prev_tdot = nd_ar[n1].prev_tdot;
                                nd_ar[n2].prev_fl_den = nd_ar[n1].prev_fl_den;
                                nd_ar[n2].prev_fl_den_dot = nd_ar[n1].prev_fl_den_dot;
                                
                                nd_ar[n2].pp_temp = nd_ar[n1].pp_temp;
                                nd_ar[n2].pp_fl_den = nd_ar[n1].pp_fl_den;
                            }
                            ct += 1;
                        }
                    }
                }

                el_ar[near_el].is_active = false;
                nd_ar[n1].initialize_disp(del_t);
                nd_ar[n1].initialize_temp(del_t);
                nd_ar[n1].initialize_fl_den(del_t);

            }
        }
    }

}