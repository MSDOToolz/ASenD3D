use crate::particle_source::ParticleSource;
use crate::constants::MAX_INT;
use crate::list_ent::QuadFloat;
use crate::element::Element;
use crate::node::Node;
use crate::nd_el_set::Set;
use crate::diff_doub::*;

impl ParticleSource {
    pub fn is_active(&self, time : f64) -> bool {
        self.active_time[0] <= time && self.active_time[1] >= time
    }

    pub fn get_coord(&self, crd : &mut [f64], time : f64) {
        let mut pt = match self.coord.front() {
            None => 0f64,
            Some(x) => x.f1,
        };
        let mut px = match self.coord.front() {
            None => 0f64,
            Some(x) => x.f2,
        };
        let mut py = match self.coord.front() {
            None => 0f64,
            Some(x) => x.f3,
        };
        let mut pz = match self.coord.front() {
            None => 0f64,
            Some(x) => x.f4,
        };
        for ent in self.coord.iter() {
            if ent.f1 > time {
                crd[0] = px + (ent.f2 - px)*(time - pt)/(ent.f1 - pt);
                crd[1] = py + (ent.f3 - py)*(time - pt)/(ent.f1 - pt);
                crd[2] = pz + (ent.f4 - pz)*(time - pt)/(ent.f1 - pt);
                return;
            }
            pt = ent.f1;
            px = ent.f2;
            py = ent.f3;
            pz = ent.f4;
        }
        crd[0] = px;
        crd[1] = py;
        crd[2] = pz;
    }

    pub fn get_mean_vel(&self, vel : &mut [f64], time : f64) {
        let mut pt = match self.mean_vel.front() {
            None => 0f64,
            Some(x) => x.f1,
        };
        let mut px = match self.mean_vel.front() {
            None => 0f64,
            Some(x) => x.f2,
        };
        let mut py = match self.mean_vel.front() {
            None => 0f64,
            Some(x) => x.f3,
        };
        let mut pz = match self.mean_vel.front() {
            None => 0f64,
            Some(x) => x.f4,
        };
        for ent in self.mean_vel.iter() {
            if ent.f1 > time {
                vel[0] = px + (ent.f2 - px)*(time - pt)/(ent.f1 - pt);
                vel[1] = py + (ent.f3 - py)*(time - pt)/(ent.f1 - pt);
                vel[2] = pz + (ent.f4 - pz)*(time - pt)/(ent.f1 - pt);
                return;
            }
            pt = ent.f1;
            px = ent.f2;
            py = ent.f3;
            pz = ent.f4;
        }
        vel[0] = px;
        vel[1] = py;
        vel[2] = pz;
    }

    pub fn get_temperature(&self, time : f64) -> f64 {
        let mut pt = match self.temp.front() {
            None => 0f64,
            Some(x) => x.f1,
        };
        let mut pv = match self.temp.front() {
            None => 0f64,
            Some(x) => x.f2,
        };
        for ent in self.temp.iter() {
            if ent.f1 > time {
                return pv + (ent.f2 - pv)*(time - pt)/(ent.f1 - pt);
            }
            pt = ent.f1;
            pv = ent.f2;
        }
        pv
    }

    pub fn get_frequency(&self, time : f64) -> f64 {
        let mut pt = match self.frequency.front() {
            None => 0f64,
            Some(x) => x.f1,
        };
        let mut pf = match self.frequency.front() {
            None => 0f64,
            Some(x) => x.f2,
        };
        for ent in self.frequency.iter() {
            if ent.f1 > time {
                return pf + (ent.f2 - pf)*(time - pt)/(ent.f1 - pt);
            }
            pt = ent.f1;
            pf = ent.f2;
        }
        pf
    }

    pub fn out_of_bounds(&self, crd : &[DiffDoub0]) -> bool {
        if crd[0].val < self.x_range[0] {
            return true;
        }
        if crd[0].val > self.x_range[1] {
            return true;
        }
        if crd[1].val < self.y_range[0] {
            return true;
        }
        if crd[1].val > self.y_range[1] {
            return true;
        }
        if crd[2].val < self.z_range[0] {
            return true;
        }
        if crd[2].val > self.z_range[1] {
            return true;
        }
        false
    }

    pub fn release_if_clear(&mut self, time : f64, del_t : f64, el_ar : &Vec<Element>, nd_ar : &mut Vec<Node>, el_sets : &Vec<Set>) {
        let mut crd = [DiffDoub0::new(); 3];
        let mut crd_ref = [DiffDoub0::new(); 3];
        let mut proj = [0f64; 3];
        let mut vel = [0f64; 3];
        let mut ndi : usize;

        self.since_release += del_t;
        if self.is_active(time) {
            let freq = self.get_frequency(time);
            if self.since_release >= freq {
                if self.ref_node_i != MAX_INT {
                    nd_ar[self.ref_node_i].get_def_crd_dfd0(&mut crd_ref);
                }
                self.get_coord(&mut proj, time);
                self.get_mean_vel(&mut vel, time);
                //--------------------------------------------------------
                // todo: come up with a way of generating random velocity
                //-------------------------------------------------------
                for el in el_sets[self.elset_pt].labels.iter() {
                    if el_ar[*el].this_type == 1 {
                        ndi = el_ar[*el].nodes[0];
                        nd_ar[ndi].get_def_crd_dfd0(&mut crd);
                        if self.out_of_bounds(&crd) {
                            for i in 0..3 {
                                nd_ar[ndi].prev_disp[i] = crd_ref[i].val - nd_ar[ndi].coord_dfd0[i].val + proj[i];
                                nd_ar[ndi].pp_disp[i] = nd_ar[ndi].prev_disp[i] - vel[i];
                                nd_ar[ndi].displacement[i] = nd_ar[ndi].prev_disp[i] + vel[i];
                                nd_ar[ndi].prev_vel[i] = vel[i];
                                nd_ar[ndi].prev_acc[i] = 0f64;
                            }
                            nd_ar[ndi].prev_temp = self.get_temperature(time);
                            nd_ar[ndi].pp_temp = nd_ar[ndi].prev_temp;
                            nd_ar[ndi].temperature = nd_ar[ndi].prev_temp;
                            nd_ar[ndi].prev_tdot = 0.0;
                            self.since_release = 0.0;
                            return;
                        }
                    }
                }
            }
        }
    }
}