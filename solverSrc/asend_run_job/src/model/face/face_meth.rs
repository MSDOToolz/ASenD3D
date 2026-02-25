use crate::model::face::*;
use crate::constants::*;
use crate::fmath::*;
use crate::model::node::*;
use crate::matrix_functions::*;


impl Face {

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
        let mut cent = [0f64; 3];
        let mut dx_ds1 = [0f64; 3];
        let mut dx_ds2 = [0f64; 3];

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
        let mut tmp = DiffDoub0::new();
        let mut shft : usize;
        let lnd = &self.loc_nodes;
        
        if self.num_nds == 4 {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd0(&nd_crd[lnd[2] + shft]);
                v1[i].sub(&nd_crd[lnd[0] + shft]);
                v2[i].set_val_dfd0(&nd_crd[lnd[3] + shft]);
                v2[i].sub(&nd_crd[lnd[1] + shft]);   
            }
        }
        else {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd0(&nd_crd[lnd[1] + shft]);
                v1[i].sub(&nd_crd[lnd[0] + shft]);
                v2[i].set_val_dfd0(&nd_crd[lnd[2] + shft]);
                v2[i].sub(&nd_crd[lnd[0] + shft]);
            }
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

    //end dup
 
//skip 
 
//DiffDoub1 versions: 
    //dup1

    pub fn get_area_normal_dfd1(&self, area : &mut DiffDoub1, norm : &mut [DiffDoub1], nd_crd : & Vec<DiffDoub1>, el_nn : usize) {
        let mut v1 = [DiffDoub1::new(); 3];
        let mut v2 = [DiffDoub1::new(); 3];
        let mut tmp = DiffDoub1::new();
        let mut shft : usize;
        let lnd = &self.loc_nodes;
        
        if self.num_nds == 4 {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd1(&nd_crd[lnd[2] + shft]);
                v1[i].sub(&nd_crd[lnd[0] + shft]);
                v2[i].set_val_dfd1(&nd_crd[lnd[3] + shft]);
                v2[i].sub(&nd_crd[lnd[1] + shft]);   
            }
        }
        else {
            for i in 0..3 {
                shft = i*el_nn;
                v1[i].set_val_dfd1(&nd_crd[lnd[1] + shft]);
                v1[i].sub(&nd_crd[lnd[0] + shft]);
                v2[i].set_val_dfd1(&nd_crd[lnd[2] + shft]);
                v2[i].sub(&nd_crd[lnd[0] + shft]);
            }
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


