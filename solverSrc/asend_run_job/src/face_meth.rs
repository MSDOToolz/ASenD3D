use crate::face::*;
use crate::diff_doub::*;
use crate::node::*;
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


