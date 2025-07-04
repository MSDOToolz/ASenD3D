use crate::mesh_face::*;
use crate::mesh_node::*;
use crate::utilities::*;
use crate::fmath::*;


impl MeshFace {
    pub fn copy_data(&mut self, f_in : &MeshFace) {
        for i1 in 0..3 {
            self.nodes[i1] = f_in.nodes[i1];
            self.norm_dir[i1] = f_in.norm_dir[i1];
        }
        self.elements[0] = f_in.elements[0];
        self.elements[1] = f_in.elements[1];
        self.proj_dist = f_in.proj_dist;
        return;
    }

    pub fn init_norm_dir(&mut self, nd_ar : &Vec<MeshNode>) {
        let mut v1 : [f64; 3] = [0f64; 3];
        let mut v2 : [f64; 3] = [0f64; 3];
        let mut cp : [f64; 3] = [0f64; 3];
        let mut mag : f64;
        let area : f64;
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        for i1 in 0..3 {
            v1[i1] = n2[i1] - n1[i1];
            v2[i1] = n3[i1] - n1[i1];
        }
        cross_prod(&mut cp, &mut  v1, &mut  v2);
        mag = sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]);
        area = 0.5 * mag;
        mag = 1.0 / mag;
        self.norm_dir[0] = mag * cp[0];
        self.norm_dir[1] = mag * cp[1];
        self.norm_dir[2] = mag * cp[2];
        self.proj_dist = 1.2408064788027997 * sqrt(area);// (2^1.5/3^0.75)*sqrt(area)
        return;
    }

    pub fn norm_dir_from_el_cent(&mut self, cent : &[f64], nd_ar : &Vec<MeshNode>) {
        let mut fc_cent : [f64; 3] = [0f64; 3];
        let mut d_vec : [f64; 3] = [0f64; 3];
        self.get_centroid(&mut fc_cent,  nd_ar);
        d_vec[0] = fc_cent[0] - cent[0];
        d_vec[1] = fc_cent[1] - cent[1];
        d_vec[2] = fc_cent[2] - cent[2];
        let dp : f64 =  d_vec[0] * self.norm_dir[0] + d_vec[1] * self.norm_dir[1] + d_vec[2] * self.norm_dir[2];
        if dp < 0.0 {
            self.norm_dir[0]  *=  -1.0;
            self.norm_dir[1]  *=  -1.0;
            self.norm_dir[2]  *=  -1.0;
        }
        return;
    }

    pub fn get_centroid(& self, cent : &mut [f64], nd_ar : &Vec<MeshNode>) {
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        cent[0] = 0.333333333333333 * (n1[0] + n2[0] + n3[0]);
        cent[1] = 0.333333333333333 * (n1[1] + n2[1] + n3[1]);
        cent[2] = 0.333333333333333 * (n1[2] + n2[2] + n3[2]);
        return;
    }

    pub fn get_longest_edge_len(& self, nd_ar : &Vec<MeshNode>) -> f64 {
        let mut d_vec : [f64; 3] = [0f64; 3];
        let mut dist : f64;
        let mut max_dist : f64 =  0.0;
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        d_vec[0] = n2[0] - n1[0];
        d_vec[1] = n2[1] - n1[1];
        d_vec[2] = n2[2] - n1[2];
        dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
        if dist > max_dist {
            max_dist = dist;
        }
        d_vec[0] = n3[0] - n2[0];
        d_vec[1] = n3[1] - n2[1];
        d_vec[2] = n3[2] - n2[2];
        dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
        if dist > max_dist {
            max_dist = dist;
        }
        d_vec[0] = n1[0] - n3[0];
        d_vec[1] = n1[1] - n3[1];
        d_vec[2] = n1[2] - n3[2];
        dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
        if dist > max_dist {
            max_dist = dist;
        }
        return  max_dist;
    }

    pub fn get_intersection(& self, out_param : &mut [f64], pt : &[f64], vec : & [f64], nd_ar : & Vec<MeshNode>) -> bool {
        let mut mat : [f64; 9] = [0f64; 9];
        let mut mat_mag : f64;
        let mut rhs : [f64; 3] = [0f64; 3];
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        mat[0] = vec[0];
        mat[3] = vec[1];
        mat[6] = vec[2];
        mat[1] = n1[0] - n2[0];
        mat[4] = n1[1] - n2[1];
        mat[7] = n1[2] - n2[2];
        mat[2] = n1[0] - n3[0];
        mat[5] = n1[1] - n3[1];
        mat[8] = n1[2] - n3[2];
        mat_mag = 0.0;
        for i1 in 0..9 {
            mat_mag  +=  mat[i1] * mat[i1];
        }
        mat_mag = sqrt(mat_mag);
        q_rfactor(&mut mat,  3,  0,  2,  0,  2,  0);
        let det : f64 =  mat[0] * mat[4] * mat[8];
        if fabs(det) < 1.0e-6*mat_mag*mat_mag*mat_mag {
            return  false;
        }
        else {
            rhs[0] = n1[0] - pt[0];
            rhs[1] = n1[1] - pt[1];
            rhs[2] = n1[2] - pt[2];
            solveq_rx_eqb(out_param, &mut  mat, &mut  rhs,  3,  0,  2,  0,  2,  0);
            if out_param[1] > 0.0000000001 && out_param[2] > 0.0000000001 {
                if (out_param[1] + out_param[2]) < 0.9999999999 {
                    return  true;
                }
            }
            return  false;
        }
    }

    pub fn edges_intersect(& self, fc : usize, dist_tol : f64, nd_ar : &Vec<MeshNode>, fc_ar : &Vec<MeshFace>) -> bool {
        let mut i3 : usize;
        let mut i4 : usize;
        let mut mat : [f64; 6] = [0f64; 6];
        let mut rhs : [f64; 3] = [0f64; 3];
        let mut soln : [f64; 2] = [0f64; 2];
        let mut det : f64;
        let mut d_vec : [f64; 3] = [0f64; 3];
        let mut dist : f64;
        let mut v1 : [f64; 3] = [0f64; 3];
        let mut v2 : [f64; 3] = [0f64; 3];
        let mut mat_mag : f64;
        
        let fc_nodes = &fc_ar[fc].nodes;
        
        for i1 in 0..3 {
            let n11 = &nd_ar[self.nodes[i1]].coord;
            i3 = i1 + 1;
            if i3 > 2 {
                i3 = 0;
            }
            let n21 = &nd_ar[self.nodes[i3]].coord;
            for i2 in 0..3 {
                let n12 = &nd_ar[fc_nodes[i2]].coord;
                i4 = i2 + 1;
                if i4 > 2 {
                    i4 = 0;
                }
                let n22 = &nd_ar[fc_nodes[i4]].coord;
                v1[0] = n21[0] - n11[0];
                v1[1] = n21[1] - n11[1];
                v1[2] = n21[2] - n11[2];
                v2[0] = n22[0] - n12[0];
                v2[1] = n22[1] - n12[1];
                v2[2] = n22[2] - n12[2];
                mat[0] = v1[0];
                mat[1] = -v2[0];
                mat[2] = v1[1];
                mat[3] = -v2[1];
                mat[4] = v1[2];
                mat[5] = -v2[2];
                rhs[0] = n12[0] - n11[0];
                rhs[1] = n12[1] - n11[1];
                rhs[2] = n12[2] - n11[2];
                mat_mag = 0.0;
                for i5 in 0..6 {
                    mat_mag  +=  mat[i5] * mat[i5];
                }
                mat_mag = sqrt(mat_mag);
                q_rfactor(&mut mat,  2,  0,  2,  0,  1,  0);
                det = mat[0] * mat[3];
                if fabs(det) > 1.0e-6*mat_mag*mat_mag {
                    solveq_rx_eqb(&mut soln, &mut  mat, &mut  rhs,  2,  0,  2,  0,  1,  0);
                    if soln[0] > 0.0000000001 && soln[0] < 0.9999999999 && soln[1] > 0.0000000001 && soln[1] < 0.9999999999 {
                        d_vec[0] = n11[0] + soln[0] * v1[0] - n12[0] - soln[1] * v2[0];
                        d_vec[1] = n11[1] + soln[0] * v1[1] - n12[1] - soln[1] * v2[1];
                        d_vec[2] = n11[2] + soln[0] * v1[2] - n12[2] - soln[1] * v2[2];
                        dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
                        if dist < dist_tol {
                            return  true;
                        }
                    }
                }
            }
        }
        return  false;
    }

    pub fn get_shared_nodes(& self, nd_pts : &mut [usize], shared : &mut [bool], fc : usize, fc_ar : & Vec<MeshFace>) -> usize {
        let fc_nds = &fc_ar[fc].nodes;
        nd_pts[0] = self.nodes[0];
        nd_pts[1] = self.nodes[1];
        nd_pts[2] = self.nodes[2];
        shared[0] = false;
        shared[1] = false;
        shared[2] = false;
        let mut num_shared : usize =  0;
        for i1 in 0..3 {
            for i2 in 0..3 {
                if nd_pts[i1] == fc_nds[i2] {
                    shared[i1] = true;
                    num_shared += 1usize;
                }
            }
        }
        return  num_shared;
    }

    pub fn print_info(&mut self, nd_ar : &mut Vec<MeshNode>) {
        
        println!("{}", "nodes:" );
        for i1 in 0..3 {
            let crd = &nd_ar[self.nodes[i1]].coord;
            println!("{}{}{}{}{}", crd[0] , ", " , crd[1] , ", " , crd[2] );
        }
        return;
    }

}


