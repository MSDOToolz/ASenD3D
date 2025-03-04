use crate::mesh_element::*;
use crate::mesh_node::*;
use crate::utilities::*;

impl MeshElement {
    pub fn get_centroid(& self, cent : &mut [f64], nd_ar : &Vec<MeshNode>) {
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        let n4 = &nd_ar[self.nodes[3]].coord;
        cent[0] = 0.25 * (n1[0] + n2[0] + n3[0] + n4[0]);
        cent[1] = 0.25 * (n1[1] + n2[1] + n3[1] + n4[1]);
        cent[2] = 0.25 * (n1[2] + n2[2] + n3[2] + n4[2]);
        return;
    }

    pub fn get_volume(& self, nd_ar : &Vec<MeshNode>) -> f64 {
        let mut mat : [f64; 9] = [0f64; 9];
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        let n4 = &nd_ar[self.nodes[3]].coord;
        mat[0] = n2[0] - n1[0];
        mat[3] = n2[1] - n1[1];
        mat[6] = n2[2] - n1[2];
        mat[1] = n3[0] - n1[0];
        mat[4] = n3[1] - n1[1];
        mat[7] = n3[2] - n1[2];
        mat[2] = n4[0] - n1[0];
        mat[5] = n4[1] - n1[1];
        mat[8] = n4[2] - n1[2];
        q_rfactor(&mut mat,  3,  0,  2,  0,  2,  0);
        let vol : f64 =  mat[0] * mat[4] * mat[8];
        return  vol;
    }

    pub fn point_in(& self, pt : & [f64], nd_ar : &Vec<MeshNode>) -> bool {
        let mut mat : [f64; 9] = [0f64; 9];
        let mut rhs : [f64; 3] = [0f64; 3];
        let mut soln : [f64; 3] = [0f64; 3];
        let sol_sum : f64;
        let n1 = &nd_ar[self.nodes[0]].coord;
        let n2 = &nd_ar[self.nodes[1]].coord;
        let n3 = &nd_ar[self.nodes[2]].coord;
        let n4 = &nd_ar[self.nodes[3]].coord;
        mat[0] = n2[0] - n1[0];
        mat[3] = n2[1] - n1[1];
        mat[6] = n2[2] - n1[2];
        mat[1] = n3[0] - n1[0];
        mat[4] = n3[1] - n1[1];
        mat[7] = n3[2] - n1[2];
        mat[2] = n4[0] - n1[0];
        mat[5] = n4[1] - n1[1];
        mat[8] = n4[2] - n1[2];
        rhs[0] = pt[0] - n1[0];
        rhs[1] = pt[1] - n1[1];
        rhs[2] = pt[2] - n1[2];
        q_rfactor(&mut mat,  3,  0,  2,  0,  2,  0);
        solveq_rx_eqb(&mut soln, &mut  mat, &mut  rhs,  3,  0,  2,  0,  2,  0);
        if soln[0] > 0.0000000001 && soln[1] > 0.0000000001 && soln[2] > 0.0000000001 {
            sol_sum = soln[0] + soln[1] + soln[2];
            if sol_sum < 0.9999999999 {
                return  true;
            }
        }
        return  false;
    }

}


