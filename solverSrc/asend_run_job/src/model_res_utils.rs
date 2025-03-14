use crate::model::*;
use crate::diff_doub::*;
use crate::matrix_functions::*;
use crate::fmath::*;

impl Model {
    pub fn get_element_stress_strain(&mut self, stress : &mut [f64], strain : &mut [f64], el_i : usize, layer : usize, ip : usize) {
        let mut dds = [DiffDoub0::new(); 6];
        let mut dde = [DiffDoub0::new(); 6];
        let mut spt = [0f64; 3];
        if ip == 0 {
            for i in 0..3 {
                spt[i] = self.elements[el_i].s_cent[i];
            }
        }
        else {
            let si = 3*(ip - 1);
            vec_to_ar(&mut spt, & self.elements[el_i].int_pts,si,si+3);
        }
        let sci = self.solve_cmd;
        let nl_g = self.job[sci].nonlinear_geom;
        self.elements[el_i].get_stress_prereq_dfd0(&mut self.d0_pre, &mut self.sections, &mut self.materials, &mut self.nodes, &self.design_vars);
        self.elements[el_i].get_stress_strain_dfd0(&mut dds, &mut dde, &mut spt, layer, nl_g, &mut self.d0_pre);
        for i in 0..6 {
            stress[i] = dds[i].val;
            strain[i] = dde[i].val;
        }
    }

    pub fn set_element_user_status(&mut self, el_i : usize, stat_str : &str) {
        self.elements[el_i].user_status = stat_str.to_string();
    }

    pub fn array_eig(vals : &mut [f64], vecs : &mut [f64], mat : & [f64]) {
        let mut t_mat = [0.0f64; 9];
        let mut t_vec : [f64; 3];
        let mut t_next = [0.0f64; 3];
        let mut dp : f64;
        let mut mag : f64;
        let mut minv : f64;
        let mut mi : usize;
        let mut ct : usize;
        for i in 0..3 {
            t_vec = [0.0f64; 3];
            t_vec[i] = 1.0;
            dp = 0.0;
            mag = 1.0;
            ct = 0;
            while fabs(dp) < 0.99999999*mag && ct < 50 {
                mi = 0;
                for j in 0..3 {
                    t_next[j] = 0.0;
                    for k in 0..3 {
                        t_next[j] += mat[mi]*t_vec[k];
                        mi += 1;
                    }
                }
                mi = 0;
                for _j in 0..i {
                    dp = t_mat[mi]*t_next[0] + t_mat[mi+1]*t_next[1] + t_mat[mi+2]*t_next[2];
                    t_next[0] -= dp*t_mat[mi];
                    t_next[1] -= dp*t_mat[mi+1];
                    t_next[2] -= dp*t_mat[mi+2];
                    mi += 3;
                }
                mag = sqrt(t_next[0]*t_next[0] + t_next[1]*t_next[1] + t_next[2]*t_next[2]);
                if mag < 1.0e-12 {
                    t_vec = [0.8824799487067698, 0.10858326094306138, 0.45764485747516953];
                    dp = 0.0;
                    mag = 1.0;
                }
                else {
                    dp = t_vec[0]*t_next[0] + t_vec[1]*t_next[1] + t_vec[2]*t_next[2];
                    minv = 1.0/mag;
                    t_vec[0] = minv*t_next[0];
                    t_vec[1] = minv*t_next[1];
                    t_vec[2] = minv*t_next[2];
                }
                ct += 1;
            }
            vals[i] = dp;
            mi = 3*i;
            t_mat[mi] = t_vec[0];
            t_mat[mi+1] = t_vec[1];
            t_mat[mi+2] = t_vec[2];
        }
        let mut order = [0usize, 1, 2];
        let mut fswap : f64;
        let mut iswap : usize;
        for _i in 0..2 {
            for j in 0..2 {
                if vals[j+1] < vals[j] {
                    fswap = vals[j];
                    vals[j] = vals[j+1];
                    vals[j+1] = fswap;
                    iswap = order[j];
                    order[j] = order[j+1];
                    order[j+1] = iswap;
                }
            }
        }
        let mut mi2 : usize;
        for i in 0..3 {
            mi = 3*i;
            mi2 = 3*order[i];
            for j in 0..3 {
                vecs[mi+j] = t_mat[mi2+j];
            }
        }
    }

    pub fn principal_stress(p_stress : &mut [f64], p_dir : &mut [f64], stress : & [f64]) {
        let mut m = [0.0f64; 9];
        m[0] = stress[0];
        m[1] = stress[3];
        m[2] = stress[4];
        m[3] = stress[3];
        m[4] = stress[1];
        m[5] = stress[5];
        m[6] = stress[4];
        m[7] = stress[5];
        m[8] = stress[2];
        let det = m[0]*(m[4]*m[8] - m[5]*m[7]) + m[1]*(m[5]*m[6] - m[3]*m[8]) + m[2]*(m[3]*m[7] - m[4]*m[6]);
        if fabs(det) < 1.0e-12f64 {
            m[0] += 1.0f64;
            m[4] += 1.0f64;
            m[8] += 1.0f64;
        }
        Model::array_eig(p_stress, p_dir, &m);
        if fabs(det) < 1.0e-12f64 {
            p_stress[0] -= 1.0f64;
            p_stress[1] -= 1.0f64;
            p_stress[2] -= 1.0f64;
        }
    }

    pub fn mises(stress : & [f64]) -> f64 {
        let a = stress[0] - stress[1];
        let b = stress[1] - stress[2];
        let c = stress[2] - stress[0];
        let d = stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5];
        sqrt(0.5f64*(a*a + b*b + c*c + 6.0f64*d))
    }

    pub fn principal_strain(p_strain : &mut [f64], p_dir : &mut [f64], strain : & [f64]) {
        let mut m = [0.0f64; 9];
        m[0] = strain[0];
        m[1] = 0.5*strain[3];
        m[2] = 0.5*strain[4];
        m[3] = 0.5*strain[3];
        m[4] = strain[1];
        m[5] = 0.5*strain[5];
        m[6] = 0.5*strain[4];
        m[7] = 0.5*strain[5];
        m[8] = strain[2];
        let det = m[0]*(m[4]*m[8] - m[5]*m[7]) + m[1]*(m[5]*m[6] - m[3]*m[8]) + m[2]*(m[3]*m[7] - m[4]*m[6]);
        if fabs(det) < 1.0e-12f64 {
            m[0] += 1.0f64;
            m[4] += 1.0f64;
            m[8] += 1.0f64;
        }
        Model::array_eig(p_strain, p_dir, &m);
        if fabs(det) < 1.0e-12f64 {
            p_strain[0] -= 1.0f64;
            p_strain[1] -= 1.0f64;
            p_strain[2] -= 1.0f64;
        }
    }

}