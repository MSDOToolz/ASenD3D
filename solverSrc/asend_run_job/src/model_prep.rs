use crate::model::*;
use crate::constants::*;
use crate::matrix_functions::*;
use crate::constraint::*;
use crate::nd_el_set::*;
use crate::node::*;
use crate::element::*;
use crate::face::*;
use crate::diff_doub::*;
use crate::cpp_str::*;
use crate::cpp_map::*;
use crate::design_var::*;
use crate::list_ent::*;
use crate::spatial_grid::*;

use std::collections::LinkedList;

impl Model {

    pub fn add_constraint_conn(nd_con : &mut Vec<Set>, con_lst : &ConstraintList, nd_sets : &Vec<Set>) {
        let mut i1 : usize;
        let mut i2 : usize;

        for this_const in con_lst.const_vec.iter() {
            for term1 in this_const.terms.iter() {
                let t1_labs = &nd_sets[term1.ns_ptr].labels;
                i1 = t1_labs.len();
                for term2 in this_const.terms.iter() {
                    let t2_labs = &nd_sets[term2.ns_ptr].labels;
                    i2 = t2_labs.len();
                    if i1 > 1 && i2 > 1 {
                        let mut iter2 = t2_labs.iter();
                        let mut i2val : usize;
                        for iter1 in t1_labs.iter() {
                            i2val = match iter2.next() {
                                None => panic!("Error: sets of mismatched size found in constraint definition"),
                                Some(x) => *x,
                            };
                            nd_con[*iter1].add_if_absent(i2val);
                            nd_con[i2val].add_if_absent(*iter1);
                        }
                    }
                    else {
                        for iter1 in t1_labs.iter() {
                            for iter2 in t2_labs.iter() {
                                nd_con[*iter1].add_if_absent(*iter2);
                                nd_con[*iter2].add_if_absent(*iter1);
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn reorder_nodes(&mut self, block_dim : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut nd1 : usize;
        let mut nd2 : usize;
        let mut min_ct : usize;
        let mut min_nd : usize;
        let mut min_dist : f64;
        let mut since_restart : usize;
        let mut dist : f64;
        let mut el_num_nds : usize;
        let mut el_dof_per_nd : usize;
        let num_nodes : usize = self.nodes.len();
        
        let mut nodal_conn = vec![Set::new(); num_nodes];
        let mut node_inserted = vec![0usize; num_nodes];
        
        // build nodal connectivity
        for el in self.elements.iter_mut() {
            el_num_nds = el.num_nds;
            el_dof_per_nd = el.dof_per_nd;
            for i1 in 0..el_num_nds {
                nd1 = el.nodes[i1];
                if el_dof_per_nd > 3 {
                    self.nodes[nd1].num_dof = el_dof_per_nd;
                }
                for i2 in i1+1..el_num_nds {
                    nd2 = el.nodes[i2];
                    nodal_conn[nd1].add_if_absent(nd2);
                    nodal_conn[nd2].add_if_absent(nd1);
                    self.nodes[nd1].add_conn_nd(nd2);
                    self.nodes[nd2].add_conn_nd(nd1);
                }
                self.nodes[nd1].add_element(el.label, i1);
                if el.this_type > 100 {
                    self.nodes[nd1].fluid = true;
                }
            }
        }
        
        Model::add_constraint_conn(&mut nodal_conn, &self.elastic_const, &self.node_sets);
        Model::add_constraint_conn(&mut nodal_conn, &self.thermal_const, &self.node_sets);
        Model::add_constraint_conn(&mut nodal_conn, &self.diff_const, &self.node_sets);
        Model::add_constraint_conn(&mut nodal_conn, &self.fluid_const, &self.node_sets);
        
        // find Node with least connectivity
        min_ct = num_nodes;
        min_nd = 0;
        for i1 in 0..num_nodes {
            i2 = nodal_conn[i1].labels.len();
            if i2 < min_ct {
                min_nd = i1;
                min_ct = i2;
            }
            node_inserted[i1] = 0;
        }
        
        // put self.nodes into an integer list in level order.
        //let mut ordered_nds : LinkedList<usize> = LinkedList::new();
        let mut ordered_nds : Vec<usize> = vec![MAX_INT; num_nodes];
        
        ordered_nds[0] = min_nd;
        let mut on_len = 1usize;
        let mut on_i = 0usize;
        nd1 = min_nd;
        //let mut this_nd = ordered_nds.iter_mut();
        node_inserted[min_nd] = 1;
        since_restart = 0;
        //let mut end_this_nd : bool = false;
        let mut end_nbrs : bool;
        nd2 = min_nd;
        
        while on_len < num_nodes {
            //match this_nd.next() {
            //    None => {end_this_nd = true;},
            //    Some(x) => {nd1 = *x; end_this_nd = false;},
            //}
            let mut neighbor_nd = nodal_conn[nd1].labels.iter();
            end_nbrs = false;
            while !end_nbrs && since_restart < block_dim {
                match neighbor_nd.next() {
                    None => {end_nbrs = true;},
                    Some(x) => {nd2 = *x; end_nbrs = false;},
                }
                if node_inserted[nd2] == 0 {
                    ordered_nds[on_len] = nd2;
                    on_len += 1usize;
                    node_inserted[nd2] = 1;
                    since_restart += 1usize;
                }
            }
            //on_i += 1usize;
            //nd1 = ordered_nds[on_i];
            if on_len < num_nodes {
                if since_restart >= block_dim || on_i == (on_len - 1) {
                    // nd1 = match ordered_nds.back() {
                    //     None => panic!("Error: failure to access the back element of ordered_nds."),
                    //     Some(x) => *x,
                    // };
                    nd1 = ordered_nds[on_len - 1];
                    min_dist = 1.0e+100;
                    for i1 in 0..num_nodes {
                        if node_inserted[i1] == 0 {
                            dist = get_dist(& self.nodes[nd1].coord, & self.nodes[i1].coord);
                            if dist < min_dist {
                                min_dist = dist;
                                min_nd = i1;
                            }
                        }
                    }
                    ordered_nds[on_len] = min_nd;
                    on_len += 1usize;
                    nd1 = min_nd;
                    node_inserted[min_nd] = 1;
                    since_restart += 1usize;
                    on_i = on_len - 1;
                    if since_restart >= block_dim {
                        since_restart = 0;
                        // nd2 = 0usize;
                        // // Move the iterator up until it's on nd1
                        // while (nd2 != nd1) {
                        //     nd2 = match this_nd.next() {
                        //         None => panic!("Error: failure to reset ordered nodes iterator after reaching end of matrix block"),
                        //         Some(x) => *x,
                        //     };
                        // }
                    }
                    // if (end_this_nd) {
                    //     this_nd = ordered_nds.iter_mut();
                    //     end_this_nd = false;
                    // }
                }
                else {
                    on_i += 1;
                    nd1 = ordered_nds[on_i];
                }
            }
        }
        
        // update the global degree of freedom indexes for the self.nodes
        
        i2 = 0;// index in elastic matrix
        i3 = 0;// sorted rank for solid nodes
        i4 = 0;// sorted rank for fluid nodes
        let mut this_nd : &mut Node;
        for ndi in ordered_nds.iter_mut() {
            this_nd = &mut self.nodes[*ndi];
            if this_nd.fluid {
                this_nd.sorted_rank = i4;
                i4 += 1;
            }
            else {
                this_nd.sorted_rank = i3;
                i3 += 1usize;
                this_nd.dof_index[0] = i2;
                i2 += 1usize;
                this_nd.dof_index[1] = i2;
                i2 += 1usize;
                this_nd.dof_index[2] = i2;
                i2 += 1usize;
                if this_nd.num_dof == 6 {
                    this_nd.dof_index[3] = i2;
                    i2 += 1usize;
                    this_nd.dof_index[4] = i2;
                    i2 += 1usize;
                    this_nd.dof_index[5] = i2;
                    i2 += 1usize;
                }
            }
        }
        self.el_mat_dim = i2;
        self.elastic_mat.set_dim(self.el_mat_dim);
        self.therm_mat.set_dim(i3);
        self.diff_mat.set_dim(i3);
        self.fl_mat_dim = 6*i4;
        self.fluid_mat.set_dim(self.fl_mat_dim);
        self.fluid_lf.set_dim(self.fl_mat_dim);
        
        let scmd = &mut self.job[self.solve_cmd];
        if scmd.solver_method.s == "iterative" && scmd.max_it == 0 {
            if self.fl_mat_dim > self.el_mat_dim {
                scmd.max_it = self.fl_mat_dim;
            }
            else {
                scmd.max_it = self.el_mat_dim;
            }
        }
        if scmd.fluid && scmd.elastic {
            self.fsi_disp_map.set_dim(3*i4);
            self.fsi_temp_map.set_dim(self.fl_mat_dim);
            self.fluid_mesh_def.set_dim(3*i4);
        }
        
        for this_el in self.elements.iter_mut() {
            i3 = this_el.num_int_dof;
            if i3 > 0 {
                this_el.int_dof_index = i2;
                i2  +=  i3;
            }
        }
        self.tot_glob_dof = i2;
        i3 = self.nodes.len();

        self.elastic_ld_vec = vec![0f64; self.el_mat_dim];
        self.elastic_sol_vec = vec![0f64; self.el_mat_dim];
        self.therm_ld_vec = vec![0f64; i3];
        self.therm_sol_vec = vec![0f64; i3];
        self.diff_ld_vec = vec![0f64; i3];
        self.diff_sol_vec = vec![0f64; i3];
        self.fluid_ld_vec = vec![0f64; self.fl_mat_dim];
        self.fluid_sol_vec = vec![0f64; self.fl_mat_dim];
        
        self.temp_v1 = vec![0f64; self.el_mat_dim];
        self.temp_v2 = vec![0f64; self.el_mat_dim];
        self.temp_v3 = vec![0f64; self.el_mat_dim];
        self.temp_d1 = vec![DiffDoub0::new(); self.el_mat_dim];
        
        self.d_ld_u = vec![0f64; self.tot_glob_dof];
        self.d_ld_v = vec![0f64; self.el_mat_dim];
        self.d_ld_a = vec![0f64; self.el_mat_dim];
        self.d_ld_t = vec![0f64; i3];
        self.d_ld_tdot = vec![0f64; i3];
        self.d_ld_con = vec![0f64; i3];
        self.d_ld_condot = vec![0f64; i3];

        self.u_adj = vec![0f64; self.el_mat_dim];
        self.v_adj = vec![0f64; self.el_mat_dim];
        self.a_adj = vec![0f64; self.el_mat_dim];
        self.t_adj = vec![0f64; i3];
        self.tdot_adj = vec![0f64; i3];
        self.con_adj = vec![0f64; i3];
        self.condot_adj = vec![0f64; i3];
        
        self.d_rud_d = vec![DiffDoub1::new(); self.el_mat_dim];
        self.d_rtd_d = vec![DiffDoub1::new(); i3];
        self.d_rdd_d = vec![DiffDoub1::new(); i3];
        
        self.el_in_d = vec![0usize; self.elements.len()];
        
        i3 = self.design_vars.len();
        if i3 > 0 {
            self.d_ld_d = vec![0f64; i3];
        }
        
        return;
    }

    pub fn build_constraint_mats(&mut self) {
        self.elastic_const.build_all_mats(&self.nodes, &self.node_sets);
        self.thermal_const.build_all_mats(&self.nodes, &self.node_sets);
        self.diff_const.build_all_mats(&self.nodes, &self.node_sets);
        self.fluid_const.build_all_mats(&self.nodes, &self.node_sets);
    }

    pub fn update_reference(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        // Set the Section pointers for all self.elements and Material pointers for all self.sections
        let mut el_set : CppStr;
        let mut mat_name : CppStr;
        i2 = 0;
        for this_sec in self.sections.iter_mut() {
            el_set = this_sec.el_set_name.clone();
            i1 = self.es_map.at(&el_set.to_string());
            for eli in self.element_sets[i1].labels.iter_mut() {
                self.elements[*eli].sect_ptr = i2;
            }
            i3 = 0;
            for this_mat in self.materials.iter_mut() {
                mat_name = this_mat.name.clone();
                if mat_name.s == this_sec.mat_name.s {
                    this_sec.mat_ptr = i3;
                }
                for this_lay in this_sec.layers.iter_mut() {
                    if mat_name.s == this_lay.mat_name.s {
                        this_lay.mat_ptr = i3;
                    }
                }
                i3 += 1usize;
            }
            i2 += 1usize;
        }
        
        // Set Node & Element Set pointers in loads, constraints, design variables and objectives
        let mut nd_set : CppStr;
        for this_load in self.elastic_loads.iter_mut() {
            nd_set = this_load.node_set.clone();
            if CppMap::key_in_map(&mut self.ns_map, &nd_set.s) {
                i1 = self.ns_map.at(&nd_set.to_string());
                this_load.nd_set_ptr = i1;
            }
            else {
                el_set = this_load.element_set.clone();
                i1 = self.es_map.at(&el_set.to_string());
                this_load.el_set_ptr = i1;
            }
        }
        for this_load in self.thermal_loads.iter_mut() {
            nd_set = this_load.node_set.clone();
            if CppMap::key_in_map(&mut self.ns_map, &nd_set.s) {
                i1 = self.ns_map.at(&nd_set.to_string());
                this_load.nd_set_ptr = i1;
            }
            else {
                el_set = this_load.element_set.clone();
                i1 = self.es_map.at(&el_set.to_string());
                this_load.el_set_ptr = i1;
            }
        }
        
        for this_const in self.elastic_const.const_vec.iter_mut() {
            for this_cterm in this_const.terms.iter_mut() {
                nd_set = this_cterm.node_set.clone();
                i1 = self.ns_map.at(&nd_set.to_string());
                this_cterm.ns_ptr = i1;
            }
        }
        for this_const in self.thermal_const.const_vec.iter_mut() {
            for this_cterm in this_const.terms.iter_mut() {
                nd_set = this_cterm.node_set.clone();
                i1 = self.ns_map.at(&nd_set.to_string());
                this_cterm.ns_ptr = i1;
            }
        }
        
        for this_dv in self.design_vars.iter_mut() {
            nd_set = this_dv.nd_set_name.clone();
            if CppMap::key_in_map(&mut self.ns_map, &nd_set.s) {
                i1 = self.ns_map.at(&nd_set.to_string());
                this_dv.nd_set_ptr = i1;
            }
            else {
                el_set = this_dv.el_set_name.clone();
                i1 = self.es_map.at(&el_set.to_string());
                this_dv.el_set_ptr = i1;
            }
        }
        
        for this_term in self.obj.terms.iter_mut() {
            el_set = this_term.el_set_name.clone();
            if CppMap::key_in_map(&mut self.es_map, &el_set.s) {
                i1 = self.es_map.at(&el_set.s);
                this_term.el_set_ptr = i1;
            }
            else {
                nd_set = this_term.nd_set_name.clone();
                i1 = self.ns_map.at(&nd_set.s);
                this_term.nd_set_ptr = i1;
            }
        }
        
        // build dv reference list for self.nodes and self.elements
        let mut coef_len : usize;
        let mut const_coef : f64;
        let mut dvi : usize = 0;
        for this_dv in self.design_vars.iter_mut() {
            coef_len = this_dv.coefs.len();
            if this_dv.el_set_ptr < MAX_INT {
                if coef_len < 2 {
                    if coef_len == 0 {
                        const_coef = 1.0;
                    } else {
                        const_coef = match this_dv.coefs.front() {
                            None => 1.0,
                            Some(x) => *x,
                        };
                    }
                    for eli in self.element_sets[this_dv.el_set_ptr].labels.iter_mut() {
                        self.elements[*eli].add_design_variable(dvi, const_coef);
                    }
                }
                else {
                    let set_labs : &LinkedList<usize> = &self.element_sets[this_dv.el_set_ptr].labels;
                    let mut coef_iter = this_dv.coefs.iter();
                    let mut c_val : f64;
                    for set_iter in set_labs.iter() {
                        c_val = match coef_iter.next() {
                            None => panic!("Error: coefficient list provided for a design variable does not match the length of the element set"),
                            Some(x) => *x,
                        };
                        self.elements[*set_iter].add_design_variable(dvi, c_val);
                    }
                }
            }
            
            if this_dv.nd_set_ptr < MAX_INT {
                if coef_len < 2 {
                    if coef_len == 0 {
                        const_coef = 1.0;
                    }
                    else {
                        const_coef = match this_dv.coefs.front() {
                            None => 1.0,
                            Some(x) => *x,
                        };
                    }
                    for ndi in self.node_sets[this_dv.nd_set_ptr].labels.iter_mut() {
                        self.nodes[*ndi].add_design_variable(dvi, const_coef);
                    }
                }
                else {
                    let set_labs : &LinkedList<usize> = &self.node_sets[this_dv.nd_set_ptr].labels;
                    let mut coef_iter = this_dv.coefs.iter();
                    let mut c_val : f64;
                    for set_iter in set_labs.iter() {
                        c_val = match coef_iter.next() {
                            None => panic!("Error: coefficient list provided for a design variable does not match the length of the node set"),
                            Some(x) => *x,
                        };
                        self.nodes[*set_iter].add_design_variable(dvi, c_val);
                    }
                }
            }
            dvi += 1usize;
        }
        
        // build comprehensive Element list for each design variable
        
        let mut el_label : usize;
        let mut el_num_nds : usize;
        let mut this_nd : & Node;
        let mut this_dv : &mut DesignVariable;
        let mut tmp_vec = vec![IDCapsule::new(); self.design_vars.len()];
        for this_el in self.elements.iter_mut() {
        //for eli in 0..self.elements.len() {
            el_label = this_el.label;
            i1 = 0usize;
            for dv in this_el.design_vars.iter() {
                tmp_vec[i1] = dv.clone();
                i1 += 1usize;
            }
            for i2 in 0..i1 {
                self.design_vars[tmp_vec[i2].int_dat].add_comp_el(el_label);
                this_el.add_comp_dvar(tmp_vec[i2].int_dat);
            }
            el_num_nds = this_el.num_nds;
            for i1 in 0..el_num_nds {
                this_nd = & self.nodes[this_el.nodes[i1]];
                for dv in this_nd.d_var_lst.iter() {
                    this_dv = &mut self.design_vars[dv.int_dat];
                    if this_dv.category.s == "nodeCoord" {
                        this_dv.add_comp_el(el_label);
                        this_el.add_comp_dvar(dv.int_dat);
                    }
                }
            }
        }
        
        return;
    }

    pub fn find_surface_faces(&mut self) {
        let mut i1 : usize;
        let mut low_nd : usize;
        let mut _added : bool;
        
        let mut fc_ct : usize =  0;
        for el in self.elements.iter_mut() {
            fc_ct  +=  el.num_faces;
        }
        
        self.faces = vec![Face::new(); fc_ct];
        
        i1 = 0;
        for el in self.elements.iter_mut() {
            el.initialize_faces(&mut self.faces, &mut i1);
        }
        
        let num_nodes : usize =  self.nodes.len();
        let mut f_larray = vec![FacePtList::new(); num_nodes];
        let mut proceed : bool;
        for this_el in self.elements.iter_mut() {
            proceed = match this_el.this_type {
                3 => false,
                41 => false,
                _ => true,
            };
            if proceed {
                for this_fc in this_el.faces.iter_mut() {
                    low_nd = self.faces[*this_fc].get_low_nd();
                    _added = f_larray[low_nd].add_if_absent(*this_fc, &mut self.faces);
                }
            }
        }

        for this_face in self.faces.iter() {
            if this_face.on_surf {
                for ndi in 0..this_face.num_nds {
                    i1 = this_face.glob_nodes[ndi];
                    self.nodes[i1].on_surf = true;
                }
            }
        }
        
        return;
    }

    pub fn prep_matrix_factorizations(&mut self) {
        let mut zero_ar : [f64; 9] = [ 0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ];
        
        //let mut scmd = &mut self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        
        if self.job[sci].thermal {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_temp();
                if self.job[sci].dynamic {
                    this_nd.update_tdot(self.job[sci].newmark_gamma,  self.job[sci].time_step);
                }
            }
            if !self.therm_lt.is_allocated() {
                self.build_thermal_soln_load(true);
                self.thermal_const.update_active_status(0.0);
                self.therm_lt.allocate_from_sparse_mat(&mut self.therm_mat,  &mut self.thermal_const, self.job[sci].solver_block_dim);
            }
            if !self.therm_scaled {
                self.therm_scaled = Model::scale_const(&mut self.thermal_const, &self.therm_mat);
            }
        }

        if self.job[sci].diffusion {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_fl_den();
                if self.job[sci].dynamic {
                    this_nd.update_fl_den_dot(self.job[sci].newmark_gamma, self.job[sci].time_step);
                }
            }
            if !self.diff_lt.is_allocated() {
                self.build_diff_soln_load(true);
                self.diff_const.update_active_status(0.0);
                self.diff_lt.allocate_from_sparse_mat(&mut self.diff_mat, &mut self.diff_const, self.job[sci].solver_block_dim);
            }
            if !self.diff_scaled {
                self.diff_scaled = Model::scale_const(&mut self.diff_const, &self.diff_mat);
            }
        }
        
        if self.job[sci].elastic {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_disp();
                if self.job[sci].dynamic {
                    this_nd.update_vel_acc(self.job[sci].newmark_beta,  self.job[sci].newmark_gamma,  self.job[sci].time_step);
                }
            }
            for this_el in self.elements.iter_mut() {
                this_el.set_int_disp(&mut zero_ar);
                this_el.set_int_prev_disp(&mut zero_ar);
            }
            if !self.elastic_lt.is_allocated() {
                self.build_elastic_soln_load(true);
                self.elastic_const.update_active_status(0.0);
                self.elastic_lt.allocate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const, 6 * self.job[sci].solver_block_dim);
            }
            if !self.elastic_scaled {
                self.elastic_scaled = Model::scale_const(&mut self.elastic_const, &self.elastic_mat);
            }
        }
        
        return;
    }

    pub fn build_lf_mat(&mut self) {
        let mut dist : f64;
        let mut coef : f64;
        let mut row : usize;
        let mut col : usize;
        self.fluid_lf.zero_all();
        let diss_lev = self.job[self.solve_cmd].dissipation;
        let diag_coef = 1.0 - diss_lev;
        for nd in self.nodes.iter() {
            if nd.fluid {
                row = 6*nd.sorted_rank;
                for ni in nd.conn_nds.iter() {
                    col = 6*self.nodes[*ni].sorted_rank;
                    dist = get_dist(&nd.coord, &self.nodes[*ni].coord);
                    coef = 1.0/dist;
                    for i in 0..6 {
                        self.fluid_lf.add_entry(row+i,col+i,coef);
                    }
                }
                for i in 0..6 {
                    self.fluid_lf.scale_row_to_sum(row+i, diss_lev);
                    self.fluid_lf.add_entry(row+i, row+i, diag_coef);
                }
            }
        }
    }

    pub fn build_mesh_def_mat(&mut self) {
        self.fluid_mesh_def.zero_all();

        let mut el10 = Element::new();
        el10.initialize_type(10);
        let mut el8 = Element::new();
        el8.initialize_type(8);
        let mut el6 = Element::new();
        el6.initialize_type(6);
        let mut el4 = Element::new();
        el4.initialize_type(4);

        let r_vec = match self.d0_scratch.front_mut() {
            None => panic!("Error: no scratch vector, build_mesh_def_mat()"),
            Some(x) => &mut x.dat,
        };
        let mut dr_du : &mut Vec<f64> = &mut Vec::new();
        let mut dr_dt : &mut Vec<f64> = &mut Vec::new();
        let mut dr_dc : &mut Vec<f64> = &mut Vec::new();
        let mut si = 0usize;
        for s in self.scratch.iter_mut() {
            if si == 0 {
                dr_du = &mut s.dat;
            }
            else if si == 1 {
                dr_dt = &mut s.dat;
            }
            else if si == 2 {
                dr_dc = &mut s.dat;
            }
            si += 1;
        }

        let mut st_el : &mut Element;
        let mut gr : usize;
        let mut gc : usize;
        let mut k : usize;
        for el in self.elements.iter() {
            if el.this_type > 100 {
                el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut self.sections, &mut self.materials, &mut self.nodes, &self.design_vars);
                st_el = match el.this_type {
                    1000 => &mut el10,
                    800 => &mut el8,
                    600 => &mut el6,
                    400 => &mut el4,
                    _ => panic!("Error: unrecognized element type, build_mesh_def_mat"),
                };
                for i in 0..el.num_nds {
                    st_el.nodes[i] = el.nodes[i];
                }
                st_el.get_ruk_dfd0(r_vec, dr_du, dr_dt, dr_dc, true, false, &mut self.d0_pre);
                k = 0;
                for i1 in el.nodes.iter() {
                    gr = 3*self.nodes[*i1].sorted_rank;
                    for _i2 in 0..3 {
                        for j1 in el.nodes.iter() {
                            gc = 3*self.nodes[*j1].sorted_rank;
                            for _j2 in 0..3 {
                                self.fluid_mesh_def.add_entry(gr,gc, dr_du[k]);
                                k += 1;
                                gc += 1;
                            }
                        }
                        gr += 1;
                    }
                }
            }
        }

        // build the surface constraints for the mesh deformation
        self.mesh_def_const.const_vec = vec![Constraint::new()];
        let new_const = &mut self.mesh_def_const.const_vec[0];
        new_const.this_type = CppStr::from("displacement");
        new_const.scale_fact = 10000f64*self.fluid_mesh_def.get_max_abs_val();
        k = 0;
        for nd in self.nodes.iter() {
            if nd.fluid && nd.on_surf {
                k += 3;
            }
        }
        new_const.mat.set_dim(k);
        k = 0;
        for nd in self.nodes.iter() {
            if nd.fluid && nd.on_surf {
                gr = 3*nd.sorted_rank;
                for _i in 0..3 {
                    new_const.mat.add_entry(k,gr, 1.0);
                    gr += 1;
                    k += 1;
                }
            }
        }

        self.mesh_def_lt.allocate_from_sparse_mat(&mut self.fluid_mesh_def, &mut self.mesh_def_const, 3);
        self.mesh_def_lt.populate_from_sparse_mat(&mut self.fluid_mesh_def, &mut self.mesh_def_const);
        self.mesh_def_lt.ldl_factor();

    }

    pub fn build_fsi_map(&mut self) {

        //put all the structrual surface faces into a spatial grid list
        let mut x_rng = [1.0e+100f64, -1.0e+100];
        let mut y_rng = [1.0e+100f64, -1.0e+100];
        let mut z_rng = [1.0e+100f64, -1.0e+100];
        let mut dist : f64;
        let mut tot_dist = 0.0f64;
        let mut hit_ct = 0usize;
        for nd in self.nodes.iter() {
            if !nd.fluid {
                for ni in nd.conn_nds.iter() {
                    dist = get_dist(&nd.coord, &self.nodes[*ni].coord);
                    tot_dist += dist;
                    hit_ct += 1;
                }
                if nd.coord[0] < x_rng[0] {
                    x_rng[0] = nd.coord[0];
                }
                if nd.coord[0] > x_rng[1] {
                    x_rng[1] = nd.coord[0];
                }
                if nd.coord[1] < y_rng[0] {
                    y_rng[0] = nd.coord[1];
                }
                if nd.coord[1] > y_rng[1] {
                    y_rng[1] = nd.coord[1];
                }
                if nd.coord[2] < z_rng[0] {
                    z_rng[0] = nd.coord[2];
                }
                if nd.coord[2] > z_rng[1] {
                    z_rng[1] = nd.coord[2];
                }
            }
        }
        let avg_dist = tot_dist/(hit_ct as f64);
        let spacing = 2.0*avg_dist;

        let mut fc_gd_lst = SpatialGrid::new();
        fc_gd_lst.initialize(&mut x_rng, spacing, &mut y_rng, spacing, &mut z_rng, spacing);

        let mut cent = [0f64; 3];
        let mut fi = 0usize;
        let mut sfct = 0usize;
        for fc in self.faces.iter() {
            if fc.on_surf && self.elements[fc.host_el].this_type < 100 {
                fc.get_centroid(&mut cent, &self.nodes);
                fc_gd_lst.add_ent(fi, &cent);
                sfct += 1;
            }
            fi += 1;
        }

        // map the fluid surface points to structural faces, build fsi matrices
        let mut max_gap = self.job[self.solve_cmd].max_fsi_gap;
        if max_gap < 0.0 {
            max_gap = avg_dist;
        }
        let mut near_fcs = vec![0usize; sfct];
        let mut scrd = [0f64, 0f64];
        let mut min_dist : f64;
        let mut min_fc = 0usize;
        let mut min_scrd = [0f64, 0f64];
        let mut lst_ln : usize;
        let mut fci : usize;
        let mut n_vec = [DiffDoub0::new(); 8];
        let mut this_fc : &Face;
        let mut row : usize;
        let mut col : usize;
        let mut num_surf_nds = 0usize;
        for nd in self.nodes.iter() {
            if nd.on_surf && nd.fluid {
                num_surf_nds += 1;
                lst_ln = fc_gd_lst.get_in_radius(&mut near_fcs, sfct, &nd.coord, spacing);
                min_dist = 1.0e+100;
                for i in 0..lst_ln {
                    fci = near_fcs[i];
                    dist = self.faces[fci].get_proj_dist(&mut scrd, &nd.coord, &self.nodes);
                    if dist < min_dist {
                        min_dist = dist;
                        min_fc = fci;
                        min_scrd[0] = scrd[0];
                        min_scrd[1] = scrd[1];
                    }
                }
                if min_dist < max_gap {
                    this_fc = &self.faces[min_fc];
                    this_fc.get_basis_dfd0(&mut n_vec, &min_scrd);
                    row = 3*nd.sorted_rank;
                    for i in 0..3 {
                        for j in 0..this_fc.num_nds {
                            col = self.nodes[this_fc.glob_nodes[j]].dof_index[i];
                            self.fsi_disp_map.add_entry(row, col, n_vec[j].val);
                        }
                        row += 1;
                    }
                    row = 6*nd.sorted_rank + 4; //temperature dof of nd in the fluid eqns
                    for j in 0..this_fc.num_nds {
                        col = self.nodes[this_fc.glob_nodes[j]].sorted_rank;
                        self.fsi_temp_map.add_entry(row, col, n_vec[j].val);
                    }
                }
            }
        }

        //Build the constraint matrices for velocity and temperature on the fluid surface nodes

        let num_fl_cn = self.fluid_const.const_vec.len();
        let mut this_con = &mut self.fluid_const.const_vec[num_fl_cn - 2];
        this_con.mat.set_dim(3*num_surf_nds);
        let mut veli = 0usize;
        for nd in self.nodes.iter() {
            if nd.on_surf && nd.fluid {
                col = 6*nd.sorted_rank;
                for i in 1..4 {
                    this_con.mat.add_entry(veli, col+i, 1.0f64);
                    veli += 1;
                }
            }
        }

        this_con = &mut self.fluid_const.const_vec[num_fl_cn - 1];
        this_con.mat.set_dim(num_surf_nds);
        let mut tempi = 0usize;
        for nd in self.nodes.iter() {
            if nd.on_surf && nd.fluid {
                col = 6*nd.sorted_rank + 4;
                this_con.mat.add_entry(tempi, col, 1.0f64);
                tempi += 1;
            }
        }

    }

    pub fn analysis_prep(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let num_nds : usize;
        let block_dim : usize;
        
        //check options for solve command
        if self.solve_cmd < MAX_INT {
            let scmd = &mut self.job[self.solve_cmd];
            if scmd.solver_method.s == "direct" {
                scmd.solver_block_dim = 2000000000;
            }
            else {
                i1 = 6;
                num_nds = self.nodes.len();
                while (i1 * i1) < num_nds {
                    i1  +=  6;
                }
                if scmd.solver_block_dim == 2000000000 {
                    scmd.solver_block_dim = i1;
                }
            }
            block_dim = scmd.solver_block_dim;
            
            if scmd.static_load_time.len() == 0 {
                scmd.static_load_time.push_back(0.0);
            }
        }
        else {
            block_dim = 2000000000;
        }
        
        //allocate layers in stress prerequisite objects
        i1 = 0;
        for sec in self.sections.iter_mut() {
            i2 = sec.layers.len();
            if i2 > i1 {
                i1 = i2;
            }
        }
        if i1 > 0 {
            self.d0_pre.allocate_layers_dfd0(i1);
            self.d1_pre.allocate_layers_dfd1(i1);
        }

        //calculate initial values for nodal coordinates as a function of design variables

        for nd in self.nodes.iter_mut() {
            nd.calc_crd_dfd0(&self.design_vars);
        }
        
        //additional preparatory functions
        self.update_reference();
        self.reorder_nodes(block_dim);
        self.build_constraint_mats();
        self.prep_matrix_factorizations();
        self.find_surface_faces();

        if self.job[self.solve_cmd].fluid {
            self.build_lf_mat();
            if self.job[self.solve_cmd].elastic {
                self.build_mesh_def_mat();
                self.build_fsi_map();
            }
            else if self.job[self.solve_cmd].thermal {
                self.build_fsi_map();
            }
        }
        
        self.an_prep_run = true;
        return;
    }

}