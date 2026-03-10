use crate::model::*;
use crate::matrix_functions::*;
use crate::cpp_str::*;
use crate::cpp_map::*;

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
        let mut i2 : usize;
        let mut i3 : usize;
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
            el_num_nds = el.num_nds();
            el_dof_per_nd = el.dof_per_nd();
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
            }
        }
        
        Model::add_constraint_conn(&mut nodal_conn, &self.elastic_const, &self.node_sets);
        Model::add_constraint_conn(&mut nodal_conn, &self.thermal_const, &self.node_sets);
        Model::add_constraint_conn(&mut nodal_conn, &self.diff_const, &self.node_sets);
        
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
        let mut this_nd : &mut Node;
        for ndi in ordered_nds.iter_mut() {
            this_nd = &mut self.nodes[*ndi];
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
        self.el_mat_dim = i2;
        self.elastic_mat.set_dim(self.el_mat_dim);
        self.therm_mat.set_dim(i3);
        self.diff_mat.set_dim(i3);
        
        let scmd = &mut self.job[self.solve_cmd];

        if scmd.explicit {
            if scmd.thermal {
                for r in 0..i3 {
                    self.therm_mat.add_entry(r, r, 0.0);
                }
            }
            if scmd.diffusion {
                for r in 0..i3 {
                    self.diff_mat.add_entry(r, r, 0.0);
                }
            }
            if scmd.elastic {
                for r in 0..i2 {
                    self.elastic_mat.add_entry(r, r, 0.0);
                }
            }
        }

        if scmd.solver_method.s == "iterative" && scmd.max_it == 0 {
            scmd.max_it = self.el_mat_dim;
        }
        
        for this_el in self.elements.iter_mut() {
            i3 = this_el.num_int_dof();
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

        if self.init_stat_file.s != "" {
            self.read_initial_state(&mut self.init_stat_file.clone());
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
            el_num_nds = this_el.num_nds();
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

         // initialize structures needed for interactions

         self.interactions.initialize(&self.nodes, &self.node_sets, &self.ns_map, &self.elements, &self.design_vars);

         // initialize references for particle sources
 
         if !self.particle_sources.is_empty() {
             for sc in self.particle_sources.iter_mut() {
                 sc.elset_pt = self.es_map.at(&sc.element_set.s);
                 for i in 0..3 {
                     i1 = self.ns_map.at(&sc.ref_nodes[i].s);
                     if i1 < MAX_INT {
                         sc.ref_nodes_i[i] = match self.node_sets[i1].labels.front() {
                             None => MAX_INT,
                             Some(x) => *x,
                         }
                     } 
                 }
                 sc.swapset_pt = self.es_map.at(&sc.swap_set.s);
                 if sc.swapset_pt < MAX_INT {
                    i1 = sc.refine_lev;
                    sc.insert_els = vec![MAX_INT; i1*i1*i1];
                 }
                 sc.deact_ob_els(&mut self.elements, &mut self.nodes, &self.element_sets);
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
            fc_ct  +=  el.num_faces();
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

        if self.job[sci].const_scale_factor < 0.0 {
            self.job[sci].const_scale_factor = match self.job[sci].explicit {
                true => 1.0,
                false => 100000.0,
            };
        }
        
        if self.job[sci].thermal {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_temp(self.job[sci].time_step);
                if self.job[sci].dynamic {
                    this_nd.update_tdot(self.job[sci].newmark_gamma,  self.job[sci].time_step, self.job[sci].explicit);
                }
            }
            if !self.therm_lt.is_allocated() {
                self.build_thermal_soln_load(true, 0.0);
                self.thermal_const.update_active_status(0.0);
                self.therm_scaled = Model::scale_const(&mut self.thermal_const, &self.therm_mat, self.job[sci].const_scale_factor);
                self.thermal_const.add_to_sparse_mat(&mut self.therm_mat);
                self.therm_lt.allocate_from_sparse_mat(&mut self.therm_mat, self.job[sci].solver_block_dim);
            }
        }

        if self.job[sci].diffusion {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_fl_den(self.job[sci].time_step);
                if self.job[sci].dynamic {
                    this_nd.update_fl_den_dot(self.job[sci].newmark_gamma, self.job[sci].time_step, self.job[sci].explicit);
                }
            }
            if !self.diff_lt.is_allocated() {
                self.build_diff_soln_load(true);
                self.diff_const.update_active_status(0.0);
                self.diff_scaled = Model::scale_const(&mut self.diff_const, &self.diff_mat, self.job[sci].const_scale_factor);
                self.diff_const.add_to_sparse_mat(&mut self.diff_mat);
                self.diff_lt.allocate_from_sparse_mat(&mut self.diff_mat, self.job[sci].solver_block_dim);
            }
        }
        
        if self.job[sci].elastic {
            for this_nd in self.nodes.iter_mut() {
                this_nd.initialize_disp(self.job[sci].time_step);
                if self.job[sci].dynamic {
                    this_nd.update_vel_acc(self.job[sci].newmark_beta,  self.job[sci].newmark_gamma,  self.job[sci].time_step, self.job[sci].explicit);
                }
            }
            for this_el in self.elements.iter_mut() {
                this_el.set_int_disp(&mut zero_ar);
                this_el.set_int_prev_disp(&mut zero_ar);
            }
            if !self.elastic_lt.is_allocated() {
                self.build_elastic_soln_load(true, 0.0);
                self.elastic_const.update_active_status(0.0);
                self.elastic_scaled = Model::scale_const(&mut self.elastic_const, &self.elastic_mat, self.job[sci].const_scale_factor);
                self.elastic_const.add_to_sparse_mat(&mut self.elastic_mat);
                self.elastic_lt.allocate_from_sparse_mat(&mut self.elastic_mat, 6*self.job[sci].solver_block_dim);
            }
        }
        
        return;
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
            
            if !self.interactions.int_vec.is_empty() && !scmd.nonlinear_geom {
                println!("Warning: the presence of active interactions inherently requires nonlinear analysis.  Switching nonlinear option to 'yes'");
                scmd.nonlinear_geom = true;
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
        
        self.an_prep_run = true;
        return;
    }

}