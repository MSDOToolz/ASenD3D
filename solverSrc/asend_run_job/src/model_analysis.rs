use crate::model::*;
use crate::constants::*;
use crate::list_ent::*;
use crate::nd_el_set::*;
use crate::node::*;
use crate::element::*;
use crate::design_var::*;
use crate::face::*;
use crate::diff_doub::*;
use crate::job::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;
use crate::fmath::*;

use std::collections::LinkedList;

impl Model {
    pub fn reorder_nodes(&mut self, block_dim : usize) {
        let mut i1 : usize;
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
                }
            }
        }
        
        for this_const in self.elastic_const.const_vec.iter() {
            for term1 in this_const.terms.iter() {
                let t1_labs = &self.node_sets[term1.ns_ptr].labels;
                i1 = t1_labs.len();
                for term2 in this_const.terms.iter() {
                    let t2_labs = &self.node_sets[term2.ns_ptr].labels;
                    i2 = t2_labs.len();
                    if i1 > 1 && i2 > 1 {
                        let mut iter2 = t2_labs.iter();
                        let mut i2val : usize;
                        for iter1 in t1_labs.iter() {
                            i2val = match iter2.next() {
                                None => panic!("Error: sets of mismatched size found in constraint definition"),
                                Some(x) => *x,
                            };
                            nodal_conn[*iter1].add_if_absent(i2val);
                            nodal_conn[i2val].add_if_absent(*iter1);
                        }
                    }
                    else {
                        for iter1 in t1_labs.iter() {
                            for iter2 in t2_labs.iter() {
                                nodal_conn[*iter1].add_if_absent(*iter2);
                                nodal_conn[*iter2].add_if_absent(*iter1);
                            }
                        }
                    }
                }
            }
        }
        
        for this_const in self.thermal_const.const_vec.iter() {
            for term1 in this_const.terms.iter() {
                let t1_labs = &self.node_sets[term1.ns_ptr].labels;
                for iter1 in t1_labs.iter() {
                    for term2 in this_const.terms.iter() {
                        let t2_labs = &self.node_sets[term2.ns_ptr].labels;
                        for iter2 in t2_labs.iter() {
                            nodal_conn[*iter1].add_if_absent(*iter2);
                            nodal_conn[*iter2].add_if_absent(*iter1);
                        }
                    }
                }
            }
        }
        
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
            on_i += 1usize;
            nd1 = ordered_nds[on_i];
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
            }
        }
        
        // update the global degree of freedom indexes for the self.nodes
        
        i2 = 0;// index in elastic matrix
        i3 = 0;// sorted rank for solid self.nodes
        let mut this_nd : &mut Node;
        for ndi in ordered_nds.iter_mut() {
            this_nd = &mut self.nodes[*ndi];
            if !this_nd.fluid {
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
        self.non_frc_el_mat.set_dim(self.el_mat_dim);
        self.therm_mat.set_dim(self.nodes.len());
        
        let scmd = &mut self.job[self.solve_cmd];
        if scmd.solver_method.s == "iterative" && scmd.max_it == 0 {
            scmd.max_it = self.el_mat_dim;
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
        
        self.temp_v1 = vec![0f64; self.el_mat_dim];
        self.temp_v2 = vec![0f64; self.el_mat_dim];
        self.temp_v3 = vec![0f64; self.el_mat_dim];
        self.temp_v4 = vec![0f64; self.el_mat_dim];
        self.temp_v5 = vec![0f64; self.el_mat_dim];
        self.temp_v6 = vec![0f64; self.el_mat_dim];
        self.temp_d1 = vec![DiffDoub0::new(); self.el_mat_dim];
        
        self.d_ld_u = vec![0f64; self.tot_glob_dof];
        self.d_ld_v = vec![0f64; self.el_mat_dim];
        self.d_ld_a = vec![0f64; self.el_mat_dim];
        self.d_ld_t = vec![0f64; i3];
        self.d_ld_tdot = vec![0f64; i3];
        self.u_adj = vec![0f64; self.el_mat_dim];
        self.v_adj = vec![0f64; self.el_mat_dim];
        self.a_adj = vec![0f64; self.el_mat_dim];
        self.t_adj = vec![0f64; i3];
        self.tdot_adj = vec![0f64; i3];
        
        self.d_rud_d = vec![DiffDoub1::new(); self.el_mat_dim];
        self.d_rtd_d = vec![DiffDoub1::new(); i3];
        
        self.el_in_d = vec![0usize; self.elements.len()];
        
        i3 = self.design_vars.len();
        if i3 > 0 {
            self.d_ld_d = vec![0f64; i3];
        }
        
        return;
    }

    pub fn build_constraint_mats(&mut self) {
        for this_const in self.elastic_const.const_vec.iter_mut() {
            this_const.build_mat(&mut self.nodes, &mut self.node_sets);
        }
        for this_const in self.thermal_const.const_vec.iter_mut() {
            this_const.build_mat(&mut self.nodes, &mut self.node_sets);
        }
        return;
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
        for this_el in self.elements.iter_mut() {
            if this_el.dof_per_nd == 3 {
                for this_fc in this_el.faces.iter_mut() {
                    low_nd = self.faces[*this_fc].get_low_nd();
                    _added = f_larray[low_nd].add_if_absent(*this_fc, &mut self.faces);
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
                self.therm_lt.allocate_from_sparse_mat(&mut self.therm_mat,  &mut  self.thermal_const,   self.job[sci].solver_block_dim);
            }
            if !self.therm_scaled {
                self.scale_thermal_const();
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
                self.build_elastic_soln_load(true,  true);
                self.elastic_const.update_active_status(0.0);
                self.scale_elastic_const();
                self.elastic_lt.allocate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const,   6 * self.job[sci].solver_block_dim);
            }
        }
        
        return;
    }

    pub fn analysis_prep(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let num_nds : usize;
        let block_dim : usize;
        
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
        
        self.update_reference();
        self.reorder_nodes(block_dim);
        self.build_constraint_mats();
        self.prep_matrix_factorizations();
        self.find_surface_faces();
        
        self.an_prep_run = true;
        return;
    }

    pub fn build_elastic_app_load(&mut self, time : f64) {
        // construct loads from input file
        let mut num_dof : usize;
        let mut dof_ind : usize;
        let mut ld_type : CppStr;
        let mut nd_dvld = [DiffDoub0::new(); 6];
        let mut this_nd : &Node;
        let mut this_el : &Element;
        let mut scmd : &JobCommand;
        
        for i1 in 0..self.el_mat_dim {
            self.temp_d1[i1].set_val(0.0);
        }
        
        //loads from Model input file
        for this_load in self.elastic_loads.iter() {
            ld_type = this_load.this_type.clone();
            if time >= this_load.active_time[0] && time <= this_load.active_time[1] {
                if ld_type.s == "nodalForce" {
                    for ndi in self.node_sets[this_load.nd_set_ptr].labels.iter_mut() {
                        this_nd = &self.nodes[*ndi];
                        num_dof = this_nd.num_dof;
                        for i1 in 0..num_dof {
                            dof_ind = this_nd.dof_index[i1];
                            self.elastic_ld_vec[dof_ind]  +=  this_load.load[i1];
                        }
                    }
                }
                else {
                    for eli in self.element_sets[this_load.el_set_ptr].labels.iter_mut() {
                        this_el = &self.elements[*eli];
                        this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                        scmd = &self.job[self.solve_cmd];
                        this_el.get_app_load_dfd0(&mut self.temp_d1, this_load, scmd.nonlinear_geom, &mut self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.sections, &mut  self.faces, &mut  self.nodes, & self.design_vars);
                    }
                }
            }
        }
        
        for i1 in 0..self.el_mat_dim {
            self.elastic_ld_vec[i1]  +=  self.temp_d1[i1].val;
        }
        
        // design variable dependent loads.
        for this_nd in self.nodes.iter_mut() {
            this_nd.get_elastic_dvload_dfd0(&mut nd_dvld, &mut self.design_vars);
            num_dof = this_nd.num_dof;
            for i1 in 0..num_dof {
                dof_ind = this_nd.dof_index[i1];
                self.elastic_ld_vec[dof_ind]  +=  nd_dvld[i1].val;
            }
        }
        
        return;
    }

    pub fn build_thermal_app_load(&mut self, time : f64) {
        // construct loads from input file
        let tot_nodes : usize;
        let mut dof_ind : usize;
        
        
        let mut ld_type : CppStr;
        let nd_load : [f64; 6] = [0f64; 6];
        let mut nd_dvld = DiffDoub0::new();
        
        let mut this_nd : &Node;
        let mut this_el : &Element;
        
        tot_nodes = self.nodes.len();
        for i1 in 0..tot_nodes {
            self.temp_d1[i1].set_val(0.0);
        }
        
        //loads from Model input file
        for this_load in self.thermal_loads.iter() {
            ld_type = this_load.this_type.clone();
            if time >= this_load.active_time[0] && time <= this_load.active_time[1] {
                if ld_type.s == "nodalHeatGen" {
                    for ndi in self.node_sets[this_load.nd_set_ptr].labels.iter_mut() {
                        this_nd = &self.nodes[*ndi];
                        dof_ind = this_nd.sorted_rank;
                        self.therm_ld_vec[dof_ind]  +=  nd_load[0];
                    }
                }
                else {
                    for eli in self.element_sets[this_load.el_set_ptr].labels.iter_mut() {
                        this_el = &self.elements[*eli];
                        this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                        this_el.get_app_therm_load_dfd0(&mut self.temp_d1, this_load, &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.sections, &mut  self.faces, &mut  self.nodes, & self.design_vars);
                    }
                }
            }
        }
        
        for i1 in 0..tot_nodes {
            self.therm_ld_vec[i1]  +=  self.temp_d1[i1].val;
        }
        
        // design variable dependent loads.
        for this_nd in self.nodes.iter_mut() {
            this_nd.get_thermal_dvload_dfd0(&mut nd_dvld, & self.design_vars);
            dof_ind = this_nd.sorted_rank;
            self.therm_ld_vec[dof_ind]  +=  nd_dvld.val;
        }
        
        return;
    }

    pub fn build_elastic_soln_load(&mut self, build_mat : bool, full_ref : bool) {
        
        for i1 in 0..self.el_mat_dim {
            self.temp_d1[i1].set_val(0.0);
        }
        if build_mat {
            self.elastic_mat.zero_all();
        }
        if full_ref {
            self.non_frc_el_mat.zero_all();
        }
        
        let scmd = &self.job[self.solve_cmd];
        for this_el in self.elements.iter_mut() {
            this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
            if this_el.this_type == 21 {
                this_el.get_ru_dfd0(&mut self.temp_d1, &mut  self.elastic_mat,  build_mat, scmd, &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.nodes);
            }
            else {
                if full_ref {
                    this_el.get_ru_dfd0(&mut self.temp_d1, &mut  self.non_frc_el_mat,  true, scmd, &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.nodes);
                }
                else {
                    this_el.get_ru_dfd0(&mut self.temp_d1, &mut  self.non_frc_el_mat,  false, scmd, &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.nodes);
                }
            }
        }
        
        if build_mat {
            self.elastic_mat.add_matrix(&mut self.non_frc_el_mat);
        }
        
        for i1 in 0..self.el_mat_dim {
            self.elastic_ld_vec[i1] -=  self.temp_d1[i1].val;
        }
        
        return;
    }

    pub fn build_thermal_soln_load(&mut self, build_mat : bool) {
        let num_nodes : usize;
        
        num_nodes = self.nodes.len();
        for i1 in 0..num_nodes {
            self.temp_d1[i1].set_val(0.0);
        }
        
        if build_mat {
            self.therm_mat.zero_all();
        }
        
        let scmd = &self.job[self.solve_cmd];
        for this_el in self.elements.iter_mut() {
            this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
            this_el.get_rt_dfd0(&mut self.temp_d1, &mut  self.therm_mat,  build_mat, scmd, &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut(), &mut  self.nodes);
        }
        
        for i1 in 0..num_nodes {
            self.therm_ld_vec[i1]  -=  self.temp_d1[i1].val;
        }
        
        return;
    }

    pub fn scale_elastic_const(&mut self) {
        let scale_fact : f64 =  100000.0*self.elastic_mat.get_max_abs_val();
        self.elastic_const.set_scale_fact(scale_fact);
        self.elastic_scaled = true;
        return;
    }

    pub fn scale_thermal_const(&mut self) {
        let scale_fact : f64 =  100000.0 * self.therm_mat.get_max_abs_val();
        self.thermal_const.set_scale_fact(scale_fact);
        self.therm_scaled = true;
        return;
    }

    pub fn build_elastic_const_load(&mut self) {
        let mut nd_dof : usize;
        let mut glob_ind : usize;
        for this_nd in self.nodes.iter() {
            nd_dof = this_nd.num_dof;
            for i1 in 0..nd_dof {
                glob_ind = this_nd.dof_index[i1];
                self.temp_v2[glob_ind] = this_nd.displacement[i1];
            }
        }
        self.elastic_const.get_total_load(&mut self.elastic_ld_vec, &mut  self.temp_v2, &mut  self.temp_v3);
        return;
    }

    pub fn build_thermal_const_load(&mut self) {
        let mut nd_temp : f64;
        let mut glob_ind : usize;
        for this_nd in self.nodes.iter_mut() {
            nd_temp = this_nd.temperature;
            glob_ind = this_nd.sorted_rank;
            self.temp_v2[glob_ind] = nd_temp;
        }
        self.thermal_const.get_total_load(&mut self.therm_ld_vec, &mut  self.temp_v2, &mut  self.temp_v3);
        return;
    }

    pub fn solve_step(&mut self, time : f64, app_ld_fact : f64, full_ref : bool) {
        let mut i2 : usize;
        let num_nodes : usize;
        let max_nlit : usize;
        let mut d_unorm : f64;
        let d_utol : f64;
        let mut absd_u : f64;
        let mut ndof : usize;
        let mut dof_ind : usize;
        let mut nd_del_disp : [f64; 6] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        
        //let mut cmd = &mut self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        if self.job[sci].thermal {
            num_nodes = self.nodes.len();
            for i1 in 0..num_nodes {
                self.therm_ld_vec[i1] = 0.0;
                self.therm_sol_vec[i1] = 0.0;
            }
            self.build_thermal_app_load(time);
            self.build_thermal_soln_load(false);
            if !self.therm_scaled {
                self.scale_thermal_const();
            }
            self.thermal_const.update_active_status(time);
            self.build_thermal_const_load();
            if self.thermal_const.any_just_activated() {
                self.therm_lt.allocate_from_sparse_mat(&mut self.therm_mat,&mut self.thermal_const, self.job[sci].solver_block_dim);
                self.therm_lt.populate_from_sparse_mat(&mut self.therm_mat, &mut self.thermal_const);
                self.therm_lt.ldl_factor();
            }
            if self.job[sci].solver_method.s == "direct" {
                self.therm_lt.ldl_solve(&mut self.therm_sol_vec, &mut  self.therm_ld_vec);
            }
            else {
                conj_grad_sparse(&mut self.therm_sol_vec, &mut  self.therm_mat, &mut  self.thermal_const, &mut  self.therm_lt, &mut  self.therm_ld_vec,  self.job[sci].conv_tol,  self.job[sci].max_it);
                //g_mres_sparse(self.temp_v2, *self.therm_mat, *self.thermal_const, *self.therm_lt, self.temp_v1, cmd->conv_tol, cmd->max_it, cmd->solver_block_dim);
            }
            for this_nd in self.nodes.iter_mut() {
                dof_ind = this_nd.sorted_rank;
                this_nd.temperature  +=  self.therm_sol_vec[dof_ind];
                if self.job[sci].dynamic {
                    this_nd.update_tdot(self.job[sci].newmark_gamma,  self.job[sci].time_step);
                }
            }
        }
        
        if self.job[sci].elastic {
            if self.job[sci].nonlinear_geom {
                max_nlit = 50;
            } else {
                max_nlit = 1;
            }
            self.elastic_const.update_active_status(time);
            d_utol = 1.0e-12;
            d_unorm = 1.0;
            i2 = 0;
            while i2 < max_nlit && d_unorm > d_utol {
                for i1 in 0..self.el_mat_dim {
                    self.elastic_ld_vec[i1] = 0.0;
                    self.elastic_sol_vec[i1] = 0.0;
                }
                self.build_elastic_app_load(time);
                for i1 in 0..self.el_mat_dim {
                    self.elastic_ld_vec[i1]  *=  app_ld_fact;
                }
                self.build_elastic_soln_load(self.job[sci].nonlinear_geom, full_ref);
                if !self.elastic_scaled {
                    self.scale_elastic_const();
                }
                self.build_elastic_const_load();
                for this_el in self.elements.iter_mut() {
                    if this_el.num_int_dof > 0 {
                        this_el.update_external(&mut self.elastic_ld_vec, 1, &mut self.nodes, &mut self.scratch.iter_mut());
                    }
                }
                if self.elastic_const.any_just_activated() {
                    self.elastic_lt.allocate_from_sparse_mat(&mut self.elastic_mat, &mut self.elastic_const, self.job[sci].solver_block_dim);
                    self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
                    self.elastic_lt.ldl_factor();
                }
                else if self.job[sci].nonlinear_geom {
                    self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
                    self.elastic_lt.ldl_factor();
                }
                if self.job[sci].solver_method.s == "direct" {
                    self.elastic_lt.ldl_solve(&mut self.elastic_sol_vec, &mut self.elastic_ld_vec);
                }
                else {
                    conj_grad_sparse(&mut self.elastic_sol_vec, &mut  self.elastic_mat, &mut  self.elastic_const, &mut  self.elastic_lt, &mut  self.elastic_ld_vec,  self.job[sci].conv_tol,  self.job[sci].max_it);
                    //g_mres_sparse(self.temp_v2, *self.elastic_mat, *self.elastic_const, *self.elastic_lt, self.temp_v1, cmd->conv_tol, cmd->max_it, 6*cmd->solver_block_dim);
                }
                for this_nd in self.nodes.iter_mut() {
                    ndof = this_nd.num_dof;
                    for i1 in 0..ndof {
                        dof_ind = this_nd.dof_index[i1];
                        nd_del_disp[i1] = self.elastic_sol_vec[dof_ind];
                    }
                    this_nd.add_to_displacement(&mut nd_del_disp);
                    if self.job[sci].dynamic {
                        this_nd.update_vel_acc(self.job[sci].newmark_beta, self.job[sci].newmark_gamma, self.job[sci].time_step);
                    }
                }
                for this_el in self.elements.iter_mut() {
                    if this_el.num_int_dof > 0 {
                        this_el.update_internal(&mut self.elastic_sol_vec, 1, &mut self.nodes, &mut self.scratch.iter_mut());
                    }
                }
                d_unorm = 0.0;
                for i1 in 0..self.el_mat_dim {
                    absd_u = fabs(self.elastic_sol_vec[i1]);
                    if absd_u > d_unorm {
                        d_unorm = absd_u;
                    }
                }
                if self.job[sci].nonlinear_geom {
                    println!("{}{}{}{}", "Nonlinear iteration: " , i2 , ", max solution step: " , d_unorm );
                }
                i2 += 1usize;
            }
        }
        
        
        return;
    }

    pub fn solve(&mut self) {
        //let mut cmd = &mut self.job[self.solve_cmd];
        let ci : usize = self.solve_cmd;
        let mut i1 : usize;
        let mut i2 : usize;
        let mut since_ref : usize;
        let mut app_ld_fact : f64;
        let mut time : f64;
        let ld_steps : usize =  self.job[ci].load_ramp_steps;
        let c1 : f64 =  1.0/((ld_steps*ld_steps) as f64);
        
        if !self.an_prep_run {
            self.analysis_prep();
        }
        
        if self.job[ci].thermal {
            self.therm_lt.populate_from_sparse_mat(&mut self.therm_mat,  &mut  self.thermal_const);
            self.therm_lt.ldl_factor();
        }
        
        if self.job[ci].elastic && !self.job[ci].nonlinear_geom {
            self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
            println!("{}", "factoring stiffness matrix" );
            self.elastic_lt.ldl_factor();
            println!("{}", "finished factoring" );
        }
        
        if self.job[ci].dynamic {
            if self.job[ci].save_soln_hist {
                self.write_time_step_soln(0);
            }
            time = 0.0;
            i1 = 1;
            since_ref = 0;
            while time < self.job[ci].sim_period {
                since_ref += 1usize;
                if since_ref == self.job[ci].full_reform {
                    self.solve_step(time,  1.0,  true);
                    since_ref = 0;
                }
                else {
                    self.solve_step(time,  1.0,  false);
                }
                for this_nd in self.nodes.iter_mut() {
                    if self.job[ci].thermal {
                        this_nd.advance_temp();
                        this_nd.update_tdot(self.job[ci].newmark_gamma, self.job[ci].time_step);
                    }
                    if self.job[ci].elastic {
                        this_nd.advance_disp();
                        this_nd.update_vel_acc(self.job[ci].newmark_beta, self.job[ci].newmark_gamma, self.job[ci].time_step);
                    }
                }
                if self.job[ci].elastic {
                    for this_el in self.elements.iter_mut() {
                        this_el.advance_int_disp();
                    }
                }
                if self.job[ci].save_soln_hist {
                    self.write_time_step_soln(i1);
                    self.time_steps_saved = i1;
                    i1 += 1usize;
                }
                time  +=  self.job[ci].time_step;
            }
        } else {
            let lt_len = self.job[ci].static_load_time.len();
            let mut lt_vec = vec![0.0f64; lt_len];
            i1 = 0;
            for lt in self.job[ci].static_load_time.iter() {
                lt_vec[i1] = *lt;
                i1 += 1usize;
            }
            let mut this_ld : f64;
            //for this_ld in self.job[ci].static_load_time.iter() {
            for i in 0..lt_len {
                this_ld = lt_vec[i];
                if self.job[ci].nonlinear_geom {
                    for i1 in 0..ld_steps {
                        i2 = ld_steps - i1 - 1;
                        app_ld_fact = 1.0 - c1 * ((i2 * i2) as f64);
                        self.solve_step(this_ld, app_ld_fact, true);
                    }
                }
                else {
                    app_ld_fact = 1.0;
                    self.solve_step(this_ld, app_ld_fact, true);
                }
                if self.job[ci].save_soln_hist {
                    for this_nd in self.nodes.iter_mut() {
                        if self.job[ci].thermal {
                            this_nd.advance_temp();
                        }
                        if self.job[ci].elastic {
                            this_nd.advance_disp();
                        }
                    }
                    if self.job[ci].elastic {
                        for this_el in self.elements.iter_mut() {
                            this_el.advance_int_disp();
                        }
                    }
                    self.write_time_step_soln(i);
                    self.time_steps_saved = i + 1;
                }
                //i3 += 1usize;
            }
        }
        
        return;
    }

    pub fn zero_solution(&mut self, fields : &mut LinkedList<CppStr>) {
        let mut zero_ar : [f64; 9] = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut n_int_dof : usize;
        
        let mut disp_inc : bool =  false;
        for this_field in fields.iter_mut() {
            if this_field.s == "displacement" {
                disp_inc = true;
            }
        }
        
        for this_nd in self.nodes.iter_mut() {
            for this_field in fields.iter_mut() {
                if this_field.s == "temperature" {
                    this_nd.temperature = 0.0;
                }
                else if this_field.s == "tdot" {
                    this_nd.temp_change_rate = 0.0;
                }
                else if this_field.s == "displacement" {
                    this_nd.set_displacement(&mut zero_ar);
                }
                else if this_field.s == "velocity" {
                    this_nd.set_velocity(&mut zero_ar);
                }
                else if this_field.s == "acceleration" {
                    this_nd.set_acceleration(&mut zero_ar);
                }
            }
        }
        
        if disp_inc {
            for this_el in self.elements.iter_mut() {
                n_int_dof = this_el.num_int_dof;
                if n_int_dof > 0 {
                    this_el.set_int_disp(&mut zero_ar);
                }
            }
        }
        
        return;
    }

    pub fn eigen_solve(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut rvec = vec![DiffDoub0::new(); 33];
        let mut rv2 = vec![0f64; 33];
        let mut d_rd_a = vec![0f64; 1089];
        let mut zeros : [f64; 9] = [0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0];
        let mut shft : f64;
        let mut v_kv : f64;
        let mut m_avg : f64;
        let m_min : f64;
        
        //let mut cmd = &mut self.job[self.modal_cmd];
        //let mut scmd = &mut self.job[self.solve_cmd];
        let ci : usize = self.modal_cmd;
        let sci : usize = self.solve_cmd;
        
        if !self.an_prep_run {
            self.analysis_prep();
        }
        
        if self.eig_vecs.len() == 0 {
            i1 = self.job[ci].num_modes;
            self.eig_vals = vec![0f64; i1];
            self.load_fact = vec![0f64; i1];
            i1  *=  self.el_mat_dim;
            self.eig_vecs = vec![0f64; i1];
            self.diag_mass = vec![0f64; self.el_mat_dim];
        }
        
        if self.job[ci].this_type.s == "buckling" {
            for i1 in 0..self.el_mat_dim {
                self.diag_mass[i1] = 1.0;
            }
        }
        else {
            for i1 in 0..self.el_mat_dim {
                self.diag_mass[i1] = 0.0;
            }
            for this_el in self.elements.iter_mut() {
                this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                i2 = this_el.num_nds * this_el.dof_per_nd;
                for i1 in 0..i2 {
                    self.d0_pre.glob_acc[i1].set_val(1.0);
                }
                this_el.get_rum_dfd0(&mut rvec, &mut  d_rd_a,  false,  true,  self.job[sci].nonlinear_geom, &mut  self.d0_pre, &mut self.d0_scratch.iter_mut());
                for i1 in 0..i2 {
                    rv2[i1] = rvec[i1].val;
                }
                this_el.add_to_glob_vec(&mut rv2, &mut  self.diag_mass,  false,  false, &mut  self.nodes);
            }
        }
        
        m_avg = 0.0;
        for i1 in 0..self.el_mat_dim {
            m_avg  +=  fabs(self.diag_mass[i1]);
        }
        m_avg  /=  self.el_mat_dim as f64;
        m_min = 1.0e-6 * m_avg;
        for i1 in 0..self.el_mat_dim {
            if self.diag_mass[i1] < m_min {
                println!("{}{}", "Warning: small or negative value in diagonal mass matrix: " , self.diag_mass[i1] );
                println!("{}{}", "Adjusting to minimum value" , m_min );
                self.diag_mass[i1] = m_min;
            }
        }
        
        self.job[ci].nonlinear_geom = self.job[sci].nonlinear_geom;
        self.job[ci].dynamic = self.job[sci].dynamic;
        self.job[sci].nonlinear_geom = true;
        self.job[sci].dynamic = false;
        println!("{}", "building stiffness matrix" );
        self.build_elastic_soln_load(true,  true);
        println!("{}", "finished building matrix" );
        if self.job[ci].this_type.s == "buckling" || self.job[ci].this_type.s == "frequency" {
            for i1 in 0..self.el_mat_dim {
                shft = -self.job[ci].tgt_eval * self.diag_mass[i1];
                self.elastic_mat.add_entry(i1,   i1,   shft);
            }
            if self.job[ci].solver_method.s == "direct" {
                self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
                println!("{}", "factoring stiffness matrix" );
                self.elastic_lt.ldl_factor();
                println!("{}", "finished factoring.  beginning eigensolve" );
                eigen_sparse_direct(&mut self.eig_vals, &mut  self.eig_vecs,  self.job[ci].num_modes, &mut  self.elastic_lt, &mut  self.diag_mass,  self.el_mat_dim);
                println!("{}", "finished eigensolve" );
            }
            else {
                
            }
            for i1 in 0..self.job[ci].num_modes {
                self.eig_vals[i1]  +=  self.job[ci].tgt_eval;
            }
        }
        else {
            //step_inc = self.time_steps_saved / cmd.num_modes;
            //if (step_inc == 0) {
            //	string er_str = "Error: not enough solution time steps have been saved to find the requested number of active eigenmodes.\n";
            //	er_str += "Be sure to Set the saveSolnHist option to true in the dynamic solve command.\n";
            //	throw invalid_argument(er_str);
            //}
            //for (i1 = 0; i1 < cmd.num_modes; i1++) {
            //	i2 = (i1 + 1) * step_inc;
            //	i3 = i1 * self.el_mat_dim;
            //	read_time_step_soln(i2);
            //	this_nd = self.nodes.get_first();
            //	while (this_nd) {
            //		this_nd->get_prev_disp(nd_disp);
            //		n_dof = this_nd->get_num_dof();
            //		for (i4 = 0; i4 < n_dof; i4++) {
            //			glob_ind = this_nd->get_dof_index(i4);
            //			self.eig_vecs[i3 + glob_ind] = nd_disp[i4];
            //		}
            //		this_nd = this_nd->get_next();
            //	}
            //}
            //
            //get_nearest_evec_subspace(self.elastic_mat, self.elastic_const, self.diag_mass, self.eig_vecs, self.eig_vals, cmd->num_modes);
        }
        
        if self.job[ci].this_type.s == "buckling" {
            self.write_time_step_soln(MAX_INT);
            for this_nd in self.nodes.iter_mut() {
                this_nd.advance_disp();
                this_nd.set_displacement(&mut zeros);
            }
            for this_el in self.elements.iter_mut() {
                this_el.advance_int_disp();
                this_el.set_int_disp(&mut zeros);
            }
            self.build_elastic_soln_load(true,  true);
            i2 = 0;
            for i1 in 0..self.job[ci].num_modes {
                for i3 in 0..self.el_mat_dim {
                    self.temp_v1[i3] = 0.0;
                }
                sub_vec(&mut self.temp_v2, &mut  self.eig_vecs,  i2,  i2 + self.el_mat_dim);
                self.elastic_mat.vector_multiply(&mut self.temp_v1, &mut  self.temp_v2,  false);
                v_kv = 0.0;
                for i3 in 0..self.el_mat_dim {
                    v_kv  +=  self.eig_vecs[i2 + i3] * self.temp_v1[i3];
                }
                if (v_kv - self.eig_vals[i1]) != 0.0 {
                    self.load_fact[i1] = v_kv / (v_kv - self.eig_vals[i1]);
                }
                else {
                    self.load_fact[i1] = 1.0e+100;
                }
                i2  +=  self.el_mat_dim;
            }
            for this_nd in self.nodes.iter_mut() {
                this_nd.backstep_disp();
            }
            for this_el in self.elements.iter_mut() {
                this_el.backstep_int_disp();
            }
            self.read_time_step_soln(MAX_INT);
        }
        else {
            for i1 in 0..self.job[ci].num_modes {
                if self.eig_vals[i1] >= 0.0 {
                    self.load_fact[i1] = 0.159154943091895335 * sqrt(self.eig_vals[i1]);
                }
                else {
                    self.load_fact[i1] = -1.0;
                }
            }
        }
        
        self.job[sci].nonlinear_geom = self.job[ci].nonlinear_geom;
        self.job[sci].dynamic = self.job[ci].dynamic;
        
        return;
    }

    pub fn set_soln_to_mode(&mut self, field : &mut CppStr, mode : usize, max_val : f64) {
        let i1 : usize;
        let mut i2 : usize;
        let mode_st : usize;
        let mut n_dof : usize;
        let mut glob_ind : usize;
        let mut disp_fields = CppStr::from("displacement velocity acceleration");
        let mut nd_dat : [f64; 6] = [0f64; 6];
        let mut scale_fact : f64;
        let mut ab_val : f64;
        
        mode_st = mode * self.el_mat_dim;
        i2 = mode_st;
        scale_fact = 0.0;
        for _i1 in 0..self.el_mat_dim {
            ab_val = fabs(self.eig_vecs[i2]);
            if ab_val > scale_fact {
                scale_fact = ab_val;
            }
            i2 += 1usize;
        }
        scale_fact = max_val / scale_fact;
        
        i1 = disp_fields.find(&field.s.as_str());
        if i1 < MAX_INT {
            for this_nd in self.nodes.iter_mut() {
                n_dof = this_nd.num_dof;
                for i2 in 0..n_dof {
                    glob_ind = this_nd.dof_index[i2];
                    nd_dat[i2] = scale_fact * self.eig_vecs[mode_st + glob_ind];
                }
                if field.s == "displacement" {
                    this_nd.set_displacement(&mut nd_dat);
                }
                else if field.s == "velocity" {
                    this_nd.set_velocity(&mut nd_dat);
                }
                else {
                    this_nd.set_acceleration(&mut nd_dat);
                }
            }
        }
        else {
            println!("Warning: {} is not a valid field to Set to a mode. Aborting setSolnToMode()", field.s);
        }
        
        return;
    }

    pub fn augmentd_ld_u(&mut self) {
        let i2 : usize;
        let mut i3 : usize;
        let mut num_nds : usize;
        let mut nd_dof : usize;
        //let mut scmd = &self.job[self.solve_cmd];
        let sci : usize = self.solve_cmd;
        let dt : f64 =  self.job[sci].time_step;
        let bet : f64 =  self.job[sci].newmark_beta;
        let gam : f64 =  self.job[sci].newmark_gamma;
        let d_tdotd_tp : f64 =  -1.0 / (gam * dt);
        let d_tdotd_tdotp : f64 =  -(1.0 - gam) / gam;
        let d_ad_up : f64 =  1.0 / (dt * dt * (bet - gam));
        let d_ad_vp : f64 =  d_ad_up * dt;
        let d_ad_ap : f64 =  d_ad_up * dt * dt * (0.5 + bet - gam);
        let d_vd_up : f64 =  dt * gam * d_ad_up;
        let d_vd_vp : f64 =  1.0 + dt * gam * d_ad_vp;
        let d_vd_ap : f64 =  dt * (1.0 - gam) + dt * gam * d_ad_ap;
        
        let c1 : f64 =  d_ad_up;
        let c2 : f64 =  d_ad_ap;
        let c3 : f64 =  dt * (1.0 - gam);
        let c4 : f64 =  1.0 / (dt*(bet - gam));
        
        let mut rvec = vec![DiffDoub0::new(); 30];
        let mut mmat = vec![0f64; 900];
        let mut dmat = vec![0f64; 900];
        let mut el_adj = vec![0f64; 30];
        let mut eld_ld_t = vec![0f64; 10];
        let mut eld_ld_tdot = vec![0f64; 10];
        let mut eld_ld_u = vec![0f64; 30];
        let mut eld_ld_v = vec![0f64; 30];
        let mut eld_ld_a = vec![0f64; 30];
        
        
        
        if self.job[sci].thermal {
            for this_el in self.elements.iter_mut() {
                this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                this_el.get_rtm_dfd0(&mut rvec, &mut  mmat,  true,  true, &mut  self.d0_pre);
                this_el.get_el_vec(&mut el_adj, &mut  self.t_adj,  true,  false, &mut  self.nodes);
                num_nds = this_el.num_nds;
                i3 = 0;
                for i1 in 0..num_nds {
                    eld_ld_t[i1] = 0.0;
                    eld_ld_tdot[i1] = 0.0;
                    for i2 in 0..num_nds {
                        eld_ld_t[i1]  -=  d_tdotd_tp * mmat[i3] * el_adj[i2];
                        eld_ld_tdot[i1]  -=  d_tdotd_tdotp * mmat[i3] * el_adj[i2];
                        i3 += 1usize;
                    }
                }
                this_el.add_to_glob_vec(&mut eld_ld_t, &mut  self.d_ld_t,  true,  false, &mut  self.nodes);
                this_el.add_to_glob_vec(&mut eld_ld_tdot, &mut  self.d_ld_tdot,  true,  false, &mut  self.nodes);
            }
            i2 = self.nodes.len();
            for i1 in 0..i2 {
                self.d_ld_t[i1]  +=  self.tdot_adj[i1];
                self.d_ld_tdot[i1]  +=  c3 * self.tdot_adj[i1];
            }
        }
        
        if self.job[sci].elastic {
            for this_el in self.elements.iter_mut() {
                this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                this_el.get_rum_dfd0(&mut rvec, &mut  mmat,  true,  true,  self.job[sci].nonlinear_geom, &mut  self.d0_pre, &mut self.d0_scratch.iter_mut());
                this_el.get_rud_dfd0(&mut rvec, &mut  dmat,  true, &self.job[sci], &mut  self.d0_pre, &mut self.scratch.iter_mut(), &mut self.d0_scratch.iter_mut());
                this_el.get_el_vec(&mut el_adj, &mut  self.u_adj,  false,  false, &mut  self.nodes);
                nd_dof = this_el.num_nds * this_el.dof_per_nd;
                i3 = 0;
                for i1 in 0..nd_dof {
                    eld_ld_u[i1] = 0.0;
                    eld_ld_v[i1] = 0.0;
                    eld_ld_a[i1] = 0.0;
                    for i2 in 0..nd_dof {
                        eld_ld_u[i1]  -=  (d_ad_up * mmat[i3] + d_vd_up * dmat[i3]) * el_adj[i2];
                        eld_ld_a[i1]  -=  (d_ad_ap * mmat[i3] + d_vd_ap * dmat[i3]) * el_adj[i2];
                        eld_ld_v[i1]  -=  (d_ad_vp * mmat[i3] + d_vd_vp * dmat[i3]) * el_adj[i2];
                        i3 += 1usize;
                    }
                }
                this_el.add_to_glob_vec(&mut eld_ld_u, &mut  self.d_ld_u,  false,  false, &mut  self.nodes);
                this_el.add_to_glob_vec(&mut eld_ld_a, &mut  self.d_ld_a,  false,  false, &mut  self.nodes);
                this_el.add_to_glob_vec(&mut eld_ld_v, &mut  self.d_ld_v,  false,  false, &mut  self.nodes);
            }
            for i1 in 0..self.el_mat_dim {
                self.d_ld_u[i1]  +=  c1 * self.a_adj[i1];
                self.d_ld_a[i1]  +=  c2 * self.a_adj[i1] + c3 * self.v_adj[i1];
                self.d_ld_v[i1]  +=  c4 * self.a_adj[i1] + self.v_adj[i1];
            }
        }
        
        return;
    }

    pub fn solve_for_adjoint(&mut self, time : f64, full_ref : bool) {
        let mut i3 : usize;
        let tot_nodes : usize;
        let mut el_num_nds : usize;
        let mut el_dof_per_nd : usize;
        let mut el_nd_dof : usize;
        let mut el_int_dof : usize;
        let c1 : f64;
        let mut c2 : f64;
        let c3 : f64;
        
        let mut rvec = vec![DiffDoub0::new(); 33];
        let mut d_rd_u = vec![0f64; 1089];
        let mut d_rd_t = vec![0f64; 330];
        let mut eld_ld_t = vec![0f64; 10];
        let mut el_adj = vec![0f64; 33];
        
        //let mut scmd = &self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        
        self.obj.calculated_ld_u(&mut self.d_ld_u, &mut  self.d_ld_v, &mut  self.d_ld_a, &mut  self.d_ld_t, &mut  self.d_ld_tdot,  time, 
            self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  self.node_sets, &mut  self.element_sets, &mut  self.sections, &mut  self.materials, & self.design_vars,
            &mut self.d0_pre, &mut self.d0_scratch.iter_mut());
        
        if self.job[sci].elastic {
            self.elastic_const.update_active_status(time);
            if self.job[sci].dynamic {
                c1 = self.job[sci].time_step * self.job[sci].newmark_gamma;
                c2 = self.job[sci].time_step;
                c2 = 1.0 / (c2 * c2 * (self.job[sci].newmark_beta - self.job[sci].newmark_gamma));
                for i1 in 0..self.el_mat_dim {
                    self.v_adj[i1] = self.d_ld_v[i1];
                    self.a_adj[i1] = self.d_ld_a[i1] + c1 * self.v_adj[i1];
                    self.d_ld_u[i1]  -=  c2 * self.a_adj[i1];
                }
            }
            if self.elastic_const.any_just_activated() {
                self.build_elastic_soln_load(true, full_ref);
                self.elastic_lt.allocate_from_sparse_mat(&mut self.elastic_mat, &mut self.elastic_const, self.job[sci].solver_block_dim);
                self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
                self.elastic_lt.ldl_factor();
            }
            else if self.job[sci].nonlinear_geom {
                self.build_elastic_soln_load(true,  full_ref);
                self.elastic_lt.populate_from_sparse_mat(&mut self.elastic_mat,  &mut  self.elastic_const);
                self.elastic_lt.ldl_factor();
            }
            for this_el in self.elements.iter_mut() {
                this_el.set_intd_ld_u(&mut self.d_ld_u);
                this_el.update_external(&mut self.d_ld_u,  0, &mut  self.nodes, &mut self.scratch.iter_mut());
            }
            if self.job[sci].solver_method.s == "direct" {
                self.elastic_lt.ldl_solve(&mut self.u_adj, &mut  self.d_ld_u);
            }
            else {
                conj_grad_sparse(&mut self.u_adj, &mut  self.elastic_mat, &mut  self.elastic_const, &mut  self.elastic_lt, &mut  self.d_ld_u,  self.job[sci].conv_tol,  self.job[sci].max_it);
                //g_mres_sparse(self.u_adj, *self.elastic_mat, *self.elastic_const, *self.elastic_lt, self.d_ld_u, self.solve_cmd->conv_tol, self.solve_cmd->max_it, 6*self.solve_cmd->solver_block_dim);
            }
            for this_el in self.elements.iter_mut() {
                this_el.update_internal(&mut self.u_adj,  0, &mut  self.nodes, &mut self.scratch.iter_mut());
            }
        }
        
        if self.job[sci].thermal {
            if self.job[sci].elastic {
                for this_el in self.elements.iter_mut() {
                    this_el.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                    this_el.get_ruk_dfd0(&mut rvec, &mut  d_rd_u, &mut  d_rd_t,  true,  self.job[sci].nonlinear_geom, &mut  self.d0_pre);
                    this_el.get_el_vec(&mut el_adj, &mut  self.u_adj,  false,  false, &mut  self.nodes);
                    el_num_nds = this_el.num_nds;
                    el_dof_per_nd = this_el.dof_per_nd;
                    el_int_dof = this_el.num_int_dof;
                    el_nd_dof = el_num_nds * el_dof_per_nd;
                    for i1 in 0..el_num_nds {
                        eld_ld_t[i1] = 0.0;
                        i3 = i1;
                        for i2 in 0..el_nd_dof {
                            //i3 = i2 * el_num_nds + i1;
                            eld_ld_t[i1]  -=  d_rd_t[i3] * el_adj[i2];
                            i3  +=  el_num_nds;
                        }
                        for i2 in 0..el_int_dof {
                            eld_ld_t[i1]  -=  d_rd_t[i3] * el_adj[i2];
                            i3  +=  el_num_nds;
                        }
                    }
                    this_el.add_to_glob_vec(&mut eld_ld_t, &mut  self.d_ld_t,  true,  false, &mut  self.nodes);
                }
            }
            tot_nodes = self.nodes.len();
            self.thermal_const.update_active_status(time);
            if self.job[sci].dynamic {
                c3 = -1.0 / (self.job[sci].time_step * self.job[sci].newmark_gamma);
                for i1 in 0..tot_nodes {
                    self.tdot_adj[i1] = c3 * self.d_ld_tdot[i1];
                    self.d_ld_t[i1]  -=  self.tdot_adj[i1];
                }
            }
            if self.thermal_const.any_just_activated() {
                self.therm_lt.allocate_from_sparse_mat(&mut self.therm_mat, &mut self.thermal_const, self.job[sci].solver_block_dim);
                self.therm_lt.populate_from_sparse_mat(&mut self.therm_mat, &mut self.thermal_const);
                self.therm_lt.ldl_factor();
            }
            if self.job[sci].solver_method.s == "direct" {
                self.therm_lt.ldl_solve(&mut self.t_adj, &mut  self.d_ld_t);
            }
            else {
                conj_grad_sparse(&mut self.t_adj, &mut  self.therm_mat, &mut  self.thermal_const, &mut  self.therm_lt, &mut  self.d_ld_t,  self.job[sci].conv_tol,  self.job[sci].max_it);
                //g_mres_sparse(self.t_adj, *self.therm_mat, *self.thermal_const, *self.therm_lt, self.d_ld_t, self.solve_cmd->conv_tol, self.solve_cmd->max_it, self.solve_cmd->solver_block_dim);
            }
        }
        
        return;
    }

    pub fn d_rthermald_d(&mut self, d_var_num : usize) {
        let tot_nodes : usize;
        let mut glob_ind : usize;
        let mut dv_val = DiffDoub0::new();
        
        let mut scmd : &JobCommand = &self.job[self.solve_cmd];
        //let mut this_dv : &mut DesignVariable = &mut self.design_vars[d_var_num];
        self.design_vars[d_var_num].get_value_dfd0(&mut dv_val);
        self.design_vars[d_var_num].diff_val.set_val_2(dv_val.val, 1.0);
        
        tot_nodes = self.nodes.len();
        for i1 in 0..tot_nodes {
            self.d_rtd_d[i1].set_val(0.0);
        }
        
        for i1 in 0..self.elements.len() {
            self.el_in_d[i1] = 0;
        }
        
        for eli in self.design_vars[d_var_num].comp_el_list.iter() {
            self.el_in_d[*eli] = 1;
        }
        
        // applied contribution
        let mut this_el : &mut Element;
        for this_ld in self.thermal_loads.iter_mut() {
            if this_ld.this_type.s != "nodalHeatGen" {
                for eli in self.element_sets[this_ld.el_set_ptr].labels.iter_mut() {
                    if self.el_in_d[*eli] == 1 {
                        this_el = &mut self.elements[*eli];
                        this_el.get_stress_prereq_dfd1(&mut self.d1_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                        this_el.get_app_therm_load_dfd1(&mut self.d_rtd_d, this_ld, &mut  self.d1_pre, &mut self.scratch.iter_mut(), &mut self.d1_scratch.iter_mut(), &mut self.sections, &mut  self.faces, &mut  self.nodes, & self.design_vars);
                    }
                }
            }
        }
        
        for i1 in 0..tot_nodes {
            self.d_rtd_d[i1].neg();
        }
        
        // solution-dependent contribution of load
        for eli in self.design_vars[d_var_num].comp_el_list.iter() {
            this_el = &mut self.elements[*eli];
            this_el.get_stress_prereq_dfd1(&mut self.d1_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
            this_el.get_rt_dfd1(&mut self.d_rtd_d, &mut self.therm_mat, false, &mut scmd, &mut self.d1_pre, &mut self.scratch.iter_mut(), &mut self.d1_scratch.iter_mut(), &mut self.nodes);
        }
        
        
        // design variable dependent contribution
        if self.design_vars[d_var_num].nd_set_ptr < MAX_INT {
            let mut nd_ld = DiffDoub1::new();
            let mut this_nd : &mut Node;
            for ndi in self.node_sets[self.design_vars[d_var_num].nd_set_ptr].labels.iter_mut() {
                this_nd = &mut self.nodes[*ndi];
                this_nd.get_thermal_dvload_dfd1(&mut nd_ld, & self.design_vars);
                glob_ind = this_nd.sorted_rank;
                self.d_rtd_d[glob_ind].sub(& nd_ld);
            }
        }
        
        self.design_vars[d_var_num].diff_val.set_val_2(dv_val.val, 0.0);
        
        return;
    }

    pub fn d_relasticd_d(&mut self, d_var_num : usize) {
        let mut num_dof : usize;
        let mut glob_ind : usize;
        let mut dv_val = DiffDoub0::new();
        
        let scmd : &JobCommand = &self.job[self.solve_cmd];
        //let mut this_dv : &mut DesignVariable = &mut self.design_vars[d_var_num];
        self.design_vars[d_var_num].get_value_dfd0(&mut dv_val);
        self.design_vars[d_var_num].diff_val.set_val_2(dv_val.val, 1.0);
        
        for i1 in 0..self.el_mat_dim {
            self.d_rud_d[i1].set_val(0.0);
        }
        
        for i1 in 0..self.elements.len() {
            self.el_in_d[i1] = 0;
        }
        
        for eli in self.design_vars[d_var_num].comp_el_list.iter_mut() {
            self.el_in_d[*eli] = 1;
        }
        
        // applied contribution
        
        for this_ld in self.elastic_loads.iter_mut() {
            if this_ld.this_type.s != "nodalForce" {
                let mut this_el : &Element;
                for eli in self.element_sets[this_ld.el_set_ptr].labels.iter_mut() {
                    if self.el_in_d[*eli] > 0 {
                        this_el = &self.elements[*eli];
                        this_el.get_stress_prereq_dfd1(&mut self.d1_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                        this_el.get_app_load_dfd1(&mut self.d_rud_d, this_ld,  scmd.nonlinear_geom, &mut  self.d1_pre, &mut self.scratch.iter_mut(), &mut self.d1_scratch.iter_mut(), &mut  self.sections, &mut  self.faces, &mut  self.nodes, & self.design_vars);
                    }
                }
            }
        }
        
        for i1 in 0..self.el_mat_dim {
            self.d_rud_d[i1].neg();
        }
        
        // solution-dependent contribution of load
        //let mut this_el : &Element;
        for eli in self.design_vars[d_var_num].comp_el_list.iter() {
            //this_el = &self.elements[*eli];
            self.elements[*eli].get_stress_prereq_dfd1(&mut self.d1_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
            self.elements[*eli].get_ru_dfd1(&mut self.d_rud_d, &mut  self.elastic_mat,  false, scmd, &mut  self.d1_pre, &mut self.scratch.iter_mut(), &mut self.d1_scratch.iter_mut(), &mut  self.nodes);
        }
        
        
        // design variable dependent contribution
        
        if self.design_vars[d_var_num].nd_set_ptr < MAX_INT {
            let mut nd_ld = [DiffDoub1::new(); 6];
            let mut this_nd : &mut  Node;
            for ndi in self.node_sets[self.design_vars[d_var_num].nd_set_ptr].labels.iter_mut() {
                this_nd = &mut self.nodes[*ndi];
                this_nd.get_elastic_dvload_dfd1(&mut nd_ld, & self.design_vars);
                num_dof = this_nd.num_dof;
                for i1 in 0..num_dof {
                    glob_ind = this_nd.dof_index[i1];
                    nd_ld[i1].neg();
                    self.d_rud_d[glob_ind].add(& nd_ld[i1]);
                }
            }
        }
        
        self.design_vars[d_var_num].diff_val.set_val_2(dv_val.val,    0.0);
        
        return;
    }

    pub fn get_objective(&mut self) {
        let mut i1 : usize;
        let mut time : f64;
        
        //let mut scmd = &self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        if !self.job[sci].save_soln_hist {
            let mut er_str = "Error: The solution history must be saved to calculate the Objective gradient.\n".to_string();
            er_str = er_str + "Save history by setting saveSolnHist=True in the solve command.\n";
            panic!("{}",er_str);
        }
        
        if self.job[sci].dynamic {
            self.obj.clear_values();
            self.read_time_step_soln(self.time_steps_saved);
            i1 = self.time_steps_saved;
            time = self.job[sci].time_step * (self.time_steps_saved as f64);
            while i1 > 0 {
                i1 -= 1usize;
                for this_nd in self.nodes.iter_mut() {
                    if self.job[sci].elastic {
                        this_nd.backstep_disp();
                    }
                    if self.job[sci].thermal {
                        this_nd.backstep_temp();
                    }
                }
                if self.job[sci].elastic {
                    for this_el in self.elements.iter_mut() {
                        this_el.backstep_int_disp();
                    }
                }
                self.read_time_step_soln(i1);
                self.obj.calculate_terms(time,  self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  self.node_sets, &mut  self.element_sets, &mut  self.sections, &mut  self.materials, & self.design_vars, &mut  self.d0_pre);
                time  -=  self.job[sci].time_step;
            }
        }
        else {
            self.obj.clear_values();
            i1 = 0;
            let lt_len = self.job[sci].static_load_time.len();
            let mut lt_vec = vec![0.0f64; lt_len];
            for lt in self.job[sci].static_load_time.iter() {
                lt_vec[i1] = *lt;
                i1 += 1usize;
            }
            let mut this_ld_tm : f64;
            //for this_ld_tm in self.job[sci].static_load_time.iter_mut() {
            for i in 0..lt_len {
                this_ld_tm = lt_vec[i];
                self.read_time_step_soln(i);
                for this_nd in self.nodes.iter_mut() {
                    if self.job[sci].elastic {
                        this_nd.backstep_disp();
                    }
                    if self.job[sci].thermal {
                        this_nd.backstep_temp();
                    }
                }
                if self.job[sci].elastic {
                    for this_el in self.elements.iter_mut() {
                        this_el.backstep_int_disp();
                    }
                }
                self.obj.calculate_terms(this_ld_tm,  self.job[sci].nonlinear_geom, &mut self.nodes, &mut self.elements, &mut self.node_sets, &mut self.element_sets, &mut self.sections, &mut self.materials, &mut self.design_vars, &mut  self.d0_pre);
                //i1 += 1usize;
            }
        }
        return;
    }

    pub fn get_obj_gradient(&mut self) {
        let mut i1 : usize;
        let mut i4 : usize;
        let time : f64;
        let mut since_ref : usize;
        let num_dv : usize =  self.design_vars.len();
        
        //let mut scmd = &self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        if !self.job[sci].save_soln_hist {
            let mut er_str = "Error: The solution history must be saved to calculate the Objective gradient.\n".to_string();
            er_str  = er_str + "Save history by setting saveSolnHist=True in the solve command.\n";
            panic!("{}", er_str);
        }
        
        self.obj.clear_values();
        for i1 in 0..num_dv {
            self.d_ld_d[i1] = 0.0;
        }
        
        if self.job[sci].dynamic {
            time = self.job[sci].time_step * (self.time_steps_saved as f64);
            i1 = self.time_steps_saved;
            since_ref = 0;
            self.read_time_step_soln(i1);
            while i1 > 0 {
                for i2 in 0..self.nodes.len() {
                    self.d_ld_t[i2] = 0.0;
                    self.d_ld_tdot[i2] = 0.0;
                }
                for i2 in 0..self.el_mat_dim {
                    self.d_ld_a[i2] = 0.0;
                    self.d_ld_v[i2] = 0.0;
                }
                for i2 in 0..self.tot_glob_dof {
                    self.d_ld_u[i2] = 0.0;
                }
                if i1 < self.time_steps_saved {
                    self.augmentd_ld_u();
                }
                for this_nd in self.nodes.iter_mut() {
                    if self.job[sci].elastic {
                        this_nd.backstep_disp();
                    }
                    if self.job[sci].thermal {
                        this_nd.backstep_temp();
                    }
                }
                for this_el in self.elements.iter_mut() {
                    this_el.backstep_int_disp();
                }
                self.read_time_step_soln(i1 - 1);
                self.obj.calculate_terms(time,  self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  self.node_sets, &mut  self.element_sets, &mut  self.sections, &mut  self.materials, & self.design_vars, &mut  self.d0_pre);
                since_ref += 1usize;
                if since_ref == self.job[sci].full_reform {
                    self.solve_for_adjoint(time,  true);
                    since_ref = 0;
                }
                else {
                    self.solve_for_adjoint(time,  false);
                }
                self.obj.calculated_ld_d(&mut self.d_ld_d,  time,  self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  &mut  self.element_sets, &mut  self.sections, &mut  self.materials, &mut self.design_vars, &mut  self.d1_pre);
                for i2 in 0..num_dv {
                    if self.job[sci].thermal {
                        self.d_rthermald_d(i2);
                        for i3 in 0..self.nodes.len() {
                            self.d_ld_d[i2] -= self.t_adj[i3] * self.d_rtd_d[i3].dval;
                        }
                    }
                    if self.job[sci].elastic {
                        self.d_relasticd_d(i2);
                        for i3 in 0..self.el_mat_dim {
                            self.d_ld_d[i2]  -=  self.u_adj[i3] * self.d_rud_d[i3].dval;
                        }
                        for this_el in self.elements.iter_mut() {
                            self.d_ld_d[i2]  -=  this_el.get_int_adjd_rd_d();
                        }
                    }
                }
                i1 -= 1usize;
            }
        }
        else {
            i4 = 0;
            let lt_len = self.job[sci].static_load_time.len();
            let mut lt_vec = vec![0.0f64; lt_len];
            for lt in self.job[sci].static_load_time.iter() {
                lt_vec[i4] = *lt;
                i4 += 1usize;
            }
            //for this_ld in self.job[sci].static_load_time.iter_mut() {
            let mut this_ld : f64;
            for i in 0..lt_len {
                this_ld = lt_vec[i];
                self.read_time_step_soln(i);
                for this_nd in self.nodes.iter_mut() {
                    if self.job[sci].elastic {
                        this_nd.backstep_disp();
                    }
                    if self.job[sci].thermal {
                        this_nd.backstep_temp();
                    }
                }
                for this_el in self.elements.iter_mut() {
                    this_el.backstep_int_disp();
                }
                self.obj.calculate_terms(this_ld,  self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  self.node_sets, &mut  self.element_sets, &mut  self.sections, &mut  self.materials, & self.design_vars, &mut  self.d0_pre);
                for i2 in 0..self.nodes.len() {
                    self.d_ld_t[i2] = 0.0;
                    self.d_ld_tdot[i2] = 0.0;
                }
                for i2 in 0..self.el_mat_dim {
                    self.d_ld_a[i2] = 0.0;
                    self.d_ld_v[i2] = 0.0;
                }
                for i2 in 0..self.tot_glob_dof {
                    self.d_ld_u[i2] = 0.0;
                }
                self.solve_for_adjoint(this_ld, true);
                self.obj.calculated_ld_d(&mut self.d_ld_d, this_ld, self.job[sci].nonlinear_geom, &mut  self.nodes, &mut  self.elements, &mut  self.element_sets, &mut  self.sections, &mut  self.materials, &mut self.design_vars, &mut  self.d1_pre);
                for i1 in 0..num_dv {
                    if self.job[sci].thermal {
                        self.d_rthermald_d(i1);
                        for i3 in 0..self.nodes.len() {
                            self.d_ld_d[i1]  -=  self.t_adj[i3] * self.d_rtd_d[i3].dval;
                        }
                    }
                    if self.job[sci].elastic {
                        self.d_relasticd_d(i1);
                        for i2 in 0..self.el_mat_dim {
                            self.d_ld_d[i1]  -=  self.u_adj[i2] * self.d_rud_d[i2].dval;
                        }
                        for this_el in self.elements.iter_mut() {
                            self.d_ld_d[i1]  -=  this_el.get_int_adjd_rd_d();
                        }
                    }
                }
                //i4 += 1usize;
            }
        }
        return;
    }

}


