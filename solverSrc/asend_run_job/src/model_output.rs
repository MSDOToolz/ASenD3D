use crate::model::*;
use crate::constants::*;
use crate::diff_doub::*;
use crate::node::*;
use crate::element::*;
use crate::nd_el_set::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;

use std::collections::LinkedList;
use std::fs::File;
use std::io::Write;

impl Model {
    pub fn write_time_step_soln(&mut self, t_step : usize) {
        let full_file = format!("{}{}{}{}",self.job[self.solve_cmd].file_name.s, "/solnTStep", t_step, ".out");

        let mut out_file = match File::create(&full_file) {
            Err(why) => panic!("couldn't open file, {}, {}", full_file, why),
            Ok(file) => file,
        };

        for nd in self.nodes.iter() {
            if self.job[self.solve_cmd].thermal {
                let _ = out_file.write(&nd.prev_temp.to_be_bytes());
                let _ = out_file.write(&nd.prev_tdot.to_be_bytes());
            }
            if self.job[self.solve_cmd].elastic {
                for i in 0..nd.num_dof {
                    let _ = out_file.write(&nd.prev_disp[i].to_be_bytes());
                }

                for i in 0..nd.num_dof {
                    let _ = out_file.write(&nd.prev_vel[i].to_be_bytes());
                }

                for i in 0..nd.num_dof {
                    let _ = out_file.write(&nd.prev_acc[i].to_be_bytes());
                }

            }
        }

        if self.job[self.solve_cmd].elastic {
            for el in self.elements.iter() {
                if el.num_int_dof > 0 {
                    for i in 0..el.num_int_dof {
                        let _ = out_file.write(&el.int_prev_disp[i].to_be_bytes());
                    }
                }
            }
        }


    }

    pub fn write_node_results(&mut self, file_name : &mut CppStr, node_set : &mut CppStr, fields : &mut LinkedList<CppStr>, time_step : usize) {
        let mut i1 : usize;
        let mut glob_ind : usize;
        let set_pt : usize;
        let mut nd_dat : [f64; 6] = [0f64; 6];
        let mut ndof : usize;
        let time : f64;
        
        //let mut scmd = &self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        
        if time_step < MAX_INT {
            // read the results from the time step file and store them in self.nodes
            self.read_time_step_soln(time_step);
            for nd_pt in self.nodes.iter_mut() {
                if self.job[sci].thermal {
                    nd_pt.backstep_temp();
                }
                if self.job[sci].elastic {
                    nd_pt.backstep_disp();
                }
            }
            if self.job[sci].dynamic {
                time = self.job[sci].time_step * (time_step as f64);
            }
            else {
                let mut lt_iter = self.job[sci].static_load_time.iter();
                let mut lst_end : bool;
                let mut this_ld_tm : f64 = 0.0f64;
                match lt_iter.next() {
                    None => {lst_end = true;},
                    Some(x) => {this_ld_tm = *x; lst_end = false;}, 
                }
                i1 = 0;
                while !lst_end && i1 < time_step {
                    match lt_iter.next() {
                        None => {lst_end = true;},
                        Some(x) => {this_ld_tm = *x; lst_end = false;}, 
                    }
                    i1 += 1usize;
                }
                time = this_ld_tm;
            }
        }
        else {
            time = -1.0;
        }
        
        if CppMap::key_in_map(&mut self.ns_map, &mut node_set.s) {
            set_pt = self.ns_map.at(&node_set.to_string());
        }
        else {
            println!("Warning: there is no Node Set named {}. Defaulting to all nodes in writeNodeResults",node_set.s);
            set_pt = self.ns_map.at(&"all".to_string());
        }
        //let mut this_set = &self.node_sets[set_pt];
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(format!("{}", "nodeResults:\n").as_bytes());
        let _ = out_file.write(format!("{}{}{}", "    time: " , time , "\n").as_bytes());
        let _ = out_file.write(format!("{}{}{}", "    nodeSet: " , node_set.s , "\n").as_bytes());
        
        for this_field in fields.iter_mut() {
            let _ = out_file.write(format!("{}{}{}", "    " , this_field.s , ":\n").as_bytes());
            // -----------------
            // calculate reaction force if necessary
            if this_field.s == "reactionForce" {
                for i1 in 0..self.el_mat_dim {
                    self.elastic_ld_vec[i1] = 0.0;
                }
                self.build_elastic_soln_load(false,  false);
                let mut nd_pt : &Node;
                for nd_label in self.node_sets[set_pt].labels.iter_mut() {
                    nd_pt = &self.nodes[*nd_label];
                    let _ = out_file.write(format!("{}{}", "        - [" , *nd_label).as_bytes());
                    ndof = nd_pt.num_dof;
                    for i1 in 0..ndof {
                        glob_ind = nd_pt.dof_index[i1];
                        let _ = out_file.write(format!("{0}{1:.12e}", ", " , -self.elastic_ld_vec[glob_ind]).as_bytes());
                    }
                    let _ = out_file.write(format!("{}", "]\n").as_bytes());
                }
            }
            else if this_field.s == "reactionHeatGen" {
                for i1 in 0..self.el_mat_dim {
                    self.therm_ld_vec[i1] = 0.0;
                }
                self.build_thermal_soln_load(false);
                let mut nd_pt : &Node;
                for nd_label in self.node_sets[set_pt].labels.iter_mut() {
                    nd_pt = &self.nodes[*nd_label];
                    let _ = out_file.write(format!("{}{}", "        - [" , *nd_label).as_bytes());
                    glob_ind = nd_pt.sorted_rank;
                    let _ = out_file.write(format!("{0}{1:.12e}{2}", ", " , -self.therm_ld_vec[glob_ind] , "\n").as_bytes());
                }
            }
            else {
                let mut nd_pt : &Node;
                for nd_label in self.node_sets[set_pt].labels.iter_mut() {
                    nd_pt = &self.nodes[*nd_label];
                    let _ = out_file.write(format!("{}{}", "        - [" , *nd_label).as_bytes());
                    if this_field.s == "displacement" {
                        ndof = nd_pt.num_dof;
                        for i1 in 0..ndof {
                            let _ = out_file.write(format!("{0}{1:.12e}", ", " , nd_pt.displacement[i1]).as_bytes());
                        }
                        let _ = out_file.write(format!("{}", "]\n").as_bytes());
                    }
                    else if this_field.s == "velocity" {
                        ndof = nd_pt.num_dof;
                        for i1 in 0..ndof {
                            let _ = out_file.write(format!("{0}{1:.12e}", ", " , nd_pt.velocity[i1]).as_bytes());
                        }
                        let _ = out_file.write(format!("{}", "]\n").as_bytes());
                    }
                    else if this_field.s == "acceleration" {
                        ndof = nd_pt.num_dof;
                        for i1 in 0..ndof {
                            let _ = out_file.write(format!("{0}{1:.12e}", ", " , nd_pt.acceleration[i1]).as_bytes());
                        }
                        let _ = out_file.write(format!("{}", "]\n").as_bytes());
                    }
                    else if this_field.s == "temperature" {
                        nd_dat[0] = nd_pt.temperature;
                        let _ = out_file.write(format!("{0}{1:.12e}{2}", ", " , nd_dat[0] , "]\n").as_bytes());
                    }
                    else if this_field.s == "tdot" {
                        nd_dat[0] = nd_pt.temp_change_rate;
                        let _ = out_file.write(format!("{0}{1:.12e}{2}", ", " , nd_dat[0] , "]\n").as_bytes());
                    }
                }
            }
        }
        
        return;
    }

    pub fn write_element_results(&mut self, file_name : &mut CppStr, el_set : &mut CppStr, fields : &mut LinkedList<CppStr>, position : &mut CppStr, time_step : usize) {
        let mut i1 : usize;
        let mut i2 : usize;
        let set_pt : usize;
        let mut this_type : usize;
        let mut num_ip : usize;
        let mut num_lay : usize;
        let mut int_pts : [f64; 24] = [0f64; 24];
        let mut strain = [DiffDoub0::new(); 6];
        let mut stress = [DiffDoub0::new(); 6];
        let mut seden : f64;
        let mut def = [DiffDoub0::new(); 9];
        let mut frc_mom = [DiffDoub0::new(); 9];
        let mut t_grad = [DiffDoub0::new(); 3];
        let mut flux = [DiffDoub0::new(); 3];
        let mut field_list : CppStr;
        let time : f64;
        
        //let mut scmd = &self.job[self.solve_cmd];
        let sci = self.solve_cmd;
        
        if time_step < MAX_INT {
            // read the results from the time step file and store them in self.nodes
            self.read_time_step_soln(time_step);
            for nd_pt in self.nodes.iter_mut() {
                if self.job[sci].thermal {
                    nd_pt.backstep_temp();
                }
                if self.job[sci].elastic {
                    nd_pt.backstep_disp();
                }
            }
            if self.job[sci].elastic {
                for el_pt in self.elements.iter_mut() {
                    el_pt.backstep_int_disp();
                }
            }
            if self.job[sci].dynamic {
                time = self.job[sci].time_step * (time_step as f64);
            }
            else {
                let mut lt_iter = self.job[sci].static_load_time.iter();
                let mut lst_end : bool;
                let mut this_ld_tm : f64 = 0.0f64;
                match lt_iter.next() {
                    None => {lst_end = true;},
                    Some(x) => {this_ld_tm = *x; lst_end = false;},
                }
                i1 = 0;
                while !lst_end && i1 < time_step {
                    i1 += 1usize;
                    match lt_iter.next() {
                        None => {lst_end = true;},
                        Some(x) => {this_ld_tm = *x; lst_end = false;},
                    }
                }
                time = this_ld_tm;
            }
        }
        else {
            time = -1.0;
        }
        
        if CppMap::key_in_map(&mut self.es_map, &el_set.s) {
            set_pt = self.es_map.at(&el_set.to_string());
        }
        else {
            println!("Warning: there is no Element Set named {}. Defaulting to all elements in writeNodeResults", el_set.s);
            set_pt = self.es_map.at(&"all".to_string());
        }
        let this_set : &Set = &self.element_sets[set_pt];
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(format!("{}", "elementResults:\n").as_bytes());
        let _ = out_file.write(format!("{}{}{}", "    time: " , time , "\n").as_bytes());
        let _ = out_file.write(format!("{}{}{}", "    elSet: " , el_set.s , "\n").as_bytes());
        
        for this_field in fields.iter_mut() {
            let _ = out_file.write(format!("{}{}{}", "    " , this_field.s , ":\n").as_bytes());
            field_list = CppStr::from("stress strain");
            i2 = field_list.find(this_field.s.as_str());
            if i2 < MAX_INT {
                let _ = out_file.write(format!("{}", "    ##  - [element label, integration pt, layer, S11, S22, S33, S12, S13, S23]\n").as_bytes());
            }
            field_list = CppStr::from("strainEnergyDen");
            i2 = field_list.find(this_field.s.as_str());
            if i2 < MAX_INT {
                let _ = out_file.write(format!("{}", "    ##  - [element label, integration pt, layer, strain energy]\n").as_bytes());
            }
            field_list = CppStr::from("sectionDef sectionFrcMom");
            i2 = field_list.find(this_field.s.as_str());
            if i2 < MAX_INT {
                let _ = out_file.write(format!("{}", "    ## for shells:  - [element label, integration pt, S11, S22, S12, K11, K22, K12]\n").as_bytes());
                let _ = out_file.write(format!("{}", "    ## for beams:  - [element label, integration pt, S11, S12, S13, K11, K12, K13]\n").as_bytes());
            }
            if this_field.s == "heatFlux" {
                let _ = out_file.write(format!("{}", "    ##  - [element label, integration pt, layer, f1, f2, f3]\n").as_bytes());
            }
            if this_field.s == "tempGradient" {
                let _ = out_file.write(format!("{}", "    ##  - [element label, integration pt, layer, dT/dx, dT/dy, dT/dz]\n").as_bytes());
            }
            
            let mut el_pt : &mut Element;
            for el_label in this_set.labels.iter() {
                el_pt = &mut self.elements[*el_label];
                this_type = el_pt.this_type;
                
                field_list = CppStr::from("stress strain strainEnergyDen");
                i2 = field_list.find(&this_field.s.as_str());
                if i2 < MAX_INT {
                    el_pt.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                    if position.s == "intPts" {
                        num_ip = el_pt.num_ip;
                        vec_to_ar(&mut int_pts, & el_pt.int_pts,  0,  3 * num_ip);
                    }
                    else {
                        num_ip = 1;
                        int_pts[0] = el_pt.s_cent[0];
                        int_pts[1] = el_pt.s_cent[1];
                        int_pts[2] = el_pt.s_cent[2];
                    }
                    for i1 in 0..num_ip {
                        if this_type == 3 || this_type == 41 {
                            num_lay = self.sections[el_pt.sect_ptr].layers.len();
                        } else {
                            num_lay = 1;
                        }
                        for i2 in 0..num_lay {
                            el_pt.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut int_pts[(3*i1)..],  i2,  self.job[sci].nonlinear_geom, &mut  self.d0_pre);
                            let _ = out_file.write(format!("{}{}{}", "        - [" , el_label , ", ").as_bytes());
                            if this_field.s == "strain" {
                                let _ = out_file.write(format!("{0}{1}{2}{3}{4:.12e}", i1 , ", " , i2 , ", " , strain[0].val).as_bytes());
                                for i3 in 1..6 {
                                    let _ = out_file.write(format!("{0}{1:.12e}", ", " , strain[i3].val).as_bytes());
                                }
                                let _ = out_file.write(format!("{}", "]\n").as_bytes());
                            } else if this_field.s == "stress" {
                                let _ = out_file.write(format!("{0}{1}{2}{3}{4:.12e}", i1 , ", " , i2 , ", " , stress[0].val).as_bytes());
                                for i3 in 1..6 {
                                    let _ = out_file.write(format!("{0}{1:.12e}", ", " , stress[i3].val).as_bytes());
                                }
                                let _ = out_file.write(format!("{}", "]\n").as_bytes());
                            } else {
                                seden = 0.0;
                                for i3 in 0..6 {
                                    seden  +=  stress[i3].val * strain[i3].val;
                                }
                                seden  *=  0.5;
                                let _ = out_file.write(format!("{0}{1}{2}{3}{4:.12e}{5}", i1 , ", " , i2 , ", " , seden , "]\n").as_bytes());
                            }
                        }
                    }
                }
                
                field_list = CppStr::from("sectionDef sectionFrcMom");
                i2 = field_list.find(&this_field.s.as_str());
                if i2 < MAX_INT && el_pt.dof_per_nd == 6 {
                    el_pt.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                    if position.s == "intPts" {
                        num_ip = el_pt.num_ip;
                        vec_to_ar(&mut int_pts, & el_pt.int_pts,  0,  3 * num_ip);
                    }
                    else {
                        num_ip = 1;
                        int_pts[0] = el_pt.s_cent[0];
                        int_pts[1] = el_pt.s_cent[1];
                        int_pts[2] = el_pt.s_cent[2];
                    }
                    for i1 in 0..num_ip {
                        el_pt.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut int_pts[(3*i1)..],  self.job[sci].nonlinear_geom, &mut self.d0_pre);
                        let _ = out_file.write(format!("{}{}{}", "        - [" , el_label , ", ").as_bytes());
                        if this_field.s == "sectionDef" {
                            let _ = out_file.write(format!("{0}{1}{2:.12e}", i1 , ", " , def[0].val).as_bytes());
                            for i3 in 1..6 {
                                let _ = out_file.write(format!("{0}{1:.12e}", ", " , def[i3].val).as_bytes());
                            }
                            let _ = out_file.write(format!("{0}", "]\n").as_bytes());
                        }
                        else if this_field.s == "sectionFrcMom" {
                            let _ = out_file.write(format!("{0}{1}{2:.12e}", i1 , ", " , frc_mom[0].val).as_bytes());
                            for i3 in 1..6 {
                                let _ = out_file.write(format!("{0}{1:.12e}", ", " , frc_mom[i3].val).as_bytes());
                            }
                            let _ = out_file.write(format!("{}", "]\n").as_bytes());
                        }
                    }
                }
                
                field_list = CppStr::from("tempGradient heatFlux");
                i2 = field_list.find(&this_field.s.as_str());
                if i2 < MAX_INT {
                    el_pt.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
                    if position.s == "intPts" {
                        num_ip = el_pt.num_ip;
                        vec_to_ar(&mut int_pts, & el_pt.int_pts,  0,  3 * num_ip);
                    }
                    else {
                        num_ip = 1;
                        int_pts[0] = el_pt.s_cent[0];
                        int_pts[1] = el_pt.s_cent[1];
                        int_pts[2] = el_pt.s_cent[2];
                    }
                    for i1 in 0..num_ip {
                        this_type = el_pt.this_type;
                        if this_type == 3 || this_type == 41 {
                            num_lay = self.sections[el_pt.sect_ptr].layers.len();
                        }
                        else {
                            num_lay = 1;
                        }
                        for i2 in 0..num_lay {
                            el_pt.get_flux_tgrad_dfd0(&mut flux, &mut t_grad, &mut int_pts[(3*i1)..],  i2, &mut self.d0_pre);
                            let _ = out_file.write(format!("{}{}{}{}{}{}{}", "        - [" , el_label , ", " , i1 , ", " , i2 , ", ").as_bytes());
                            if this_field.s == "tempGradient" {
                                let _ = out_file.write(format!("{0:.12e}{1}{2:.12e}{3}{4:.12e}{5}", t_grad[0].val , ", " , t_grad[1].val , ", " , t_grad[2].val , "]\n").as_bytes());
                            }
                            else if this_field.s == "heatFlux" {
                                let _ = out_file.write(format!("{0:.12e}{1}{2:.12e}{3}{4:.12e}{5}", flux[0].val , ", " , flux[1].val , ", " , flux[2].val , "]\n").as_bytes());
                            }
                        }
                    }
                }
                
            }
        }
        
        return;
    }

    pub fn write_modal_results(&mut self, file_name : &mut CppStr) {
        let mut i3 : usize;
        let mut nd : usize;
        let mut dof_per_nd : usize;
        let mut glob_ind : usize;
        let mcmd = &self.job[self.modal_cmd];
        let n_mds : usize =  mcmd.num_modes;
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(format!("{}", "modalResults:\n").as_bytes());
        let _ = out_file.write(format!("{}", "    eigenValues:\n").as_bytes());
        for i1 in 0..n_mds {
            let _ = out_file.write(format!("{0}{1:.12e}{2}", "      - " , self.eig_vals[i1] , "\n").as_bytes());
        }
        if mcmd.this_type.s == "buckling" {
            let _ = out_file.write(format!("{}", "    loadFactors:\n").as_bytes());
        }
        else {
            let _ = out_file.write(format!("{}", "    frequencies:\n").as_bytes());
        }
        for i1 in 0..n_mds {
            let _ = out_file.write(format!("{0}{1:.12e}{2}", "      - " , self.load_fact[i1] , "\n").as_bytes());
        }
        
        if mcmd.write_modes {
            let _ = out_file.write(format!("{}", "    modes:\n").as_bytes());
            for i1 in 0..n_mds {
                let _ = out_file.write(format!("{}{}{}", "      - mode: " , i1 , "\n").as_bytes());
                let _ = out_file.write(format!("{}", "        displacement:\n").as_bytes());
                for this_nd in self.nodes.iter_mut() {
                    nd = this_nd.label;
                    let _ = out_file.write(format!("{}{}", "          - [" , nd).as_bytes());
                    dof_per_nd = this_nd.num_dof;
                    for i2 in 0..dof_per_nd {
                        glob_ind = this_nd.dof_index[i2];
                        i3 = i1 * self.el_mat_dim + glob_ind;
                        let _ = out_file.write(format!("{0}{1:.12e}", ", " , self.eig_vecs[i3]).as_bytes());
                    }
                    let _ = out_file.write(format!("{}", "]\n").as_bytes());
                }
            }
        }
        
        return;
    }

    pub fn write_objective(&mut self, file_name : &mut CppStr, include_fields : &mut LinkedList<CppStr>, write_grad : bool) {
        let num_dv : usize;
        let mut tot_obj : f64;
        let mut this_val : f64;
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(format!("{}", "objective:\n").as_bytes());
        let _ = out_file.write(format!("{}", "    terms:\n").as_bytes());
        tot_obj = 0.0;
        for this_term in self.obj.terms.iter_mut() {
            this_val = this_term.value;
            let _ = out_file.write(format!("{0}{1:.12e}{2}", "        - value: " , this_val , "\n").as_bytes());
            for fld_str in include_fields.iter_mut() {
                if fld_str.s == "category" {
                    let _ = out_file.write(format!("{}{}{}", "          category: " , this_term.category.s , "\n").as_bytes());
                }
                else if fld_str.s == "operator" {
                    let _ = out_file.write(format!("{}{}{}", "          operator: " , this_term.optr.s , "\n").as_bytes());
                }
                else if fld_str.s == "component" {
                    let _ = out_file.write(format!("{}{}{}", "          component: " , this_term.component , "\n").as_bytes());
                }
                else if fld_str.s == "layer" {
                    let _ = out_file.write(format!("{}{}{}", "          Layer: " , this_term.layer , "\n").as_bytes());
                }
                else if fld_str.s == "coefficient" {
                    let _ = out_file.write(format!("{0}{1:.12e}{2}", "          coefficient: " , this_term.coef , "\n").as_bytes());
                }
                else if fld_str.s == "exponent" {
                    let _ = out_file.write(format!("{0}{1:.12e}{2}", "          exponent: " , this_term.expnt , "\n").as_bytes());
                }
                else if fld_str.s == "elementSet" {
                    let _ = out_file.write(format!("{}{}{}", "          elementSet: " , this_term.el_set_name.s , "\n").as_bytes());
                }
                else if fld_str.s == "nodeSet" {
                    let _ = out_file.write(format!("{}{}{}", "          nodeSet: " , this_term.nd_set_name.s , "\n").as_bytes());
                }
                else if fld_str.s == "activeTime" {
                    let _ = out_file.write(format!("{0}{1:.12e}{2}{3:.12e}{4}", "          activeTime: [" , this_term.active_time[0] , ", " , this_term.active_time[1] , "]\n").as_bytes());
                }
            }
            tot_obj  +=  this_val;
        }
        let _ = out_file.write(format!("{0}{1:.12e}{2}", "    totalValue: " , tot_obj , "\n").as_bytes());
        if write_grad {
            let _ = out_file.write(format!("{}", "objectiveGradient:\n").as_bytes());
            num_dv = self.design_vars.len();
            for i1 in 0..num_dv {
                let _ = out_file.write(format!("{0}{1}{2}{3:.12e}{4}", "    - [" , i1 , ", " , self.d_ld_d[i1] , "]\n").as_bytes());
            }
        }
        
        return;
    }

}


