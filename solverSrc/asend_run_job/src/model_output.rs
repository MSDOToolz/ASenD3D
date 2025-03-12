use crate::model::*;
use crate::constants::*;
use crate::diff_doub::*;
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

    pub fn write_node_results(&mut self, file_name : &CppStr, node_set : &CppStr, fields : &mut LinkedList<CppStr>, time_step : usize) {
        let mut i1 : usize;
        let mut glob_ind : usize;
        let set_pt : usize;
        let mut nd_dat : [f64; 6] = [0f64; 6];
        let time : f64;
        
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
        
        if CppMap::key_in_map(&mut self.ns_map, &node_set.s) {
            set_pt = self.ns_map.at(&node_set.s);
        }
        else {
            println!("Warning: there is no Node Set named {}. Defaulting to all nodes in writeNodeResults",node_set.s);
            set_pt = self.ns_map.at(&"all".to_string());
        }
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };

        let _ = out_file.write(b"node,");
        for this_field in fields.iter() {
            if this_field.s == "displacement" {
                let _ = out_file.write(b"U1,U2,U3,R1,R2,R3,");
            }
            else if this_field.s == "velocity" {
                let _ = out_file.write(b"V1,V2,V3,RV1,RV2,RV3,");
            }
            else if this_field.s == "acceleration" {
                let _ = out_file.write(b"A1,A2,A3,RA1,RA2,RA3,");
            }
            else if this_field.s == "temperature" {
                let _ = out_file.write(b"T,");
            }
            else if this_field.s == "tdot" {
                let _ = out_file.write(b"TDOT,");
            }
            else if this_field.s == "reactionForce" {
                let _ = out_file.write(b"RF1,RF2,RF3,RM1,RM2,RM3,");
            }
            else if this_field.s == "reactionHeatGen" {
                let _ = out_file.write(b"RHG,");
            }
        }
        let _ = out_file.write(b"time\n");

        let mut lab_vec = vec![0usize; self.node_sets[set_pt].labels.len()];
        let mut i2 = 0usize;
        for ni in self.node_sets[set_pt].labels.iter() {
            lab_vec[i2] = *ni;
            i2 += 1;
        }
        
        for nd_label in lab_vec.iter() {
            let _ = out_file.write(format!("{},", *nd_label).as_bytes());
            for this_field in fields.iter_mut() {
                // -----------------
                // calculate reaction force if necessary
                if this_field.s == "reactionForce" {
                    for i1 in 0..self.el_mat_dim {
                        self.elastic_ld_vec[i1] = 0.0;
                    }
                    self.build_elastic_soln_load(false,  false);
                    let nd_pt = &self.nodes[*nd_label];
                    for i1 in 0..6 {
                        glob_ind = nd_pt.dof_index[i1];
                        let _ = out_file.write(format!("{0:.12e},", -self.elastic_ld_vec[glob_ind]).as_bytes());
                    }
                }
                else if this_field.s == "reactionHeatGen" {
                    for i1 in 0..self.el_mat_dim {
                        self.therm_ld_vec[i1] = 0.0;
                    }
                    self.build_thermal_soln_load(false);
                    let nd_pt = &self.nodes[*nd_label];
                    glob_ind = nd_pt.sorted_rank;
                    let _ = out_file.write(format!("{0:.12e},", -self.therm_ld_vec[glob_ind]).as_bytes());
                }
                else {
                    let nd_pt = &self.nodes[*nd_label];
                    if this_field.s == "displacement" {
                        for i1 in 0..6 {
                            let _ = out_file.write(format!("{0:.12e},", nd_pt.displacement[i1]).as_bytes());
                        }
                    }
                    else if this_field.s == "velocity" {
                        for i1 in 0..6 {
                            let _ = out_file.write(format!("{0:.12e},", nd_pt.velocity[i1]).as_bytes());
                        }
                    }
                    else if this_field.s == "acceleration" {
                        for i1 in 0..6 {
                            let _ = out_file.write(format!("{0:.12e},", nd_pt.acceleration[i1]).as_bytes());
                        }
                    }
                    else if this_field.s == "temperature" {
                        nd_dat[0] = nd_pt.temperature;
                        let _ = out_file.write(format!("{0:.12e},", nd_dat[0]).as_bytes());
                    }
                    else if this_field.s == "tdot" {
                        nd_dat[0] = nd_pt.temp_change_rate;
                        let _ = out_file.write(format!("{0:.12e},", nd_dat[0]).as_bytes());
                    }
                }
            }
            let _ = out_file.write(format!("{0:.12e}\n", time).as_bytes());
        }

    }

    pub fn write_element_results(&mut self, file_name : &mut CppStr, el_set : &mut CppStr, fields : &mut LinkedList<CppStr>, position : &mut CppStr, time_step : usize) {
        let mut i1 : usize;
        let mut str_ind : usize;
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

        let _ = out_file.write(b"element,int_pt,layer,");
        for this_field in fields.iter() {
            if this_field.s == "stress" {
                let _ = out_file.write(b"S11,S22,S33,S12,S13,S23,");
            }
            else if this_field.s == "strain" {
                let _ = out_file.write(b"E11,E22,E33,E12,E13,E23,");
            }
            else if this_field.s == "strainEnergDen" {
                let _ = out_file.write(b"SE,");
            }
            else if this_field.s == "sectionFrcMom" {
                let _ = out_file.write(b"SECT_F1,SECT_F2,SECT_F3,SECT_M1,SECT_M2,SECT_M3,");
            }
            else if this_field.s == "sectionDef" {
                let _ = out_file.write(b"SECT_E1,SECT_E2,SECT_E3,SECT_K1,SECT_K2,SECT_K3,");
            }
            else if this_field.s == "heatFlux" {
                let _ = out_file.write(b"HFLX1,HFLX2,HFLX3,");
            }
            else if this_field.s == "tempGradient" {
                let _ = out_file.write(b"TGRAD1,TGRAD2,TGRAD3,");
            }
        }
        let _ = out_file.write(b"time\n");

        for el_label in this_set.labels.iter() {
            let el_pt = &mut self.elements[*el_label];
            this_type = el_pt.this_type;
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
            el_pt.get_stress_prereq_dfd0(&mut self.d0_pre, &mut  self.sections, &mut  self.materials, &mut  self.nodes, & self.design_vars);
            for i1 in 0..num_ip {
                if this_type == 3 || this_type == 41 {
                    num_lay = self.sections[el_pt.sect_ptr].layers.len();
                } else {
                    num_lay = 1;
                }
                for i2 in 0..num_lay {
                    let _ = out_file.write(format!("{},{},{},", *el_label, i1, i2).as_bytes());
                    for this_field in fields.iter_mut() {
                                                
                        field_list = CppStr::from("stress strain strainEnergyDen");
                        str_ind = field_list.find(&this_field.s.as_str());
                        if str_ind < MAX_INT {
                            el_pt.get_stress_strain_dfd0(&mut stress, &mut  strain, &mut int_pts[(3*i1)..],  i2,  self.job[sci].nonlinear_geom, &mut  self.d0_pre);
                            if this_field.s == "strain" {
                                for i3 in 0..6 {
                                    let _ = out_file.write(format!("{0:.12e},", strain[i3].val).as_bytes());
                                }
                            } else if this_field.s == "stress" {
                                for i3 in 0..6 {
                                    let _ = out_file.write(format!("{0:.12e},", stress[i3].val).as_bytes());
                                }
                            } else {
                                seden = 0.0;
                                for i3 in 0..6 {
                                    seden  +=  stress[i3].val * strain[i3].val;
                                }
                                seden  *=  0.5;
                                let _ = out_file.write(format!("{0:.12e},", seden).as_bytes());
                            }
                        }
                        
                        field_list = CppStr::from("sectionDef sectionFrcMom");
                        str_ind = field_list.find(&this_field.s.as_str());
                        if str_ind < MAX_INT && el_pt.dof_per_nd == 6 {
                            el_pt.get_def_frc_mom_dfd0(&mut def, &mut  frc_mom, &mut int_pts[(3*i1)..],  self.job[sci].nonlinear_geom, &mut self.d0_pre);
                            if this_field.s == "sectionDef" {
                                for i3 in 0..6 {
                                    let _ = out_file.write(format!("{0:.12e},", def[i3].val).as_bytes());
                                }
                            }
                            else if this_field.s == "sectionFrcMom" {
                                for i3 in 0..6 {
                                    let _ = out_file.write(format!("{0:.12e},", frc_mom[i3].val).as_bytes());
                                }
                            }
                        }
                        
                        field_list = CppStr::from("tempGradient heatFlux");
                        str_ind = field_list.find(&this_field.s.as_str());
                        if str_ind < MAX_INT {
                            el_pt.get_flux_tgrad_dfd0(&mut flux, &mut t_grad, &mut int_pts[(3*i1)..],  i2, &mut self.d0_pre);
                            if this_field.s == "tempGradient" {
                                let _ = out_file.write(format!("{0:.12e},{1:.12e},{2:.12e},", t_grad[0].val, t_grad[1].val, t_grad[2].val).as_bytes());
                            }
                            else if this_field.s == "heatFlux" {
                                let _ = out_file.write(format!("{0:.12e},{1:.12e},{2:.12e},", flux[0].val, flux[1].val, flux[2].val).as_bytes());
                            }
                        }
                    }
                    let _ = out_file.write(format!("{0:.12e}\n", time).as_bytes());
                }
            }
        }
        
        
        return;
    }

    pub fn write_modal_results(&mut self, file_name : &mut CppStr) {
        let mut i3 : usize;
        let mut glob_ind : usize;
        let mcmd = &self.job[self.modal_cmd];
        let n_mds : usize =  mcmd.num_modes;
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(b"mode,eigenvalue,");
        if mcmd.this_type.s == "buckling" {
            let _ = out_file.write(b"load_factor\n");
        }
        else {
            let _ = out_file.write(b"frequency\n");
        }

        for i1 in 0..n_mds {
            let _ = out_file.write(format!("{0},{1:.12e},{2:.12e}\n", i1, self.eig_vals[i1], self.load_fact[i1]).as_bytes());
        }

        if mcmd.write_modes {
            let dot_ind = file_name.find(".");
            let exten : CppStr;
            let rt_file : CppStr;
            if dot_ind < MAX_INT {
                rt_file = file_name.substr(0, dot_ind);
                exten = file_name.substr(dot_ind,MAX_INT);
            }
            else {
                rt_file = file_name.clone();
                exten = CppStr::from("");
            }
            for i1 in 0..n_mds {
                let full_file = format!("{}_mode{}{}", rt_file.s, i1, exten.s);
                let mut out_file = match File::create(&full_file) {
                    Err(_why) => panic!("could not open file {}", full_file),
                    Ok(file) => file,
                };
                let _ = out_file.write(b"Node,U1,U2,U3,R1,R2,R3\n");
                for this_nd in self.nodes.iter() {
                    let _ = out_file.write(format!("{}", this_nd.label).as_bytes());
                    for i2 in 0..6 {
                        glob_ind = this_nd.dof_index[i2];
                        i3 = i1 * self.el_mat_dim + glob_ind;
                        let _ = out_file.write(format!(",{0:.12e}", self.eig_vecs[i3]).as_bytes());
                    }
                    let _ = out_file.write(b"\n");
                }
            }
        }
        
    }

    pub fn write_objective(&mut self, file_name : &mut CppStr, include_fields : &mut LinkedList<CppStr>, write_grad : bool) {
        let num_dv : usize;
        let mut this_val : f64;
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(b"objectiveTerm,value");
        for fld_str in include_fields.iter() {
            let _ = out_file.write(format!(",{}",fld_str.s).as_bytes());
        }
        let _ = out_file.write(b"\n");
        let mut ti = 0usize;
        for this_term in self.obj.terms.iter() {
            this_val = this_term.value;
            let _ = out_file.write(format!("{0},{1:.12e}", ti, this_val).as_bytes());
            for fld_str in include_fields.iter_mut() {
                if fld_str.s == "category" {
                    let _ = out_file.write(format!(",{}", this_term.category.s).as_bytes());
                }
                else if fld_str.s == "operator" {
                    let _ = out_file.write(format!(",{}", this_term.optr.s).as_bytes());
                }
                else if fld_str.s == "component" {
                    let _ = out_file.write(format!(",{}", this_term.component).as_bytes());
                }
                else if fld_str.s == "layer" {
                    let _ = out_file.write(format!(",{}", this_term.layer).as_bytes());
                }
                else if fld_str.s == "coefficient" {
                    let _ = out_file.write(format!(",{0:.12e}", this_term.coef).as_bytes());
                }
                else if fld_str.s == "exponent" {
                    let _ = out_file.write(format!(",{0:.12e}", this_term.expnt).as_bytes());
                }
                else if fld_str.s == "elementSet" {
                    let _ = out_file.write(format!(",{}", this_term.el_set_name.s).as_bytes());
                }
                else if fld_str.s == "nodeSet" {
                    let _ = out_file.write(format!(",{}", this_term.nd_set_name.s).as_bytes());
                }
                else if fld_str.s == "activeTime" {
                    let _ = out_file.write(format!(",{0:.12e}", this_term.active_time[0]).as_bytes());
                }
            }
            let _ = out_file.write(b"\n");
            ti += 1;
        }
        
        if write_grad {
            let dot_ind = file_name.find(".");
            let exten : CppStr;
            let rt_file : CppStr;
            if dot_ind < MAX_INT {
                rt_file = file_name.substr(0, dot_ind);
                exten = file_name.substr(dot_ind,MAX_INT);
            }
            else {
                rt_file = file_name.clone();
                exten = CppStr::from("");
            }

            let full_file = format!("{}_grad{}", rt_file.s, exten.s);
            let mut out_file = match File::create(&full_file) {
                Err(_why) => panic!("Error: could not open file, {}", full_file),
                Ok(file) => file,
            };

            let _ = out_file.write(b"designVariable,objGrad\n");
            num_dv = self.design_vars.len();
            for i1 in 0..num_dv {
                let _ = out_file.write(format!("{0},{1:.12e}\n", i1, self.d_ld_d[i1]).as_bytes());
            }
        }
        
    }   

}
