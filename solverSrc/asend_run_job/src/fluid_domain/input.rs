use std::ops::Sub;

use crate::constants::MAX_INT;
use crate::fluid_domain::*;
use crate::cpp_str::*;
use crate::file_util::*;
use crate::list_ent::*;

pub fn increment_ct(ct : usize) -> usize {
    match ct {
        MAX_INT => 0,
        _ => ct + 1,
    }
}

impl FluidDomain {

    pub fn read_domain_input(&mut self, file_name : &CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;

        let mut hd_updated : bool;

        let mut nd_ct = 0usize;
        let mut cl_ct = 0usize;
        let mut ns_ct = 0usize;
        let mut cs_ct = 0usize;
        let mut sd_ct = 0usize;
        let mut fl_ct = 0usize;
        let mut cnst_ct = 0usize;
        let mut ld_ct = 0usize;

        let mut isfnd = false;


        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                hd_updated = read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "nodes" => {if data_len == 4 {
                                    nd_ct += 1;
                                }},
                    "cells" => {if data_len == 5 {
                                    cl_ct += 1;
                                }},
                    "sets" => {if headings[2].s != "" && hd_updated {
                                   match headings[1].s.as_str() {
                                       "node" => ns_ct += 1,
                                       "cell" => cs_ct += 1,
                                       &_ => (),
                                   }
                               }},
                    "subDomains" => {if headings[1].s == "cellSet" && data_len == 1 {
                                         sd_ct += 1;
                                     }},
                    "fluids" => {if headings[1].s != "" && headings[2].s == "" && hd_updated {
                                     fl_ct += 1;
                                 }},
                    "constraints" => self.const_loop1(&mut cnst_ct, &headings, data_len),
                    "loads" => self.load_loop1(&mut ld_ct, &headings, data_len),
                    "initialState" => {if !isfnd {
                                           self.init_stat_file = file_name.clone();
                                           isfnd = true;
                                       }},
                    &_ => (),
                }
            }
        }

        self.nodes = vec![Node::new(); nd_ct];
        self.cells = vec![Cell::new(); cl_ct];
        self.node_sets = vec![Set::new(); ns_ct + nd_ct + 1];
        self.cell_sets = vec![Set::new(); cs_ct + cl_ct + 1];
        self.sub_domains = vec![SubDomain::new(); sd_ct];
        self.fluids = vec![Fluid::new(); fl_ct];
        self.constraints.const_vec = vec![Constraint::new(); cnst_ct];
        self.loads = vec![Load::new(); ld_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ns_ct = MAX_INT;
            cs_ct = MAX_INT;
            sd_ct = MAX_INT;
            fl_ct = MAX_INT;
            cnst_ct = MAX_INT;
            ld_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                hd_updated = read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "nodes" => {if data_len == 4 {
                                    nd_ct = data[0].stoi();
                                    self.nodes[nd_ct].label = nd_ct;
                                    for i in 0..3 {
                                        self.nodes[nd_ct].coord[i] = data[i+1].stod();
                                    }
                                }},
                    "cells" => {if data_len == 5 {
                                    cl_ct = data[0].stoi();
                                    self.cells[cl_ct].label = cl_ct;
                                    for i in 0..4 {
                                        self.cells[cl_ct].nodes[i] = data[i+1].stoi();
                                    }
                                }},
                    "sets" => {if headings[2].s != "" {
                                   if hd_updated {
                                       match headings[1].s.as_str() {
                                           "node" => {ns_ct = increment_ct(ns_ct);
                                                      self.node_sets[ns_ct].name = headings[2].clone();},
                                           "cell" => {cs_ct = increment_ct(cs_ct);
                                                      self.cell_sets[cs_ct].name = headings[2].clone();},
                                           &_ => (),
                                       }
                                   }
                                   if data_len == 1 {
                                       match headings[1].s.as_str() {
                                           "node" => self.node_sets[ns_ct].labels.push_back(data[0].stoi()),
                                           "cell" => self.cell_sets[cs_ct].labels.push_back(data[0].stoi()),
                                           &_ => (),
                                       };
                                   }
                               }},
                    "subDomains" => {if data_len == 1 {
                                         match headings[1].s.as_str() {
                                            "cellSet" => {sd_ct = increment_ct(sd_ct);
                                                          self.sub_domains[sd_ct].cell_set_name = data[0].clone();},
                                            "fluid" => {self.sub_domains[sd_ct].fluid_name = data[0].clone();},
                                            &_ => (),
                                         }
                                     }},
                    "fluids" => {if headings[1].s != "" && headings[2].s == "" && hd_updated {
                                     fl_ct = increment_ct(fl_ct);
                                     self.fluids[fl_ct].name = headings[1].clone();
                                 }
                                 if headings[2].s != "" {
                                    if data_len == 1 {
                                        match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                                            "viscosity" => self.fluids[fl_ct].viscosity = data[0].stod(),
                                            "thermalconductivity" => self.fluids[fl_ct].therm_cond = data[0].stod(),
                                            "thermalexpansion" => self.fluids[fl_ct].expansion = data[0].stod(),
                                            "thermalspecHeat" => self.fluids[fl_ct].spec_heat = data[0].stod(),
                                            "idealGasConst" => self.fluids[fl_ct].ideal_gas = data[0].stod(),
                                            "bulkModulus" => self.fluids[fl_ct].bulk_modulus = data[0].stod(),
                                            "refTemp" => self.fluids[fl_ct].ref_temp = data[0].stod(),
                                            "refPres" => self.fluids[fl_ct].ref_pres = data[0].stod(),
                                            "refDen" => self.fluids[fl_ct].ref_den = data[0].stod(),
                                            "refEnth" => self.fluids[fl_ct].ref_enth = data[0].stod(),
                                            "tempVisCoef" => self.fluids[fl_ct].temp_vis_coef = data[0].stod(),
                                            "turbVisCoef" => self.fluids[fl_ct].turb_vis_coef = data[0].stod(),
                                            "gradTurbCoef" => self.fluids[fl_ct].grad_turb_coef = data[0].stod(),
                                            "dissTurbCoef" => self.fluids[fl_ct].diss_turb_coef = data[0].stod(),
                                            "compressible" => self.fluids[fl_ct].compressible = data[0].s.contains("yes"),
                                            &_ => (),
                                        }
                                    }
                                 }},
                    "constraints" => self.const_loop2(&mut cnst_ct, &headings, &mut data, data_len),
                    "loads" => self.load_loop2(&mut ld_ct, &headings, &mut data, data_len),
                    &_ => (),
                }
            }
        }

    }

    pub fn const_loop1(&self, ct : &mut usize, headings : &Vec<CppStr>, data_len : usize) {
        if headings[1].s == "type" && data_len == 1 {
            *ct += 1;
        }
    }

    pub fn get_curr_term(&mut self, ct : usize) -> &mut ConstraintTerm {
        match self.constraints.const_vec[ct].terms.back_mut() {
            None => panic!("Error: tried to access the back constraint term of an empty list"),
            Some(x) => x,
        }
    }

    pub fn const_loop2(&mut self, ct : &mut usize, headings : &Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize) {
        if headings[1].s == "type" && data_len == 1 {
            *ct = increment_ct(*ct);
            self.constraints.const_vec[*ct].this_type = data[0].clone();
        }
        if headings[2].s != "" && data_len == 1 {
            match headings[2].s.as_str() {
                "nodeSet" => {let mut new_tm = ConstraintTerm::new();
                              new_tm.node_set = data[0].clone();
                              self.constraints.const_vec[*ct].terms.push_back(new_tm);},
                "dof" => self.get_curr_term(*ct).dof = data[0].stoi(),
                "coef" => self.get_curr_term(*ct).coef = data[0].stod(),
                &_ => (),
            }
        }
        if data_len == 2 {
            match headings[1].s.as_str() {
                "rhs" => self.constraints.const_vec[*ct].rhs.push_back(ConstTimePt { time: data[0].stod(), value: data[1].stod()}),
                "activeTime" => {self.constraints.const_vec[*ct].active_time[0] = data[0].stod();
                                 self.constraints.const_vec[*ct].active_time[1] = data[1].stod();},
                &_ => (),
            }
        }
    }

    pub fn read_constraint_input(&mut self, file_name : &CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;

        let mut cnst_ct = 0usize;

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "constraints" => self.const_loop1(&mut cnst_ct, &headings, data_len),
                    &_ => (),
                }
            }
        }

        self.constraints.const_vec = vec![Constraint::new(); cnst_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            cnst_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "constraints" => self.const_loop2(&mut cnst_ct, &headings, &mut data, data_len),
                    &_ => (),
                }
            }
        }

    }

    pub fn load_loop1(&self, ct : &mut usize, headings : &Vec<CppStr>, data_len : usize) {
        if headings[1].s == "type" && data_len == 1 {
            *ct += 1;
        }
    }

    pub fn load_loop2(&mut self, ct : &mut usize, headings : &Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize) {
        if data_len == 1 {
            match headings[1].s.as_str() {
                "type" => {*ct = increment_ct(*ct);
                           self.loads[*ct].this_type = data[0].clone();},
                "cellSet" => self.loads[*ct].cell_set = data[0].clone(),
                "angularVel" => self.loads[*ct].angular_vel = data[0].stod(),
                &_ => (),
            }
        }
        if data_len == 2 && headings[1].s == "activeTime" {
            self.loads[*ct].active_time[0] = data[0].stod();
            self.loads[*ct].active_time[1] = data[1].stod();
        }
        if data_len == 3 {
            match headings[1].s.as_str() {
                "center" => {self.loads[*ct].center[0] = data[0].stod();
                             self.loads[*ct].center[1] = data[1].stod();
                             self.loads[*ct].center[2] = data[2].stod();},
                "axis" => {self.loads[*ct].axis[0] = data[0].stod();
                           self.loads[*ct].axis[1] = data[1].stod();
                           self.loads[*ct].axis[2] = data[2].stod();},
                &_ => (),
            }
        }
        if data_len == 4 && headings[1].s == "load" {
            let newpt = QuadFloat { f1: data[0].stod(), f2: data[1].stod(), f3: data[2].stod(), f4: data[3].stod() };
            self.loads[*ct].load.push_back(newpt);
        }
    }

    pub fn read_load_input(&mut self, file_name : &CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;

        let mut ld_ct = 0usize;

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "loads" => self.load_loop1(&mut ld_ct, &headings, data_len),
                    &_ => (),
                }
            }
        }

        self.loads = vec![Load::new(); ld_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ld_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "loads" => self.load_loop2(&mut ld_ct, &headings, &mut data, data_len),
                    &_ => (),
                }
            }
        }

    }

    pub fn read_initial_state(&mut self) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;

        let mut nsi : usize;
        let mut in_dat = [0f64; 6];
        let mut this_nd : &mut Node;

        if let Ok(lines) = read_lines(self.init_stat_file.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "initialState" && data_len == 7 {
                    nsi = match self.ns_map.get(&data[0].s) {
                        None => panic!("Error: no node set named {}, in reading of initial state", data[0].s),
                        Some(x) => *x,
                    };
                    for i in 0..6 {
                        in_dat[i] = data[i+1].stod();
                    }
                    match headings[1].s.as_str() {
                        "flow" => {for ndi in self.node_sets[nsi].labels.iter() {
                                       this_nd = &mut self.nodes[*ndi];
                                       this_nd.initial_fl_den = in_dat[0];
                                       for i in 0..3 {
                                           this_nd.initial_fl_vel[i] = in_dat[i+1];
                                       }
                                       this_nd.initial_temperature = in_dat[4];
                                       this_nd.initial_turb_e = in_dat[5];
                                   }},
                        "flowDot" => {for ndi in self.node_sets[nsi].labels.iter() {
                                          this_nd = &mut self.nodes[*ndi];
                                          this_nd.initial_fl_den_dot = in_dat[0];
                                          for i in 0..3 {
                                              this_nd.initial_fl_vel_dot[i] = in_dat[i+1];
                                          }
                                          this_nd.initial_temp_dot = in_dat[4];
                                          this_nd.initial_turb_e_dot = in_dat[5];
                                      }},
                        &_ => (),
                    }
                }
            }
        }
    }

    pub fn read_des_var_input(&mut self, file_name : &CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;

        let mut dv_ct = 0usize;

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[1].s == "category" && data_len == 1 {
                    dv_ct += 1;
                }
            }
        }

        self.design_vars = vec![DesignVariable::new(); dv_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            dv_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if data_len == 1 {
                    match headings[1].s.as_str() {
                        "category" => {dv_ct = increment_ct(dv_ct);
                                       self.design_vars[dv_ct].category = data[0].clone()},
                        "component" => self.design_vars[dv_ct].component = data[0].stoi(),
                        "cellSet" => self.design_vars[dv_ct].cell_set_name = data[0].clone(),
                        "nodeSet" => self.design_vars[dv_ct].nd_set_name = data[0].clone(),
                        "coefficients" => self.design_vars[dv_ct].coefs.push_back(data[0].stod()),
                        &_ => (),
                    }
                }
                else if data_len == 2 && headings[1].s == "activeTime" {
                    self.design_vars[dv_ct].active_time[0] = data[0].stod();
                    self.design_vars[dv_ct].active_time[1] = data[1].stod();
                }
            }
        }

    }

}