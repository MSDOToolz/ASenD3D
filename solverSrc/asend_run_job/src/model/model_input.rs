use crate::model::*;
use crate::cpp_str::CppStr;
use crate::list_ent::{DualFloat, QuadFloat};

use crate::file_util::*;
use std::fs::File;
use std::io::{self, Read, BufRead};
use std::path::Path;
use std::collections::LinkedList;
use std::usize::MAX;

fn increment_ct(ct : usize) -> usize {
    match ct {
        MAX_INT => 0,
        _ => ct + 1,
    }
}

impl Model {

    pub fn read_job(&mut self, file_name : &mut CppStr) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        let mut cmd_ct : usize = 0usize;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            cmd_ct = 0;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[1].s == "command" && data_len == 1 {
                    cmd_ct += 1usize;
                }
            }
        }
        
        self.job = vec![JobCommand::new(); cmd_ct];

        let all_last = String::from("all_last");
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            cmd_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if data_len == 1 {
                    match headings[1].s.as_str() {
                        "command" => {cmd_ct = increment_ct(cmd_ct);
                                      self.job[cmd_ct].cmd_string = data[0].clone();},
                        "fileName" => self.job[cmd_ct].file_name = data[0].clone(),
                        "dynamic" => self.job[cmd_ct].dynamic = data[0].s.contains("yes"),
                        "explicit" => self.job[cmd_ct].explicit = data[0].s.contains("yes"),
                        "elastic" => self.job[cmd_ct].elastic = data[0].s.contains("yes"),
                        "loadRampSteps" => self.job[cmd_ct].load_ramp_steps = data[0].stoi(),
                        "newmarkBeta" => self.job[cmd_ct].newmark_beta = data[0].stod(),
                        "newmarkGamma" => self.job[cmd_ct].newmark_gamma = data[0].stod(),
                        "nonlinearGeom" => self.job[cmd_ct].nonlinear_geom = data[0].s.contains("yes"),
                        "maxNonlinIterations" => self.job[cmd_ct].max_nl_it = data[0].stoi(),
                        "nonlinConvTol" => self.job[cmd_ct].nl_conv_tol = data[0].stod(),
                        "abortNonlinDiv" => self.job[cmd_ct].abort_nl = data[0].s.contains("yes"),
                        "constScaleFactor" => self.job[cmd_ct].const_scale_factor = data[0].stod(),
                        "enforceMaxCon" => self.job[cmd_ct].enforce_max_c = data[0].s.contains("yes"),
                        "saveSolnHist" => self.job[cmd_ct].save_soln_hist = data[0].s.contains("yes"),
                        "solnHistFreq" => self.job[cmd_ct].soln_hist_freq = data[0].stoi(),
                        "solnHistDir" => self.job[cmd_ct].file_name = data[0].clone(),
                        "lumpMass" => self.job[cmd_ct].lump_mass = data[0].s.contains("yes"),
                        "simPeriod" => self.job[cmd_ct].sim_period = data[0].stod(),
                        "solverBandwidth" => self.job[cmd_ct].solver_bandwidth = data[0].stoi(),
                        "solverBlockDim" => self.job[cmd_ct].solver_block_dim = data[0].stoi(),
                        "solverMethod" => self.job[cmd_ct].solver_method = data[0].clone(),
                        "maxIterations" => self.job[cmd_ct].max_it = data[0].stoi(),
                        "convergenceTol" => self.job[cmd_ct].conv_tol = data[0].stod(),
                        "staticLoadTime" => self.job[cmd_ct].static_load_time.push_back(data[0].stod()),
                        "thermal" => self.job[cmd_ct].thermal = data[0].s.contains("yes"),
                        "diffusion" => self.job[cmd_ct].diffusion = data[0].s.contains("yes"),
                        "timeStep" => self.job[cmd_ct].time_step = data[0].stod(),
                        "userUpdate" => self.job[cmd_ct].run_user_update = data[0].s.contains("yes"),
                        "type" => self.job[cmd_ct].this_type = data[0].clone(),
                        "numModes" => self.job[cmd_ct].num_modes = data[0].stoi(),
                        "targetEigenvalue" => self.job[cmd_ct].tgt_eval = data[0].stod(),
                        "solnField" => self.job[cmd_ct].soln_field = data[0].clone(),
                        "mode" => self.job[cmd_ct].mode = data[0].stoi(),
                        "maxAmplitude" => self.job[cmd_ct].max_amplitude = data[0].stod(),
                        "nodeSet" => self.job[cmd_ct].node_set = data[0].clone(),
                        "fields" => self.job[cmd_ct].fields.push_back(data[0].clone()),
                        "timeSteps" => {if all_last.contains(data[0].s.as_str()) {
                                            self.job[cmd_ct].time_step_tag = data[0].clone();
                                        }
                                        else {
                                            match data[0].is_int() {
                                                true => self.job[cmd_ct].time_steps.push_back(data[0].stoi()),
                                                false => println!("Warning: possible invalid entry for {} timeSteps in job file {}", headings[0].s, file_name.s),
                                            }
                                        }},
                        "elementSet" => self.job[cmd_ct].element_set = data[0].clone(),
                        "position" => self.job[cmd_ct].position = data[0].clone(),
                        "writeModes" => self.job[cmd_ct].write_modes = data[0].s.contains("yes"),
                        "properties" => self.job[cmd_ct].properties.push_back(data[0].clone()),
                        "include" => self.job[cmd_ct].obj_include.push_back(data[0].clone()),
                        "writeGradient" => self.job[cmd_ct].write_gradient = data[0].s.contains("yes"),
                        &_ => (),
                    }
                }
                else if data_len > 1 {
                    match headings[1].s.as_str() {
                        "timeSteps" => {i2 = data[0].stoi();
                                        i3 = data[1].stoi();
                                        i4 = match data_len {
                                            3 => data[2].stoi(),
                                            _ => 1,
                                        };
                                        i1 = i2;
                                        while i1 < i3 {
                                            self.job[cmd_ct].time_steps.push_back(i1);
                                            i1 += i4;
                                        }},
                        &_ => (),
                    }
                }
            }
        }
        else {
            panic!("Error: could not open job file {}", file_line.s);
        }
        
    }

    pub fn read_model_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut hd_updated : bool;
        let mut i1 : usize;
        let mut i2 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut doub_inp : [f64; 10] = [0.0; 10];
        let mut int_inp : [usize; 10] = [0,0,0,0,0,0,0,0,0,0];
        let mut el_type : usize = 0usize;
        let mut lay_pt : &mut Layer;
        
        let mut nd_ct : usize =  0;
        let mut el_ct : usize =  0;
        let mut ns_ct : usize =  0;
        let mut es_ct : usize =  0;
        let mut sec_ct : usize =  0;
        let mut mat_ct : usize =  0;
        let mut const_ct = [0usize; 4];
        let mut load_ct = [0usize; 4];
        let mut int_ct : usize = 0;
        let mut ps_ct : usize = 0;

        let mut lst_ar = vec![CppStr::new(); 4];
        lst_ar[0] = CppStr::from("nodalForce bodyForce gravitational centrifugal surfacePressure surfaceTraction");
        lst_ar[1] = CppStr::from("nodalHeatGen bodyHeatGen surfaceFlux");
        lst_ar[2] = CppStr::from("nodalMassGen massGen massFlux");
        let mut load_type : CppStr = CppStr::new();
        let mut ld_ind = 0usize;

        let mut const_ind = 0usize;
        let mut con_type : CppStr = CppStr::new();
        let mut all_types = CppStr::from("displacement temperature concentration fluid");
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                hd_updated = read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                match format!("{}{}", headings[0].s, headings[1].s).as_str() {
                    "nodes" => {if data_len == 4 {
                                    nd_ct += 1;
                                }},
                    "elementsconnectivity" => {if data_len > 1 {
                                                   el_ct += 1;
                                               }},
                    "setsnode" => {if headings[2].s != "" && hd_updated {
                                       ns_ct += 1;
                                   }},
                    "setselement" => {if headings[2].s != "" && hd_updated {
                                          es_ct += 1;
                                      }},
                    "sectionstype" => {if data_len == 1 {
                                           sec_ct += 1;
                                       }},
                    &_ => {if headings[0].s == "materials" {
                               if headings[1].s != "" && headings[2].s == "" && hd_updated {
                                   mat_ct += 1;
                               }
                           }
                           else {
                               self.const_loop1(&mut headings, &mut data, data_len, &mut const_ct);
                               self.load_loop1(&mut headings, &mut data, data_len, &mut load_ct, &mut lst_ar);
                               self.interaction_loop1(&mut headings, data_len, &mut int_ct);
                               self.ps_loop1(&headings, data_len, &mut ps_ct);
                               if headings[0].s == "initialState" {
                                   self.init_stat_file = file_name.clone();
                               }
                           }},
                }
            }
        }
        
        self.nodes = vec![Node::new(); nd_ct];
        self.elements = vec![Element::new(); el_ct];
        self.node_sets = vec![Set::new(); ns_ct + nd_ct + 1];
        self.element_sets = vec![Set::new(); es_ct + el_ct + 1];
        self.sections = vec![Section::new(); sec_ct];
        self.materials = vec![Material::new(); mat_ct];

        self.elastic_const.const_vec = vec![Constraint::new(); const_ct[0]];
        self.thermal_const.const_vec = vec![Constraint::new(); const_ct[1]];
        self.diff_const.const_vec = vec![Constraint::new(); const_ct[2]];

        self.elastic_loads = vec![Load::new(); load_ct[0]];
        self.thermal_loads = vec![Load::new(); load_ct[1]];
        self.diff_loads = vec![Load::new(); load_ct[2]];

        self.interactions.int_vec = vec![Interaction::new(); int_ct];
        self.particle_sources = vec![ParticleSource::new(); ps_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ns_ct = MAX_INT;
            es_ct = MAX_INT;
            sec_ct = MAX_INT;
            mat_ct = MAX_INT;
            for i in 0..4 {
                const_ct[i] = MAX_INT;
                load_ct[i] = MAX_INT;
            }
            int_ct = MAX_INT;
            ps_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                hd_updated = read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                match headings[0].s.as_str() {
                    "nodes" => {if data_len == 4 {
                                    nd_ct = data[0].stoi();
                                    self.nodes[nd_ct].label = nd_ct;
                                    doub_inp[0] = data[1].stod();
                                    doub_inp[1] = data[2].stod();
                                    doub_inp[2] = data[3].stod();
                                    self.nodes[nd_ct].set_crd(&mut doub_inp);
                                }},
                    "elements" => {if headings[1].s == "type" && data_len == 1 {
                                       el_type = match data[0].s.as_str() {
                                           "tet4" => 4,
                                           "wedge6" => 6,
                                           "brick8" => 8,
                                           "brickIM" => 81,
                                           "tet10" => 10,
                                           "shell3" => 3,
                                           "shell4" => 41,
                                           "beam2" => 2,
                                           "frcFld" => 21,
                                           "mass" => 1,
                                           &_ => panic!("Error: unrecognized element type: {}", data[0].s),
                                       };
                                   }
                                   else if headings[1].s == "connectivity" && data_len > 1 {
                                       el_ct = data[0].stoi();
                                       let new_el = &mut self.elements[el_ct];
                                       new_el.initialize_type(el_type);
                                       new_el.label = el_ct;
                                       for i1 in 1..data_len {
                                           int_inp[i1-1] = data[i1].stoi();
                                       }
                                       new_el.set_nodes(&mut int_inp);
                                   }},
                    "sets" => {if headings[2].s != "" {
                                   if hd_updated {
                                       match headings[1].s.as_str() {
                                           "node" => {ns_ct = increment_ct(ns_ct);
                                                      self.node_sets[ns_ct].name = headings[2].clone();},
                                           "element" => {es_ct = increment_ct(es_ct);
                                                         self.element_sets[es_ct].name = headings[2].clone();},
                                           &_ => {},
                                       }
                                   }
                                   else if data_len == 1 {
                                       match headings[1].s.as_str() {
                                           "node" => self.node_sets[ns_ct].labels.push_back(data[0].stoi()),
                                           "element" => self.element_sets[es_ct].labels.push_back(data[0].stoi()),
                                           &_ => {},
                                       }
                                   }
                               }},
                    "sections" => {if data_len == 1 {
                                       match format!("{}{}{}", headings[1].s, headings[2].s, headings[3].s).as_str() {
                                           "type" => {sec_ct = increment_ct(sec_ct);
                                                      self.sections[sec_ct].this_type = data[0].clone();},
                                           "material" => self.sections[sec_ct].mat_name = data[0].clone(),
                                           "layupzOffset" => self.sections[sec_ct].z_offset = data[0].stod(),
                                           "layuplayersmaterial" => {let mut new_lay = Layer::new();
                                                                     new_lay.mat_name = data[0].clone();
                                                                     self.sections[sec_ct].layers.push_back(new_lay);},
                                           "layuplayersthickness" => {lay_pt = match self.sections[sec_ct].layers.back_mut() {
                                                                          None => panic!("could not access back of section layer list in input reader"),
                                                                          Some(x) => x,
                                                                      };
                                                                      lay_pt.thickness = data[0].stod();},
                                           "layuplayersangle" => {lay_pt = match self.sections[sec_ct].layers.back_mut() {
                                                                      None => panic!("could not access back of section layer list in input reader"),
                                                                      Some(x) => x,
                                                                  };
                                                                  lay_pt.angle = data[0].stod();},
                                           "beamPropsarea" => self.sections[sec_ct].area = data[0].stod(),
                                           "beamPropsJ" => self.sections[sec_ct].polar_moment = data[0].stod(),
                                           "beamPropsconductivity" => self.sections[sec_ct].conductivity = data[0].stod(),
                                           "beamPropsspecHeat" => self.sections[sec_ct].spec_heat = data[0].stod(),
                                           "potFieldcoef" => self.sections[sec_ct].pot_coef = data[0].stod(),
                                           "potFieldexp" => self.sections[sec_ct].pot_exp = data[0].stod(),
                                           "dampFieldcoef" => self.sections[sec_ct].damp_coef = data[0].stod(),
                                           "dampFielddistExp" => self.sections[sec_ct].damp_dist_exp = data[0].stod(),
                                           "dampFieldvelExp" => self.sections[sec_ct].damp_vel_exp = data[0].stod(),
                                           "magFieldcoef" => self.sections[sec_ct].mag_coef = data[0].stod(),
                                           "magFielddistExp" => self.sections[sec_ct].mag_dist_exp = data[0].stod(),
                                           "magFieldvelExp" => self.sections[sec_ct].mag_vel_exp = data[0].stod(),
                                           "thermFieldcondCoef" => self.sections[sec_ct].cond_coef = data[0].stod(),
                                           "thermFieldradCoef" => self.sections[sec_ct].rad_coef = data[0].stod(),
                                           "thermFieldrefTemp" => self.sections[sec_ct].ref_temp = data[0].stod(),
                                           "massPerEl" => self.sections[sec_ct].mass_per_el = data[0].stod(),
                                           "specHeat" => self.sections[sec_ct].spec_heat = data[0].stod(),
                                           "elementSet" => self.sections[sec_ct].el_set_name = data[0].clone(),
                                           &_ => (),
                                       }
                                   }
                                   else if data_len == 3 {
                                       match format!("{}{}", headings[1].s, headings[2].s).as_str() {
                                           "beamPropsstiffness" => {int_inp[0] = data[0].stoi() - 1;
                                                                    int_inp[1] = data[1].stoi() - 1;
                                                                    doub_inp[0] = data[2].stod();
                                                                    self.sections[sec_ct].set_stiffness(int_inp[0],  int_inp[1],  doub_inp[0]);},
                                           "beamPropsmass" => {int_inp[0] = data[0].stoi() - 1;
                                                               int_inp[1] = data[1].stoi() - 1;
                                                               doub_inp[0] = data[2].stod();
                                                               self.sections[sec_ct].set_mass(int_inp[0], int_inp[1], doub_inp[0]);},
                                           "beamPropsdamping" => {int_inp[0] = data[0].stoi() - 1;
                                                                  int_inp[1] = data[1].stoi() - 1;
                                                                  doub_inp[0] = data[2].stod();
                                                                  self.sections[sec_ct].set_damping(int_inp[0],   int_inp[1],   doub_inp[0]);},
                                           &_ => (),
                                       }
                                   }
                                   else if data_len == 5 {
                                       match format!("{}{}", headings[1].s, headings[2].s).as_str() {
                                           "beamPropsI" => {for i1 in 0..5 {
                                                                doub_inp[i1] = data[i1].stod();
                                                            }
                                                            self.sections[sec_ct].set_area_moment(&mut doub_inp);},
                                           &_ => (),
                                       }
                                   }
                                   else if data_len == 6 {
                                       match format!("{}{}", headings[1].s, headings[2].s).as_str() {
                                           "orientation" => {for i1 in 0..6 {
                                                                 doub_inp[i1] = data[i1].stod();
                                                             }
                                                             self.sections[sec_ct].set_orientation(&mut doub_inp);},
                                           "beamPropsexpLoadCoef" => {for i1 in 0..6 {
                                                                          doub_inp[i1] = data[i1].stod();
                                                                      }
                                                                      self.sections[sec_ct].set_exp_ld(&mut doub_inp);},
                                           &_ => (),
                                       }
                                   }},
                    "materials" => {if headings[1].s != "" {
                                        if headings[2].s == "" && hd_updated {
                                            mat_ct = increment_ct(mat_ct);
                                            self.materials[mat_ct].name = headings[1].clone();
                                        }
                                        if data_len == 1 {
                                            match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                                                "density" => self.materials[mat_ct].density = data[0].stod(),
                                                "elasticE" => {let this_mat = &mut self.materials[mat_ct];
                                                               this_mat.modulus[0] = data[0].stod();
                                                               this_mat.modulus[1] = this_mat.modulus[0];
                                                               this_mat.modulus[2] = this_mat.modulus[0];},
                                                "elasticnu" => {let this_mat = &mut self.materials[mat_ct];
                                                                this_mat.poisson_ratio[0] = data[0].stod();
                                                                this_mat.poisson_ratio[1] = this_mat.poisson_ratio[0];
                                                                this_mat.poisson_ratio[2] = this_mat.poisson_ratio[0];},
                                                "elasticG" => {let this_mat = &mut self.materials[mat_ct];
                                                               this_mat.shear_mod[0] = data[0].stod();
                                                               this_mat.shear_mod[1] = this_mat.shear_mod[0];
                                                               this_mat.shear_mod[2] = this_mat.shear_mod[0];},
                                                "thermalspecHeat" => self.materials[mat_ct].spec_heat = data[0].stod(),
                                                "diffusionmaxConcentration" => self.materials[mat_ct].max_concentration = data[0].stod(),
                                                &_ => (),
                                            }
                                        }
                                        else if data_len == 3 {
                                            match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                                                "elasticE" => {let this_mat = &mut self.materials[mat_ct];
                                                               this_mat.modulus[0] = data[0].stod();
                                                               this_mat.modulus[1] = data[1].stod();
                                                               this_mat.modulus[2] = data[2].stod();},
                                                "elasticnu" => {let this_mat = &mut self.materials[mat_ct];
                                                                this_mat.poisson_ratio[0] = data[0].stod();
                                                                this_mat.poisson_ratio[1] = data[1].stod();
                                                                this_mat.poisson_ratio[2] = data[2].stod();},
                                                "elasticG" => {let this_mat = &mut self.materials[mat_ct];
                                                               this_mat.shear_mod[0] = data[0].stod();
                                                               this_mat.shear_mod[1] = data[1].stod();
                                                               this_mat.shear_mod[2] = data[2].stod();},
                                                "elasticstiffness" => {int_inp[0] = data[0].stoi() - 1;
                                                                       int_inp[1] = data[1].stoi() - 1;
                                                                       doub_inp[0] = data[2].stod();
                                                                       self.materials[mat_ct].set_stiffness(int_inp[0], int_inp[1], doub_inp[0]);},
                                                "damping" => {int_inp[0] = data[0].stoi() - 1;
                                                              int_inp[1] = data[1].stoi() - 1;
                                                              doub_inp[0] = data[2].stod();
                                                              self.materials[mat_ct].set_damping(int_inp[0], int_inp[1], doub_inp[0]);},
                                                &_ => (),
                                            }
                                        }
                                        else if data_len == 6 {
                                            match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                                                "thermalconductivity" => {let this_mat = &mut self.materials[mat_ct];
                                                                          for i1 in 0..6 {
                                                                              this_mat.conductivity[i1] = data[i1].stod();
                                                                          }},
                                                "thermalexpansion" => {let this_mat = &mut self.materials[mat_ct];
                                                                       for i1 in 0..6 {
                                                                           this_mat.expansion[i1] = data[i1].stod();
                                                                       }},
                                                "diffusiondiffusivity" => {let this_mat = &mut self.materials[mat_ct];
                                                                           for i1 in 0..6 {
                                                                               this_mat.diffusivity[i1] = data[i1].stod();
                                                                           }},
                                                "diffusionexpansion" => {let this_mat = &mut self.materials[mat_ct];
                                                                         for i1 in 0..6 {
                                                                             this_mat.diff_exp[i1] = data[i1].stod();
                                                                         }},
                                                &_ => (),
                                            }
                                        }
                                        if headings[2].s == "custom" && headings[3].s != "" && data_len > 0 {
                                            let mut dat_vec = vec![0.0f64; data_len];
                                            for di in 0..data_len {
                                                dat_vec[di] = match data[di].s.parse::<f64>() {
                                                    Err(_why) => panic!("Error: could not parse input for custom material property '{}' to float type", data[di].s),
                                                    Ok(x) => x,
                                                };
                                            }
                                            self.materials[mat_ct].custom.insert(headings[3].s.clone(),dat_vec);
                                        }
                                    }},
                    &_ => {self.const_loop2(&mut headings, &mut data, data_len, &mut const_ct, &mut const_ind, &mut con_type, &mut all_types);
                           self.load_loop2(&mut headings, &mut data, data_len, &mut load_ct, &mut lst_ar, &mut ld_ind, &mut load_type);
                           self.interaction_loop2(&mut headings, &mut data, data_len, &mut int_ct);
                           self.ps_loop2(&headings, &mut data, data_len, &mut ps_ct);},
                }
            }
        } 
        else {
            panic!("Error: could not open Model input file: {}", file_name.s);
        }
        
        // create all and individual nd/element sets
        
        i1 = self.nodes.len();
        for i2 in 0..i1 {
            if ns_ct == MAX_INT {
                ns_ct = 0usize;
            }
            else {
                ns_ct += 1usize;
            }
            self.node_sets[ns_ct].name.s = i2.to_string();
            self.node_sets[ns_ct].labels.push_back(i2);
        }
        
        ns_ct += 1usize;
        self.node_sets[ns_ct].name = CppStr::from("all");
        let labs : &mut LinkedList<usize> = &mut self.node_sets[ns_ct].labels;
        for i2 in 0..i1 {
            labs.push_back(i2);
        }
        
        i1 = self.elements.len();
        for i2 in 0..i1 {
            if es_ct == MAX_INT {
                es_ct = 0usize;
            }
            else {
                es_ct += 1usize;
            }
            self.element_sets[es_ct].name.s = i2.to_string();
            self.element_sets[es_ct].labels.push_back(i2);
        }
        
        es_ct += 1usize;
        let new_set : &mut Set = &mut self.element_sets[es_ct];
        new_set.name.s = "all".to_string();
        for i2 in 0..i1 {
            new_set.labels.push_back(i2);
        }
        
        // populate Set arrays
        
        i2 = 0;
        for ns in self.node_sets.iter_mut() {
            self.ns_map.insert(ns.name.to_string(), i2);
            i2 += 1usize;
        }
        
        i2 = 0;
        for es in self.element_sets.iter_mut() {
            self.es_map.insert(es.name.to_string(), i2);
            i2 += 1usize;
        }
        
        return;
    }

    pub fn get_curr_constraint(&mut self, curr_type : &CppStr, curr_ct : usize) -> &mut Constraint {
        match curr_type.s.as_str() {
            "displacement" => &mut self.elastic_const.const_vec[curr_ct],
            "temperature" => &mut self.thermal_const.const_vec[curr_ct],
            "concentration" => &mut self.diff_const.const_vec[curr_ct],
            &_ => panic!("Error, unrecognized constraint type, {}", curr_type.s),
        }
    }

    pub fn const_loop1(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, ct_ar : &mut [usize]) {
        if headings[0].s == "constraints" {
            if headings[1].s == "type" && data_len == 1 {
                match data[0].s.as_str() {
                    "displacement" => ct_ar[0] += 1,
                    "temperature" => ct_ar[1] += 1,
                    "concentration" => ct_ar[2] += 1,
                    &_ => (),
                }
            }
        }
    }

    pub fn const_loop2(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, ct_ar : &mut [usize], ct_ind : &mut usize, curr_type : &mut CppStr, all_types : &mut CppStr) {
        let mut flt_in = [0f64; 2];
        
        if headings[0].s == "constraints" {
            if data_len == 1 {
                match format!("{}{}", headings[1].s, headings[2].s).as_str() {
                    "type" => {*ct_ind = match data[0].s.as_str() {
                               "displacement" => 0,
                               "temperature" => 1,
                               "concentration" => 2,
                               &_ => panic!("Error: {} is not a valid constraint type. Allowable values are {}", data[0].s, all_types.s),
                               };
                               ct_ar[*ct_ind] = increment_ct(ct_ar[*ct_ind]);
                               curr_type.s = data[0].s.clone();
                               self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).this_type = curr_type.clone();},
                    "termsnodeSet" => {let mut new_cn = ConstraintTerm::new();
                                       new_cn.node_set = data[0].clone();
                                       self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).terms.push_back(new_cn);},
                    "termsdof" => {match self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).terms.back_mut() {
                                       None => {panic!("failed to access back of constraint terms list");},
                                       Some(x) => {x.dof = data[0].stoi();},
                                   }},
                    "termscoef" => {match self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).terms.back_mut() {
                                        None => {panic!("failed to access back of constraint terms list");},
                                        Some(x) => {x.coef = data[0].stod();},
                                    }},
                    "rhs" => {self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).rhs.push_back(ConstTimePt {time : 0.0, value : data[0].stod()});
                              self.get_curr_constraint(curr_type, ct_ar[*ct_ind]).rhs.push_back(ConstTimePt {time : 1.0e+100, value : data[0].stod()});},
                    &_ => (),
                }
            }
            else if data_len == 2 {
                match format!("{}{}", headings[1].s, headings[2].s).as_str() {
                    "rhs" => self.get_curr_constraint(&curr_type, ct_ar[*ct_ind]).rhs.push_back(ConstTimePt {time : data[0].stod(), value : data[1].stod()}),
                    "activeTime" => {flt_in[0] = data[0].stod();
                                     flt_in[1] = data[1].stod();
                                     self.get_curr_constraint(&curr_type, ct_ar[*ct_ind]).set_act_time(&flt_in);},
                    &_ => (),
                }
            }
        }
    }

    pub fn read_constraint_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut ct_ar = [0usize; 4];
        let mut ct_ind = 0usize;
        let mut curr_type : CppStr = CppStr::new();
        let mut all_types = CppStr::from("displacement temperature concentration fluid");
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.const_loop1(&mut headings, &mut data, data_len, &mut ct_ar);
            }
        }

        if ct_ar[0] + ct_ar[1] + ct_ar[2] + ct_ar[3] == 0 {
            return;
        }
        
        self.elastic_const.const_vec = vec![Constraint::new(); ct_ar[0]];
        self.thermal_const.const_vec = vec![Constraint::new(); ct_ar[1]];
        self.diff_const.const_vec = vec![Constraint::new(); ct_ar[2]];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ct_ar[0] = MAX_INT;
            ct_ar[1] = MAX_INT;
            ct_ar[2] = MAX_INT;
            ct_ar[3] = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                self.const_loop2(&mut headings, &mut data, data_len, &mut ct_ar, &mut ct_ind, &mut curr_type, &mut all_types);
            }
        } else {
            panic!("Error: could not open Constraint input file: {}",file_name.s);
        }
        
        return;
    }

    pub fn get_curr_ld(&mut self, curr_type : &CppStr, curr_ct : usize) -> &mut Load {
        match curr_type.s.as_str() {
            "elastic" => &mut self.elastic_loads[curr_ct],
            "thermal" => &mut self.thermal_loads[curr_ct],
            "diffusion" => &mut self.diff_loads[curr_ct],
            &_ => panic!("Error: unrecognized load type '{}' in get_curr_ld()", curr_type.s),
        }
    }

    pub fn load_loop1(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, ct_ar : &mut [usize], lst_ar : &mut [CppStr]) {
        let mut i1 : usize;
        if headings[0].s == "loads" {
            if headings[1].s == "type" && data_len == 1 {
                if lst_ar[0].s.contains(data[0].s.as_str()) {
                    ct_ar[0] += 1usize;
                }
                if lst_ar[1].s.contains(data[0].s.as_str()) {
                    ct_ar[1] += 1usize;
                }
                if lst_ar[2].s.contains(data[0].s.as_str()) {
                    ct_ar[2] += 1usize;
                }
            }
        }
    }

    pub fn load_loop2(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, ct_ar : &mut [usize], lst_ar : &mut [CppStr], ct_ind : &mut usize, curr_type : &mut CppStr) {
        let mut doub_inp = [0f64; 3];
        let mut i1 : usize;
        
        if headings[0].s == "loads" {
            if data_len == 1 {
                match headings[1].s.as_str() {
                    "type" => {if lst_ar[0].s.contains(data[0].s.as_str()) {
                                   ct_ar[0] = increment_ct(ct_ar[0]);
                                   *curr_type = CppStr::from("elastic");
                                   *ct_ind = 0;
                                   self.get_curr_ld(curr_type,ct_ar[*ct_ind]).this_type = data[0].clone();
                               }
                               if lst_ar[1].s.contains(data[0].s.as_str()) {
                                   ct_ar[1] = increment_ct(ct_ar[1]);
                                   *curr_type = CppStr::from("thermal");
                                   *ct_ind = 1;
                                   self.get_curr_ld(curr_type,ct_ar[*ct_ind]).this_type = data[0].clone();
                               }
                               if lst_ar[2].s.contains(data[0].s.as_str()) {
                                   ct_ar[2] = increment_ct(ct_ar[2]);
                                   *curr_type = CppStr::from("diffusion");
                                   *ct_ind = 2;
                                   self.get_curr_ld(curr_type, ct_ar[*ct_ind]).this_type = data[0].clone();
                               }},
                    "activeTime" => {doub_inp[0] = data[0].stod();
                                     doub_inp[1] = 1.0e+100;
                                     self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_act_time(&mut doub_inp);},
                    "nodeSet" => self.get_curr_ld(curr_type,ct_ar[*ct_ind]).node_set = data[0].clone(),
                    "elementSet" => self.get_curr_ld(curr_type,ct_ar[*ct_ind]).element_set = data[0].clone(),
                    "normTolerance" => self.get_curr_ld(curr_type,ct_ar[*ct_ind]).norm_tol = data[0].stod(),
                    "angularVelocity" => self.get_curr_ld(curr_type,ct_ar[*ct_ind]).angular_vel = data[0].stod(),
                    &_ => (),
                }
            }
            else if data_len == 2 {
                match headings[1].s.as_str() {
                    "activeTime" => {doub_inp[0] = data[0].stod();
                                     doub_inp[1] = data[1].stod();
                                     self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_act_time(&mut doub_inp);},
                    "load" => {let mut new_ld = LoadTimePt::new();
                               new_ld.time = data[0].stod();
                               new_ld.value[0] = data[1].stod();
                               self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_load(new_ld);},
                    &_ => (),
                }
            }
            else if data_len == 3 {
                match headings[1].s.as_str() {
                    "normDir" => {doub_inp[0] = data[0].stod();
                                  doub_inp[1] = data[1].stod();
                                  doub_inp[2] = data[2].stod();
                                  self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_norm_dir(&mut doub_inp);},
                    "center" => {doub_inp[0] = data[0].stod();
                                 doub_inp[1] = data[1].stod();
                                 doub_inp[2] = data[2].stod();
                                 self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_center(&mut doub_inp);},
                    "axis" => {doub_inp[0] = data[0].stod();
                               doub_inp[1] = data[1].stod();
                               doub_inp[2] = data[2].stod();
                               self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_axis(&mut doub_inp);}
                    &_ => (),
                }
            }
            else if data_len == 4 || data_len == 7 {
                match headings[1].s.as_str() {
                    "load" => {let mut new_ld = LoadTimePt::new();
                               new_ld.time = data[0].stod();
                               for i1 in 1..data_len {
                                   new_ld.value[i1-1] = data[i1].stod();
                               }
                               self.get_curr_ld(curr_type,ct_ar[*ct_ind]).set_load(new_ld);},
                    &_ => (),
                }
            }
        }
    }

    pub fn read_load_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut lst_ar = vec![CppStr::new(); 4];
        lst_ar[0] = CppStr::from("nodalForce bodyForce gravitational centrifugal surfacePressure surfaceTraction");
        lst_ar[1] = CppStr::from("nodalHeatGen bodyHeatGen surfaceFlux");
        lst_ar[2] = CppStr::from("nodalMassGen massGen massFlux");
        
        let mut ct_ar = [0usize; 4];
        let mut curr_type : CppStr = CppStr::new();
        let mut ct_ind = 0usize;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.load_loop1(&mut headings, &mut data, data_len, &mut ct_ar, &mut lst_ar);
            }
        }

        if ct_ar[0] + ct_ar[1] + ct_ar[2] + ct_ar[3] == 0 {
            return;
        }
        
        self.elastic_loads = vec![Load::new(); ct_ar[0]];
        self.thermal_loads = vec![Load::new(); ct_ar[1]];
        self.diff_loads = vec![Load::new(); ct_ar[2]];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ct_ar[0] = MAX_INT;
            ct_ar[1] = MAX_INT;
            ct_ar[2] = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                self.load_loop2(&mut headings, &mut data, data_len, &mut ct_ar, &mut lst_ar, &mut ct_ind, &mut curr_type);
            }
        } else {
            panic!("Error: could not open Load input file: {}", file_name.s);
        }
        
        return;
    }

    pub fn interaction_loop1(&mut self, headings : &mut Vec<CppStr>, data_len : usize, int_ct : &mut usize) {
        if headings[0].s == "interactions" && headings[1].s != "" {
            if headings[2].s == "nodeSet1" && data_len == 1 {
                *int_ct += 1;
            }
        }
    }

    pub fn interaction_loop2(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, int_ct : &mut usize) {
        if headings[0].s == "interactions" && headings[1].s != "" {
            if data_len == 1 {
                if headings[2].s == "nodeSet1" {
                    if *int_ct == MAX_INT {
                        *int_ct = 0;
                    }
                    else {
                        *int_ct += 1;
                    }
                    self.interactions.int_vec[*int_ct].name = headings[1].clone();
                }
                match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                    "nodeSet1" => self.interactions.int_vec[*int_ct].node_set1 = data[0].clone(),
                    "nodeSet2" => self.interactions.int_vec[*int_ct].node_set2 = data[0].clone(),
                    "potFieldexp" => self.interactions.int_vec[*int_ct].pot_exp = data[0].stod(),
                    "dampFielddistExp" => self.interactions.int_vec[*int_ct].damp_dist_exp = data[0].stod(),
                    "dampFieldvelExp" => self.interactions.int_vec[*int_ct].damp_vel_exp = data[0].stod(),
                    "magFielddistExp" => self.interactions.int_vec[*int_ct].mag_dist_exp = data[0].stod(),
                    "magFieldvelExp" => self.interactions.int_vec[*int_ct].mag_vel_exp = data[0].stod(),
                    "thermFieldcondCoef" => self.interactions.int_vec[*int_ct].cond_coef = data[0].stod(),
                    "thermFieldradCoef" => self.interactions.int_vec[*int_ct].rad_coef = data[0].stod(),
                    "thermFieldrefTemp" => self.interactions.int_vec[*int_ct].ref_temp = data[0].stod(),
                    "maxDistance" => self.interactions.int_vec[*int_ct].max_dist = data[0].stod(),
                    "maxNeighbors" => self.interactions.int_vec[*int_ct].max_nbrs = data[0].stoi(),
                    "maxDistRatio" => self.interactions.int_vec[*int_ct].max_ratio = data[0].stod(),
                    "idealGasConstant" => self.interactions.int_vec[*int_ct].ideal_gas = data[0].stod(),
                    "bulkModulus" => self.interactions.int_vec[*int_ct].bulk_mod = data[0].stod(),
                    "expansion" => self.interactions.int_vec[*int_ct].therm_exp = data[0].stod(),
                    "refDen" => self.interactions.int_vec[*int_ct].ref_den = data[0].stod(),
                    "refPres" => self.interactions.int_vec[*int_ct].ref_pres = data[0].stod(),
                    &_ => {},
                }
            }
            else if data_len == 2 {
                match format!("{}{}", headings[2].s, headings[3].s).as_str() {
                    "potFieldcoef" => {let new_ent = DualFloat{f1 : data[0].stod(), f2 : data[1].stod()};
                                       self.interactions.int_vec[*int_ct].pot_coef.push_back(new_ent);},
                    "dampFieldcoef" => {let new_ent = DualFloat{f1 : data[0].stod(), f2 : data[1].stod()};
                                        self.interactions.int_vec[*int_ct].damp_coef.push_back(new_ent);},
                    "magFieldcoef" => {let new_ent = DualFloat{f1 : data[0].stod(), f2 : data[1].stod()};
                                        self.interactions.int_vec[*int_ct].mag_coef.push_back(new_ent);},
                    "activeTime" => {self.interactions.int_vec[*int_ct].active_time[0] = data[0].stod();
                                     self.interactions.int_vec[*int_ct].active_time[1] = data[1].stod();},
                    &_ => {},
                }
            }
            
            
        }
    }

    pub fn read_interaction_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        let mut int_ct : usize = 0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            int_ct = 0;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.interaction_loop1(&mut headings, data_len, &mut int_ct);
            }
        }

        if int_ct == 0 {
            return;
        }

        self.interactions.int_vec = vec![Interaction::new(); int_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            int_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.interaction_loop2(&mut headings, &mut data, data_len, &mut int_ct);
            }
        }
    }

    pub fn ps_loop1(&self, headings : &Vec<CppStr>, data_len : usize, ps_ct : &mut usize) {
        if headings[0].s == "particleSources" {
            if headings[1].s == "elementSet" && data_len == 1 {
                *ps_ct += 1;
            }
        }
    }

    pub fn ps_loop2(&mut self, headings : &Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, ps_ct : &mut usize) {
        if headings[0].s == "particleSources" {
            if data_len == 1 {
                match headings[1].s.as_str() {
                    "elementSet" => {*ps_ct = increment_ct(*ps_ct);
                                     self.particle_sources[*ps_ct].element_set = data[0].clone();},
                    "swapSet" => self.particle_sources[*ps_ct].swap_set = data[0].clone(),
                    "swapDistance" => self.particle_sources[*ps_ct].swap_dist = data[0].stod(),
                    "swapRefLevel" => self.particle_sources[*ps_ct].refine_lev = data[0].stoi(),
                    "swapSpacing" => self.particle_sources[*ps_ct].spacing = data[0].stod(),
                    "randomVel" => self.particle_sources[*ps_ct].random_vel = data[0].stod(),
                    "velInLocal" => self.particle_sources[*ps_ct].vel_in_local = data[0].s.contains("yes"),
                    &_ => (),
                }
            }
            else if data_len == 2 {
                match headings[1].s.as_str() {
                    "temperature" => self.particle_sources[*ps_ct].temp.push_back(DualFloat {f1: data[0].stod(), f2: data[1].stod()}),
                    "frequency" => self.particle_sources[*ps_ct].frequency.push_back(DualFloat {f1: data[0].stod(), f2: data[1].stod()}),
                    "boundXRange" => {self.particle_sources[*ps_ct].x_range[0] = data[0].stod();
                                      self.particle_sources[*ps_ct].x_range[1] = data[1].stod();},
                    "boundYRange" => {self.particle_sources[*ps_ct].y_range[0] = data[0].stod();
                                      self.particle_sources[*ps_ct].y_range[1] = data[1].stod();},
                    "boundZRange" => {self.particle_sources[*ps_ct].z_range[0] = data[0].stod();
                                      self.particle_sources[*ps_ct].z_range[1] = data[1].stod();},
                    "activeTime" => {self.particle_sources[*ps_ct].active_time[0] = data[0].stod();
                                      self.particle_sources[*ps_ct].active_time[1] = data[1].stod();},
                    &_ => (),
                }
            }
            else if data_len == 3 {
                match headings[1].s.as_str() {
                    "refNodes" => {self.particle_sources[*ps_ct].ref_nodes[0] = data[0].clone();
                                   self.particle_sources[*ps_ct].ref_nodes[1] = data[1].clone();
                                   self.particle_sources[*ps_ct].ref_nodes[2] = data[2].clone();},
                    &_ => (),
                }
            }
            else if data_len == 4 {
                let mut new_pt = QuadFloat::new();
                new_pt.f1 = data[0].stod();
                new_pt.f2 = data[1].stod();
                new_pt.f3 = data[2].stod();
                new_pt.f4 = data[3].stod();
                match headings[1].s.as_str() {
                    "coordinates" => self.particle_sources[*ps_ct].coord.push_back(new_pt),
                    "meanVel" => self.particle_sources[*ps_ct].mean_vel.push_back(new_pt),
                    &_ => (),
                }
            }
        }
    }

    pub fn read_part_src_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        let mut ps_ct : usize = 0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ps_ct = 0;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.ps_loop1(&headings, data_len, &mut ps_ct);
            }
        }

        if ps_ct == 0 {
            return;
        }

        self.particle_sources = vec![ParticleSource::new(); ps_ct];

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ps_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                self.ps_loop2(&headings, &mut data, data_len, &mut ps_ct);
            }
        }
    }

    pub fn init_state_loop(&mut self, headings : &mut Vec<CppStr>, data : &mut Vec<CppStr>, data_len : usize, disp_hdings : &mut CppStr, fl_hdings : &mut CppStr) {
        let mut doub_inp : [f64; 10] = [ 0.0; 10];
        let mut i2 : usize;
        let mut i3 : usize;
        let mut seti : usize;
        let mut this_nd : &mut Node;

        if headings[0].s == "initialState" {
            if disp_hdings.s.contains(headings[1].s.as_str()) && data_len > 3 {
                seti = self.ns_map.at(&data[0].to_string());
                for ndi in self.node_sets[seti].labels.iter_mut() {
                    this_nd = &mut self.nodes[*ndi];
                    i2 = 1;
                    for i1 in 0..6 {
                        if i2 < data_len {
                            doub_inp[i1] = data[i2].stod();
                        }
                        else {
                            doub_inp[i1] = 0.0;
                        }
                        i2 += 1usize;
                    }
                    match headings[1].s.as_str() {
                        "displacement" => this_nd.set_initial_disp(&mut doub_inp),
                        "velocity" => this_nd.set_initial_vel(&mut doub_inp),
                        "acceleration" => this_nd.set_initial_acc(&mut doub_inp),
                        &_ => (),
                    }
                }
            }
            if data_len == 2 {
                match headings[1].s.as_str() {
                    "temperature" => {seti = self.ns_map.at(&data[0].to_string());
                                      for ndi in self.node_sets[seti].labels.iter_mut() {
                                          self.nodes[*ndi].initial_temp = data[1].stod();
                                      }},
                    "tdot" => {seti = self.ns_map.at(&data[0].to_string());
                               for ndi in self.node_sets[seti].labels.iter_mut() {
                                   self.nodes[*ndi].initial_tdot = data[1].stod();
                               }},
                    "concentration" => {seti = self.ns_map.at(&data[0].to_string());
                                        for ndi in self.node_sets[seti].labels.iter_mut() {
                                            self.nodes[*ndi].initial_fl_den = data[1].stod();
                                        }},
                    "cdot" => {seti = self.ns_map.at(&data[0].to_string());
                               for ndi in self.node_sets[seti].labels.iter_mut() {
                                   self.nodes[*ndi].initial_fl_den_dot = data[1].stod();
                               }},
                    &_ => (),
                }
            }
        }
    }

    pub fn read_initial_state(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut disp_hdings = CppStr::from(" displacement velocity acceleration");
        let mut fl_hdings = CppStr::from("flow flowdot");
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                self.init_state_loop(&mut headings, &mut data, data_len, &mut disp_hdings, &mut fl_hdings);
            }
        } else {
            panic!("Error: could not open initial state input file: {}", file_name.s);
        }
        
        return;
    }

    pub fn read_des_var_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut i1 : usize;
        let mut doub_inp : [f64; 10] = [0.0; 10];
        let mut int_inp : [usize; 10] = [0,0,0,0,0,0,0,0,0,0];
        
        let mut dv_ct : usize =  0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "designVariables" {
                    if headings[1].s == "category" && data_len == 1 {
                        dv_ct += 1usize;
                    }
                }
            }
        }
        
        self.design_vars = vec![DesignVariable::new(); dv_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            dv_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "designVariables" {
                    if data_len == 1 {
                        match headings[1].s.as_str() {
                            "category" => {dv_ct = increment_ct(dv_ct);
                                           self.design_vars[dv_ct].category = data[0].clone();},
                            "elementSet" => self.design_vars[dv_ct].el_set_name = data[0].clone(),
                            "nodeSet" => self.design_vars[dv_ct].nd_set_name = data[0].clone(),
                            "interaction" => self.design_vars[dv_ct].int_name = data[0].clone(),
                            "activeTime" => {doub_inp[0] = data[0].stod();
                                             doub_inp[1] = 1.0e+100;
                                             self.design_vars[dv_ct].set_active_time(&mut doub_inp);},
                            "component" => self.design_vars[dv_ct].component = data[0].stoi(),
                            "layer" => self.design_vars[dv_ct].layer = data[0].stoi(),
                            "coefficients" => self.design_vars[dv_ct].coefs.push_back(data[0].stod()),
                            &_ => (),
                        }
                    }
                    else if data_len == 2 {
                        match headings[1].s.as_str() {
                            "activeTime" => {doub_inp[0] = data[0].stod();
                                             doub_inp[1] = data[1].stod();
                                             self.design_vars[dv_ct].set_active_time(&mut doub_inp);},
                            "component" => {int_inp[0] = data[0].stoi() - 1;
                                            int_inp[1] = data[1].stoi() - 1;
                                            if int_inp[0] >= int_inp[1] {
                                                i1 = 6*int_inp[1] + int_inp[0];
                                            } else {
                                                i1 = 6*int_inp[0] + int_inp[1];
                                            }
                                            self.design_vars[dv_ct].component = i1;},
                            &_ => (),
                        }
                    }
                }
            }
        } else {
            panic!("Error: could not open design variable input file: {}",file_name.s);
        }
        
        return;
    }

    pub fn read_objective_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut doub_inp : [f64; 10] = [0.0; 10];
        
        let mut ob_ct : usize =  0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "objectiveTerms" {
                    if headings[1].s == "category" && data_len == 1 {
                        ob_ct += 1usize;
                    }
                }
            }
        }
        
        self.obj.terms = vec![ObjectiveTerm::new(); ob_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ob_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "objectiveTerms" {
                    if data_len == 1 {
                        match headings[1].s.as_str() {
                            "category" => {ob_ct = increment_ct(ob_ct);
                                           self.obj.terms[ob_ct].category = data[0].clone();},
                            "operator" => self.obj.terms[ob_ct].optr = data[0].clone(),
                            "activeTime" => {doub_inp[0] = data[0].stod();
                                             doub_inp[1] = 1.0e+100;
                                             self.obj.terms[ob_ct].set_active_time(&mut doub_inp);},
                            "component" => self.obj.terms[ob_ct].component = data[0].stoi(),
                            "layer" => self.obj.terms[ob_ct].layer = data[0].stoi(),
                            "coefficient" => self.obj.terms[ob_ct].coef = data[0].stod(),
                            "exponent" => self.obj.terms[ob_ct].expnt = data[0].stod(),
                            "elementSet" => self.obj.terms[ob_ct].el_set_name = data[0].clone(),
                            "nodeSet" => self.obj.terms[ob_ct].nd_set_name = data[0].clone(),
                            "targetValue" => {if data[0].is_doub() {
                                                  doub_inp[0] = data[0].stod();
                                                  self.obj.terms[ob_ct].tgt_vals.push_back(doub_inp[0]);
                                              }
                                              else {
                                                  self.obj.terms[ob_ct].tgt_tag = data[0].clone();
                                              }},
                            &_ => (),
                        }
                    }
                    else if data_len == 2 {
                        match headings[1].s.as_str() {
                            "activeTime" => {doub_inp[0] = data[0].stod();
                                             doub_inp[1] = data[1].stod();
                                             self.obj.terms[ob_ct].set_active_time(&mut doub_inp);},
                            &_ => (),
                        }
                    }
                }
            }
        } else {
            panic!("Error: could not open Objective input file: {}",file_name.s);
        }
        
        return;
    }

    pub fn read_des_var_values(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        let mut label : usize;
        let mut value : f64;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if data_len == 2 {
                    label = data[0].stoi();
                    value = data[1].stod();
                    self.design_vars[label].value.set_val(value);
                }
            }
        } else {
            panic!("Error: could not open design variable value input file: {}", file_name.s);
        }
        
        return;
    }

    pub fn read_node_results(&mut self, file_name : &mut CppStr) {
        let mut col_hd : Vec<&str> = vec![" "; 27];
        let mut ln_dat : Vec<&str>;
        let mut dat : f64;
        let mut ndi : usize;
        let mut ln_cpy : String;

        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                if line.contains("node") {
                    ln_cpy = line.clone();
                    col_hd = ln_cpy.split(',').collect();
                }
                else if line.contains(',') {
                    ln_dat = line.split(',').collect();
                    ndi = match ln_dat[0].parse::<usize>() {
                        Err(_why) => panic!("Error: problem reading the node results file, {}", file_name.s),
                        Ok(x) => x,
                    };
                    for i in 1..ln_dat.len() {
                        dat = match ln_dat[i].parse::<f64>() {
                            Err(_why) => panic!("Error: problem reading the node results file, {}", file_name.s),
                            Ok(x) => x,
                        };
                        match col_hd[i] {
                            "U1" => {self.nodes[ndi].displacement[0] = dat;},
                            "U2" => {self.nodes[ndi].displacement[1] = dat;},
                            "U3" => {self.nodes[ndi].displacement[2] = dat;},
                            "R1" => {self.nodes[ndi].displacement[3] = dat;},
                            "R2" => {self.nodes[ndi].displacement[4] = dat;},
                            "R3" => {self.nodes[ndi].displacement[5] = dat;},
                            "V1" => {self.nodes[ndi].velocity[0] = dat;},
                            "V2" => {self.nodes[ndi].velocity[1] = dat;},
                            "V3" => {self.nodes[ndi].velocity[2] = dat;},
                            "RV1" => {self.nodes[ndi].velocity[3] = dat;},
                            "RV2" => {self.nodes[ndi].velocity[4] = dat;},
                            "RV3" => {self.nodes[ndi].velocity[5] = dat;},
                            "A1" => {self.nodes[ndi].acceleration[0] = dat;},
                            "A2" => {self.nodes[ndi].acceleration[1] = dat;},
                            "A3" => {self.nodes[ndi].acceleration[2] = dat;},
                            "RA1" => {self.nodes[ndi].acceleration[3] = dat;},
                            "RA2" => {self.nodes[ndi].acceleration[4] = dat;},
                            "RA3" => {self.nodes[ndi].acceleration[5] = dat;},
                            "T" => {self.nodes[ndi].temperature = dat;},
                            "TDOT" => {self.nodes[ndi].temp_change_rate = dat;},
                            "C" => self.nodes[ndi].fl_den = dat,
                            "CDOT" => self.nodes[ndi].fl_den_dot = dat,
                            "DEN" => self.nodes[ndi].fl_den = dat,
                            "DENDOT" => self.nodes[ndi].fl_den_dot = dat,
                            &_ => {},
                        }
                    }
                }
            }
        }
    }

    pub fn read_time_step_soln(&mut self, t_step : usize) {
        let full_file = format!("{}{}{}{}",self.job[self.solve_cmd].file_name.s, "/solnTStep", t_step, ".out");
        let path = Path::new(full_file.as_str());

        let in_file = match File::open(&path) {
            Err(why) => panic!("couldn't open file, {}, {}", full_file, why),
            Ok(file) => file,
        };

        let mut reader = io::BufReader::new(in_file);

        let mut buf8 = [0u8; 8];
        let mut _b_read = 0usize;

        for nd in self.nodes.iter_mut() {
            if self.job[self.solve_cmd].thermal {
                _b_read = match reader.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_temp = f64::from_be_bytes(buf8);
                
                _b_read = match reader.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_tdot = f64::from_be_bytes(buf8);
            }

            if self.job[self.solve_cmd].diffusion {
                _b_read = match reader.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_fl_den = f64::from_be_bytes(buf8);
                
                _b_read = match reader.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_fl_den_dot = f64::from_be_bytes(buf8);
            }

            if self.job[self.solve_cmd].elastic {
                for i in 0..nd.num_dof {
                    _b_read = match reader.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_disp[i] = f64::from_be_bytes(buf8);
                }

                for i in 0..nd.num_dof {
                    _b_read = match reader.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_vel[i] = f64::from_be_bytes(buf8);
                }

                for i in 0..nd.num_dof {
                    _b_read = match reader.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_acc[i] = f64::from_be_bytes(buf8);
                }

            }
            
        }

        if self.job[self.solve_cmd].elastic {
            for el in self.elements.iter_mut() {
                if el.num_int_dof() > 0 {
                    for i in 0..el.num_int_dof() {
                        _b_read = match reader.read(&mut buf8) {
                            Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                            Ok(n) => n,
                        };
                        el.int_prev_disp[i] = f64::from_be_bytes(buf8);
                    }
                }
            }
        }


    }

}


