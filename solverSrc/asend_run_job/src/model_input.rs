use crate::model::*;
use crate::constants::*;
use crate::job::*;
use crate::node::*;
use crate::element::*;
use crate::nd_el_set::*;
use crate::section::*;
use crate::constraint::*;
use crate::load::*;
use crate::design_var::*;
use crate::objective::*;
use crate::cpp_str::CppStr;

use std::fs::File;
use std::io::{self, Read, BufRead};
use std::path::Path;
use std::collections::LinkedList;
//use std::fs::read_to_string;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

impl Model {
    pub fn read_input_line(& self, file_line : &mut CppStr, headings : &mut Vec<CppStr>, hd_ld_space : &mut [usize], data : &mut Vec<CppStr>, data_len : &mut usize) {
        let mut i1 : usize;
        let i2 : usize;
        let mut ln_len : usize;
        let wrd_len : usize;
        i1 = file_line.find("#");
        if i1 < MAX_INT {
            *file_line = file_line.substr(0,i1);
        }
        file_line.s = file_line.s.clone() + " ";
        ln_len = file_line.len();
        i1 = file_line.find(":");
        *data_len = 0;
        if i1 < MAX_INT {
            i2 = file_line.find_first_not_of(" -\n\t");
            wrd_len = i1 - i2;
            if headings[0].s == "" || hd_ld_space[0] == i2 {
                headings[0] = file_line.substr(i2,wrd_len);
                hd_ld_space[0] = i2;
                headings[1] = CppStr::from("");
                hd_ld_space[1] = 0;
                headings[2] = CppStr::from("");
                hd_ld_space[2] = 0;
                headings[3] = CppStr::from("");
                hd_ld_space[3] = 0;
            } else if headings[1].s == "" || hd_ld_space[1] == i2 {
                headings[1] = file_line.substr(i2,wrd_len);
                hd_ld_space[1] = i2;
                headings[2] = CppStr::from("");
                hd_ld_space[2] = 0;
                headings[3] = CppStr::from("");
                hd_ld_space[3] = 0;
            } else if headings[2].s == "" || hd_ld_space[2] == i2 {
                headings[2] = file_line.substr(i2,wrd_len);
                hd_ld_space[2] = i2;
                headings[3] = CppStr::from("");
                hd_ld_space[3] = 0;
            } else {
                headings[3] = file_line.substr(i2,wrd_len);
                hd_ld_space[3] = i2;
            }
            i1 += 1usize;
            while i1 < ln_len {
                *file_line = file_line.substr(i1, MAX_INT);
                i1 = file_line.find_first_not_of(" ,[]\t\n");
                if i1 < MAX_INT {
                    *file_line = file_line.substr(i1, MAX_INT);
                    ln_len = file_line.len();
                    i1 = file_line.find_first_of(" ,[]\t\n");
                    if i1 < MAX_INT {
                        data[*data_len] = file_line.substr(0,i1);
                        *data_len += 1usize;
                    } else {
                        i1 = ln_len;
                    }
                } else {
                    i1 = ln_len;
                }
            }
        } else {
            i1 = file_line.find("- ");
            if i1 < MAX_INT {
                i1 += 1usize;
                while i1 < ln_len {
                    *file_line = file_line.substr(i1, MAX_INT);
                    i1 = file_line.find_first_not_of(" ,[]\t\n");
                    if i1 < MAX_INT {
                        *file_line = file_line.substr(i1, MAX_INT);
                        ln_len = file_line.len();
                        i1 = file_line.find_first_of(" ,[]\t\n");
                        if i1 < MAX_INT {
                            data[*data_len] = file_line.substr(0,i1);
                            *data_len += 1usize;
                        } else {
                            i1 = ln_len;
                        }
                    } else {
                        i1 = ln_len;
                    }
                }
            }
        }
        
        return;
    }

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
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[1].s == "command" && data_len == 1 {
                    cmd_ct += 1usize;
                }
            }
        }
        
        self.job = vec![JobCommand::new(); cmd_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            cmd_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[1].s == "command" && data_len == 1 {
                    if cmd_ct == MAX_INT {
                        cmd_ct = 0;
                    }
                    else {
                        cmd_ct += 1usize;
                    }
                    self.job[cmd_ct].cmd_string = data[0].clone();
                } else if headings[1].s == "fileName" && data_len == 1 {
                    self.job[cmd_ct].file_name = data[0].clone();
                } else if headings[1].s == "dynamic" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].dynamic = true;
                    } else {
                        self.job[cmd_ct].dynamic = false;
                    }
                }
                else if headings[1].s == "elastic" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].elastic = true;
                    } else {
                        self.job[cmd_ct].elastic = false;
                    }
                } else if headings[1].s == "loadRampSteps" && data_len == 1 {
                    self.job[cmd_ct].load_ramp_steps = CppStr::stoi(&mut data[0]);
                } else if headings[1].s == "newmarkBeta" && data_len == 1 {
                    self.job[cmd_ct].newmark_beta = CppStr::stod(&mut data[0]);
                } else if headings[1].s == "newmarkGamma" && data_len == 1 {
                    self.job[cmd_ct].newmark_gamma = CppStr::stod(&mut data[0]);
                } else if headings[1].s == "nonlinearGeom" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].nonlinear_geom = true;
                    } else {
                        self.job[cmd_ct].nonlinear_geom = false;
                    }
                } else if headings[1].s == "saveSolnHist" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].save_soln_hist = true;
                    } else {
                        self.job[cmd_ct].save_soln_hist = false;
                    }
                }
                else if headings[1].s == "solnHistDir" && data_len == 1 {
                    self.job[cmd_ct].file_name = data[0].clone();
                }
                else if headings[1].s == "lumpMass" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].lump_mass = true;
                    }
                    else {
                        self.job[cmd_ct].lump_mass = false;
                    }
                }
                else if headings[1].s == "fullReformFreq" {
                    self.job[cmd_ct].full_reform = CppStr::stoi(&mut data[0]);
                }
                else if headings[1].s == "simPeriod" && data_len == 1 {
                    self.job[cmd_ct].sim_period = CppStr::stod(&mut data[0]);
                } else if headings[1].s == "solverBandwidth" && data_len == 1 {
                    self.job[cmd_ct].solver_bandwidth = CppStr::stoi(&mut data[0]);
                } else if headings[1].s == "solverBlockDim" && data_len == 1 {
                    self.job[cmd_ct].solver_block_dim = CppStr::stoi(&mut data[0]);
                } else if headings[1].s == "solverMethod" && data_len == 1 {
                    self.job[cmd_ct].solver_method = data[0].clone();
                }
                else if headings[1].s == "maxIterations" && data_len == 1 {
                    self.job[cmd_ct].max_it = CppStr::stoi(&mut data[0]);
                }
                else if headings[1].s == "convergenceTol" && data_len == 1 {
                    self.job[cmd_ct].conv_tol = CppStr::stod(&mut data[0]);
                }
                else if headings[1].s == "staticLoadTime" && data_len == 1 {
                    self.job[cmd_ct].static_load_time.push_back(CppStr::stod(&mut data[0]));
                    //new_cmd->static_load_time = stod(data[0]);
                } else if headings[1].s == "thermal" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].thermal = true;
                    } else {
                        self.job[cmd_ct].thermal = false;
                    }
                } else if headings[1].s == "timeStep" && data_len == 1 {
                    self.job[cmd_ct].time_step = CppStr::stod(&mut data[0]);
                } else if headings[1].s == "type" && data_len == 1 {
                    self.job[cmd_ct].this_type = data[0].clone();
                } else if headings[1].s == "numModes" && data_len == 1 {
                    self.job[cmd_ct].num_modes = CppStr::stoi(&mut data[0]);
                }
                else if headings[1].s == "targetEigenvalue" && data_len == 1 {
                    self.job[cmd_ct].tgt_eval = CppStr::stod(&mut data[0]);
                }
                else if headings[1].s == "solnField" && data_len == 1 {
                    self.job[cmd_ct].soln_field = data[0].clone();
                } else if headings[1].s == "mode" && data_len == 1 {
                    self.job[cmd_ct].mode = CppStr::stoi(&mut data[0]);
                } else if headings[1].s == "maxAmplitude" && data_len == 1 {
                    self.job[cmd_ct].max_amplitude = CppStr::stod(&mut data[0]);
                } else if headings[1].s == "nodeSet" && data_len == 1 {
                    self.job[cmd_ct].node_set = data[0].clone();
                } else if headings[1].s == "fields" && data_len == 1 {
                    self.job[cmd_ct].fields.push_back(data[0].clone());
                } else if headings[1].s == "timeSteps" {
                    if data_len == 1 {
                        if data[0].s == "all" {
                            self.job[cmd_ct].time_step_tag = CppStr::from("all");
                        } else if data[0].s == "last" {
                            self.job[cmd_ct].time_step_tag = CppStr::from("last");
                        } else {
                            if CppStr::is_int(&mut data[0]) {
                                i1 = CppStr::stoi(&mut data[0]);
                                self.job[cmd_ct].time_steps.push_back(i1);
                            }
                            else {
                                println!("Warning: possible invalid entry for {} timeSteps in job file {}", headings[0].s, file_name.s);
                            }
                        }
                    } else if data_len > 1 {
                        i2 = CppStr::stoi(&mut data[0]);
                        i3 = CppStr::stoi(&mut data[1]);
                        if data_len == 3 {
                            i4 = CppStr::stoi(&mut data[2]);
                        } else {
                            i4 = 1;
                        }
                        i1 = i2;
                        while i1 < i3 {
                            self.job[cmd_ct].time_steps.push_back(i1);
                            i1 += i4;
                        }
                    }
                } else if headings[1].s == "elementSet" && data_len == 1 {
                    self.job[cmd_ct].element_set = data[0].clone();
                }
                else if headings[1].s == "position" && data_len == 1 {
                    self.job[cmd_ct].position = data[0].clone();
                }
                else if headings[1].s == "writeModes" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].write_modes = true;
                    } else {
                        self.job[cmd_ct].write_modes = false;
                    }
                } else if headings[1].s == "properties" && data_len == 1 {
                    self.job[cmd_ct].properties.push_back(data[0].clone());
                } else if headings[1].s == "include" && data_len == 1 {
                    self.job[cmd_ct].obj_include.push_back(data[0].clone());
                } else if headings[1].s == "writeGradient" && data_len == 1 {
                    if data[0].s == "yes" {
                        self.job[cmd_ct].write_gradient = true;
                    } else {
                        self.job[cmd_ct].write_gradient = false;
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
        let mut fl_ct : usize =  0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "nodes" && data_len == 4 {
                    nd_ct += 1usize;
                }
                else if headings[0].s == "elements" {
                    if headings[1].s == "connectivity" && data_len > 1 {
                        el_ct += 1usize;
                    }
                }
                else if headings[0].s == "sets" {
                    if headings[1].s == "node" {
                        if headings[2].s == "name" && data_len == 1 {
                            ns_ct += 1usize;
                        }
                    }
                    else if headings[1].s == "element" {
                        if headings[2].s == "name" && data_len == 1 {
                            es_ct += 1usize;
                        }
                    }
                }
                else if headings[0].s == "sections" {
                    if headings[1].s == "type" && data_len == 1 {
                        sec_ct += 1usize;
                    }
                }
                else if headings[0].s == "materials" {
                    if headings[1].s == "name" && data_len == 1 {
                        mat_ct += 1usize;
                    }
                }
                else if headings[0].s == "fluids" {
                    if headings[1].s == "name" && data_len == 1 {
                        fl_ct += 1usize;
                    }
                }
            }
        }
        
        self.nodes = vec![Node::new(); nd_ct];
        self.elements = vec![Element::new(); el_ct];
        self.node_sets = vec![Set::new(); ns_ct + nd_ct + 1];
        self.element_sets = vec![Set::new(); es_ct + el_ct + 1];
        self.sections = vec![Section::new(); sec_ct];
        self.materials = vec![Material::new(); mat_ct];
        self.fluids = vec![Fluid::new(); fl_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ns_ct = MAX_INT;
            es_ct = MAX_INT;
            sec_ct = MAX_INT;
            mat_ct = MAX_INT;
            fl_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "nodes" && data_len == 4 {
                    nd_ct = CppStr::stoi(&mut data[0]);
                    self.nodes[nd_ct].label = nd_ct;
                    doub_inp[0] = CppStr::stod(&mut data[1]);
                    doub_inp[1] = CppStr::stod(&mut data[2]);
                    doub_inp[2] = CppStr::stod(&mut data[3]);
                    self.nodes[nd_ct].set_crd(&mut doub_inp);
                } else if headings[0].s == "elements" {
                    if headings[1].s == "type" && data_len == 1 {
                        if data[0].s == "tet4" {
                            el_type = 4;
                        } else if data[0].s == "wedge6" {
                            el_type = 6;
                        } else if data[0].s == "brick8" {
                            el_type = 8;
                        } else if data[0].s == "brickIM" {
                            el_type = 81;
                        }
                        else if data[0].s == "tet10" {
                            el_type = 10;
                        }
                        else if data[0].s == "shell3" {
                            el_type = 3;
                        } else if data[0].s == "shell4" {
                            el_type = 41;
                        } else if data[0].s == "beam2" {
                            el_type = 2;
                        }
                        else if data[0].s == "frcFld" {
                            el_type = 21;
                        }
                        else if data[0].s == "mass" {
                            el_type = 1;
                        }
                        else if data[0].s == "fl4" {
                            el_type = 400;
                        }
                        else if data[0].s == "fl6" {
                            el_type = 600;
                        }
                        else if data[0].s == "fl8" {
                            el_type = 800;
                        }
                        else if data[0].s == "fl10" {
                            el_type = 1000;
                        }
                        else {
                            panic!("Error: unrecognized element type: {}", data[0].s);
                        }
                    } else if headings[1].s == "connectivity" && data_len > 1 {
                        el_ct = CppStr::stoi(&mut data[0]);
                        let new_el = &mut self.elements[el_ct];
                        new_el.initialize_type(el_type);
                        new_el.label = el_ct;
                        for i1 in 1..data_len {
                            int_inp[i1-1] = CppStr::stoi(&mut data[i1]);
                        }
                        new_el.set_nodes(&mut int_inp);
                    }
                } else if headings[0].s == "sets" {
                    if headings[1].s == "node" {
                        if headings[2].s == "name" && data_len == 1 {
                            if ns_ct == MAX_INT {
                                ns_ct = 0;
                            }
                            else {
                                ns_ct += 1usize;
                            }
                            let new_set : &mut Set = &mut self.node_sets[ns_ct];
                            new_set.name = data[0].clone();
                        } else if headings[2].s == "labels" {
                            if data_len == 1 {
                                self.node_sets[ns_ct].labels.push_back(CppStr::stoi(&mut data[0]));
                            } else if data_len > 1 {
                                i2 = CppStr::stoi(&mut data[0]);
                                i3 = CppStr::stoi(&mut data[1]);
                                if data_len == 3 {
                                    i4 = CppStr::stoi(&mut data[2]);
                                } else {
                                    i4 = 1;
                                }
                                i1 = i2;
                                while i1 < i3 {
                                    self.node_sets[ns_ct].labels.push_back(i1);
                                    i1 += i4;
                                }
                            }
                        }
                    } else if headings[1].s == "element" {
                        if headings[2].s == "name" && data_len == 1 {
                            if es_ct == MAX_INT {
                                es_ct = 0;
                            }
                            else {
                                es_ct += 1usize;
                            }
                            let new_set = &mut self.element_sets[es_ct];
                            new_set.name = data[0].clone();
                        } else if headings[2].s == "labels" {
                            if data_len == 1 {
                                self.element_sets[es_ct].labels.push_back(CppStr::stoi(&mut data[0]));
                            } else if data_len > 1 {
                                i2 = CppStr::stoi(&mut data[0]);
                                i3 = CppStr::stoi(&mut data[1]);
                                if data_len == 3 {
                                    i4 = CppStr::stoi(&mut data[2]);
                                } else {
                                    i4 = 1;
                                }
                                i1 = i2;
                                while i1 < i3 {
                                    self.element_sets[es_ct].labels.push_back(i1);
                                    i1 += i4;
                                }
                            }
                        }
                    }
                } else if headings[0].s == "sections" {
                    if headings[1].s == "type" && data_len == 1 {
                        if sec_ct == MAX_INT {
                            sec_ct = 0;
                        }
                        else {
                            sec_ct += 1usize;
                        }
                        self.sections[sec_ct].this_type = data[0].clone();
                    } else if headings[1].s == "material" && data_len == 1 {
                        self.sections[sec_ct].mat_name = data[0].clone();
                    }
                    else if headings[1].s == "fluid" && data_len == 1 {
                        self.sections[sec_ct].fl_name = data[0].clone();
                    }
                    else if headings[1].s == "orientation" && data_len == 6 {
                        for i1 in 0..6 {
                            doub_inp[i1] = CppStr::stod(&mut data[i1]);
                        }
                        self.sections[sec_ct].set_orientation(&mut doub_inp);
                    } else if headings[1].s == "layup" {
                        if headings[2].s == "zOffset" && data_len == 1 {
                            self.sections[sec_ct].z_offset = CppStr::stod(&mut data[0]);
                        } else if headings[2].s == "layers" {
                            if headings[3].s == "material" && data_len == 1 {
                                let mut new_lay = Layer::new();
                                new_lay.mat_name = data[0].clone();
                                self.sections[sec_ct].layers.push_back(new_lay);
                            } else if headings[3].s == "thickness" && data_len == 1 {
                                lay_pt = match self.sections[sec_ct].layers.back_mut() {
                                    None => panic!("could not access back of section layer list in input reader"),
                                    Some(x) => x,
                                };
                                lay_pt.thickness = CppStr::stod(&mut data[0]);
                            } else if headings[3].s == "angle" && data_len == 1 {
                                lay_pt = match self.sections[sec_ct].layers.back_mut() {
                                    None => panic!("could not access back of section layer list in input reader"),
                                    Some(x) => x,
                                };
                                lay_pt.angle = CppStr::stod(&mut data[0]);
                            }
                        }
                    } else if headings[1].s == "beamProps" {
                        if headings[2].s == "area" && data_len == 1 {
                            self.sections[sec_ct].area = CppStr::stod(&mut data[0]);
                        } else if headings[2].s == "I" && data_len == 5 {
                            for i1 in 0..5 {
                                doub_inp[i1] = CppStr::stod(&mut data[i1]);
                            }
                            self.sections[sec_ct].set_area_moment(&mut doub_inp);
                        } else if headings[2].s == "J" && data_len == 1 {
                            self.sections[sec_ct].polar_moment = CppStr::stod(&mut data[1]);
                        } else if headings[2].s == "stiffness" && data_len == 3 {
                            int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                            int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                            doub_inp[0] = CppStr::stod(&mut data[2]);
                            self.sections[sec_ct].set_stiffness(int_inp[0],  int_inp[1],  doub_inp[0]);
                        } else if headings[2].s == "mass" && data_len == 3 {
                            int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                            int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                            doub_inp[0] = CppStr::stod(&mut data[2]);
                            self.sections[sec_ct].set_mass(int_inp[0], int_inp[1], doub_inp[0]);
                        }
                        else if headings[2].s == "damping" && data_len == 3 {
                            int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                            int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                            doub_inp[0] = CppStr::stod(&mut data[2]);
                            self.sections[sec_ct].set_damping(int_inp[0],   int_inp[1],   doub_inp[0]);
                        }
                        else if headings[2].s == "expLoadCoef" && data_len == 6 {
                            for i1 in 0..6 {
                                doub_inp[i1] = CppStr::stod(&mut data[i1]);
                            }
                            self.sections[sec_ct].set_exp_ld(&mut doub_inp);
                        } else if headings[2].s == "conductivity" && data_len == 1 {
                            self.sections[sec_ct].conductivity = CppStr::stod(&mut data[0]);
                        } else if headings[2].s == "specHeat" && data_len == 1 {
                            self.sections[sec_ct].spec_heat = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "potField" {
                        if headings[2].s == "coef" && data_len == 1 {
                            self.sections[sec_ct].pot_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "exp" && data_len == 1 {
                            self.sections[sec_ct].pot_exp = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "dampField" {
                        if headings[2].s == "coef" && data_len == 1 {
                            self.sections[sec_ct].damp_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "exp" && data_len == 1 {
                            self.sections[sec_ct].damp_exp = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "thermField" {
                        if headings[2].s == "condCoef" && data_len == 1 {
                            self.sections[sec_ct].cond_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "radCoef" && data_len == 1 {
                            self.sections[sec_ct].rad_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "refTemp" && data_len == 1 {
                            self.sections[sec_ct].ref_temp = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "massPerEl" && data_len == 1 {
                        self.sections[sec_ct].mass_per_el = CppStr::stod(&mut data[0]);
                    }
                    else if headings[1].s == "specHeat" && data_len == 1 {
                        self.sections[sec_ct].spec_heat = CppStr::stod(&mut data[0]);
                    }
                    else if headings[1].s == "fluidParams" {
                        if headings[2].s == "denVisCoef" && data_len == 1 {
                            self.sections[sec_ct].den_vis_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "tempVisCoef" && data_len == 1 {
                            self.sections[sec_ct].temp_vis_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "turbVisCoef" && data_len == 1 {
                            self.sections[sec_ct].turb_vis_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "gradVTurbCoef" && data_len == 1 {
                            self.sections[sec_ct].grad_vturb_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "dissTurbCoef" && data_len == 1 {
                            self.sections[sec_ct].diss_turb_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "enthCoef" && data_len == 1 {
                            self.sections[sec_ct].enth_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "enthExp" && data_len == 1 {
                            self.sections[sec_ct].enth_exp = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "presCoef" && data_len == 1 {
                            self.sections[sec_ct].pres_coef = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "presExp" && data_len == 1 {
                            self.sections[sec_ct].pres_exp = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "refTemp" && data_len == 1 {
                            self.sections[sec_ct].ref_temp = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "refDen" && data_len == 1 {
                            self.sections[sec_ct].ref_den = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "refTurbE" && data_len == 1 {
                            self.sections[sec_ct].ref_turb_e = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "refEnth" && data_len == 1 {
                            self.sections[sec_ct].ref_enth = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "elementSet" && data_len == 1 {
                        self.sections[sec_ct].el_set_name = data[0].clone();
                    }
                } else if headings[0].s == "materials" {
                    if headings[1].s == "name" && data_len == 1 {
                        if mat_ct == MAX_INT {
                            mat_ct = 0;
                        }
                        else {
                            mat_ct += 1usize;
                        }
                        self.materials[mat_ct].name = data[0].clone();
                    } else if headings[1].s == "density" && data_len == 1 {
                        self.materials[mat_ct].density = CppStr::stod(&mut data[0]);
                    } else if headings[1].s == "elastic" {
                        if headings[2].s == "E" && data_len > 0 {
                            let this_mat = &mut self.materials[mat_ct];
                            this_mat.modulus[0] = CppStr::stod(&mut data[0]);
                            if data_len == 3 {
                                this_mat.modulus[1] = CppStr::stod(&mut data[1]);
                                this_mat.modulus[2] = CppStr::stod(&mut data[2]);
                            } else {
                                this_mat.modulus[1] = this_mat.modulus[0];
                                this_mat.modulus[2] = this_mat.modulus[0];
                            }
                        } else if headings[2].s == "nu" && data_len > 0 {
                            let this_mat = &mut self.materials[mat_ct];
                            this_mat.poisson_ratio[0] = CppStr::stod(&mut data[0]);
                            if data_len == 3 {
                                this_mat.poisson_ratio[1] = CppStr::stod(&mut data[1]);
                                this_mat.poisson_ratio[2] = CppStr::stod(&mut data[2]);
                            } else {
                                this_mat.poisson_ratio[1] = this_mat.poisson_ratio[0];
                                this_mat.poisson_ratio[2] = this_mat.poisson_ratio[0];
                            }
                        } else if headings[2].s == "G" && data_len > 0 {
                            let this_mat = &mut self.materials[mat_ct];
                            this_mat.shear_mod[0] = CppStr::stod(&mut data[0]);
                            if data_len == 3 {
                                this_mat.shear_mod[1] = CppStr::stod(&mut data[1]);
                                this_mat.shear_mod[2] = CppStr::stod(&mut data[2]);
                            } else {
                                this_mat.shear_mod[1] = this_mat.shear_mod[0];
                                this_mat.shear_mod[2] = this_mat.shear_mod[0];
                            }
                        } else if headings[2].s == "stiffness" && data_len == 3 {
                            int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                            int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                            doub_inp[0] = CppStr::stod(&mut data[2]);
                            self.materials[mat_ct].set_stiffness(int_inp[0],  int_inp[1],  doub_inp[0]);
                        }
                    } else if headings[1].s == "thermal" {
                        if headings[2].s == "conductivity" && data_len == 6 {
                            let this_mat = &mut self.materials[mat_ct];
                            for i1 in 0..6 {
                                this_mat.conductivity[i1] = CppStr::stod(&mut data[i1]);
                            }
                        } else if headings[2].s == "expansion" && data_len == 6 {
                            let this_mat = &mut self.materials[mat_ct];
                            for i1 in 0..6 {
                                this_mat.expansion[i1] = CppStr::stod(&mut data[i1]);
                            }
                        } else if headings[2].s == "specHeat" && data_len == 1 {
                            self.materials[mat_ct].spec_heat = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "damping" && data_len == 3 {
                        int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                        int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                        doub_inp[0] = CppStr::stod(&mut data[2]);
                        self.materials[mat_ct].set_damping(int_inp[0],   int_inp[1],   doub_inp[0]);
                    }
                    else if headings[1].s == "failure" {
                        if headings[2].s == "maxStress" {
                            if headings[3].s == "tensile" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_ten_stress[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_ten_stress[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_ten_stress[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_ten_stress[1] = n_mat.max_ten_stress[0];
                                    n_mat.max_ten_stress[2] = n_mat.max_ten_stress[0];
                                }
                            } else if headings[3].s == "compressive" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_comp_stress[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_comp_stress[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_comp_stress[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_comp_stress[1] = n_mat.max_comp_stress[0];
                                    n_mat.max_comp_stress[2] = n_mat.max_comp_stress[0];
                                }
                            } else if headings[3].s == "shear" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_shear_stress[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_shear_stress[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_shear_stress[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_shear_stress[1] = n_mat.max_shear_stress[0];
                                    n_mat.max_shear_stress[2] = n_mat.max_shear_stress[0];
                                }
                            }
                        } else if headings[2].s == "maxStrain" {
                            if headings[3].s == "tensile" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_ten_strain[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_ten_strain[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_ten_strain[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_ten_strain[1] = n_mat.max_ten_strain[0];
                                    n_mat.max_ten_strain[2] = n_mat.max_ten_strain[0];
                                }
                            } else if headings[3].s == "compressive" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_comp_strain[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_comp_strain[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_comp_strain[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_comp_strain[1] = n_mat.max_comp_strain[0];
                                    n_mat.max_comp_strain[2] = n_mat.max_comp_strain[0];
                                }
                            } else if headings[3].s == "shear" && data_len > 0 {
                                let n_mat = &mut self.materials[mat_ct];
                                n_mat.max_shear_strain[0] = CppStr::stod(&mut data[0]);
                                if data_len == 3 {
                                    n_mat.max_shear_strain[1] = CppStr::stod(&mut data[1]);
                                    n_mat.max_shear_strain[2] = CppStr::stod(&mut data[2]);
                                } else {
                                    n_mat.max_shear_strain[1] = n_mat.max_shear_strain[0];
                                    n_mat.max_shear_strain[2] = n_mat.max_shear_strain[0];
                                }
                            }
                        } else if headings[2].s == "maxStrainEnergy" && data_len == 1 {
                            self.materials[mat_ct].max_str_eng = CppStr::stod(&mut data[0]);
                        } else if headings[2].s == "maxMises" && data_len == 1 {
                            self.materials[mat_ct].max_mises = CppStr::stod(&mut data[0]);
                        }
                    }
                }
                else if headings[0].s == "fluids" {
                    if headings[1].s == "name" && data_len == 1 {
                        if fl_ct == MAX_INT {
                            fl_ct = 0;
                        }
                        else {
                            fl_ct += 1usize;
                        }
                        self.fluids[fl_ct].name = data[0].clone();
                    }
                    else if headings[1].s == "viscosity" && data_len == 1 {
                        self.fluids[fl_ct].viscosity = CppStr::stod(&mut data[0]);
                    }
                    else if headings[1].s == "thermal" {
                        if headings[2].s == "conductivity" && data_len == 1 {
                            self.fluids[fl_ct].therm_cond = CppStr::stod(&mut data[0]);
                        }
                        else if headings[2].s == "specHeat" && data_len == 1 {
                            self.fluids[fl_ct].spec_heat = CppStr::stod(&mut data[0]);
                        }
                    }
                    else if headings[1].s == "idealGasConst" && data_len == 1 {
                        self.fluids[fl_ct].ideal_gas = CppStr::stod(&mut data[0]);
                    }
                }
            }
        } else {
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

    pub fn read_constraint_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let all_types = CppStr::from("displacement temperature");
        
        let mut ec_ct : usize =  0;
        let mut tc_ct : usize =  0;
        let mut new_con = Constraint::new();
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "constraints" {
                    if headings[1].s == "type" && data_len == 1 {
                        if data[0].s == "displacement" {
                            ec_ct += 1usize;
                        }
                        else if data[0].s == "temperature" {
                            tc_ct += 1usize;
                        }
                    }
                }
            }
        }

        if ec_ct == 0 && tc_ct == 0 {
            return;
        }
        
        self.elastic_const.const_vec = vec![Constraint::new(); ec_ct];
        self.thermal_const.const_vec = vec![Constraint::new(); tc_ct];
        
        let mut tm_pt : &mut ConstraintTerm;
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            ec_ct = MAX_INT;
            tc_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "constraints" {
                    if headings[1].s == "type" && data_len == 1 {
                        if new_con.this_type.s == "displacement" {
                            self.elastic_const.const_vec[ec_ct] = new_con.clone();
                        }
                        else if new_con.this_type.s == "temperature" {
                            self.thermal_const.const_vec[tc_ct] = new_con.clone();
                        }
                        if data[0].s == "displacement" {
                            if ec_ct == MAX_INT {
                                ec_ct = 0;
                            }
                            else {
                                ec_ct += 1usize;
                            }
                            new_con = Constraint::new();
                            new_con.this_type = CppStr::from("displacement");
                        }
                        else if data[0].s == "temperature" {
                            if tc_ct == MAX_INT {
                                tc_ct = 0;
                            }
                            else {
                                tc_ct += 1usize;
                            }
                            new_con = Constraint::new();
                            new_con.this_type = CppStr::from("temperature");
                        }
                        else {
                            panic!("Error: {} is not a valid constraint type. Allowable values are {}", data[0].s, all_types.s);
                        }
                    }
                    else if headings[1].s == "terms" {
                        if headings[2].s == "nodeSet" && data_len == 1 {
                            let mut new_cn = ConstraintTerm::new();
                            new_cn.node_set = data[0].clone();
                            new_con.terms.push_back(new_cn);
                        } else if headings[2].s == "dof" && data_len == 1 {
                            tm_pt = match new_con.terms.back_mut() {
                                None => panic!("failed to access back of constraint terms list"),
                                Some(x) => x,
                            };
                            tm_pt.dof = CppStr::stoi(&mut data[0]);
                        } else if headings[2].s == "coef" && data_len == 1 {
                            tm_pt = match new_con.terms.back_mut() {
                                None => panic!("failed to acces back of constraint terms list"),
                                Some(x) => x,
                            };
                            tm_pt.coef = CppStr::stod(&mut data[0]);
                        }
                    } else if headings[1].s == "rhs" && data_len == 1 {
                        new_con.rhs = CppStr::stod(&mut data[0]);
                    }
                }
            }
        } else {
            panic!("Error: could not open Constraint input file: {}",file_name.s);
        }

        if new_con.this_type.s == "displacement" {
            self.elastic_const.const_vec[ec_ct] = new_con.clone();
        }
        else if new_con.this_type.s == "temperature" {
            self.thermal_const.const_vec[tc_ct] = new_con.clone();
        }
        
        return;
    }

    pub fn read_load_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut i1 : usize;
        let mut doub_inp : [f64; 10] = [0.0; 10 ];
        let mut elastic_list = CppStr::from("nodalForce bodyForce gravitational centrifugal surfacePressure surfaceTraction");
        let mut thermal_list = CppStr::from("bodyHeadGen surfaceFlux");
        let mut new_ld = Load::new();
        
        let mut e_ld_ct : usize =  0;
        let mut t_ld_ct : usize =  0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "loads" {
                    if headings[1].s == "type" && data_len == 1 {
                        i1 = elastic_list.find(data[0].s.as_str());
                        if i1 < MAX_INT {
                            e_ld_ct += 1usize;
                        }
                        i1 = thermal_list.find(data[0].s.as_str());
                        if i1 < MAX_INT {
                            t_ld_ct += 1usize;
                        }
                    }
                }
            }
        }

        if e_ld_ct == 0 && t_ld_ct == 0 {
            return;
        }
        
        self.elastic_loads = vec![Load::new(); e_ld_ct];
        self.thermal_loads = vec![Load::new(); t_ld_ct];
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            e_ld_ct = MAX_INT;
            t_ld_ct = MAX_INT;
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "loads" {
                    if headings[1].s == "type" && data_len == 1 {
                        i1 = elastic_list.find(new_ld.this_type.s.as_str());
                        if i1 < MAX_INT {
                            self.elastic_loads[e_ld_ct] = new_ld.clone();
                        }
                        i1 = thermal_list.find(new_ld.this_type.s.as_str());
                        if i1 < MAX_INT {
                            self.thermal_loads[t_ld_ct] = new_ld.clone();
                        }
                        i1 = elastic_list.find(data[0].s.as_str());
                        if i1 < MAX_INT {
                            if e_ld_ct == MAX_INT {
                                e_ld_ct = 0;
                            }
                            else {
                                e_ld_ct += 1usize;
                            }
                            //self.elastic_loads[e_ld_ct].this_type = data[0].clone();
                            //new_ld = &mut self.elastic_loads[e_ld_ct];
                            new_ld = Load::new();
                            new_ld.this_type = data[0].clone();
                        }
                        i1 = thermal_list.find(data[0].s.as_str());
                        if i1 < MAX_INT {
                            if t_ld_ct == MAX_INT {
                                t_ld_ct = 0;
                            }
                            else {
                                t_ld_ct += 1usize;
                            }
                            //self.thermal_loads[t_ld_ct].this_type = data[0].clone();
                            //new_ld = &mut self.thermal_loads[t_ld_ct];
                            new_ld = Load::new();
                            new_ld.this_type = data[0].clone();
                        }
                    } else if headings[1].s == "activeTime" && data_len > 0 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        if data_len == 2 {
                            doub_inp[1] = CppStr::stod(&mut data[1]);
                        } else {
                            doub_inp[1] = 1.0e+100;
                        }
                        new_ld.set_act_time(&mut doub_inp);
                    } else if headings[1].s == "load" && data_len > 0 {
                        for i1 in 0..6 {
                            if i1 < data_len {
                                doub_inp[i1] = CppStr::stod(&mut data[i1]);
                            } else {
                                doub_inp[i1] = 0.0;
                            }
                        }
                        new_ld.set_load(&mut doub_inp);
                    } else if headings[1].s == "nodeSet" && data_len == 1 {
                        new_ld.node_set = data[0].clone();
                    } else if headings[1].s == "elementSet" && data_len == 1 {
                        new_ld.element_set = data[0].clone();
                    } else if headings[1].s == "normDir" && data_len == 3 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        doub_inp[1] = CppStr::stod(&mut data[1]);
                        doub_inp[2] = CppStr::stod(&mut data[2]);
                        new_ld.set_norm_dir(&mut doub_inp);
                    } else if headings[1].s == "normTolerance" && data_len == 1 {
                        new_ld.norm_tol = CppStr::stod(&mut data[0]);
                    } else if headings[1].s == "center" && data_len == 3 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        doub_inp[1] = CppStr::stod(&mut data[1]);
                        doub_inp[2] = CppStr::stod(&mut data[2]);
                        new_ld.set_center(&mut doub_inp);
                    } else if headings[1].s == "axis" && data_len == 3 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        doub_inp[1] = CppStr::stod(&mut data[1]);
                        doub_inp[2] = CppStr::stod(&mut data[2]);
                        new_ld.set_axis(&mut doub_inp);
                    } else if headings[1].s == "angularVelocity" {
                        new_ld.angular_vel = CppStr::stod(&mut data[0]);
                    }
                }
            }
        } else {
            panic!("Error: could not open Load input file: {}", file_name.s);
        }

        i1 = elastic_list.find(new_ld.this_type.s.as_str());
        if i1 < MAX_INT {
            self.elastic_loads[e_ld_ct] = new_ld.clone();
        }
        i1 = thermal_list.find(new_ld.this_type.s.as_str());
        if i1 < MAX_INT {
            self.thermal_loads[t_ld_ct] = new_ld.clone();
        }
        
        return;
    }

    pub fn read_initial_state(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0,0,0,0];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        
        let mut i2 : usize;
        let mut i3 : usize;
        let mut seti : usize;
        let mut doub_inp : [f64; 10] = [ 0.0; 10];
        let mut disp_hdings = CppStr::from(" displacement velocity acceleration");
        let mut this_nd : &mut Node;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "initialState" {
                    i3 = disp_hdings.find(headings[1].s.as_str());
                    if i3 < MAX_INT && data_len > 3 {
                        seti = self.ns_map.at(&data[0].to_string());
                        for ndi in self.node_sets[seti].labels.iter_mut() {
                            this_nd = &mut self.nodes[*ndi];
                            i2 = 1;
                            for i1 in 0..6 {
                                if i2 < data_len {
                                    doub_inp[i1] = CppStr::stod(&mut data[i2]);
                                }
                                else {
                                    doub_inp[i1] = 0.0;
                                }
                                i2 += 1usize;
                            }
                            if headings[1].s == "displacement" {
                                this_nd.set_initial_disp(&mut doub_inp);
                            }
                            else if headings[1].s == "velocity" {
                                this_nd.set_initial_vel(&mut doub_inp);
                            }
                            else if headings[1].s == "acceleration" {
                                this_nd.set_initial_acc(&mut doub_inp);
                            }
                        }
                    } else if headings[1].s == "temperature" && data_len == 2 {
                        seti = self.ns_map.at(&data[0].to_string());
                        for ndi in self.node_sets[seti].labels.iter_mut() {
                            self.nodes[*ndi].initial_temp = CppStr::stod(&mut data[1]);
                        }
                    } else if headings[1].s == "tdot" && data_len == 2 {
                        seti = self.ns_map.at(&data[0].to_string());
                        for ndi in self.node_sets[seti].labels.iter_mut() {
                            self.nodes[*ndi].initial_tdot = CppStr::stod(&mut data[1]);
                        }
                    }
                }
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
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
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
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "designVariables" {
                    if headings[1].s == "category" && data_len == 1 {
                        if dv_ct == MAX_INT {
                            dv_ct = 0;
                        }
                        else {
                            dv_ct += 1usize;
                        }
                        self.design_vars[dv_ct].category = data[0].clone();
                    } else if headings[1].s == "elementSet" && data_len == 1 {
                        self.design_vars[dv_ct].el_set_name = data[0].clone();
                    } else if headings[1].s == "nodeSet" && data_len == 1 {
                        self.design_vars[dv_ct].nd_set_name = data[0].clone();
                    } else if headings[1].s == "activeTime" && data_len > 0 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        if data_len == 2 {
                            doub_inp[1] = CppStr::stod(&mut data[1]);
                        } else {
                            doub_inp[1] = 1.0e+100;
                        }
                        self.design_vars[dv_ct].set_active_time(&mut doub_inp);
                    } else if headings[1].s == "component" {
                        if data_len == 1 {
                            self.design_vars[dv_ct].component = CppStr::stoi(&mut data[0]);
                        } else if data_len == 2 {
                            int_inp[0] = CppStr::stoi(&mut data[0]) - 1;
                            int_inp[1] = CppStr::stoi(&mut data[1]) - 1;
                            if int_inp[0] >= int_inp[1] {
                                i1 = 6*int_inp[1] + int_inp[0];
                            } else {
                                i1 = 6*int_inp[0] + int_inp[1];
                            }
                            self.design_vars[dv_ct].component = i1;
                        }
                    } else if headings[1].s == "layer" && data_len == 1 {
                        self.design_vars[dv_ct].layer = CppStr::stoi(&mut data[0]);
                    } else if headings[1].s == "coefficients" && data_len == 1 {
                        self.design_vars[dv_ct].coefs.push_back(CppStr::stod(&mut data[0]));
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
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
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
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if headings[0].s == "objectiveTerms" {
                    if headings[1].s == "category" && data_len == 1 {
                        if ob_ct == MAX_INT {
                            ob_ct = 0;
                        }
                        else {
                            ob_ct += 1usize;
                        }
                        self.obj.terms[ob_ct].category = data[0].clone();
                    } else if headings[1].s == "operator" && data_len == 1 {
                        self.obj.terms[ob_ct].optr = data[0].clone();
                    } else if headings[1].s == "activeTime" && data_len > 0 {
                        doub_inp[0] = CppStr::stod(&mut data[0]);
                        if data_len == 2 {
                            doub_inp[1] = CppStr::stod(&mut data[1]);
                        } else {
                            doub_inp[1] = 1.0e+100;
                        }
                        self.obj.terms[ob_ct].set_active_time(&mut doub_inp);
                    } else if headings[1].s == "component" && data_len == 1 {
                        self.obj.terms[ob_ct].component = CppStr::stoi(&mut data[0]);
                    } else if headings[1].s == "layer" && data_len == 1 {
                        self.obj.terms[ob_ct].layer = CppStr::stoi(&mut data[0]);
                    } else if headings[1].s == "coefficient" && data_len == 1 {
                        self.obj.terms[ob_ct].coef = CppStr::stod(&mut data[0]);
                    } else if headings[1].s == "exponent" && data_len == 1 {
                        self.obj.terms[ob_ct].expnt = CppStr::stod(&mut data[0]);
                    } else if headings[1].s == "elementSet" && data_len == 1 {
                        self.obj.terms[ob_ct].el_set_name = data[0].clone();
                    } else if headings[1].s == "nodeSet" && data_len == 1 {
                        self.obj.terms[ob_ct].nd_set_name = data[0].clone();
                    } else if headings[1].s == "targetValue" && data_len == 1 {
                        if CppStr::is_doub(&mut data[0]) {
                            doub_inp[0] = CppStr::stod(&mut data[0]);
                            self.obj.terms[ob_ct].tgt_vals.push_back(doub_inp[0]);
                        }
                        else {
                            self.obj.terms[ob_ct].tgt_tag = data[0].clone();
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
                self.read_input_line(&mut file_line, &mut headings, &mut hd_ld_space, &mut data, &mut data_len);
                if data_len == 2 {
                    label = CppStr::stoi(&mut data[0]);
                    value = CppStr::stod(&mut data[1]);
                    self.design_vars[label].value.set_val(value);
                }
            }
        } else {
            panic!("Error: could not open design variable value input file: {}", file_name.s);
        }
        
        return;
    }

    pub fn read_node_results(&mut self, file_name : &mut CppStr) {
        let mut i2 : usize;
        let mut nd : usize;
        let mut nd_dat : [f64; 6] = [0f64; 6];
        
        let mut file_line = CppStr::new();
        let mut headings = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [ 0,0,0,0 ];
        let mut data = vec![CppStr::new(); 11];
        let mut data_len : usize = 0usize;
        let mut disp_fields = CppStr::from("displacement velocity acceleration");
        let mut thrm_fields = CppStr::from("temperature tdot");
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                self.read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "nodeResults" && data_len > 1 {
                    nd = CppStr::stoi(&mut data[0]);
                    let this_nd = &mut self.nodes[nd];
                    i2 = disp_fields.find(headings[1].s.as_str());
                    if i2 < MAX_INT {
                        for i1 in 0..this_nd.num_dof {
                            nd_dat[i1] = CppStr::stod(&mut data[i1 + 1]);
                        }
                        if headings[1].s == "displacement" {
                            this_nd.set_displacement(&mut nd_dat);
                        }
                        else if headings[1].s == "velocity" {
                            this_nd.set_velocity(&mut nd_dat);
                        }
                        else {
                            this_nd.set_acceleration(&mut nd_dat);
                        }
                    }
                    i2 = thrm_fields.find(headings[1].s.as_str());
                    if i2 < MAX_INT {
                        nd_dat[0] = CppStr::stod(&mut data[1]);
                        if headings[1].s == "temperature" {
                            this_nd.temperature = nd_dat[0];
                        }
                        else {
                            this_nd.temp_change_rate = nd_dat[0];
                        }
                    }
                }
            }
        }
        else {
            panic!("Error: could not open Node result input file: {}",file_name.s);
        }
        
        return;
    }

    pub fn read_time_step_soln(&mut self, t_step : usize) {
        let full_file = format!("{}{}{}{}",self.job[self.solve_cmd].file_name.s, "/solnTStep", t_step, ".out");
        let path = Path::new(full_file.as_str());

        let mut in_file = match File::open(&path) {
            Err(why) => panic!("couldn't open file, {}, {}", full_file, why),
            Ok(file) => file,
        };

        let mut buf8 = [0u8; 8];
        let mut _b_read = 0usize;

        for nd in self.nodes.iter_mut() {
            if self.job[self.solve_cmd].thermal {
                _b_read = match in_file.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_temp = f64::from_be_bytes(buf8);
                _b_read = match in_file.read(&mut buf8) {
                    Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                    Ok(n) => n,
                };
                nd.prev_tdot = f64::from_be_bytes(buf8);
            }
            if self.job[self.solve_cmd].elastic {
                for i in 0..nd.num_dof {
                    _b_read = match in_file.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_disp[i] = f64::from_be_bytes(buf8);
                }

                for i in 0..nd.num_dof {
                    _b_read = match in_file.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_vel[i] = f64::from_be_bytes(buf8);
                }

                for i in 0..nd.num_dof {
                    _b_read = match in_file.read(&mut buf8) {
                        Err(why) => panic!("problem reading file, {}, {}", full_file, why),
                        Ok(n) => n,
                    };
                    nd.prev_acc[i] = f64::from_be_bytes(buf8);
                }

            }
        }

        if self.job[self.solve_cmd].elastic {
            for el in self.elements.iter_mut() {
                if el.num_int_dof > 0 {
                    for i in 0..el.num_int_dof {
                        _b_read = match in_file.read(&mut buf8) {
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


