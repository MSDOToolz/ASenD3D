use crate::model::*;
use crate::constants::*;
use crate::list_ent::*;
use crate::node::*;
use crate::element::*;
use crate::nd_el_set::*;
use crate::section::*;
use crate::constraint::*;
use crate::design_var::*;
use crate::objective::*;
use crate::matrix_functions::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


impl Model {
    
    pub fn execute_job(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;
        let mut num_tsteps : usize;
        let mut time : f64;
        let mut cmd_str = CppStr::new();
        let mut file_name = CppStr::new();
        let mut exten = CppStr::new();
        let mut full_fname = CppStr::new();
        
        //for this_cmd in self.job.iter_mut() {
        for ci in 0..self.job.len() {
            cmd_str = self.job[ci].cmd_string.clone();
            file_name = self.job[ci].file_name.clone();
            if(cmd_str.s == "readModelInput") {
                println!("{}{}", "reading Model input: ", file_name.s);
                self.read_model_input(&mut file_name);
                self.read_constraint_input(&mut file_name);
                self.read_load_input(&mut file_name);
                self.read_initial_state(&mut file_name);
            } else if(cmd_str.s == "readConstraints") {
                println!("{}{}", "reading constraints: ", file_name.s);
                self.read_constraint_input(&mut file_name);
            } else if(cmd_str.s == "readLoads") {
                println!("{}{}", "reading loads: " , file_name.s);
                self.read_load_input(&mut file_name);
            } else if(cmd_str.s == "readInitialState") {
                println!("{}{}", "reading initial state: " , file_name.s);
                self.read_initial_state(&mut file_name);
            } else if(cmd_str.s == "readDesignVarInput") {
                println!("{}{}", "reading design variable input: " , file_name.s);
                self.read_des_var_input(&mut file_name);
            } else if(cmd_str.s == "readObjectiveInput") {
                println!("{}{}", "reading Objective input: " , file_name.s);
                self.read_objective_input(&mut file_name);
            }
            else if (cmd_str.s == "readDesignVarValues") {
                self.read_des_var_values(&mut file_name);
            }
            else if (cmd_str.s == "readNodeResults") {
                self.read_node_results(&mut file_name);
            }
            else if (cmd_str.s == "solvePrep") {
                self.solve_cmd = ci;
                self.analysis_prep();
            }
            else if (cmd_str.s == "solve") {
                println!("{}", "running main analysis " );
                self.solve_cmd = ci;
                self.solve();
            }
            else if (cmd_str.s == "zeroSolution") {
                let mut fld_cln = self.job[ci].fields.clone();
                self.zero_solution(&mut fld_cln);
            }
            else if (cmd_str.s == "modalAnalysis") {
                println!("{}{}", "running modal analysis " , file_name.s);
                self.modal_cmd = ci;
                self.eigen_solve();
            }
            else if (cmd_str.s == "setSolnToMode") {
                let mut sol_fld_cln = self.job[ci].soln_field.clone();
                let mut md = self.job[ci].mode;
                let mut mx_amp = self.job[ci].max_amplitude;
                self.set_soln_to_mode(&mut sol_fld_cln, md,  mx_amp);
            }
            else if (cmd_str.s == "calcObjective") {
                println!("{}{}", "calculating Objective function " , file_name.s);
                self.get_objective();
            }
            else if (cmd_str.s == "calcObjGradient") {
                println!("{}{}", "calculating Objective gradient " , file_name.s);
                self.get_obj_gradient();
            } else if (cmd_str.s == "writeNodeResults" || cmd_str.s == "writeElementResults") {
                if (cmd_str.s == "writeNodeResults") {
                    println!("writing Node results");
                }
                else {
                    println!("writing Element results");
                }
                num_tsteps = self.job[ci].time_steps.len();
                let mut ns_cln = self.job[ci].node_set.clone();
                let mut es_cln = self.job[ci].element_set.clone();
                let mut fld_cln = self.job[ci].fields.clone();
                let mut pos_cln = self.job[ci].position.clone();
                let mut ts_cln = self.job[ci].time_steps.clone();
                if (num_tsteps > 0) {
                    i2 = self.job[ci].file_name.find(".");
                    if (i2 < MAX_INT) {
                        exten = self.job[ci].file_name.substr(i2, MAX_INT);
                        file_name = self.job[ci].file_name.substr(0, i2);
                    }
                    else {
                        exten = CppStr::from("");
                        file_name = self.job[ci].file_name.clone();
                    }
                    for ts in ts_cln.iter() {
                        full_fname.s = format!("{}_timestep{}{}", file_name.s, ts, exten.s);
                        if (cmd_str.s == "writeNodeResults") {
                            self.write_node_results(&mut full_fname, &mut ns_cln, &mut fld_cln, *ts);
                        }
                        else {
                            self.write_element_results(&mut full_fname, &mut es_cln, &mut fld_cln, &mut pos_cln, *ts);
                        }
                    }
                }
                else {
                    file_name = self.job[ci].file_name.clone();
                    if (cmd_str.s == "writeNodeResults") {
                        self.write_node_results(&mut file_name, &mut ns_cln, &mut fld_cln,  MAX_INT);
                    }
                    else {
                        self.write_element_results(&mut file_name, &mut es_cln, &mut fld_cln, &mut pos_cln,  MAX_INT);
                    }
                }
            }
            else if (cmd_str.s == "writeModalResults") {
                println!("writing modal results");
                file_name = self.job[ci].file_name.clone();
                let mut wm = self.job[ci].write_modes;
                self.write_modal_results(&mut file_name, wm);
            }
            else if (cmd_str.s == "writeObjective") {
                println!("writing Objective results");
                file_name = self.job[ci].file_name.clone();
                let mut ob_in = self.job[ci].obj_include.clone();
                let mut wg = self.job[ci].write_gradient;
                self.write_objective(&mut file_name, &mut ob_in, wg);
            }
            
        }
        
        return;
    }

    pub fn print_contents(&self) {
        let mut ct = 0usize;
        println!("nodes:");
        for nd in self.nodes.iter() {
            println!("{}, {}, {}, {}", ct, nd.coord[0], nd.coord[1], nd.coord[2]);
            ct += 1;
        }
        println!("elements:");
        ct = 0;
        for el in self.elements.iter() {
            print!("{}",ct);
            for i in 0..el.num_nds {
                print!(", {}",el.nodes[i]);
            }
            println!("");
            ct += 1;
        }
        println!("node sets:");
        for ns in self.node_sets.iter() {
            println!("name: {}", ns.name.s);
            for nd in ns.labels.iter() {
                println!("{},",*nd);
            }
        }
        println!("element sets:");
        for es in self.element_sets.iter() {
            println!("name: {}", es.name.s);
            for el in es.labels.iter() {
                println!("{},",*el);
            }
        }
        println!("sections:");
        for sec in self.sections.iter() {
            println!("element set: {}", sec.el_set_name.s);
            println!("material: {}", sec.mat_name.s);
        }
        println!("materials:");
        for mat in self.materials.iter() {
            println!("name: {}", mat.name.s);
            println!("modulus: {}, {}, {}", mat.modulus[0], mat.modulus[1], mat.modulus[2]);
            println!("poisson: {}, {}, {}", mat.poisson_ratio[0], mat.poisson_ratio[1], mat.poisson_ratio[2]);
            println!("shear modulus: {}, {}, {}", mat.shear_mod[0], mat.shear_mod[1], mat.shear_mod[2]);
        }
        println!("constraints:");
        for cn in self.elastic_const.const_vec.iter() {
            println!("terms:");
            for tm in cn.terms.iter() {
                println!("node_set: {}", tm.node_set.s);
                println!("dof: {}", tm.dof);
                println!("coef: {}", tm.coef);
            }
        }
        println!("loads:");
        for ld in self.elastic_loads.iter() {
            println!("type: {}", ld.this_type.s);
            println!("node_set: {}", ld.node_set.s);
            print!("load: ");
            for i in 0..6 {
                print!("{}, ", ld.load[i]);
            }
            println!("");
        }
    }

    pub fn print_vec(in_vec : &Vec<f64>) {
        println!("Vector contents:");
        for f in in_vec.iter() {
            println!("{},",*f);
        }
    }

}


