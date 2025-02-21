use crate::constants::MAX_INT;
use crate::cpp_str::CppStr;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct JobCommand {
    pub cmd_string : CppStr,
    pub file_name : CppStr,
    pub dynamic : bool,
    pub elastic : bool,
    pub load_ramp_steps : usize,
    pub newmark_beta : f64,
    pub newmark_gamma : f64,
    pub ray_damp_ck : f64,
    pub ray_damp_cm : f64,
    pub nonlinear_geom : bool,
    pub save_soln_hist : bool,
    pub sim_period : f64,
    pub lump_mass : bool,
    pub full_reform : usize,
    pub solver_bandwidth : usize,
    pub solver_block_dim : usize,
    pub solver_method : CppStr,
    pub max_it : usize,
    pub conv_tol : f64,
    pub static_load_time : LinkedList<f64>,
    pub thermal : bool,
    pub time_step : f64,
    pub this_type : CppStr,
    pub num_modes : usize,
    pub tgt_eval : f64,
    pub soln_field : CppStr,
    pub mode : usize,
    pub max_amplitude : f64,
    pub node_set : CppStr,
    pub fields : LinkedList<CppStr>,
    pub time_steps : LinkedList<usize>,
    pub time_step_tag : CppStr,
    pub element_set : CppStr,
    pub position : CppStr,
    pub write_modes : bool,
    pub properties : LinkedList<CppStr>,
    pub obj_include : LinkedList<CppStr>,
    pub write_gradient : bool,
}

impl JobCommand {
    pub fn new() -> JobCommand {
        JobCommand {
            // solve
            cmd_string : CppStr::new(),
            file_name : CppStr::new(),
            dynamic : false,
            elastic : true,
            load_ramp_steps : 1usize,
            newmark_beta : 0.25f64,
            newmark_gamma : 0.5f64,
            ray_damp_ck : 0f64,
            ray_damp_cm : 0f64,
            nonlinear_geom : false,
            save_soln_hist : true,
            sim_period : 1f64,
            lump_mass : false,
            full_reform : 1usize,
            solver_bandwidth : MAX_INT,
            solver_block_dim : MAX_INT,
            solver_method : CppStr::new(),
            max_it : 0usize,
            conv_tol : 1.0e-12f64,
            static_load_time : LinkedList::new(),
            thermal : false,
            time_step : 1f64,
            
            // modal
            this_type : CppStr::from("buckling"),
            num_modes : 10usize, 
            tgt_eval : 0f64,
            
            // set_soln_to_mode
            soln_field : CppStr::from("displacement"),
            mode : 1usize,
            max_amplitude : 1f64,

            // write_node_results
            node_set : CppStr::from("all"),
            fields : LinkedList::new(),
            time_steps : LinkedList::new(),
            time_step_tag : CppStr::from(""),

            // write_element_results
            element_set : CppStr::from("all"),
            position : CppStr::from("centroid"),

            // write_modal_results
            write_modes : true,

            // write_objective
            properties : LinkedList::new(),
            obj_include : LinkedList::new(),
            write_gradient : true,
        }
    }
}


