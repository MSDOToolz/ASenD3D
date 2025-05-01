pub mod constants;
pub mod constraint;
pub mod design_var;
pub mod diff_doub;
pub mod element;
pub mod face;
pub mod job;
pub mod list_ent;
pub mod load;
pub mod lower_tri_mat;
pub mod lu_mat;
pub mod matrix_functions;
pub mod model;
pub mod node;
pub mod objective;
pub mod section;
pub mod nd_el_set;
pub mod const_meth;
pub mod d_var_meth;
pub mod diff_doub_meth;
pub mod element_meth;
pub mod element_equations;
pub mod element_properties;
pub mod element_soln_fields;
pub mod face_meth;
pub mod job_meth;
pub mod list_ent_meth;
pub mod spatial_grid;
pub mod grid_meth;
pub mod load_meth;
pub mod lower_tri_mat_meth;
pub mod lu_mat_meth;
pub mod model_analysis;
pub mod model_meth;
pub mod model_input;
pub mod model_output;
pub mod model_res_utils;
pub mod user;
pub mod node_meth;
pub mod objective_meth;
pub mod section_meth;
pub mod nd_el_set_meth;
pub mod cpp_str;
pub mod cpp_map;
pub mod fmath;
pub mod scratch;

use crate::model::Model;
use crate::cpp_str::CppStr;
use std::env;

fn main() {
    let args : Vec<String> = env::args().collect();
    
    // let args = vec!["".to_string(),"C:/Users/evaande/ASenDHome/ASenD3D/examples/testCases/singleHex/staticElastic/staticElasticJob.yaml".to_string()];
    
    let job_file = match args.get(1) {
        None => panic!("Error: job input file not specified in call to analysis solver."),
        Some(x) => x,
    };
    
    println!("Beginning Job: {}",job_file);
    let mut new_model = Model::new();
    let mut cpp_jf = CppStr{s : job_file.clone()};
    new_model.read_job(&mut cpp_jf);
    new_model.execute_job();
}
