pub mod constants;
pub mod constraint;
pub mod design_var;
pub mod diff_doub;
pub mod element;
pub mod face;
pub mod fluid_cell;
pub mod fluid_node;
pub mod fluid_face;
pub mod interaction;
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
pub mod particle_source;
pub mod spatial_grid;
pub mod user;
pub mod cpp_str;
pub mod cpp_map;
pub mod fmath;
pub mod scratch;

use crate::model::Model;
use crate::cpp_str::CppStr;
use std::env;

fn main() {
    let args : Vec<String> = env::args().collect();
    
    // let args = vec!["".to_string(),"C:/Users/evaande/BladeReliability/modeling/triax_coupon/x_compression/job.yaml".to_string()];
    
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
