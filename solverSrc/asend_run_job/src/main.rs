pub mod constants;
pub mod diff_doub;
pub mod file_util;
pub mod list_ent;
pub mod lower_tri_mat;
pub mod matrix_functions;
pub mod model;
pub mod nd_el_set;
pub mod spatial_grid;
pub mod cpp_str;
pub mod cpp_map;
pub mod fmath;

use crate::model::Model;
use crate::cpp_str::CppStr;
use std::env;

fn main() {
    let args : Vec<String> = env::args().collect();
    
    //let args = vec!["".to_string(),"C:/Users/evaande/ASenDHome/ASenD3D/examples/testCases/shellBeam/transverseTipLoading/job.yaml".to_string()];
    //let args = vec!["".to_string(),"C:/Users/evaande/ASenDHome/single_shell_debug/staticElastic/staticElasticJob.yaml".to_string()];
    //let args = vec!["".to_string(),"C:/Users/evans/ASenDHome/ASenD3D/examples/modalAnalysis/bouncingBall/job.yaml".to_string()];

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
