pub mod constants;
pub mod cpp_str;
pub mod mesh_element;
pub mod mesher;
pub mod mesh_face;
pub mod mesh_node;
pub mod spatial_grid;
pub mod utilities;
pub mod element_meth;
pub mod mesher_meth;
pub mod face_meth;
pub mod node_meth;
pub mod grid_meth;
pub mod fmath;

use crate::cpp_str::*;
use crate::mesher::*;
use std::env;

fn main() {
    let args : Vec<String> = env::args().collect();
    let n_arg = args.len();
    let mut input_file : CppStr;
    let mut output_file : CppStr;
    
    if (n_arg > 2) {
        input_file = CppStr::from(args[1].as_str());
        output_file = CppStr::from(args[2].as_str());
    }
    else if(n_arg > 1) {
        input_file = CppStr::from(args[1].as_str());
        output_file = CppStr::from("mesh_out.yaml");
    }
    else {
        panic!("Error: no input file provided for unstructured tet mesh generator.");
    }

    let mut mesher = Mesher::new();
    mesher.read_input(&mut input_file);
    mesher.prep();
    println!("generating mesh");
    let resolved = mesher.generate_mesh();
    println!("finished mesh");
    if(resolved) {
        println!("distributing nodes");
        mesher.distribute_nodes();
        mesher.generate_mesh();
    }
    else {
        println!("Warning: mesh reached maximum number of elements unresolved");
    }

    mesher.write_output(&mut output_file);

}
