use crate::mesh_node::*;
use crate::cpp_str::CppStr;


#[derive(Clone)]
pub struct MeshElement {
    pub nodes : [usize; 4],
}

impl MeshElement {
    pub fn new() -> MeshElement {
        MeshElement {
            nodes : [0usize; 4],
        }
    }
}
