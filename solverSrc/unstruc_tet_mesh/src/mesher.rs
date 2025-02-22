use crate::mesh_node::*;
use crate::mesh_element::*;
use crate::mesh_face::*;
use crate::spatial_grid::*;
use crate::cpp_str::CppStr;


#[derive(Clone)]
pub struct Mesher {
    pub nodes : Vec<MeshNode>,
    pub nd_ct : usize,
    pub nd_cap : usize,
    pub node_grid : SpatialGrid,
    pub elements : Vec<MeshElement>,
    pub el_ct : usize,
    pub el_cap : usize,
    pub element_grid : SpatialGrid,
    pub faces : Vec<MeshFace>,
    pub fc_ct : usize,
    pub fc_cap : usize,
    pub face_grid : SpatialGrid,
    pub grid_out1 : Vec<usize>,
    pub grid_out2 : Vec<usize>,
    pub g_olen : usize,
    pub num_bound_nds : usize,
    pub avg_proj : f64,
    pub max_proj : f64,
    pub max_edge_len : f64,
    pub glob_proj_wt : f64,
    pub max_num_els : usize,
    pub new_el : usize,
    pub new_el_fcs : Vec<MeshFace>,
    pub new_nd : usize,
}

impl Mesher {
    pub fn new() -> Mesher {
        Mesher {
            nodes : Vec::new(),
            nd_ct : 0,
            nd_cap : 0,
            node_grid : SpatialGrid::new(),
            elements : Vec::new(),
            el_ct : 0,
            el_cap : 0,
            element_grid : SpatialGrid::new(),
            faces : Vec::new(),
            fc_ct : 0,
            fc_cap : 0,
            face_grid : SpatialGrid::new(),
            grid_out1 : Vec::new(),
            grid_out2 : Vec::new(),
            g_olen : 0,
            num_bound_nds : 0,
            avg_proj : 0.0,
            max_proj : 0.0,
            max_edge_len : 0.0,
            glob_proj_wt : 0.0,
            max_num_els : 0,
            new_el : 0,
            new_el_fcs : Vec::new(),
            new_nd : 0,
        }
    }
}


