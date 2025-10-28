use crate::constants::MAX_INT;
use crate::list_ent::*;
use crate::cpp_str::CppStr;
use std::collections::LinkedList;

#[derive(Clone)]
pub struct ParticleSource {
    pub element_set : CppStr,
    pub elset_pt : usize,
    pub coord : LinkedList<QuadFloat>,
    pub ref_nodes : Vec<CppStr>,
    pub ref_nodes_i : [usize; 3],
    pub mean_vel : LinkedList<QuadFloat>,
    pub random_vel : f64,
    pub vel_in_local : bool,
    pub temp : LinkedList<DualFloat>,
    pub frequency : LinkedList<DualFloat>,
    pub x_range : [f64; 2],
    pub y_range : [f64; 2],
    pub z_range : [f64; 2],
    pub since_release : f64,
    pub active_time : [f64; 2],
}

impl ParticleSource {
    pub fn new() -> ParticleSource {
        ParticleSource {
            element_set : CppStr::new(),
            elset_pt : MAX_INT,
            coord : LinkedList::new(),
            ref_nodes : vec![CppStr::new(); 3],
            ref_nodes_i : [MAX_INT; 3],
            mean_vel : LinkedList::new(),
            random_vel : 0f64,
            vel_in_local : true,
            temp : LinkedList::new(),
            frequency : LinkedList::new(),
            x_range : [0f64; 2],
            y_range : [0f64; 2],
            z_range : [0f64; 2],
            since_release : 0f64,
            active_time : [0.0, 1.0e+100],
        }
    }
}