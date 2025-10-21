use crate::constants::MAX_INT;
use crate::list_ent::*;
use crate::cpp_str::CppStr;
use std::collections::LinkedList;

#[derive(Clone)]
pub struct ParticleSource {
    pub element_set : CppStr,
    pub elset_pt : usize,
    pub coord : LinkedList<QuadFloat>,
    pub ref_node : CppStr,
    pub ref_node_i : usize,
    pub mean_vel : LinkedList<QuadFloat>,
    pub random_vel : f64,
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
            ref_node : CppStr::new(),
            ref_node_i : MAX_INT,
            mean_vel : LinkedList::new(),
            random_vel : 0f64,
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