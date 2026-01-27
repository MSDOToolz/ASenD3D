use crate::constants::MAX_INT;
use crate::list_ent::*;
use crate::diff_doub::*;
use crate::spatial_grid::SpatialGrid;
use crate::cpp_str::CppStr;
use std::collections::LinkedList;

#[derive(Clone)]
pub struct Interaction {
    pub name : CppStr,
    pub node_set1 : CppStr,
    pub set_pt1 : usize,
    pub node_set2 : CppStr,
    pub set_pt2 : usize,
    pub pot_coef : LinkedList<DualFloat>,
    pub pot_exp : f64,
    pub damp_coef : LinkedList<DualFloat>,
    pub damp_dist_exp : f64,
    pub damp_vel_exp : f64,
    pub mag_coef : LinkedList<DualFloat>,
    pub mag_dist_exp : f64,
    pub mag_vel_exp : f64,
    pub cond_coef : f64,
    pub rad_coef : f64,
    pub ref_temp : f64,
    pub max_dist : f64,
    pub max_nbrs : usize,
    pub max_ratio : f64,
    pub ideal_gas : f64,
    pub bulk_mod : f64,
    pub therm_exp : f64,
    pub ref_den : f64,
    pub ref_pres : f64,
    pub active_time : [f64; 2],
    pub dvars : LinkedList<IDCapsule>,
}

impl Interaction {
    pub fn new() -> Interaction {
        Interaction {
            name : CppStr::new(),
            node_set1 : CppStr::new(),
            set_pt1 : MAX_INT,
            node_set2 : CppStr::new(),
            set_pt2 : MAX_INT,
            pot_coef : LinkedList::new(),
            pot_exp : 1.0,
            damp_coef : LinkedList::new(),
            damp_dist_exp : 1.0,
            damp_vel_exp : 1.0,
            mag_coef : LinkedList::new(),
            mag_dist_exp : 2.0,
            mag_vel_exp : 1.0,
            cond_coef : 0.0,
            rad_coef : 0.0,
            ref_temp : 0.0,
            max_dist : 0.0,
            max_nbrs : MAX_INT,
            max_ratio : 1.0e+100,
            ideal_gas : -1.0,
            bulk_mod : -1.0,
            therm_exp : 0.0,
            ref_den : 1.0,
            ref_pres : 0.0,
            active_time : [0.0, 1.0e+100],
            dvars : LinkedList::new(),
        }
    }
}

#[derive(Clone)]
pub struct InteractionList {
    pub int_vec : Vec<Interaction>,
    pub nd_interaction : Vec<bool>,
    pub nd_in_set : Vec<bool>,
    pub nd_active : Vec<bool>,
    pub nd_mass_dfd0 : Vec<DiffDoub0>,
    pub nd_mass_dfd1 : Vec<DiffDoub1>,
    pub interact_grid : SpatialGrid,
    pub grid_out : Vec<usize>,
    pub nearest : Vec<usize>,
    pub near_dist : Vec<f64>,
}

impl InteractionList {
    pub fn new() -> InteractionList {
        InteractionList {
            int_vec : Vec::new(),
            nd_interaction : Vec::new(),
            nd_in_set : Vec::new(),
            nd_active : Vec::new(),
            nd_mass_dfd0 : Vec::new(),
            nd_mass_dfd1 : Vec::new(),
            interact_grid : SpatialGrid::new(),
            grid_out : Vec::new(),
            nearest : Vec::new(),
            near_dist : Vec::new(),
        }
    }
}

pub mod interaction_meth;