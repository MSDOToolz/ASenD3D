use crate::constants::*;
use crate::list_ent::*;
use crate::diff_doub::*;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Node {
    pub label : usize,
    pub sorted_rank : usize,
    pub on_surf : bool,
    pub num_cells : usize,
    pub coord : [f64; 3],
    pub coord_dfd1 : [DiffDoub1; 3],
    pub displacement : [DiffDoub1; 3],
    pub prev_disp : [f64; 3],
    pub velocity : [DiffDoub1; 3],
    pub prev_vel : [f64; 3],
    
    pub fl_den : DiffDoub1,
    pub fl_den_dot : DiffDoub1,
    pub fl_vel : [DiffDoub1; 3],
    pub fl_vel_dot : [DiffDoub1; 3],
    pub temperature : DiffDoub1,
    pub temp_dot : DiffDoub1,
    pub turb_e : DiffDoub1,
    pub turb_e_dot : DiffDoub1,
    
    pub prev_fl_den : f64,
    pub prev_fl_den_dot : f64,
    pub prev_fl_vel : [f64; 3],
    pub prev_fl_vel_dot : [f64; 3],
    pub prev_temperature : f64,
    pub prev_temp_dot : f64,
    pub prev_turb_e : f64,
    pub prev_turb_e_dot : f64,
    
    pub lf_fl_den : f64,
    pub lf_fl_den_dot : f64,
    pub lf_fl_vel : [f64; 3],
    pub lf_fl_vel_dot : [f64; 3],
    pub lf_temperature : f64,
    pub lf_temp_dot : f64,
    pub lf_turb_e : f64,
    pub lf_turb_e_dot : f64,

    pub initial_fl_den : f64,
    pub initial_fl_den_dot : f64,
    pub initial_fl_vel : [f64; 3],
    pub initial_fl_vel_dot : [f64; 3],
    pub initial_temperature : f64,
    pub initial_temp_dot : f64,
    pub initial_turb_e : f64,
    pub initial_turb_e_dot : f64,
    
    pub d_var_lst : LinkedList<IDCapsule>,
    pub cell_lst : LinkedList<usize>,
    pub nd_lst : LinkedList<usize>,
}

impl Node {
    pub fn new() -> Node {
        Node {
            label : MAX_INT,
            sorted_rank : MAX_INT,
            on_surf : false,
            num_cells : MAX_INT,
            coord : [0f64; 3],
            coord_dfd1 : [DiffDoub1::new(); 3],
            displacement : [DiffDoub1::new(); 3],
            prev_disp : [0f64; 3],
            velocity : [DiffDoub1::new(); 3],
            prev_vel : [0f64; 3],
            
            fl_den : DiffDoub1::new(),
            fl_den_dot : DiffDoub1::new(),
            fl_vel : [DiffDoub1::new(); 3],
            fl_vel_dot : [DiffDoub1::new(); 3],
            temperature : DiffDoub1::new(),
            temp_dot : DiffDoub1::new(),
            turb_e : DiffDoub1::new(),
            turb_e_dot : DiffDoub1::new(),
            
            prev_fl_den : 0f64,
            prev_fl_den_dot : 0f64,
            prev_fl_vel : [0f64; 3],
            prev_fl_vel_dot : [0f64; 3],
            prev_temperature : 0f64,
            prev_temp_dot : 0f64,
            prev_turb_e : 0f64,
            prev_turb_e_dot : 0f64,
            
            lf_fl_den : 0f64,
            lf_fl_den_dot : 0f64,
            lf_fl_vel : [0f64; 3],
            lf_fl_vel_dot : [0f64; 3],
            lf_temperature : 0f64,
            lf_temp_dot : 0f64,
            lf_turb_e : 0f64,
            lf_turb_e_dot : 0f64,

            initial_fl_den : 0f64,
            initial_fl_den_dot : 0f64,
            initial_fl_vel : [0f64; 3],
            initial_fl_vel_dot : [0f64; 3],
            initial_temperature : 0f64,
            initial_temp_dot : 0f64,
            initial_turb_e : 0f64,
            initial_turb_e_dot : 0f64,
            
            d_var_lst : LinkedList::new(),
            cell_lst : LinkedList::new(),
            nd_lst : LinkedList::new(),
        }
    }
}

pub mod node_meth;