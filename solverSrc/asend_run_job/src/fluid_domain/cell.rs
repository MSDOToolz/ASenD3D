use crate::constants::MAX_INT;
use crate::diff_doub::*;
use crate::fluid_domain::design_var::*;
use crate::fluid_domain::sub_domain::*;
use crate::list_ent::*;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Cell {
    pub label : usize,
    pub nodes : [usize; 4],
    pub faces : [usize; 4],
    pub volume : DiffDoub1,
    pub dvars : LinkedList<IDCapsule>,
    pub loads : LinkedList<usize>,
    pub sub_dom_pt : usize,
}

impl Cell {
    pub fn new() -> Cell {
        Cell {
            label : MAX_INT,
            nodes : [MAX_INT; 4],
            faces : [MAX_INT; 4],
            volume : DiffDoub1::new(),
            dvars : LinkedList::new(),
            loads : LinkedList::new(),
            sub_dom_pt : MAX_INT,
        }
    }
}

pub mod cell_meth;

pub struct EqnPrereq {
    pub def_coord : [DiffDoub1; 12],

    pub fl_den : [DiffDoub1; 4],
    pub fl_den_dot : [DiffDoub1; 4],
    pub fl_vel : [DiffDoub1; 12],
    pub fl_vel_dot : [DiffDoub1; 12],
    pub temp : [DiffDoub1; 4],
    pub temp_dot : [DiffDoub1; 4],
    pub turb : [DiffDoub1; 4],
    pub turb_dot : [DiffDoub1; 4],
    pub vrel : [DiffDoub1; 12],

    pub viscosity : DiffDoub1,
    pub conductivity : DiffDoub1,
    pub expansion : DiffDoub1,
    pub spec_heat : DiffDoub1,
    pub ideal_gas : DiffDoub1,
    pub bulk_mod : DiffDoub1,
    pub ref_temp : DiffDoub1,
    pub ref_pres : DiffDoub1,
    pub ref_den : DiffDoub1,
    pub ref_enth : DiffDoub1,
    pub temp_vis_coef : DiffDoub1,
    pub turb_vis_coef : DiffDoub1,
    pub grad_turb_coef : DiffDoub1,
    pub diss_turb_coef : DiffDoub1,

    pub compressible : bool,
}