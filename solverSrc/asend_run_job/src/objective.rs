use crate::constants::*;
use crate::list_ent::*;
use crate::cpp_str::CppStr;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct ObjectiveTerm {
    pub category : CppStr,
    pub optr : CppStr,
    pub active_time : [f64; 2],
    pub component : usize,
    pub layer : usize,
    pub coef : f64,
    pub expnt : f64,
    pub el_set_name : CppStr,
    pub el_set_ptr : usize,
    pub nd_set_name : CppStr,
    pub nd_set_ptr : usize,
    pub tgt_tag : CppStr,
    pub tgt_vals : LinkedList<f64>,
    pub value : f64,
    pub q_vec : Vec<f64>,
    pub el_vol_vec : Vec<f64>,
    pub tgt_vec : Vec<f64>,
    pub err_norm_vec : Vec<f64>,
    pub q_len : usize,
    pub d_qd_u : SparseMat,
    pub d_qd_v : SparseMat,
    pub d_qd_a : SparseMat,
    pub d_qd_t : SparseMat,
    pub d_qd_tdot : SparseMat,
    pub d_qd_d : SparseMat,
    pub d_vd_d : SparseMat,
}

impl ObjectiveTerm {
    pub fn new() -> ObjectiveTerm {
        ObjectiveTerm {
            category : CppStr::new(),
            optr : CppStr::new(),
            active_time : [0f64, 1.0e+100f64],
            component : 0usize,
            layer : 0usize,
            coef : 0f64,
            expnt : 0f64,
            el_set_name : CppStr::new(),
            el_set_ptr : MAX_INT,
            nd_set_name : CppStr::new(),
            nd_set_ptr : MAX_INT,
            tgt_tag : CppStr::new(),
            tgt_vals : LinkedList::new(),
            value : 0f64,
            q_vec : Vec::new(),
            el_vol_vec : Vec::new(),
            tgt_vec : Vec::new(),
            err_norm_vec : Vec::new(),
            q_len : 0usize,
            d_qd_u : SparseMat::new(),
            d_qd_v : SparseMat::new(),
            d_qd_a : SparseMat::new(),
            d_qd_t : SparseMat::new(),
            d_qd_tdot : SparseMat::new(),
            d_qd_d : SparseMat::new(),
            d_vd_d : SparseMat::new(),
        }
    }
}

#[derive(Clone)]
pub struct Objective {
    pub terms : Vec<ObjectiveTerm>,
}

impl Objective {
    pub fn new() -> Objective {
        Objective {
            terms : Vec::new(),
        }
    }
}
