use crate::list_ent::*;
use crate::cpp_str::CppStr;
use crate::constants::*;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct ConstraintTerm {
    pub node_set : CppStr,
    pub ns_ptr : usize,
    pub dof : usize,
    pub coef : f64,
}

impl ConstraintTerm {
    pub fn new() -> ConstraintTerm {
        ConstraintTerm {
            node_set : CppStr::from("none"),
            ns_ptr : MAX_INT,
            dof : 0usize,
            coef : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct Constraint {
    pub this_type : CppStr,
    pub terms : LinkedList<ConstraintTerm>,
    pub rhs : f64,
    pub rhs_vec : Vec<f64>,
    pub active_time : [f64; 2],
    pub is_active : bool,
    pub was_active : bool,
    pub mat : SparseMat,
    pub scale_fact : f64,
}

impl Constraint {
    pub fn new() -> Constraint {
        Constraint {
            this_type : CppStr::from("none"),
            terms : LinkedList::new(),
            rhs : 0f64,
            rhs_vec : Vec::new(),
            active_time : [0.0, 1.0e+100],
            is_active : false,
            was_active : false,
            mat : SparseMat::new(),
            scale_fact : 0f64,
        }
    }

}

#[derive(Clone)]
pub struct ConstraintList {
    pub const_vec : Vec<Constraint>,
}

impl ConstraintList {
    pub fn new() -> ConstraintList {
        ConstraintList {
            const_vec : Vec::new(),
        }
    }
}

