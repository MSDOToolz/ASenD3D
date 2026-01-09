mod cell;
mod node;
mod face;
mod design_var;
mod sub_domain;
mod constraint;
mod load;

use crate::fluid_domain::node::*;
use crate::fluid_domain::cell::*;
use crate::fluid_domain::face::*;
use crate::fluid_domain::sub_domain::*;
use crate::fluid_domain::constraint::*;
use crate::fluid_domain::load::*;
use crate::fluid_domain::design_var::*;

use crate::nd_el_set::*;
use crate::cpp_str::*;

use std::collections::BTreeMap;

pub struct FluidDomain {
    pub nodes : Vec<Node>,
    pub faces : Vec<Face>,
    pub cells : Vec<Cell>,
    pub node_sets : Vec<Set>,
    pub ns_map : BTreeMap<String, usize>,
    pub cell_sets : Vec<Set>,
    pub cs_map : BTreeMap<String, usize>,
    pub sub_domains : Vec<SubDomain>,
    pub fluids : Vec<Fluid>,
    pub constraints : ConstraintList,
    pub loads : Vec<Load>,
    pub design_vars : Vec<DesignVariable>,

    pub init_stat_file : CppStr,
}

impl FluidDomain {
    pub fn new() -> FluidDomain {
        FluidDomain {
            nodes : Vec::new(),
            faces : Vec::new(),
            cells : Vec::new(),
            node_sets : Vec::new(),
            ns_map : BTreeMap::new(),
            cell_sets : Vec::new(),
            cs_map : BTreeMap::new(),
            sub_domains : Vec::new(),
            fluids : Vec::new(),
            constraints : ConstraintList::new(),
            loads : Vec::new(),
            design_vars : Vec::new(),

            init_stat_file : CppStr::new(),
        }
    }
}

mod input;
mod prep;