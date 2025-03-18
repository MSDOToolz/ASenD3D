use crate::constants::*;
use crate::list_ent::*;
use crate::lower_tri_mat::*;
use crate::job::*;
use crate::node::*;
use crate::element::*;
use crate::face::*;
use crate::nd_el_set::*;
use crate::section::*;
use crate::load::*;
use crate::constraint::*;
use crate::design_var::*;
use crate::objective::*;
use crate::diff_doub::*;
use crate::scratch::*;
use crate::cpp_map::CppMap;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Model {
    pub nodes : Vec<Node>,
    pub elements : Vec<Element>,
    pub faces : Vec<Face>,
    pub node_sets : Vec<Set>,
    pub ns_map : CppMap,
    pub element_sets : Vec<Set>,
    pub es_map : CppMap,
    pub sections : Vec<Section>,
    pub materials : Vec<Material>,
    pub fluids : Vec<Fluid>,
    pub elastic_const : ConstraintList,
    pub thermal_const : ConstraintList,
    pub elastic_loads : Vec<Load>,
    pub thermal_loads : Vec<Load>,
    pub design_vars : Vec<DesignVariable>,
    pub obj : Objective,
    pub job : Vec<JobCommand>,
    pub el_mat_dim : usize,
    pub tot_glob_dof : usize,
    pub an_prep_run : bool,
    pub time_steps_saved : usize,
    pub solve_cmd : usize,
    pub modal_cmd : usize,
    pub d0_pre : DiffDoub0StressPrereq,
    pub d1_pre : DiffDoub1StressPrereq,
    pub scratch : LinkedList<FltScr>,
    pub d0_scratch : LinkedList<DiffDoub0Scr>,
    pub d1_scratch : LinkedList<DiffDoub1Scr>,
    pub elastic_mat : SparseMat,
    pub elastic_lt : LowerTriMat,
    pub elastic_ld_vec : Vec<f64>,
    pub elastic_sol_vec : Vec<f64>,
    pub elastic_scaled : bool,
    pub therm_mat : SparseMat,
    pub therm_lt : LowerTriMat,
    pub therm_ld_vec : Vec<f64>,
    pub therm_sol_vec : Vec<f64>,
    pub therm_scaled : bool,
    pub eig_vecs : Vec<f64>,
    pub eig_vals : Vec<f64>,
    pub diag_mass : Vec<f64>,
    pub load_fact : Vec<f64>,
    pub temp_v1 : Vec<f64>,
    pub temp_v2 : Vec<f64>,
    pub temp_v3 : Vec<f64>,
    pub temp_v4 : Vec<f64>,
    pub temp_v5 : Vec<f64>,
    pub temp_v6 : Vec<f64>,
    pub temp_d1 : Vec<DiffDoub0>,
    pub d_ld_u : Vec<f64>,
    pub d_ld_v : Vec<f64>,
    pub d_ld_a : Vec<f64>,
    pub d_ld_t : Vec<f64>,
    pub d_ld_tdot : Vec<f64>,
    pub u_adj : Vec<f64>,
    pub v_adj : Vec<f64>,
    pub a_adj : Vec<f64>,
    pub t_adj : Vec<f64>,
    pub tdot_adj : Vec<f64>,
    pub d_rud_d : Vec<DiffDoub1>,
    pub d_rtd_d : Vec<DiffDoub1>,
    pub el_in_d : Vec<usize>,
    pub d_ld_d : Vec<f64>,
}

impl Model {
    pub fn new() -> Model {
        let mut new_mod = Model {
            nodes : Vec::new(),
            elements : Vec::new(),
            faces : Vec::new(),
            node_sets : Vec::new(),
            ns_map : CppMap::new(),
            element_sets : Vec::new(),
            es_map : CppMap::new(),
            sections : Vec::new(),
            materials : Vec::new(),
            fluids : Vec::new(),
            elastic_const : ConstraintList::new(),
            thermal_const : ConstraintList::new(),
            elastic_loads : Vec::new(),
            thermal_loads : Vec::new(),
            design_vars : Vec::new(),
            obj : Objective::new(),
            job : Vec::new(),
            el_mat_dim : 0usize,
            tot_glob_dof : 0usize,
            an_prep_run : false,
            time_steps_saved : 0usize,
            solve_cmd : MAX_INT,
            modal_cmd : MAX_INT,
            d0_pre : DiffDoub0StressPrereq::new(),
            d1_pre : DiffDoub1StressPrereq::new(),
            scratch : LinkedList::new(),
            d0_scratch : LinkedList::new(),
            d1_scratch : LinkedList::new(),
            elastic_mat : SparseMat::new(),
            elastic_lt : LowerTriMat::new(),
            elastic_ld_vec : Vec::new(),
            elastic_sol_vec : Vec::new(),
            elastic_scaled : false,
            therm_mat : SparseMat::new(),
            therm_lt : LowerTriMat::new(),
            therm_ld_vec : Vec::new(),
            therm_sol_vec : Vec::new(),
            therm_scaled : false,
            eig_vecs : Vec::new(),
            eig_vals : Vec::new(),
            diag_mass : Vec::new(),
            load_fact : Vec::new(),
            temp_v1 : Vec::new(),
            temp_v2 : Vec::new(),
            temp_v3 : Vec::new(),
            temp_v4 : Vec::new(),
            temp_v5 : Vec::new(),
            temp_v6 : Vec::new(),
            temp_d1 : Vec::new(),
            d_ld_u : Vec::new(),
            d_ld_v : Vec::new(),
            d_ld_a : Vec::new(),
            d_ld_t : Vec::new(),
            d_ld_tdot : Vec::new(),
            u_adj : Vec::new(),
            v_adj : Vec::new(),
            a_adj : Vec::new(),
            t_adj : Vec::new(),
            tdot_adj : Vec::new(),
            d_rud_d : Vec::new(),
            d_rtd_d : Vec::new(),
            el_in_d : Vec::new(),
            d_ld_d : Vec::new(),
        };
        for _i in 0..10 {
            new_mod.d0_scratch.push_back(DiffDoub0Scr::new());
            new_mod.d1_scratch.push_back(DiffDoub1Scr::new());
        }
        for _i in 0..5 {
            new_mod.scratch.push_back(FltScr::new());
        }
        new_mod
    }
}
