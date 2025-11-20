use crate::constants::*;
use crate::cpp_str::CppStr;

#[derive(Clone)]
pub struct Fluid {
    pub name : CppStr,
    pub viscosity : f64,
    pub ideal_gas : f64,
    pub therm_cond : f64,
    pub expansion : f64,
    pub spec_heat : f64,
    pub bulk_modulus : f64,
    pub compressible : bool,
    pub ref_temp : f64,
    pub ref_pres : f64,
    pub ref_den : f64,
    pub ref_enth : f64,
    pub temp_vis_coef : f64,
    pub turb_vis_coef : f64,
    pub grad_turb_coef : f64,
    pub diss_turb_coef : f64,
}

impl Fluid {
    pub fn new() -> Fluid {
        Fluid {
            name : CppStr::new(),
            viscosity : 0f64,
            ideal_gas : 0f64,
            therm_cond : 0f64,
            expansion : 0f64,
            spec_heat : 0f64,
            bulk_modulus : 0f64,
            compressible : true,
            ref_temp : 0f64,
            ref_pres : 0f64,
            ref_den : 0f64,
            ref_enth : 0f64,
            temp_vis_coef : 0f64,
            turb_vis_coef : 0f64,
            grad_turb_coef : 0f64,
            diss_turb_coef : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct SubDomain {
    pub cell_set_name : CppStr,
    pub fluid_name : CppStr,
    pub fluid_ptr : usize,
}

impl SubDomain {
    pub fn new() -> SubDomain {
        SubDomain {
            cell_set_name : CppStr::new(),
            fluid_name : CppStr::new(),
            fluid_ptr : 0usize,
        }
    }
}