use crate::constants::*;
use crate::cpp_str::CppStr;

use std::collections::LinkedList;
use std::collections::BTreeMap;

#[derive(Clone)]
pub struct Material {
    pub name : CppStr,
    pub density : f64,
    pub modulus : [f64; 3],
    pub poisson_ratio : [f64; 3],
    pub shear_mod : [f64; 3],
    pub stiffness : [f64; 36],
    pub conductivity : [f64; 6],
    pub expansion : [f64; 6],
    pub spec_heat : f64,
    pub damping : [f64; 36],
    pub custom : BTreeMap<String, Vec<f64>>,
}

impl Material {
    pub fn new() -> Material {
        Material {
            name : CppStr::new(),
            density : 0f64,
            modulus : [0f64; 3],
            poisson_ratio : [0f64; 3],
            shear_mod : [0f64; 3],
            stiffness : [0f64; 36],
            conductivity : [0f64; 6],
            expansion : [0f64; 6],
            spec_heat : 0f64,
            damping : [0f64; 36],
            custom : BTreeMap::new(),
        }
    }
}

#[derive(Clone)]
pub struct Fluid {
    pub name : CppStr,
    pub viscosity : f64,
    pub ideal_gas : f64,
    pub therm_cond : f64,
    pub spec_heat : f64,
}

impl Fluid {
    pub fn new() -> Fluid {
        Fluid {
            name : CppStr::new(),
            viscosity : 0f64,
            ideal_gas : 0f64,
            therm_cond : 0f64,
            spec_heat : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct Layer {
    pub mat_name : CppStr,
    pub mat_ptr : usize,
    pub thickness : f64,
    pub angle : f64,
}

impl Layer {
    pub fn new() -> Layer {
        Layer {
            mat_name : CppStr::new(),
            mat_ptr : MAX_INT,
            thickness : 0f64,
            angle : 0f64,
        }
    }
}

#[derive(Clone)]
pub struct Section {
    pub this_type : CppStr,
    pub el_set_name : CppStr,
    pub mat_name : CppStr,
    pub fl_name : CppStr,
    pub mat_ptr : usize,
    pub fl_ptr : usize,
    pub orientation : [f64; 9],
    pub z_offset : f64,
    pub layers : LinkedList<Layer>,
    pub area : f64,
    pub area_moment : [f64; 5],
    pub polar_moment : f64,
    pub stiffness : [f64; 36],
    pub mass : [f64; 36],
    pub damping : [f64; 36],
    pub exp_load_coef : [f64; 6],
    pub conductivity : f64,
    pub spec_heat : f64,
    pub mass_per_el : f64,
    pub pot_coef : f64,
    pub pot_exp : f64,
    pub damp_coef : f64,
    pub damp_exp : f64,
    pub cond_coef : f64,
    pub rad_coef : f64,
    pub den_vis_coef : f64,
    pub temp_vis_coef : f64,
    pub turb_vis_coef : f64,
    pub grad_vturb_coef : f64,
    pub diss_turb_coef : f64,
    pub enth_coef : f64,
    pub enth_exp : f64,
    pub pres_coef : f64,
    pub pres_exp : f64,
    pub ref_temp : f64,
    pub ref_den : f64,
    pub ref_turb_e : f64,
    pub ref_enth : f64,
}

impl Section {
    pub fn new() -> Section {
        Section {
            this_type : CppStr::new(),
            el_set_name : CppStr::new(),
            mat_name : CppStr::new(),
            fl_name : CppStr::new(),
            mat_ptr : MAX_INT,
            fl_ptr : MAX_INT,
            orientation : [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            z_offset : 0f64,
            layers : LinkedList::new(),
            area : 0f64,
            area_moment : [0f64; 5],
            polar_moment : 0f64,
            stiffness : [0f64; 36],
            mass : [0f64; 36],
            damping : [0f64; 36],
            exp_load_coef : [0f64; 6],
            conductivity : 0f64,
            spec_heat : 0f64,
            mass_per_el : 0f64,
            pot_coef : 0f64,
            pot_exp : 0f64,
            damp_coef : 0f64,
            damp_exp : 0f64,
            cond_coef : 0f64,
            rad_coef : 0f64,
            den_vis_coef : 0f64,
            temp_vis_coef : 0f64,
            turb_vis_coef : 0f64,
            grad_vturb_coef : 0f64,
            diss_turb_coef : 0f64,
            enth_coef : 0f64,
            enth_exp : 0f64,
            pres_coef : 0f64,
            pres_exp : 0f64,
            ref_temp : 0f64,
            ref_den : 0f64,
            ref_turb_e : 0f64,
            ref_enth : 0f64,
        }
    
    }
}
