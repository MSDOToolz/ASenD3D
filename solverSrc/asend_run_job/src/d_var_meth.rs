use crate::design_var::*;
use crate::constants::*;
use crate::list_ent::*;
use crate::diff_doub::*;
use crate::nd_el_set::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


impl DesignVariable {
    pub fn set_active_time(&mut self, new_at : &mut [f64]) {
        self.active_time[0] = new_at[0];
        self.active_time[1] = new_at[1];
    }

    pub fn get_value_dfd0(&self, inp : &mut DiffDoub0) {
        inp.set_val_dfd0(& self.value);
        return;
    }

    pub fn get_value_dfd1(&self, inp : &mut DiffDoub1) {
        inp.set_val_dfd1(& self.diff_val);
        return;
    }

    pub fn add_comp_el(&mut self, eli : usize) {
        for el in self.comp_el_list.iter_mut() {
            if (*el == eli) {
                return;
            }
        }
        self.comp_el_list.push_back(eli);
        return;
    }

}


