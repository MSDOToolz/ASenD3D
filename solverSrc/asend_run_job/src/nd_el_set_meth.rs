use crate::nd_el_set::*;
use crate::list_ent::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;


impl Set {
    pub fn add_if_absent(&mut self, new_i : usize) -> bool {
        for i in self.labels.iter_mut() {
            if (*i == new_i) {
                return  false;
            }
        }
        self.labels.push_back(new_i);
        return  true;
    }

}


