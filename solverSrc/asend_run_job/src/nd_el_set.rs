use crate::list_ent::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;

use std::collections::LinkedList;

#[derive(Clone)]
pub struct Set {
    pub name : CppStr,
    pub labels : LinkedList<usize>,
}

impl Set {
    pub fn new() -> Set {
        Set {
            name : CppStr::new(),
            labels : LinkedList::new(),
        }
    }
}
