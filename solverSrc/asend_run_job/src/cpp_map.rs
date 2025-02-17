use std::collections::BTreeMap;
use crate::constants::max_int;

#[derive(Clone)]
pub struct CppMap {
    m : BTreeMap<String, usize>,
}

impl CppMap {
    pub fn new() -> CppMap {
        CppMap {
            m : BTreeMap::new(),
        }
    }

    pub fn insert(&mut self, key : String, val : usize) {
        self.m.insert(key,val);
    }

    pub fn key_in_map(&mut self, key : & String) -> bool {
        match self.m.get(key) {
            None => false,
            Some(x) => true,
        }
    }

    pub fn at(&mut self, key : & String) -> usize {
        match self.m.get(key) {
            None => max_int,
            Some(x) => *x,
        }
    }
}