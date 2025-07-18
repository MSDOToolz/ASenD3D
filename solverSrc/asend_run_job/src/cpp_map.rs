use std::collections::BTreeMap;
use crate::constants::MAX_INT;

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

    pub fn key_in_map(&self, key : & String) -> bool {
        match self.m.get(key) {
            None => false,
            Some(_x) => true,
        }
    }

    pub fn at(&self, key : & String) -> usize {
        match self.m.get(key) {
            None => MAX_INT,
            Some(x) => *x,
        }
    }
}