use crate::constants::MAX_INT;

#[derive(Clone)]
pub struct FluidCell {
    pub nodes : [usize; 4],
    pub faces : [usize; 4],
}

impl FluidCell {
    pub fn new() -> FluidCell {
        FluidCell {
            nodes : [MAX_INT; 4],
            faces : [MAX_INT; 4],
        }
    }
}

pub mod fluid_cell_meth;