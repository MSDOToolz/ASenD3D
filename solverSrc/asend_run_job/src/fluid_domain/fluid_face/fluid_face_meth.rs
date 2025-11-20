use crate::fluid_face::*;

impl FluidFace {
    pub fn set_node(&mut self, place : usize, glob_nd : usize) {
        self.glob_nodes[place] = glob_nd;
    }

    pub fn sorted_nodes(&mut self, srt_nds : &mut [usize]) {
        let mut i4 : usize;
        let mut swap : usize;
        for i1 in 0..3 {
            srt_nds[i1] = self.glob_nodes[i1];
        }
        for _i1 in 0..2 {
            for i2 in 0..2 {
                i4 = i2 + 1;
                if srt_nds[i4] < srt_nds[i2] {
                    swap = srt_nds[i2];
                    srt_nds[i2] = srt_nds[i4];
                    srt_nds[i4] = swap;
                }
            }
        }
    }

    pub fn get_low_nd(&mut self) -> usize {
        let mut low_nd : usize = self.glob_nodes[0];
        for i1 in 1..3 {
            if self.glob_nodes[i1] < low_nd {
                low_nd = self.glob_nodes[i1];
            }
        }
        return  low_nd;
    }

}