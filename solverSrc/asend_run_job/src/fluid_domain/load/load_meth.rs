use crate::fluid_domain::load::*;
use crate::list_ent::*;

impl Load {
    pub fn get_current_ld(&self, ld : &mut [f64], time : f64) {
        let mut prev : &QuadFloat = match self.load.front() {
            None => panic!("Error: empty time point list in load"),
            Some(x) => x,
        };
        for pt in self.load.iter() {
            if pt.f1 > time {
                ld[0] = prev.f2 + (pt.f2 - prev.f2)*(time - prev.f1)/(pt.f1 - prev.f1);
                ld[1] = prev.f3 + (pt.f3 - prev.f3)*(time - prev.f1)/(pt.f1 - prev.f1);
                ld[2] = prev.f4 + (pt.f4 - prev.f4)*(time - prev.f1)/(pt.f1 - prev.f1);
                return;
            }
            prev = pt;
        }
        ld[0] = prev.f2;
        ld[1] = prev.f3;
        ld[2] = prev.f4;
    }
}