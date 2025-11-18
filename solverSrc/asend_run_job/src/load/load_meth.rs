use crate::load::*;
use crate::fmath::*;


impl Load {
    pub fn set_act_time(&mut self, new_at : &mut [f64]) {
        self.active_time[0] = new_at[0];
        self.active_time[1] = new_at[1];
        return;
    }

    pub fn set_load(&mut self, new_ld : LoadTimePt) {
        self.load.push_back(new_ld);
        return;
    }

    pub fn set_norm_dir(&mut self, new_ndir : &mut [f64]) {
        self.normal_dir[0] = new_ndir[0];
        self.normal_dir[1] = new_ndir[1];
        self.normal_dir[2] = new_ndir[2];
        let mut mag : f64 =  self.normal_dir[0]*self.normal_dir[0] + self.normal_dir[1]*self.normal_dir[1] + self.normal_dir[2]*self.normal_dir[2];
        mag = 1.0/sqrt(mag);
        self.normal_dir[0] = mag*self.normal_dir[0];
        self.normal_dir[1] = mag*self.normal_dir[1];
        self.normal_dir[2] = mag*self.normal_dir[2];
        return;
    }

    pub fn set_center(&mut self, new_cent : &mut [f64]) {
        self.center[0] = new_cent[0];
        self.center[1] = new_cent[1];
        self.center[2] = new_cent[2];
        return;
    }

    pub fn set_axis(&mut self, new_axis : &mut [f64]) {
        self.axis[0] = new_axis[0];
        self.axis[1] = new_axis[1];
        self.axis[2] = new_axis[2];
        let mut mag : f64 =  self.axis[0]*self.axis[0] + self.axis[1]*self.axis[1] + self.axis[2]*self.axis[2];
        mag = 1.0/sqrt(mag);
        self.axis[0] = mag*self.axis[0];
        self.axis[1] = mag*self.axis[1];
        self.axis[2] = mag*self.axis[2];
        return;
    }

    pub fn get_load(&self, ld : &mut [f64], time : f64) {
        let mut prev_pt = match self.load.front() {
            None => panic!("Error: time series for load is and empty list"),
            Some(x) => x,
        };
        for pt in self.load.iter() {
            if time >= prev_pt.time && time < pt.time {
                let dt = time - prev_pt.time;
                let mut slope : f64;
                for i in 0..6 {
                    slope = (pt.value[i] - prev_pt.value[i])/(pt.time - prev_pt.time);
                    ld[i] = prev_pt.value[i] + slope*dt;
                }
                return;
            }
            prev_pt = pt;
        }
    }

}


