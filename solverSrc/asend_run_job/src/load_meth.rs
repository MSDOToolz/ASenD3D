use crate::load::*;
use crate::constants::*;
use crate::nd_el_set::*;
use crate::cpp_str::CppStr;
use crate::cpp_map::CppMap;
use crate::fmath::*;


impl Load {
    pub fn set_act_time(&mut self, new_at : &mut [f64]) {
        self.active_time[0] = new_at[0];
        self.active_time[1] = new_at[1];
        return;
    }

    pub fn set_load(&mut self, new_ld : &mut [f64]) {
        self.load[0] = new_ld[0];
        self.load[1] = new_ld[1];
        self.load[2] = new_ld[2];
        self.load[3] = new_ld[3];
        self.load[4] = new_ld[4];
        self.load[5] = new_ld[5];
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

}


