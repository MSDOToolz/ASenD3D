use crate::fluid_domain::node::*;
use crate::fluid_domain::design_var::*;

impl Node {
    pub fn add_cell(&mut self, cell_i : usize) {
        for nc in self.cell_lst.iter() {
            if *nc == cell_i {
                return;
            }
        }
        self.cell_lst.push_back(cell_i);
    }

    pub fn add_conn_nd(&mut self, nd_i : usize) {
        for nd in self.nd_lst.iter() {
            if *nd == nd_i {
                return;
            }
        }
        self.nd_lst.push_back(nd_i);
    }

    pub fn calc_coord(&mut self, dvars : &Vec<DesignVariable>) {
        let mut dv : &DesignVariable;
        let mut comp : usize;
        let mut coef = DiffDoub1::new();

        self.coord_dfd1[0].set_val(self.coord[0]);
        self.coord_dfd1[1].set_val(self.coord[1]);
        self.coord_dfd1[2].set_val(self.coord[2]);

        for dvi in self.d_var_lst.iter() {
            dv = &dvars[dvi.int_dat];
            if dv.category.s == "nodeCoord" {
                comp = dv.component - 1;
                coef.set_val(dvi.doub_dat);
                coef.mult(&dv.diff_val);
                self.coord_dfd1[comp].add(&coef);
            }
        }
    }

    pub fn get_def_crd(&self, coord : &mut [DiffDoub1]) {
        for i in 0..3 {
            coord[i].set_val_dfd1(&self.coord_dfd1[i]);
            coord[i].add(&self.displacement[i]);
        }
    }

    pub fn get_vrel(&self, vrel : &mut [DiffDoub1]) {
        for i in 0..3 {
            vrel[i].set_val_dfd1(&self.fl_vel[i]);
            vrel[i].sub(&self.velocity[i]);
        }
    }

}