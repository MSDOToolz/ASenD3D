use crate::constants::*;
use crate::fluid_domain::face::*;
use crate::fluid_domain::node::*;
use crate::fluid_domain::cell::EqnPrereq;
use crate::matrix_functions::*;

impl Face {
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

    pub fn update_area_norm(&self, fdat : &mut FaceData, nd_ar : &Vec<Node>) {
        let mut v1 = [DiffDoub1::new(); 3];
        let mut v2 = [DiffDoub1::new(); 3];
        let mut nd1 = [DiffDoub1::new(); 3];
        let mut nd2 = [DiffDoub1::new(); 3];
        let mut nd3 = [DiffDoub1::new(); 3];
        nd_ar[self.glob_nodes[0]].get_def_crd(&mut nd1);
        nd_ar[self.glob_nodes[1]].get_def_crd(&mut nd2);
        nd_ar[self.glob_nodes[2]].get_def_crd(&mut nd3);
        let mut mag = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();

        for i in 0..3 {
            v1[i].set_val_dfd1(&nd2[i]);
            v1[i].sub(&nd1[i]);
            v2[i].set_val_dfd1(&nd3[i]);
            v2[i].sub(&nd1[i]);
        }

        let mut cp = [DiffDoub1::new(); 3];
        
        cross_prod_dfd1(&mut cp, &v1, &v2);

        for i in 0..3 {
            tmp.set_val_dfd1(&cp[i]);
            tmp.sqr();
            mag.add(&tmp);
        }
        mag.sqt();
        fdat.area.set_val(0.5);
        fdat.area.mult(&mag);

        tmp.set_val(1.0);
        tmp.dvd(&mag);

        for i in 0..3 {
            fdat.normal[i].set_val_dfd1(&cp[i]);
            fdat.normal[i].mult(&tmp);
        }

    }

    pub fn set_data_to_cell(fdat : &mut FaceData, cell : &CellData) {
        fdat.den.set_val_dfd1(&cell.den);
        fdat.temp.set_val_dfd1(&cell.temp);
        fdat.turb.set_val_dfd1(&cell.turb);
        for i in 0..3 {
            fdat.vel[i].set_val_dfd1(&cell.vel[i]);
            fdat.v_rel[i].set_val_dfd1(&cell.v_rel[i]);
            fdat.t_grad[i].set_val_dfd1(&cell.t_grad[i]);
        }
        for i in 0..9 {
            fdat.v_grad[i].set_val_dfd1(&cell.v_grad[i]);
        }
    }

    pub fn set_single_to_avg(res : &mut DiffDoub1, d1 : &DiffDoub1, d2 : &DiffDoub1, wt1 : &DiffDoub1, wt2 : &DiffDoub1) {
        let mut tmp = DiffDoub1::new();

        tmp.set_val_dfd1(&wt1);
        tmp.mult(d1);
        res.set_val_dfd1(&tmp);
        tmp.set_val_dfd1(&wt2);
        tmp.mult(d2);
        res.add(&tmp);
    }

    pub fn set_data_to_avg(fdat : &mut FaceData, cell1 : &CellData, cell2 : &CellData, wt1 : &DiffDoub1, wt2 : &DiffDoub1) {

        Face::set_single_to_avg(&mut fdat.den, &cell1.den, &cell2.den, wt1, wt2);
        Face::set_single_to_avg(&mut fdat.temp, &cell1.temp, &cell2.temp, wt1, wt2);
        Face::set_single_to_avg(&mut fdat.turb, &cell1.turb, &cell2.turb, wt1, wt2);
        for i in 0..3 {
            Face::set_single_to_avg(&mut fdat.vel[i], &cell1.vel[i], &cell2.vel[i], wt1, wt2);
            Face::set_single_to_avg(&mut fdat.v_rel[i], &cell1.v_rel[i], &cell2.v_rel[i], wt1, wt2);
            Face::set_single_to_avg(&mut fdat.t_grad[i], &cell1.t_grad[i], &cell2.t_grad[i], wt1, wt2);
        }
        for i in 0..9 {
            Face::set_single_to_avg(&mut fdat.v_grad[i], &cell1.v_grad[i], &cell2.v_grad[i], wt1, wt2);
        }
    }

    pub fn get_incomp_pressure(pressure : &mut DiffDoub1, den : &DiffDoub1, temp : &DiffDoub1, pre : &EqnPrereq) {
        let mut tmp = DiffDoub1::new();

        pressure.set_val_dfd1(&pre.ref_pres);
            
        tmp.set_val_dfd1(den);
        tmp.sub(&pre.ref_den);
        tmp.mult(&pre.bulk_mod);
        tmp.dvd(&pre.ref_den);
        pressure.add(&tmp);
        
        tmp.set_val(3.0);
        tmp.mult(&pre.bulk_mod);
        tmp.mult(&pre.expansion);
        tmp.mult(temp);
        pressure.add(&tmp);

        // if pressure.val < 0.0 {
        //     pressure.set_val(0.0);
        // }
    }

    pub fn get_flux_vars(&self, fdat : &mut FaceData, pre : &EqnPrereq, fc_ar : &Vec<Face>, cdat : &Vec<CellData>) {
        let mut tmp = DiffDoub1::new();
        
        let hci = self.host_cell;
        let thci = match self.twin_id {
            MAX_INT => MAX_INT,
            _ => fc_ar[self.twin_id].host_cell,
        };

        if thci == MAX_INT {
            Face::set_data_to_cell(fdat, &cdat[hci]);
        }
        else {
            let mut t_avg = DiffDoub1::new();
            let mut dp = DiffDoub1::new();
            let mut sp_snd = DiffDoub1::new();
            let mut mach = DiffDoub1::new();
            let mut wt1 = DiffDoub1::new();
            let mut wt2 = DiffDoub1::new();
            let c1 = &cdat[hci];
            let c2 = &cdat[thci];

            wt1.set_val(0.5);
            Face::set_single_to_avg(&mut t_avg, &c1.temp, &c2.temp, &wt1, &wt1);
            for i in 0..3 {
                Face::set_single_to_avg(&mut tmp, &c1.v_rel[i], &c2.v_rel[i], &wt1, &wt1);
                tmp.mult(&fdat.normal[i]);
                dp.add(&tmp);
            }

            if pre.compressible {
                sp_snd.set_val(1.0);
                tmp.set_val_dfd1(&pre.ideal_gas);
                tmp.dvd(&pre.spec_heat);
                sp_snd.add(&tmp); // gamma
                sp_snd.mult(&pre.ideal_gas); // gamma*R
                tmp.set_val_dfd1(&pre.ref_temp);
                tmp.add(&t_avg);
                sp_snd.mult(&tmp); // gamma*R*T
                sp_snd.sqt();
            }
            else {
                sp_snd.set_val_dfd1(&pre.bulk_mod);
                sp_snd.dvd(&pre.ref_den);
                sp_snd.sqt();
            }

            mach.set_val_dfd1(&dp);
            mach.dvd(&sp_snd);

            if mach.val >= 1.0 {
                wt1.set_val(1.0);
                wt2.set_val(0.0);
            }
            else if mach.val > -1.0 {
                wt1.set_val(0.5);
                wt1.mult(&mach);
                tmp.set_val(0.5);
                wt1.add(&tmp);
                wt2.set_val(1.0);
                wt2.sub(&wt1);
            }
            else {
                wt1.set_val(0.0);
                wt2.set_val(1.0);
            }

            Face::set_data_to_avg(fdat, c1, c2, &wt1, &wt2);

        }
            
    }

    pub fn convective(&self, r_out : &mut [DiffDoub1], pre : &EqnPrereq, fdat : &FaceData) {
        let mut enth = DiffDoub1::new();
        let mut kin_e = DiffDoub1::new();
        let mut comn = DiffDoub1::new();

        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut k : usize;

        //common integrand
        comn.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&fdat.v_rel[i]);
            tmp2.mult(&fdat.normal[i]);
            comn.add(&tmp2);
        }
        comn.mult(&fdat.den);
        comn.mult(&fdat.area);

        //mass

        r_out[0].add(&comn);

        // momentum

        for i in 0..3 {
            k = i+1;
            tmp.set_val_dfd1(&comn);
            tmp.mult(&fdat.vel[i]);
            r_out[k].add(&tmp);
        }

        // energy

        tmp2.set_val_dfd1(&pre.spec_heat);
        tmp2.mult(&fdat.temp);
        enth.set_val_dfd1(&pre.ref_enth);
        enth.add(&tmp2);

        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&fdat.vel[i]);
            tmp2.sqr();
            tmp.add(&tmp2);
        }
        kin_e.set_val(0.5);
        kin_e.mult(&tmp);

        tmp.set_val_dfd1(&kin_e);
        tmp.add(&enth);
        tmp.add(&fdat.turb);

        tmp2.set_val_dfd1(&comn);
        tmp2.mult(&tmp);
        r_out[4].add(&tmp2);

        // turbulence

        tmp.set_val_dfd1(&comn);
        tmp.mult(&fdat.turb);
        r_out[5].add(&tmp);


    }

    pub fn viscous(&self, r_out : &mut [DiffDoub1], pre : &EqnPrereq, fdat : &FaceData) {
        let mut vis = DiffDoub1::new();
        let mut comn = [DiffDoub1::new(); 3];

        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut k : usize;

        // calculate viscosity
        if pre.compressible {
            vis.set_val_dfd1(&pre.viscosity);
            tmp.set_val_dfd1(&pre.turb_vis_coef);
            tmp.mult(&fdat.turb);
            vis.add(&tmp);
            vis.mult(&fdat.den);
            vis.dvd(&pre.ref_den);
        }
        else {
            vis.set_val_dfd1(&pre.viscosity);
            tmp.set_val_dfd1(&pre.temp_vis_coef);
            tmp.mult(&fdat.temp);
            vis.add(&tmp);
            tmp.set_val_dfd1(&pre.turb_vis_coef);
            tmp.mult(&fdat.turb);
            vis.add(&tmp);
        }

        //calculate common integrand
        mat_mul_ar_dfd1(&mut comn, &fdat.v_grad, &fdat.normal, 3, 3, 1);
        for i in 0..3 {
            comn[i].mult(&vis);
            comn[i].mult(&fdat.area);
        }

        // momentum

        for i in 0..3 {
            k = i + 1;
            tmp.set_val_dfd1(&comn[i]);
            tmp.neg();
            r_out[k].add(&tmp);
        }

        // energy

        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&fdat.vel[i]);
            tmp2.mult(&comn[i]);
            tmp.add(&tmp2);
        }

        tmp2.set_val_dfd1(&tmp);
        tmp2.neg();
        r_out[4].add(&tmp2);

    }

    pub fn pressure(&self, r_out : &mut [DiffDoub1], pre : &EqnPrereq, fdat : &FaceData) {
        let mut pressure = DiffDoub1::new();
        let mut comn = [DiffDoub1::new(); 3];

        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();
        let mut k : usize;

        if pre.compressible {
            pressure.set_val_dfd1(&pre.ref_temp);
            pressure.add(&fdat.temp);
            pressure.mult(&fdat.den);
            pressure.mult(&pre.ideal_gas);
        }
        else {
            Face::get_incomp_pressure(&mut pressure, &fdat.den, &fdat.temp, pre);

            if pressure.val == 0.0 {
                return;
            }
        }

        for i in 0..3 {
            comn[i].set_val_dfd1(&fdat.normal[i]);
            comn[i].mult(&pressure);
            comn[i].mult(&fdat.area);
        }

        //momentum

        for i in 0..3 {
            k = i + 1;
            r_out[k].add(&comn[i]);
        }

        //energy
        tmp.set_val(0.0);
        for i in 0..3 {
            tmp2.set_val_dfd1(&fdat.vel[i]);
            tmp2.mult(&comn[i]);
            tmp.add(&tmp2);
        }

        r_out[4].add(&tmp);
    }

    pub fn heat_flux(&self, r_out : &mut [DiffDoub1], pre : &EqnPrereq, fdat : &FaceData) {
        let mut tmp = DiffDoub1::new();
        let mut tmp2 = DiffDoub1::new();

        for i in 0..3 {
            tmp2.set_val_dfd1(&fdat.t_grad[i]);
            tmp2.mult(&fdat.normal[i]);
            tmp.add(&tmp2);
        }
        tmp.mult(&pre.conductivity);
        tmp.mult(&fdat.area);
        tmp.neg();

        r_out[4].add(&tmp);

    }

}