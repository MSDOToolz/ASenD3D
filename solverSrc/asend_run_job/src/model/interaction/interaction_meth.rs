use crate::model::interaction::*;
use crate::constants::MAX_INT;
use crate::model::node::Node;
use crate::model::element::*;
use crate::model::section::*;
use crate::model::design_var::*;
use crate::model::job::*;
use crate::matrix_functions::*;
use crate::spatial_grid::SpatialGrid;
use crate::nd_el_set::Set;
use crate::cpp_map::*;
use crate::cpp_str::*;

impl Interaction {
    pub fn is_active(&self, time : f64) -> bool {
        time >= self.active_time[0] && time <= self.active_time[1]
    }

    pub fn get_pot_coef(&self, time : f64) -> f64 {
        let mut pt = match self.pot_coef.front() {
            None => panic!("Error: empty potential coefficients in interaction"),
            Some(x) => x.f1,
        };
        let mut pv = match self.pot_coef.front() {
            None => panic!("Error: empty potential coefficients in interaction"),
            Some(x) => x.f2,
        };
        for ent in self.pot_coef.iter() {
            if ent.f1 > time {
                return pv + (ent.f2 - pv)*(time - pt)/(ent.f1 - pt);
            }
            pt = ent.f1;
            pv = ent.f2;
        }
        return pv;
    }

    pub fn get_damp_coef(&self, time : f64) -> f64 {
        let mut pt = match self.damp_coef.front() {
            None => panic!("Error: empty potential coefficients in interaction"),
            Some(x) => x.f1,
        };
        let mut pv = match self.damp_coef.front() {
            None => panic!("Error: empty potential coefficients in interaction"),
            Some(x) => x.f2,
        };
        for ent in self.damp_coef.iter() {
            if ent.f1 > time {
                return pv + (ent.f2 - pv)*(time - pt)/(ent.f1 - pt);
            }
            pt = ent.f1;
            pv = ent.f2;
        }
        return pv;
    }

//dup1

    pub fn get_particle_area_vol_dfd0(&self, area : &mut DiffDoub0, vol : &mut DiffDoub0, nearest : &Vec<usize>, near_dist : &Vec<f64>) {
        let mut tmp = DiffDoub0::new();
        area.set_val(0.0);
        vol.set_val(0.0);
        let mut ct = 0;
        for i in 0..self.max_nbrs {
            if nearest[i] < MAX_INT && near_dist[i]/near_dist[0] < self.max_ratio  {
                tmp.set_val(near_dist[i].powf(2.0));
                area.add(&tmp);
                tmp.set_val(near_dist[i].powf(3.0));
                vol.add(&tmp);
                ct += 1;
            }
        }
        if ct == 0 {
            area.set_val(1.0e+100);
            vol.set_val(1.0e+100);
            return;
        }
        tmp.set_val(3.14159265358979/(ct as f64));
        area.mult(&tmp);
        tmp.set_val(0.523598775598298/(ct as f64));
        vol.mult(&tmp)
    }

    pub fn get_ig_pressure_dfd0(pressure : &mut DiffDoub0, pre : &DiffDoub0StressPrereq, volume : &DiffDoub0) {
        pressure.set_val_dfd0(&pre.ref_temp);
        pressure.add(&pre.glob_temp[0]);
        pressure.mult(&pre.ideal_gas);
        pressure.mult(&pre.mass_per_el);
        pressure.dvd(volume);
    }

    pub fn get_incomp_pressure_dfd0(pressure : &mut DiffDoub0, pre : &DiffDoub0StressPrereq, volume : &DiffDoub0) {
        let mut tmp = DiffDoub0::new();
        
        pressure.set_val_dfd0(&pre.ref_pres);
        
        tmp.set_val_dfd0(&pre.mass_per_el);
        tmp.dvd(volume);
        tmp.sub(&pre.ref_den);
        tmp.dvd(&pre.ref_den);
        tmp.mult(&pre.bulk_mod);
        pressure.add(&tmp);

        tmp.set_val(3.0);
        tmp.mult(&pre.bulk_mod);
        tmp.mult(&pre.therm_exp[0]);
        tmp.mult(&pre.glob_temp[0]);
        pressure.add(&tmp);

        if pressure.val < 0.0 {
            pressure.set_val(0.0);
        }
    }

    pub fn get_pres_frc_coef_dfd0(pre : &mut DiffDoub0StressPrereq, pressure : &DiffDoub0, area : &DiffDoub0, dist : f64) {
        pre.frc_fld_coef[0].set_val(0.083333333333333*dist.powf(pre.frc_fld_exp[0].val)); // 1/12 * d^(exp)
        pre.frc_fld_coef[0].mult(pressure);
        pre.frc_fld_coef[0].mult(area);
        pre.frc_fld_coef[0].neg();
    }

    pub fn get_global_r_dfd0(&self, glob_r : &mut Vec<DiffDoub0>, dr_du : &mut SparseMat, discipline : usize, nd_in_set : &mut Vec<bool>, nd_mass : &Vec<DiffDoub0>, 
        g_list : &SpatialGrid, g_out : &mut Vec<usize>, nearest : &mut Vec<usize>, near_dist : &mut Vec<f64>, dummy_el : &mut Element, pre : &mut DiffDoub0StressPrereq, 
        time : f64, get_matrix : bool, cmd : &JobCommand, n_sets : &Vec<Set>, nodes : &Vec<Node>, dv_ar : &Vec<DesignVariable>) {
        //discipline = 0: elastic, 1: thermal

        let mut crd1 = [DiffDoub0::new(); 3];
        let mut fcrd1 = [0f64; 3];
        let mut crd2 = [DiffDoub0::new(); 3];
        let mut fcrd2 = [0f64; 3];
        let mut dist : f64;
        let mut lst_len : usize;
        let mut inserted : bool;
        let mut part_area = DiffDoub0::new();
        let mut part_vol = DiffDoub0::new();
        let mut pressure = DiffDoub0::new();
        let mut tmp = DiffDoub0::new();
        let mut i1 : usize;
        
        for nd in n_sets[self.set_pt2].labels.iter() {
            nd_in_set[*nd] = true;
        }

        dummy_el.design_vars = self.dvars.clone();


        pre.frc_fld_coef[0].set_val(self.get_pot_coef(time));
        dummy_el.get_gen_prop_dfd0(&mut pre.frc_fld_coef[0],&mut CppStr::from("potFldCoef"), dv_ar);
        
        pre.frc_fld_coef[1].set_val(self.get_damp_coef(time));
        dummy_el.get_gen_prop_dfd0(&mut pre.frc_fld_coef[1],&mut CppStr::from("dampFldCoef"), dv_ar);
        
        pre.frc_fld_exp[0].set_val(self.pot_exp);
        
        pre.frc_fld_exp[1].set_val(self.damp_exp);
        
        pre.thrm_fld_coef[0].set_val(self.cond_coef);
        dummy_el.get_gen_prop_dfd0(&mut pre.thrm_fld_coef[0],&mut CppStr::from("condCoef"), dv_ar);
        
        pre.thrm_fld_coef[1].set_val(self.rad_coef);
        dummy_el.get_gen_prop_dfd0(&mut pre.thrm_fld_coef[1],&mut CppStr::from("radCoef"), dv_ar);
        
        pre.ref_temp.set_val(self.ref_temp);

        pre.ideal_gas.set_val(self.ideal_gas);
        dummy_el.get_gen_prop_dfd0(&mut pre.ideal_gas, &mut CppStr::from("idealGasConstant"), dv_ar);

        pre.bulk_mod.set_val(self.bulk_mod);
        dummy_el.get_gen_prop_dfd0(&mut pre.bulk_mod, &mut CppStr::from("bulkModulus"), dv_ar);

        pre.therm_exp[0].set_val(self.therm_exp);
        dummy_el.get_gen_prop_dfd0(&mut pre.therm_exp[0], &mut CppStr::from("thermalExp"), dv_ar);

        pre.ref_den.set_val(self.ref_den);

        pre.ref_pres.set_val(self.ref_pres);

        
        for nd in n_sets[self.set_pt1].labels.iter() {
            if self.max_nbrs < MAX_INT {
                for i2 in 0..self.max_nbrs {
                    nearest[i2] = MAX_INT;
                    near_dist[i2] = 1.0e+100;
                }
            }
            nodes[*nd].get_def_crd_dfd0(&mut crd1);
            fcrd1[0] = crd1[0].val;
            fcrd1[1] = crd1[1].val;
            fcrd1[2] = crd1[2].val;
            lst_len = g_list.get_in_radius(g_out, g_out.len(), &fcrd1, self.max_dist);
            //if self.max_nbrs == MAX_INT {
            for nb in 0..lst_len {
                if nd_in_set[g_out[nb]] && g_out[nb] != *nd {
                    nodes[g_out[nb]].get_def_crd_dfd0(&mut crd2);
                    fcrd2[0] = crd2[0].val;
                    fcrd2[1] = crd2[1].val;
                    fcrd2[2] = crd2[2].val;
                    dist = get_dist(&fcrd1, &fcrd2);
                    if dist <= self.max_dist {
                        if self.max_nbrs == MAX_INT {
                            dummy_el.nodes[0] = *nd;
                            dummy_el.nodes[1] = g_out[nb];
                            dummy_el.get_all_nd_var_dfd0(pre, nodes);
                            // if self.ideal_gas > 0.0 {
                            //     self.get_ig_frc_coef_dfd0(pre, &nd_mass[*nd], &nd_mass[g_out[nb]], dist);
                            // }
                            match discipline {
                                0 => dummy_el.put_ru_frc_fld_dfd0(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                                1 => dummy_el.put_rt_frc_fld_dfd0(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                                _ => (),
                            }
                        }
                        else {
                            inserted = false;
                            i1 = 0;
                            while !inserted && i1 < self.max_nbrs {
                                if nearest[i1] == MAX_INT {
                                    nearest[i1] = g_out[nb];
                                    near_dist[i1] = dist;
                                    inserted = true;
                                }
                                else if dist < near_dist[i1] {
                                    for i2 in (i1..(self.max_nbrs - 1)).rev() {
                                        nearest[i2+1] = nearest[i2];
                                        near_dist[i2+1] = near_dist[i2];
                                    }
                                    nearest[i1] = g_out[nb];
                                    near_dist[i1] = dist;
                                    inserted = true;
                                }
                                i1 += 1;
                            }
                        }
                    }
                }
            }
            //}
            if self.max_nbrs < MAX_INT {
                if self.ideal_gas > 0.0 || self.bulk_mod > 0.0 {
                    //part_vol.set_val(self.get_particle_volume(nearest, near_dist));
                    self.get_particle_area_vol_dfd0(&mut part_area, &mut part_vol, nearest, near_dist);
                    pre.glob_temp[0].set_val(nodes[*nd].temperature);
                    pre.mass_per_el.set_val_dfd0(&nd_mass[*nd]);
                    if self.ideal_gas > 0.0 {
                        Interaction::get_ig_pressure_dfd0(&mut pressure, pre, &part_vol);
                    }
                    else {
                        Interaction::get_incomp_pressure_dfd0(&mut pressure, pre, &part_vol);
                    }
                    if self.set_pt1 == self.set_pt2 {
                        tmp.set_val(0.5);
                        pressure.mult(&tmp);
                    }
                }
                for i2 in 0..self.max_nbrs {
                    if nearest[i2] < MAX_INT && near_dist[i2]/near_dist[0] < self.max_ratio {
                        dummy_el.nodes[0] = *nd;
                        dummy_el.nodes[1] = nearest[i2];
                        dummy_el.get_all_nd_var_dfd0(pre, nodes);
                        if self.ideal_gas > 0.0 || self.bulk_mod > 0.0 {
                            Interaction::get_pres_frc_coef_dfd0(pre, &pressure, &part_area, near_dist[i2]);
                        }
                        match discipline {
                            0 => dummy_el.put_ru_frc_fld_dfd0(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                            1 => dummy_el.put_rt_frc_fld_dfd0(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                            _ => (),
                        }
                    }
                }
            }
        }

        for nd in n_sets[self.set_pt2].labels.iter() {
            nd_in_set[*nd] = false;
        }
    }

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

    pub fn get_particle_area_vol_dfd1(&self, area : &mut DiffDoub1, vol : &mut DiffDoub1, nearest : &Vec<usize>, near_dist : &Vec<f64>) {
        let mut tmp = DiffDoub1::new();
        area.set_val(0.0);
        vol.set_val(0.0);
        let mut ct = 0;
        for i in 0..self.max_nbrs {
            if nearest[i] < MAX_INT && near_dist[i]/near_dist[0] < self.max_ratio  {
                tmp.set_val(near_dist[i].powf(2.0));
                area.add(&tmp);
                tmp.set_val(near_dist[i].powf(3.0));
                vol.add(&tmp);
                ct += 1;
            }
        }
        if ct == 0 {
            area.set_val(1.0e+100);
            vol.set_val(1.0e+100);
            return;
        }
        tmp.set_val(3.14159265358979/(ct as f64));
        area.mult(&tmp);
        tmp.set_val(0.523598775598298/(ct as f64));
        vol.mult(&tmp)
    }

    pub fn get_ig_pressure_dfd1(pressure : &mut DiffDoub1, pre : &DiffDoub1StressPrereq, volume : &DiffDoub1) {
        pressure.set_val_dfd1(&pre.ref_temp);
        pressure.add(&pre.glob_temp[0]);
        pressure.mult(&pre.ideal_gas);
        pressure.mult(&pre.mass_per_el);
        pressure.dvd(volume);
    }

    pub fn get_incomp_pressure_dfd1(pressure : &mut DiffDoub1, pre : &DiffDoub1StressPrereq, volume : &DiffDoub1) {
        let mut tmp = DiffDoub1::new();
        
        pressure.set_val_dfd1(&pre.ref_pres);
        
        tmp.set_val_dfd1(&pre.mass_per_el);
        tmp.dvd(volume);
        tmp.sub(&pre.ref_den);
        tmp.dvd(&pre.ref_den);
        tmp.mult(&pre.bulk_mod);
        pressure.add(&tmp);

        tmp.set_val(3.0);
        tmp.mult(&pre.bulk_mod);
        tmp.mult(&pre.therm_exp[0]);
        tmp.mult(&pre.glob_temp[0]);
        pressure.add(&tmp);

        if pressure.val < 0.0 {
            pressure.set_val(0.0);
        }
    }

    pub fn get_pres_frc_coef_dfd1(pre : &mut DiffDoub1StressPrereq, pressure : &DiffDoub1, area : &DiffDoub1, dist : f64) {
        pre.frc_fld_coef[0].set_val(0.083333333333333*dist.powf(pre.frc_fld_exp[0].val)); // 1/12 * d^(exp)
        pre.frc_fld_coef[0].mult(pressure);
        pre.frc_fld_coef[0].mult(area);
        pre.frc_fld_coef[0].neg();
    }

    pub fn get_global_r_dfd1(&self, glob_r : &mut Vec<DiffDoub1>, dr_du : &mut SparseMat, discipline : usize, nd_in_set : &mut Vec<bool>, nd_mass : &Vec<DiffDoub1>, 
        g_list : &SpatialGrid, g_out : &mut Vec<usize>, nearest : &mut Vec<usize>, near_dist : &mut Vec<f64>, dummy_el : &mut Element, pre : &mut DiffDoub1StressPrereq, 
        time : f64, get_matrix : bool, cmd : &JobCommand, n_sets : &Vec<Set>, nodes : &Vec<Node>, dv_ar : &Vec<DesignVariable>) {
        //discipline = 0: elastic, 1: thermal

        let mut crd1 = [DiffDoub1::new(); 3];
        let mut fcrd1 = [0f64; 3];
        let mut crd2 = [DiffDoub1::new(); 3];
        let mut fcrd2 = [0f64; 3];
        let mut dist : f64;
        let mut lst_len : usize;
        let mut inserted : bool;
        let mut part_area = DiffDoub1::new();
        let mut part_vol = DiffDoub1::new();
        let mut pressure = DiffDoub1::new();
        let mut tmp = DiffDoub1::new();
        let mut i1 : usize;
        
        for nd in n_sets[self.set_pt2].labels.iter() {
            nd_in_set[*nd] = true;
        }

        dummy_el.design_vars = self.dvars.clone();


        pre.frc_fld_coef[0].set_val(self.get_pot_coef(time));
        dummy_el.get_gen_prop_dfd1(&mut pre.frc_fld_coef[0],&mut CppStr::from("potFldCoef"), dv_ar);
        
        pre.frc_fld_coef[1].set_val(self.get_damp_coef(time));
        dummy_el.get_gen_prop_dfd1(&mut pre.frc_fld_coef[1],&mut CppStr::from("dampFldCoef"), dv_ar);
        
        pre.frc_fld_exp[0].set_val(self.pot_exp);
        
        pre.frc_fld_exp[1].set_val(self.damp_exp);
        
        pre.thrm_fld_coef[0].set_val(self.cond_coef);
        dummy_el.get_gen_prop_dfd1(&mut pre.thrm_fld_coef[0],&mut CppStr::from("condCoef"), dv_ar);
        
        pre.thrm_fld_coef[1].set_val(self.rad_coef);
        dummy_el.get_gen_prop_dfd1(&mut pre.thrm_fld_coef[1],&mut CppStr::from("radCoef"), dv_ar);
        
        pre.ref_temp.set_val(self.ref_temp);

        pre.ideal_gas.set_val(self.ideal_gas);
        dummy_el.get_gen_prop_dfd1(&mut pre.ideal_gas, &mut CppStr::from("idealGasConstant"), dv_ar);

        pre.bulk_mod.set_val(self.bulk_mod);
        dummy_el.get_gen_prop_dfd1(&mut pre.bulk_mod, &mut CppStr::from("bulkModulus"), dv_ar);

        pre.therm_exp[0].set_val(self.therm_exp);
        dummy_el.get_gen_prop_dfd1(&mut pre.therm_exp[0], &mut CppStr::from("thermalExp"), dv_ar);

        pre.ref_den.set_val(self.ref_den);

        pre.ref_pres.set_val(self.ref_pres);

        
        for nd in n_sets[self.set_pt1].labels.iter() {
            if self.max_nbrs < MAX_INT {
                for i2 in 0..self.max_nbrs {
                    nearest[i2] = MAX_INT;
                    near_dist[i2] = 1.0e+100;
                }
            }
            nodes[*nd].get_def_crd_dfd1(&mut crd1);
            fcrd1[0] = crd1[0].val;
            fcrd1[1] = crd1[1].val;
            fcrd1[2] = crd1[2].val;
            lst_len = g_list.get_in_radius(g_out, g_out.len(), &fcrd1, self.max_dist);
            //if self.max_nbrs == MAX_INT {
            for nb in 0..lst_len {
                if nd_in_set[g_out[nb]] && g_out[nb] != *nd {
                    nodes[g_out[nb]].get_def_crd_dfd1(&mut crd2);
                    fcrd2[0] = crd2[0].val;
                    fcrd2[1] = crd2[1].val;
                    fcrd2[2] = crd2[2].val;
                    dist = get_dist(&fcrd1, &fcrd2);
                    if dist <= self.max_dist {
                        if self.max_nbrs == MAX_INT {
                            dummy_el.nodes[0] = *nd;
                            dummy_el.nodes[1] = g_out[nb];
                            dummy_el.get_all_nd_var_dfd1(pre, nodes);
                            // if self.ideal_gas > 0.0 {
                            //     self.get_ig_frc_coef_dfd1(pre, &nd_mass[*nd], &nd_mass[g_out[nb]], dist);
                            // }
                            match discipline {
                                0 => dummy_el.put_ru_frc_fld_dfd1(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                                1 => dummy_el.put_rt_frc_fld_dfd1(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                                _ => (),
                            }
                        }
                        else {
                            inserted = false;
                            i1 = 0;
                            while !inserted && i1 < self.max_nbrs {
                                if nearest[i1] == MAX_INT {
                                    nearest[i1] = g_out[nb];
                                    near_dist[i1] = dist;
                                    inserted = true;
                                }
                                else if dist < near_dist[i1] {
                                    for i2 in (i1..(self.max_nbrs - 1)).rev() {
                                        nearest[i2+1] = nearest[i2];
                                        near_dist[i2+1] = near_dist[i2];
                                    }
                                    nearest[i1] = g_out[nb];
                                    near_dist[i1] = dist;
                                    inserted = true;
                                }
                                i1 += 1;
                            }
                        }
                    }
                }
            }
            //}
            if self.max_nbrs < MAX_INT {
                if self.ideal_gas > 0.0 || self.bulk_mod > 0.0 {
                    //part_vol.set_val(self.get_particle_volume(nearest, near_dist));
                    self.get_particle_area_vol_dfd1(&mut part_area, &mut part_vol, nearest, near_dist);
                    pre.glob_temp[0].set_val(nodes[*nd].temperature);
                    pre.mass_per_el.set_val_dfd1(&nd_mass[*nd]);
                    if self.ideal_gas > 0.0 {
                        Interaction::get_ig_pressure_dfd1(&mut pressure, pre, &part_vol);
                    }
                    else {
                        Interaction::get_incomp_pressure_dfd1(&mut pressure, pre, &part_vol);
                    }
                    if self.set_pt1 == self.set_pt2 {
                        tmp.set_val(0.5);
                        pressure.mult(&tmp);
                    }
                }
                for i2 in 0..self.max_nbrs {
                    if nearest[i2] < MAX_INT && near_dist[i2]/near_dist[0] < self.max_ratio {
                        dummy_el.nodes[0] = *nd;
                        dummy_el.nodes[1] = nearest[i2];
                        dummy_el.get_all_nd_var_dfd1(pre, nodes);
                        if self.ideal_gas > 0.0 || self.bulk_mod > 0.0 {
                            Interaction::get_pres_frc_coef_dfd1(pre, &pressure, &part_area, near_dist[i2]);
                        }
                        match discipline {
                            0 => dummy_el.put_ru_frc_fld_dfd1(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                            1 => dummy_el.put_rt_frc_fld_dfd1(glob_r, dr_du, get_matrix, cmd, pre, nodes),
                            _ => (),
                        }
                    }
                }
            }
        }

        for nd in n_sets[self.set_pt2].labels.iter() {
            nd_in_set[*nd] = false;
        }
    }

//end dup
 
//end skip 
 
 
 
 
}

impl InteractionList {

//dup1

    pub fn update_nd_mass_dfd0(&mut self, el_ar : &Vec<Element>, sec_ar : &Vec<Section>, dv_ar : &Vec<DesignVariable>) {
        for nd in self.nd_mass_dfd0.iter_mut() {
            nd.set_val(-1.0);
        }

        let mut el_nd : usize;
        for el in el_ar.iter() {
            if el.this_type == 1 {
                el_nd = el.nodes[0];
                el.get_mass_per_el_dfd0(&mut self.nd_mass_dfd0[el_nd], sec_ar, dv_ar);
            }
        }
    }

    pub fn get_global_r_dfd0(&mut self, glob_r : &mut Vec<DiffDoub0>, dr_du : &mut SparseMat, discipline : usize, pre : &mut DiffDoub0StressPrereq, 
        time : f64, get_matrix : bool, cmd : &JobCommand, n_sets : &Vec<Set>, nodes : &Vec<Node>, dv_ar : &Vec<DesignVariable>) {

        let mut dummy_el = Element::new();
        dummy_el.initialize_type(21);

        for i in self.int_vec.iter() {
            if i.is_active(time) {
                i.get_global_r_dfd0(glob_r, dr_du, discipline, &mut self.nd_in_set, &self.nd_mass_dfd0, &self.interact_grid, &mut self.grid_out, 
                    &mut self.nearest, &mut self.near_dist, &mut dummy_el, pre, time, get_matrix, cmd, n_sets, nodes, dv_ar);
            }
        }
    }

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

    pub fn update_nd_mass_dfd1(&mut self, el_ar : &Vec<Element>, sec_ar : &Vec<Section>, dv_ar : &Vec<DesignVariable>) {
        for nd in self.nd_mass_dfd1.iter_mut() {
            nd.set_val(-1.0);
        }

        let mut el_nd : usize;
        for el in el_ar.iter() {
            if el.this_type == 1 {
                el_nd = el.nodes[0];
                el.get_mass_per_el_dfd1(&mut self.nd_mass_dfd1[el_nd], sec_ar, dv_ar);
            }
        }
    }

    pub fn get_global_r_dfd1(&mut self, glob_r : &mut Vec<DiffDoub1>, dr_du : &mut SparseMat, discipline : usize, pre : &mut DiffDoub1StressPrereq, 
        time : f64, get_matrix : bool, cmd : &JobCommand, n_sets : &Vec<Set>, nodes : &Vec<Node>, dv_ar : &Vec<DesignVariable>) {

        let mut dummy_el = Element::new();
        dummy_el.initialize_type(21);

        for i in self.int_vec.iter() {
            if i.is_active(time) {
                i.get_global_r_dfd1(glob_r, dr_du, discipline, &mut self.nd_in_set, &self.nd_mass_dfd1, &self.interact_grid, &mut self.grid_out, 
                    &mut self.nearest, &mut self.near_dist, &mut dummy_el, pre, time, get_matrix, cmd, n_sets, nodes, dv_ar);
            }
        }
    }

//end dup
 
//end skip 
 
 
 
 

    pub fn initialize(&mut self, nodes : &Vec<Node>, node_sets : &Vec<Set>, ns_map : &CppMap, el_ar : &Vec<Element>, dv_ar : &Vec<DesignVariable>) {

        if self.int_vec.len() > 0 {
            let num_nds = nodes.len();
            self.nd_interaction = vec![false; num_nds];
            self.nd_in_set = vec![false; num_nds];
            self.nd_mass_dfd0 = vec![DiffDoub0::new(); num_nds];
            self.nd_mass_dfd1 = vec![DiffDoub1::new(); num_nds];
            self.grid_out = vec![0usize; num_nds];

            let mut near_len = 0usize;
            for i in self.int_vec.iter_mut() {
                i.set_pt1 = ns_map.at(&i.node_set1.s);
                for nd in node_sets[i.set_pt1].labels.iter() {
                    self.nd_interaction[*nd] = true;
                }
                i.set_pt2 = ns_map.at(&i.node_set2.s);
                for nd in node_sets[i.set_pt2].labels.iter() {
                    self.nd_interaction[*nd] = true;
                }
                if i.max_nbrs < MAX_INT && i.max_nbrs > near_len {
                    near_len = i.max_nbrs;
                }
            }
            self.nearest = vec![MAX_INT; near_len];
            self.near_dist = vec![1.0e+100; near_len];
            
            let mut x_r = [1.0e+100, -1.0e+100];
            let mut y_r = [1.0e+100, -1.0e+100];
            let mut z_r = [1.0e+100, -1.0e+100];
            
            for nd in nodes.iter() {
                if nd.coord[0] < x_r[0] {
                    x_r[0] = nd.coord[0];
                }
                if nd.coord[0] > x_r[1] {
                    x_r[1] = nd.coord[0];
                }
                if nd.coord[1] < y_r[0] {
                    y_r[0] = nd.coord[1];
                }
                if nd.coord[1] > y_r[1] {
                    y_r[1] = nd.coord[1];
                }
                if nd.coord[2] < z_r[0] {
                    z_r[0] = nd.coord[2];
                }
                if nd.coord[2] > z_r[1] {
                    z_r[1] = nd.coord[2];
                }
            }

            let mut tot_dist = 0f64;
            let mut num_hit = 0usize;
            for el in el_ar.iter() {
                for n1 in 0..el.num_nds() {
                    for n2 in 0..el.num_nds() {
                        if n1 != n2 {
                            tot_dist += get_dist(&nodes[el.nodes[n1]].coord, &nodes[el.nodes[n2]].coord);
                            num_hit += 1;
                        }
                    }
                }
            }
            let spacing = 3.0*tot_dist/(num_hit as f64);

            self.interact_grid.initialize(&mut x_r, spacing, &mut y_r, spacing, &mut z_r, spacing, num_nds);

            let num_int = self.int_vec.len();
            let mut inserted : bool;
            let mut i1 : usize;
            let mut dv_i = 0usize;
            for dv in dv_ar.iter() {
                if dv.int_name.s != "" {
                    inserted = false;
                    i1 = 0;
                    while i1 < num_int && !inserted {
                        if self.int_vec[i1].name.s == dv.int_name.s {
                            match dv.coefs.front() {
                                None => self.int_vec[i1].dvars.push_back(IDCapsule {int_dat: dv_i, doub_dat: 1.0}),
                                Some(x) => self.int_vec[i1].dvars.push_back(IDCapsule {int_dat: dv_i, doub_dat: *x}),
                            }
                            inserted = true;
                        }
                        i1 += 1;
                    }
                    if !inserted {
                        println!("Warning: no interaction named {} was found for design variable assignment. The variable will have no effect.", dv.int_name.s);
                    }
                }
                dv_i += 1;
            }
        }
    }

}