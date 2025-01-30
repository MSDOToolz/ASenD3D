#include <cmath>
#include "ObjectiveClass.h"
#include "constants.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

using namespace std;


ObjectiveTerm::ObjectiveTerm() {
	category = "";
	el_set_ptr = max_int;
	nd_set_ptr = max_int;
	q_vec.clear();
	el_vol_vec.clear();
	tgt_vec.clear();
	err_norm_vec.clear();
	active_time[0] = -1.0;
	active_time[1] = 1.0e+100;
	q_len = 0;
}

void ObjectiveTerm::set_active_time(double new_at[]) {
	active_time[0] = new_at[0];
	active_time[1] = new_at[1];
	return;
}

void ObjectiveTerm::allocate_obj(vector<Set>& nd_sets, vector<Set>& el_sets) {
	int i1;
	if (q_len == 0) {
		if (el_set_ptr < max_int) {
			q_len = el_sets[el_set_ptr].labels.size();
		}
		else if (nd_set_ptr < max_int) {
			q_len = nd_sets[nd_set_ptr].labels.size();
		}
		if (q_len > 0) {
			q_vec = vector<double>(q_len);
			if (optr == "powerNorm") {
				tgt_vec = vector<double>(q_len);
			} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
				el_vol_vec = vector<double>(q_len);
				tgt_vec = vector<double>(1);
			}
		}
	}

	for (i1 = 0; i1 < q_len; i1++) {
		q_vec[i1] = 0.0;
		if (optr == "powerNorm") {
			tgt_vec[i1] = 0.0;
		}
		else {
			el_vol_vec[i1] = 0.0;
		}
	}
	if (optr == "volumeIntegral" || optr == "volumeAverage") {
		tgt_vec[0] = 0.0;
	}

	return;
}

void ObjectiveTerm::allocate_obj_grad() {
	int i1;
	if (q_len > 0 && d_qd_u.dim == 0) {
		d_qd_u.set_dim(q_len);
		d_qd_v.set_dim(q_len);
		d_qd_a.set_dim(q_len);
		d_qd_t.set_dim(q_len);
		d_qd_tdot.set_dim(q_len);
		d_qd_d.set_dim(q_len);
		d_vd_d.set_dim(q_len);
		if (optr == "powerNorm") {
			err_norm_vec = vector<double>(q_len);
		}
	}
	d_qd_u.zero_all();
	d_qd_v.zero_all();
	d_qd_a.zero_all();
	d_qd_t.zero_all();
	d_qd_tdot.zero_all();
	d_qd_d.zero_all();
	d_vd_d.zero_all();
	for (i1 = 0; i1 < q_len; i1++) {
		err_norm_vec[i1] = 0.0;
	}
	return;
}

double ObjectiveTerm::get_power_norm() {
	int i1;
	double q_err;
	double p_sum = 0.0;
	for (i1 = 0; i1 < q_len; i1++) {
		q_err = q_vec[i1] - tgt_vec[i1];
		p_sum += pow(q_err, expnt);
	}
	return coef * p_sum;
}

void ObjectiveTerm::d_power_normd_u(vector<double>& d_ld_u, vector<double>& d_ld_v, vector<double>& d_ld_a, vector<double>& d_ld_t, vector<double>& d_ld_tdot) {
	d_qd_u.vector_multiply(d_ld_u, err_norm_vec, true);
	d_qd_v.vector_multiply(d_ld_v, err_norm_vec, true);
	d_qd_a.vector_multiply(d_ld_a, err_norm_vec, true);
	d_qd_t.vector_multiply(d_ld_t, err_norm_vec, true);
	d_qd_tdot.vector_multiply(d_ld_tdot, err_norm_vec, true);

	return;
}

void ObjectiveTerm::d_power_normd_d(vector<double>& d_ld_d) {
	d_qd_d.vector_multiply(d_ld_d, err_norm_vec, true);

	return;
}

double ObjectiveTerm::get_vol_integral() {
	int i1;
	double v_int = 0.0;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
	}
	v_int -= tgt_vec[0];
	return coef * pow(v_int, expnt);
}

void ObjectiveTerm::d_vol_integrald_u(vector<double>& d_ld_u, vector<double>& d_ld_v, vector<double>& d_ld_a, vector<double>& d_ld_t, vector<double>& d_ld_tdot) {
	int i1;
	double v_int = 0.0;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
	}
	v_int -= tgt_vec[0];
	v_int = coef*expnt*pow(v_int, expnt - 1.0);
	
	for (i1 = 0; i1 < q_len; i1++) {
		 el_vol_vec[i1]*= v_int;
	}

	d_qd_u.vector_multiply(d_ld_u, el_vol_vec, true);
	d_qd_v.vector_multiply(d_ld_v, el_vol_vec, true);
	d_qd_a.vector_multiply(d_ld_a, el_vol_vec, true);
	d_qd_t.vector_multiply(d_ld_t, el_vol_vec, true);
	d_qd_tdot.vector_multiply(d_ld_tdot, el_vol_vec, true);

	v_int = 1.0 / v_int;
	for (i1 = 0; i1 < q_len; i1++) {
		el_vol_vec[i1] *= v_int;
	}

	return;
}

void ObjectiveTerm::d_vol_integrald_d(vector<double>& d_ld_d) {
	int i1;
	double v_int = 0.0;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
	}
	v_int -= tgt_vec[0];
	v_int = coef * expnt * pow(v_int, expnt - 1.0);

	for (i1 = 0; i1 < q_len; i1++) {
		q_vec[i1] *= v_int;
		el_vol_vec[i1] *= v_int;
	}

	d_qd_d.vector_multiply(d_ld_d, el_vol_vec, true);
	d_vd_d.vector_multiply(d_ld_d, q_vec, true);

	v_int = 1.0 / v_int;
	for (i1 = 0; i1 < q_len; i1++) {
		q_vec[i1] *= v_int;
		el_vol_vec[i1] *= v_int;
	}

	return;
}

double ObjectiveTerm::get_vol_average() {
	int i1;
	double v_int = 0.0;
	double tot_vol = 0.0;
	double vol_avg;
	double va_err;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
		tot_vol += el_vol_vec[i1];
	}
	vol_avg = v_int/tot_vol;
	va_err = vol_avg - tgt_vec[0];
	return coef * pow(va_err, expnt);
}

void ObjectiveTerm::d_vol_averaged_u(vector<double>& d_ld_u, vector<double>& d_ld_v, vector<double>& d_ld_a, vector<double>& d_ld_t, vector<double>& d_ld_tdot) {
	int i1;
	double v_int = 0.0;
	double tot_vol = 0.0;
	double vol_avg;
	double va_err;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
		tot_vol += el_vol_vec[i1];
	}
	vol_avg = v_int / tot_vol;
	va_err = vol_avg - tgt_vec[0];
	va_err = coef * expnt * pow(va_err, expnt - 1.0)/tot_vol;

	for (i1 = 0; i1 < q_len; i1++) {
		el_vol_vec[i1] *= va_err;
	}

	d_qd_u.vector_multiply(d_ld_u, el_vol_vec, true);
	d_qd_v.vector_multiply(d_ld_v, el_vol_vec, true);
	d_qd_a.vector_multiply(d_ld_a, el_vol_vec, true);
	d_qd_t.vector_multiply(d_ld_t, el_vol_vec, true);
	d_qd_tdot.vector_multiply(d_ld_tdot, el_vol_vec, true);

	va_err = 1.0 / va_err;
	for (i1 = 0; i1 < q_len; i1++) {
		el_vol_vec[i1] *= va_err;
	}

	return;
}

void ObjectiveTerm::d_vol_averaged_d(vector<double>& d_ld_d) {
	int i1;
	int col;
	double v_int = 0.0;
	double tot_vol = 0.0;
	double vol_avg;
	double va_err;
	for (i1 = 0; i1 < q_len; i1++) {
		v_int += q_vec[i1] * el_vol_vec[i1];
		tot_vol += el_vol_vec[i1];
	}
	vol_avg = v_int / tot_vol;
	va_err = vol_avg - tgt_vec[0];
	va_err = coef * expnt * pow(va_err, expnt - 1.0) / tot_vol;

	for (i1 = 0; i1 < q_len; i1++) {
		el_vol_vec[i1] *= va_err;
		q_vec[i1] *= va_err;
	}

	d_qd_d.vector_multiply(d_ld_d, el_vol_vec, true);
	d_vd_d.vector_multiply(d_ld_d, q_vec, true);

	va_err = 1.0 / va_err;
	for (i1 = 0; i1 < q_len; i1++) {
		el_vol_vec[i1] *= va_err;
		q_vec[i1] *= va_err;
	}
	va_err = 1.0 / va_err;

	va_err = va_err * vol_avg;

	for (i1 = 0; i1 < q_len; i1++) {
		MatrixRow& this_row = d_vd_d.matrix[i1];
		for (auto& me : this_row.row_vec) {
			col = me.col;
			d_ld_d[col] -= va_err * me.value;
		}
	}

	return;
}

void ObjectiveTerm::get_obj_val(double time, bool n_lgeom, vector<Node>& nd_ar, vector<Element>& el_ar, vector<Set>& nd_sets, vector<Set>& el_sets, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre) {
	if (time < active_time[0] || time > active_time[1]) {
		return;
	}

	int i1;
	int fi;
	int num_lay;
	int tgt_len;
	double tgt_val;
	allocate_obj(nd_sets,el_sets);
	int q_ind;
	double nd_data[6];
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	double se_den;
	DiffDoub0 def[9];
	DiffDoub0 frc_mom[9];
	DiffDoub0 flux[3];
	DiffDoub0 t_grad[3];
	DiffDoub0 e_vol;
	DiffDoub0 e_den;
	DiffDoub0 tmp;

	string cat_list = "displacement velocity acceleration temperature tdot";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (nd_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Node Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& ndi : nd_sets[nd_set_ptr].labels) {
			if (category == "displacement") {
				q_vec[q_ind] = nd_ar[ndi].displacement[component - 1];
			} else if (category == "velocity") {
				q_vec[q_ind] = nd_ar[ndi].velocity[component - 1];
			} else if (category == "acceleration") {
				q_vec[q_ind] = nd_ar[ndi].acceleration[component - 1];
			} else if (category == "temperature") {
				q_vec[q_ind] = nd_ar[ndi].temperature;
			} else if (category == "tdot") {
				q_vec[q_ind] = nd_ar[ndi].temp_change_rate;
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			tgt_len = tgt_vals.size();
			if (tgt_len == 0) {
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = 0.0;
				}
			} else if (tgt_len == 1) {
				tgt_val = *tgt_vals.begin();
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = tgt_val;
				}
			} else {
				q_ind = 0;
				for (auto& tv : tgt_vals) {
					tgt_vec[q_ind] = tv;
					q_ind++;
				}
			}
			value+= get_power_norm();
			return;
		} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
			for (i1 = 0; i1 < q_len; i1++) {
				el_vol_vec[i1] = 1.0;
			}
			tgt_len = tgt_vals.size();
			if (tgt_len == 0) {
				tgt_vec[0] = 0.0;
			} else {
				tgt_vec[0] = *tgt_vals.begin();
			}
			if (optr == "volumeIntegral") {
				value+= get_vol_integral();
				return;
			} else {
				value+= get_vol_average();
				return;
			}
		}
	}
	cat_list = "stress strain strainEnergyDen";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			this_el.get_stress_strain_dfd0(stress, strain, this_el.s_cent, layer, n_lgeom, st_pre);
			if (category == "stress") {
				q_vec[q_ind] = stress[component - 1].val;
			} else if (category == "strain") {
				q_vec[q_ind] = strain[component - 1].val;
			} else {
				se_den = 0.0;
				for (i1 = 0; i1 < 6; i1++) {
					se_den += stress[i1].val * strain[i1].val;
				}
				se_den *= 0.5;
				q_vec[q_ind] = se_den;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				this_el.get_volume_dfd0(e_vol, st_pre, layer, sec_ar, dv_ar);
				el_vol_vec[q_ind] = e_vol.val;
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			tgt_len = tgt_vals.size();
			if (tgt_len == 0) {
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = 0.0;
				}
			}
			else if (tgt_len == 1) {
				tgt_val = *tgt_vals.begin();
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = tgt_val;
				}
			}
			else {
				q_ind = 0;
				for (auto& tv : tgt_vals) {
					tgt_vec[q_ind] = tv;
					q_ind++;
				}
			}
			value+= get_power_norm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgt_vals.size() == 0) {
				tgt_vec[0] = 0.0;
			} else {
				tgt_vec[0] = *tgt_vals.begin();
			}
			if (optr == "volumeIntegral") {
				value+= get_vol_integral();
				return;
			} else {
				value+= get_vol_average();
				return;
			}
		}
	}

	cat_list = "sectionDef sectionFrcMom";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			//this_el->get_stress_strain_dfd0(stress, strain, spt, layer, n_lgeom, st_pre);
			this_el.get_def_frc_mom_dfd0(def, frc_mom, this_el.s_cent, n_lgeom, st_pre);
			if (category == "sectionDef") {
				q_vec[q_ind] = def[component - 1].val;
			}
			else if (category == "sectionFrcMom") {
				q_vec[q_ind] = frc_mom[component - 1].val;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				num_lay = sec_ar[this_el.sect_ptr].layers.size();
				if (num_lay > 1) {
					e_vol.set_val(0.0);
					for (i1 = 0; i1 < num_lay; i1++) {
						this_el.get_volume_dfd0(tmp, st_pre, i1, sec_ar, dv_ar);
						e_vol.add(tmp);
					}
				}
				else {
					this_el.get_volume_dfd0(e_vol, st_pre, 0, sec_ar, dv_ar);
				}
				el_vol_vec[q_ind] = e_vol.val;
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			tgt_len = tgt_vals.size();
			if (tgt_len == 0) {
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = 0.0;
				}
			}
			else if (tgt_len == 1) {
				tgt_val = *tgt_vals.begin();
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = tgt_val;
				}
			}
			else {
				q_ind = 0;
				for (auto& tv : tgt_vals) {
					tgt_vec[q_ind] = tv;
					q_ind++;
				}
			}
			value += get_power_norm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgt_vals.size() == 0) {
				tgt_vec[0] = 0.0;
			}
			else {
				tgt_vec[0] = *tgt_vals.begin();
			}
			if (optr == "volumeIntegral") {
				value += get_vol_integral();
				return;
			}
			else {
				value += get_vol_average();
				return;
			}
		}
	}
	
	cat_list = "flux tempGradient";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			this_el.get_flux_tgrad_dfd0(flux, t_grad, this_el.s_cent, layer, st_pre);
			if (category == "flux") {
				q_vec[q_ind] = flux[component - 1].val;
			}
			else if (category == "tempGradient") {
				q_vec[q_ind] = t_grad[component - 1].val;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				this_el.get_volume_dfd0(e_vol, st_pre, layer, sec_ar, dv_ar);
				el_vol_vec[q_ind] = e_vol.val;
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			tgt_len = tgt_vals.size();
			if (tgt_len == 0) {
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = 0.0;
				}
			}
			else if (tgt_len == 1) {
				tgt_val = *tgt_vals.begin();
				for (i1 = 0; i1 < q_len; i1++) {
					tgt_vec[i1] = tgt_val;
				}
			}
			else {
				q_ind = 0;
				for (auto& tv : tgt_vals) {
					tgt_vec[q_ind] = tv;
					q_ind++;
				}
			}
			value += get_power_norm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgt_vals.size() == 0) {
				tgt_vec[0] = 0.0;
			}
			else {
				tgt_vec[0] = *tgt_vals.begin();
			}
			if (optr == "volumeIntegral") {
				value += get_vol_integral();
				return;
			}
			else {
				value += get_vol_average();
				return;
			}
		}
	}

	cat_list = "mass volume";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			this_el.get_volume_dfd0(e_vol,st_pre,layer,sec_ar,dv_ar);
			el_vol_vec[q_ind] = e_vol.val;
			if (category == "volume") {
				q_vec[q_ind] = 1.0;
			} else {
				this_el.get_density_dfd0(e_den, layer, sec_ar, mat_ar, dv_ar);
				q_vec[q_ind] = e_den.val;
			}
			q_ind++;
		}
		if (tgt_vals.size() == 0) {
			tgt_vec[0] = 0.0;
		} else {
			tgt_vec[0] = *tgt_vals.begin();
		}
		value+= get_vol_integral();
		return;
	}

	return;
}

void ObjectiveTerm::getd_ld_u(vector<double>& d_ld_u, vector<double>& d_ld_v, vector<double>& d_ld_a, vector<double>& d_ld_t, vector<double>& d_ld_tdot, double time, bool n_lgeom, vector<Node>& nd_ar, vector<Element>& el_ar, vector<Set>& nd_sets, vector<Set>& el_sets, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre) {
	if (time < active_time[0] || time > active_time[1]) {
		return;
	}

	int i1;
	int i2;
	int i3;
	int fi;
	allocate_obj_grad();
	int q_ind;
	int dof_ind;
	int curr_rank;
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	vector<DiffDoub0> dsd_u(288);
	vector<DiffDoub0> ded_u(288);
	vector<DiffDoub0> dsd_t(90);
	DiffDoub0 dse_dend_u[33];
	DiffDoub0 dse_dend_t[10];
	DiffDoub0 def[9];
	DiffDoub0 frc_mom[9];
	DiffDoub0 flux[3];
	DiffDoub0 t_grad[3];
	vector<DiffDoub0> d_fd_t(30);
	vector<DiffDoub0> d_tgd_t(30);
	int el_num_nds;
	int el_dof_per_nd;
	int el_num_int_dof;
	int el_tot_dof;
	DiffDoub0 e_vol;
	DiffDoub0 e_den;

	string cat_list = "displacement velocity acceleration temperature tdot";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (nd_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Node Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Node Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& ndi : nd_sets[nd_set_ptr].labels) {
			dof_ind = nd_ar[ndi].dof_index[component - 1];
			curr_rank = nd_ar[ndi].sorted_rank;
			if (category == "displacement") {
				d_qd_u.add_entry(q_ind, dof_ind, 1.0);
			}
			else if (category == "velocity") {
				d_qd_v.add_entry(q_ind, dof_ind, 1.0);
			}
			else if (category == "acceleration") {
				d_qd_a.add_entry(q_ind, dof_ind, 1.0);
			}
			else if (category == "temperature") {
				d_qd_t.add_entry(q_ind, curr_rank, 1.0);
			}
			else if (category == "tdot") {
				d_qd_tdot.add_entry(q_ind, curr_rank, 1.0);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
			return;
		}
		else if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
			else {
				d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
		}
	}

	cat_list = "stress strain strainEnergyDen";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			this_el.get_stress_strain_dfd0(stress, strain, this_el.s_cent, layer, n_lgeom, st_pre);
			this_el.d_stress_straind_u_dfd0(dsd_u, ded_u, dsd_t, this_el.s_cent, layer, n_lgeom, st_pre);
			el_num_nds = this_el.num_nds;
			el_dof_per_nd = this_el.dof_per_nd;
			el_num_int_dof = this_el.num_int_dof;
			el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
			i1 = el_tot_dof * (component - 1);
			i2 = el_num_nds * (component - 1);
			if (category == "stress") {
				sub_vec_dfd0(st_pre.scr_vec1, dsd_u, i1, i1 + el_tot_dof);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_u, st_pre.scr_vec1, false, q_ind, nd_ar);
				sub_vec_dfd0(st_pre.scr_vec1, dsd_t, i2, i2 + el_num_nds);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_t, st_pre.scr_vec1, true, q_ind, nd_ar);
			}
			else if (category == "strain") {
				sub_vec_dfd0(st_pre.scr_vec1, ded_u, i1, i1 + el_tot_dof);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_u, st_pre.scr_vec1, false, q_ind, nd_ar);
			}
			else {
				for (i2 = 0; i2 < el_tot_dof; i2++) {
					dse_dend_u[i2].set_val(0.0);
					i3 = i2;
					for (i1 = 0; i1 < 6; i1++) {
						dse_dend_u[i2].val += stress[i1].val * ded_u[i3].val + dsd_u[i3].val * strain[i1].val;
						i3 += el_tot_dof;
					}
					dse_dend_u[i2].val *= 0.5;
				}
				for (i2 = 0; i2 < el_num_nds; i2++) {
					dse_dend_t[i2].set_val(0.0);
					i3 = i2;
					for (i1 = 0; i1 < 6; i1++) {
						dse_dend_t[i2].val += dsd_t[i3].val * strain[i1].val;
						i3 += el_num_nds;
					}
					dse_dend_t[i2].val *= 0.5;
				}
				ar_to_vec_dfd0(dse_dend_u, st_pre.scr_vec1, 0, 33);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_u, st_pre.scr_vec1, false, q_ind, nd_ar);
				ar_to_vec_dfd0(dse_dend_t, st_pre.scr_vec1, 0, 10);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_t, st_pre.scr_vec1, true, q_ind, nd_ar);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
			else {
				d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
		}
	}
	
	cat_list = "sectionDef sectionFrcMom";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			//this_el->get_stress_strain_dfd0(stress, strain, spt, layer, n_lgeom, st_pre);
			//this_el->d_stress_straind_u_dfd0(dsd_u, ded_u, dsd_t, spt, layer, n_lgeom, st_pre);
			this_el.get_def_frc_mom_dfd0(def, frc_mom, this_el.s_cent, n_lgeom, st_pre);
			this_el.d_def_frc_momd_u_dfd0(ded_u, dsd_u, dsd_t, this_el.s_cent, n_lgeom, st_pre);
			el_num_nds = this_el.num_nds;
			el_dof_per_nd = this_el.dof_per_nd;
			el_num_int_dof = this_el.num_int_dof;
			el_tot_dof = el_num_nds * el_dof_per_nd + el_num_int_dof;
			i1 = el_tot_dof * (component - 1);
			i2 = el_num_nds * (component - 1);
			if (category == "sectionFrcMom") {
				sub_vec_dfd0(st_pre.scr_vec1, dsd_u, i1, i1 + el_tot_dof);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_u, st_pre.scr_vec1, false, q_ind, nd_ar);
				sub_vec_dfd0(st_pre.scr_vec1, dsd_t, i2, i2 + el_num_nds);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_t, st_pre.scr_vec1, true, q_ind, nd_ar);
			}
			else if (category == "sectionDef") {
				sub_vec_dfd0(st_pre.scr_vec1, ded_u, i1, i1 + el_tot_dof);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_u, st_pre.scr_vec1, false, q_ind, nd_ar);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
			else {
				d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
		}
	}

	cat_list = "flux tempGradient";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			this_el.get_stress_prereq_dfd0(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
			this_el.get_flux_tgrad_dfd0(flux, t_grad, this_el.s_cent, layer, st_pre);
			this_el.d_flux_tgradd_t_dfd0(d_fd_t, d_tgd_t, this_el.s_cent, layer, st_pre);
			el_num_nds = this_el.num_nds;
			i1 = el_num_nds * (component - 1);
			if (category == "flux") {
				sub_vec_dfd0(st_pre.scr_vec1, d_fd_t, i1, i1 + el_num_nds);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_t, st_pre.scr_vec1, true, q_ind, nd_ar);
			}
			else if (category == "tempGradient") {
				sub_vec_dfd0(st_pre.scr_vec1, d_tgd_t, i1, i1 + el_num_nds);
				this_el.put_vec_to_glob_mat_dfd0(d_qd_t, st_pre.scr_vec1, true, q_ind, nd_ar);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
			else {
				d_vol_averaged_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot);
				return;
			}
		}
	}

	return;
}

void ObjectiveTerm::getd_ld_d(vector<double>& d_ld_d, double time, bool n_lgeom, vector<Node>& nd_ar,  vector<Element>& el_ar, vector<Set>& nd_sets, vector<Set>& el_sets, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub1StressPrereq& st_pre) {
	if (time < active_time[0] || time > active_time[1]) {
		return;
	}

	int i1;
	int num_lay;
	int fi;
	int q_ind;
	DiffDoub0 dv_val;
	DiffDoub1 strain[6];
	DiffDoub1 stress[6];
	double se_den;
	DiffDoub1 def[9];
	DiffDoub1 frc_mom[9];
	DiffDoub1 flux[3];
	DiffDoub1 t_grad[3];
	DiffDoub1 e_vol;
	DiffDoub1 e_den;
	DiffDoub1 tmp;

	string cat_list = "stress strain strainEnergyDen";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			for (auto& dvi : this_el.comp_dvars) {
				DesignVariable& this_dv = dv_ar[dvi];
				this_dv.get_value_dfd0(dv_val);
				this_dv.diff_val.set_val_2(dv_val.val, 1.0);
				this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
				this_el.get_stress_strain_dfd1(stress, strain, this_el.s_cent, layer, n_lgeom, st_pre);
				if (category == "stress") {
					d_qd_d.add_entry(q_ind, dvi, stress[component - 1].dval);
				}
				else if (category == "strain") {
					d_qd_d.add_entry(q_ind, dvi, strain[component - 1].dval);
				}
				else {
					se_den = 0.0;
					for (i1 = 0; i1 < 6; i1++) {
						se_den += stress[i1].val * strain[i1].dval + stress[i1].dval * strain[i1].val;
					}
					se_den *= 0.5;
					d_qd_d.add_entry(q_ind, dvi, se_den);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					this_el.get_volume_dfd1(e_vol, st_pre, layer, sec_ar, dv_ar);
					d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
				}
				this_dv.diff_val.set_val_2(dv_val.val, 0.0);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_d(d_ld_d);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_d(d_ld_d);
				return;
			}
			else {
				d_vol_averaged_d(d_ld_d);
				return;
			}
		}
	}

	cat_list = "sectionDef sectionFrcMom";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			for (auto& dvi : this_el.comp_dvars) {
				DesignVariable& this_dv = dv_ar[dvi];
				this_dv.get_value_dfd0(dv_val);
				this_dv.diff_val.set_val_2(dv_val.val, 1.0);
				this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
				//this_el->get_stress_strain_dfd0(stress, strain, spt, layer, n_lgeom, st_pre);
				this_el.get_def_frc_mom_dfd1(def, frc_mom, this_el.s_cent, n_lgeom, st_pre);
				if (category == "sectionFrcMom") {
					d_qd_d.add_entry(q_ind, dvi, frc_mom[component - 1].dval);
				}
				else if (category == "sectionDef") {
					d_qd_d.add_entry(q_ind, dvi, def[component - 1].dval);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					num_lay = sec_ar[this_el.sect_ptr].layers.size();
					if (num_lay > 1) {
						e_vol.set_val(0.0);
						for (i1 = 0; i1 < num_lay; i1++) {
							this_el.get_volume_dfd1(tmp, st_pre, i1, sec_ar, dv_ar);
							e_vol.add(tmp);
						}
					}
					else {
						this_el.get_volume_dfd1(e_vol, st_pre, 0, sec_ar, dv_ar);
					}
					d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
				}
				this_dv.diff_val.set_val_2(dv_val.val, 0.0);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_d(d_ld_d);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_d(d_ld_d);
				return;
			}
			else {
				d_vol_averaged_d(d_ld_d);
				return;
			}
		}
	}
	
	cat_list = "flux tempGradient";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			for (auto& dvi : this_el.comp_dvars) {
				DesignVariable& this_dv = dv_ar[dvi];
				this_dv.get_value_dfd0(dv_val);
				this_dv.diff_val.set_val_2(dv_val.val, 1.0);
				this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
				this_el.get_flux_tgrad_dfd1(flux, t_grad, this_el.s_cent, layer, st_pre);
				if (category == "flux") {
					d_qd_d.add_entry(q_ind, dvi, flux[component - 1].dval);
				}
				else if (category == "tempGradient") {
					d_qd_d.add_entry(q_ind, dvi, t_grad[component - 1].dval);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					this_el.get_volume_dfd1(e_vol, st_pre, layer, sec_ar, dv_ar);
					d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
				}
				this_dv.diff_val.set_val_2(dv_val.val, 0.0);
			}
			q_ind++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < q_len; i1++) {
				err_norm_vec[i1] = coef * expnt * pow((q_vec[i1] - tgt_vec[i1]), (expnt - 1.0));
			}
			d_power_normd_d(d_ld_d);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				d_vol_integrald_d(d_ld_d);
				return;
			}
			else {
				d_vol_averaged_d(d_ld_d);
				return;
			}
		}
	}

	cat_list = "mass volume";
	fi = cat_list.find(category);
	if (fi > -1) {
		if (el_set_ptr == max_int) {
			string err_str = "Error: Objective terms of category '" + category + "' must have a valid Element Set specified.\n";
			err_str = err_str + "Check the Objective input file to make sure the Element Set name is correct and defined in the Model input file.";
			throw runtime_error(err_str);
		}
		q_ind = 0;
		for (auto& eli : el_sets[el_set_ptr].labels) {
			Element& this_el = el_ar[eli];
			for (auto& dvi : this_el.comp_dvars) {
				DesignVariable& this_dv = dv_ar[dvi];
				this_dv.get_value_dfd0(dv_val);
				this_dv.diff_val.set_val_2(dv_val.val, 1.0);
				this_el.get_stress_prereq_dfd1(st_pre, sec_ar, mat_ar, nd_ar, dv_ar);
				this_el.get_volume_dfd1(e_vol, st_pre, layer, sec_ar, dv_ar);
				d_vd_d.add_entry(q_ind, dvi, e_vol.dval);
				if(category == "mass") {
					this_el.get_density_dfd1(e_den, layer, sec_ar, mat_ar, dv_ar);
					d_qd_d.add_entry(q_ind, dvi, e_den.dval);
				}
				this_dv.diff_val.set_val_2(dv_val.val, 0.0);
			}
			q_ind++;
		}
		d_vol_integrald_d(d_ld_d);
		return;
	}

	return;
}

Objective::Objective() {
	terms.clear();
}

void Objective::clear_values() {
	for (auto& tm : terms) {
		tm.value = 0.0;
	}
	return;
}

void Objective::calculate_terms(double time, bool n_lgeom, vector<Node>& nd_ar,  vector<Element>& el_ar, vector<Set>& nd_sets, vector<Set>& el_sets, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre) {
	for (auto& tm : terms) {
		tm.get_obj_val(time, n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
	}
	return;
}

void Objective::calculated_ld_u(vector<double>& d_ld_u, vector<double>& d_ld_v, vector<double>& d_ld_a, vector<double>& d_ld_t, vector<double>& d_ld_tdot, double time, bool n_lgeom, vector<Node>& nd_ar, vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre) {
	for (auto& tm : terms) {
		tm.getd_ld_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot, time, n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
	}
	return;
}

void Objective::calculated_ld_d(vector<double>& d_ld_d, double time, bool n_lgeom, vector<Node>& nd_ar,  vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, vector<DesignVariable>& dv_ar, DiffDoub1StressPrereq& st_pre) {
	for (auto& tm : terms) {
		tm.getd_ld_d(d_ld_d, time, n_lgeom, nd_ar, el_ar, nd_sets, el_sets, sec_ar, mat_ar, dv_ar, st_pre);
	}
	return;
}