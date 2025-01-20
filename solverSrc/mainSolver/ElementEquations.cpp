#include <string>
#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "JobClass.h"
#include "matrixFunctions.h"

using namespace std;

const double r_2o3 = 0.666666666666666667;
const double r_1o3 = 0.333333333333333333;
const double r_1o6 = 0.166666666666666667;
const double r_pi = 3.14159265358979324;
const double r_pio180 = 0.0174532925199432958;
const double r_1ort3 = 0.577350269189625765;
const int max_int = 2000000000;

void Element::condense_mat(vector<double>& mat, vector<double>& scr1, vector<double>& scr2) {
	int st_row;
	int end_row;
	int end_col;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	
	st_row = num_nds*dof_per_nd;
	end_row = st_row + num_int_dof - 1;
	end_col = num_int_dof - 1;
	q_rfactor(internal_mat, num_int_dof, st_row, end_row, 0, end_col, 0);
	
	for (i1 = 0; i1 < st_row; i1++) {
		scr1[i1] = 0.0;
	}
	
	for (i1 = 0; i1 < st_row; i1++) {
		i3 = st_row;
		i4 = i1*num_int_dof;
		for (i2 = 0; i2 < num_int_dof; i2++) {
			scr1[i3] = internal_mat[i4];
			i3++;
			i4++;
		}
		solveq_rx_eqb(scr2,internal_mat,scr1,num_int_dof,st_row,end_row,0,end_col,0);
		i4 = i1;
		for (i2 = 0; i2 < st_row; i2++) {
			i5 = i2*num_int_dof;
			for (i3 = 0; i3 < num_int_dof; i3++) {
				mat[i4]-= internal_mat[i5]*scr2[i3];
				i5++;
			}
			i4+= (end_row+1);
		}
	}
	
	return;
}

void Element::update_external(vector<double>& ext_vec, int for_soln, vector<Node>& nd_ar, vector<double>& scr1, vector<double>& scr2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int st_row;
	int end_row;
	int end_col;
	int nd;
	int dof;
	int glb_ind;

	if (num_int_dof == 0) {
		return;
	}

    st_row = num_nds*dof_per_nd;
    end_row = st_row + num_int_dof - 1;
	end_col = num_int_dof - 1;
	for (i1 = 0; i1 < st_row; i1++) {
		scr1[i1] = 0.0;
	}
	
	if(for_soln == 1) {
		i2 = st_row;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			scr1[i2] = internal_ru[i1].val;
			i2++;
		}
	} else {
		i2 = st_row;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			scr1[i2] = -internald_ldu[i1];
			i2++;
		}
	}
	
	solveq_rx_eqb(scr2,internal_mat,scr1,num_int_dof,st_row,end_row,0,end_col,0);
	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < st_row; i1++) {
		nd = nodes[dof_table[i4]];
		dof = dof_table[i4+1];
		glb_ind = nd_ar[nd].dof_index[dof];
		for (i2 = 0; i2 < num_int_dof; i2++) {
			ext_vec[glb_ind]+= internal_mat[i3]*scr2[i2];
			i3++;
		}
		i4 += 2;
	}		
	
	return;
}

void Element::update_internal(vector<double>& ext_vec, int for_soln, vector<Node>& nd_ar, vector<double>& scr1, vector<double>& scr2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int st_row;
	int end_row;
	int end_col;
	int nd;
	int dof;
	int glb_ind;
	double x_vec[33];
	double b_vec[33];

	if (num_int_dof == 0) {
		return;
	}

    st_row = num_nds*dof_per_nd;
    end_row = st_row + num_int_dof - 1;
	end_col = num_int_dof - 1;
	for (i1 = 0; i1 < st_row; i1++) {
		scr1[i1] = 0.0;
	}
	
	if(for_soln == 1) {
		i2 = st_row;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			scr1[i2] = -internal_ru[i1].val;
			i2++;
		}
	} else {
		i2 = st_row;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			scr1[i2] = internald_ldu[i1];
			i2++;
		}
	}
	
	for (i1 = 0; i1 < st_row; i1++) {
		nd = nodes[dof_table[2*i1]];
		dof = dof_table[2*i1+1];
		glb_ind = nd_ar[nd].dof_index[dof];
		i3 = i1*num_int_dof;
		i4 = st_row;
		for (i2 = 0; i2 < num_int_dof; i2++) {
			scr1[i4] -= internal_mat[i3]*ext_vec[glb_ind];
			i3++;
			i4++;
		}
	}
	
	solveq_rx_eqb(scr2,internal_mat,scr1,num_int_dof,st_row,end_row,0,end_col,0);
	
	if(for_soln == 1) {
		for (i1 = 0; i1 < num_int_dof; i1++) {
			internal_disp[i1] += scr2[i1];
		}
	} else {
		for (i1 = 0; i1 < num_int_dof; i1++) {
			internal_adj[i1] = scr2[i1];
		}
	}
	
	return;
}

double Element::get_int_adjd_rd_d() {
	int i1;
	double prod = 0.0;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		prod += internal_adj[i1] * internal_ru[i1].dval;
	}
	return prod;
}

//dup1

void Element::get_ruk(vector<DiffDoub0>& rvec, vector<double>& d_rdu, vector<double>& d_rd_t, bool get_matrix, bool n_lgeom, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int tot_dof;
	
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 d_jwt;
	
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	DiffDoub0 thrm_stn[6];
	DiffDoub0 ip_temp;
	DiffDoub0 cte[6];
	DiffDoub0 cten[90];
	DiffDoub0 ux[9];
	DiffDoub0 sec_def[9];
	DiffDoub0 sec_fc_mom[9];
	
	DiffDoub0 tmp;
	DiffDoub0 tmp_c[81];
	DiffDoub0 tmp_gd[60];
	DiffDoub0 tmp_te[6];

	vec_to_ar(tmp_c, pre.cmat, 0, 81);
	vec_to_ar(tmp_gd, pre.glob_disp, 0, 60);
	vec_to_ar(tmp_te, pre.therm_exp, 0, 6);
	
	tot_dof = num_nds*dof_per_nd + num_int_dof;

	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < tot_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < tot_dof; i2++) {
				d_rdu[i3] = 0.0;
				i3++;
			}
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_t[i4] = 0.0;
				i4++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (dof_per_nd == 6) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}
	}
	
	for (i1 = 0; i1 < num_ip; i1++) {
		get_ip_data(n_vec,d_ndx,det_j,pre.loc_nds,&int_pts[3*i1]);
		d_jwt.set_val(det_j);
		tmp.set_val(ip_wt[i1]);
		d_jwt.mult(tmp);
		ip_temp.set_val(0.0);
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(pre.glob_temp[i2]);
			tmp.mult(n_vec[i2]);
			ip_temp.add(tmp);
		}
		if(dof_per_nd == 3) {
			mat_mul(ux,tmp_gd,d_ndx,3,n_dim,3);
			get_solid_strain(strain,ux,d_ndx,pre.loc_ori,max_int,max_int,n_lgeom);
			for (i2 = 0; i2 < 6; i2++) {
				thrm_stn[i2].set_val(pre.therm_exp[i2]);
				thrm_stn[i2].mult(ip_temp);
				strain[i2].sub(thrm_stn[i2]);
				strain[i2].sub(pre.einit[i2]);
			}
			mat_mul(stress,tmp_c,strain,6,6,1);
			if (get_matrix) {
				mat_mul(cte, tmp_c, tmp_te, 6, 6, 1);
				mat_mul(cten, cte, n_vec, 6, 1, num_nds);
			}
			for (i2 = 0; i2 < tot_dof; i2++) {
				get_solid_strain(strain,ux,d_ndx,pre.loc_ori,i2,max_int,n_lgeom);
				i4 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					tmp.set_val(stress[i3]);
					tmp.mult(strain[i3]);
					tmp.mult(d_jwt);
					rvec[i2].add(tmp);
					if(get_matrix) {
						pre.bmat[i4].set_val(strain[i3]);
						i4+= tot_dof;
					}
				}
				if(n_lgeom && get_matrix) {
					i5 = (tot_dof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < tot_dof; i3++) {
						get_solid_strain(strain,ux,d_ndx,pre.loc_ori,i2,i3,n_lgeom);
						for (i4 = 0; i4 < 6; i4++) {
							d_rdu[i5]+= stress[i4].val*strain[i4].val*d_jwt.val;
						}
						d_rdu[i6] = d_rdu[i5];
						i5++;
						i6+= tot_dof;
					}
				}
			}
		} else {
			get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,max_int,max_int);
			mat_mul(sec_fc_mom,tmp_c,sec_def,def_dim,def_dim,1);
			for (i2 = 0; i2 < 6; i2++) {
				tmp.set_val(pre.therm_exp[i2]);
				tmp.mult(ip_temp);
				tmp.add(pre.einit[i2]);
				sec_fc_mom[i2].sub(tmp);
			}
			if (get_matrix) {
				mat_mul(cten, tmp_te, n_vec, 6, 1, num_nds);
			}
			for (i2 = 0; i2 < tot_dof; i2++) {
				get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,i2,max_int);
				i4 = i2;
				for (i3 = 0; i3 < def_dim; i3++) {
					tmp.set_val(sec_fc_mom[i3]);
					tmp.mult(sec_def[i3]);
					tmp.mult(d_jwt);
					rvec[i2].add(tmp);
					if(get_matrix) {
						pre.bmat[i4].set_val(sec_def[i3]);
						i4+= tot_dof;
					}
				}
				if(n_lgeom && get_matrix) {
					i5 = (tot_dof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < tot_dof; i3++) {
						get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,i2,i3);
						for (i4 = 0; i4 < def_dim; i4++) {
							d_rdu[i5]+= sec_fc_mom[i4].val*sec_def[i4].val*d_jwt.val;
						}
						d_rdu[i6] = d_rdu[i5];
						i5++;
						i6+= tot_dof;
					}					
				}
			}
		}
		if(get_matrix) {
			mat_mul(pre.cbmat,pre.cmat,pre.bmat,def_dim,def_dim,tot_dof);
			i3 = def_dim*tot_dof;
			for (i2 = 0; i2 < i3; i2++) {
				pre.cbmat[i2].mult(d_jwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < tot_dof; i2++) {
				for (i3 = 0; i3 < tot_dof; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < def_dim; i4++) {
						d_rdu[i5]+= pre.bmat[i6].val*pre.cbmat[i7].val;
						i6+= tot_dof;
						i7+= tot_dof;
					}
					i5++;
				}
			}
			for (i2 = 0; i2 < 6 * num_nds; i2++) {
				cten[i2].mult(d_jwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < tot_dof; i2++) {
				for (i3 = 0; i3 < num_nds; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < 6; i4++) {
						d_rd_t[i5] -= pre.bmat[i6].val * cten[i7].val;
						i6 += tot_dof;
						i7 += num_nds;
					}
					i5++;
				}
			}
		}		
	}
	
	return;
}

void Element::get_rum(vector<DiffDoub0>& rvec, vector<double>& d_rd_a, bool get_matrix, bool actual_props, bool n_lgeom, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int nd1;
	int dof1;
	int nd2;
	int dof2;
	int nd_dof = num_nds * dof_per_nd;

	DiffDoub0 inst_disp[60];

	double tmp_s[3];
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 d_jwt;

	DiffDoub0 tmp;
	DiffDoub0 tmp61[6];
	DiffDoub0 tmp62[6];
	DiffDoub0 det_jwt;

	DiffDoub0 rtmp[30];

	DiffDoub0 save_m[36];


	i3 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < nd_dof; i2++) {
				d_rd_a[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1) {
		i2 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			rvec[i1].set_val(pre.glob_acc[i1]);
			rvec[i1].mult(pre.mass_per_el);
			if (get_matrix) {
				d_rd_a[i2] = pre.mass_per_el.val;
				i2 += 4;
			}
		}
		return;
	}
	else if (type == 21) {
		return;
	}

	if (dof_per_nd == 6) {
		if (n_lgeom) {
		    get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 1);
		}
		if (!actual_props) {
			for (i1 = 0; i1 < 36; i1++) {
				save_m[i1].set_val(pre.mmat[i1]);
				pre.mmat[i1].set_val(0.0);
			}
			for (i1 = 0; i1 < 36; i1 += 7) {
				pre.mmat[i1].set_val(1.0);
			}
		}
	}
	else {
		if (!actual_props) {
			save_m[0].set_val(pre.mmat[0]);
			pre.mmat[0].set_val(1.0);
		}
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		det_jwt.set_val(ip_wt[i1]);
		det_jwt.mult(det_j);
		if (dof_per_nd == 6) {
			// build matrix [d_u^i/d_u^g] * {n}
			i7 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd1 = dof_table[i7];
				dof1 = dof_table[i7 + 1];
				get_inst_disp(inst_disp, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, n_lgeom, i2, max_int);
				i6 = 0;
				i5 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					pre.bmat[i5].set_val(0.0);
					for (i4 = 0; i4 < num_nds; i4++) {
						tmp.set_val(n_vec[i4]);
						tmp.mult(inst_disp[i6]);
						pre.bmat[i5].add(tmp);
						i6++;
					}
					i5 += nd_dof;
					i6 += (n_dim - num_nds);
				}
				i7 += 2;
			}
			// tmp61 = bmat*{acc}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				i4 = 0;
				tmp61[i2].set_val(0.0);
				for (i3 = 0; i3 < nd_dof; i3++) {
					nd1 = dof_table[i4];
					dof1 = dof_table[i4 + 1];
					i6 = dof1 * num_nds + nd1;
					tmp.set_val(pre.bmat[i5]);
					tmp.mult(pre.glob_acc[i6]);
					tmp61[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			// pre.scr_vec5 = [m][b]{a}, originally tmp62
			ar_to_vec(tmp61, pre.scr_vec4, 0, 6);
			mat_mul(pre.scr_vec5, pre.mmat, pre.scr_vec4, 6, 6, 1);
			mat_mul(pre.scr_vec4, pre.scr_vec5, pre.bmat, 1, 6, nd_dof);
			// update rvec
			for (i2 = 0; i2 < nd_dof; i2++) {
				tmp.set_val(pre.scr_vec4[i2]);
				tmp.mult(det_jwt);
				rvec[i2].add(tmp);
			}
			if (get_matrix) {
				mat_mul(pre.cbmat, pre.mmat, pre.bmat, 6, 6, nd_dof);
				i4 = 6 * nd_dof;
				for (i2 = 0; i2 < i4; i2++) {
					pre.cbmat[i2].mult(det_jwt);
				}
				i5 = 0;
				for (i2 = 0; i2 < nd_dof; i2++) {
					i7 = 0;
					for (i3 = 0; i3 < nd_dof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							d_rd_a[i5] += pre.cbmat[i6].val * pre.bmat[i7].val;
							i6 += nd_dof;
							i7 += nd_dof;
						}
						i5++;
					}
				}
			}
		}
		else {
			ar_to_vec(n_vec, pre.scr_vec4, 0, 11);
			mat_mul(pre.bmat, pre.glob_acc, pre.scr_vec4, 3, num_nds, 1);
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd1 = dof_table[2 * i2];
				dof1 = dof_table[2 * i2 + 1];
				tmp.set_val(n_vec[nd1]);
				tmp.mult(pre.mmat[0]);
				tmp.mult(pre.bmat[dof1]);
				tmp.mult(det_jwt);
				rvec[i2].add(tmp);
				if (get_matrix) {
					for (i3 = 0; i3 < nd_dof; i3++) {
						nd2 = dof_table[2 * i3];
						dof2 = dof_table[2 * i3 + 1];
						if (dof2 == dof1) {
							tmp.set_val(n_vec[nd1]);
							tmp.mult(n_vec[nd2]);
							tmp.mult(pre.mmat[0]);
							tmp.mult(det_jwt);
							i4 = i2 * nd_dof + i3;
							d_rd_a[i4] += tmp.val;
						}
					}
				}
			}
		}
	}

	if (dof_per_nd == 6) {
		if (!actual_props) {
			for (i1 = 0; i1 < 36; i1++) {
				pre.mmat[i1].set_val(save_m[i1]);
			}
		}
	}
	else {
		if (!actual_props) {
			pre.mmat[0].set_val(save_m[0]);
		}
	}

	return;
}

void Element::get_rud(vector<DiffDoub0>& rvec, vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd_dof = num_nds * dof_per_nd;
	DiffDoub0 tmp;
	vector<DiffDoub0>& rtmp = pre.scr_vec3;
	//double d_rtmp[1089];
	vector<double>& d_rtmp = pre.scr_mat4;
	//double d_rd_t[330];
	vector<double>& d_rd_t = pre.scr_mat5;
	bool d_non_zero;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 wt_det_j;
	DiffDoub0 strain[6];
	DiffDoub0 sec_def[9];
	DiffDoub0 ux[9];
	DiffDoub0 bvel[9];
	DiffDoub0 dbvel[9];

	double tmp_s[3];

	i3 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < nd_dof; i2++) {
				d_rd_v[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (cmd.ray_damp_cm > 0.0) {
		tmp.set_val(cmd.ray_damp_cm);
		for (i1 = 0; i1 < 36; i1++) {
			pre.mmat[i1].mult(tmp);
		}
		for (i1 = 0; i1 < nd_dof; i1++) {
			pre.glob_acc[i1].set_val(pre.glob_vel[i1]);
		}
		get_rum(rtmp, d_rtmp, get_matrix, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
		for (i1 = 0; i1 < nd_dof; i1++) {
			rvec[i1].add(rtmp[i1]);
		}
		if (get_matrix) {
			i3 = 0;
			for (i1 = 0; i1 < nd_dof; i1++) {
				for (i2 = 0; i2 < nd_dof; i2++) {
					d_rd_v[i3] += d_rtmp[i3];
					i3++;
				}
			}
		}
	}
	if (cmd.ray_damp_ck > 0.0) {
		tmp.set_val(cmd.ray_damp_ck);
		for (i1 = 0; i1 < 81; i1++) {
			pre.cmat[i1].mult(tmp);
		}
		i3 = 0;
		i4 = 0;
		for (i1 = 0; i1 < dof_per_nd; i1++) {
			for (i2 = 0; i2 < n_dim; i2++) {
				if (i2 >= num_nds) {
					pre.glob_disp[i3].set_val(0.0);
				}
				else {
					pre.glob_disp[i3].set_val(pre.glob_vel[i4]);
					i4++;
				}
				i3++;
			}
		}
		get_ruk(rtmp, d_rtmp, d_rd_t, get_matrix, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
		for (i1 = 0; i1 < nd_dof; i1++) {
			rvec[i1].add(rtmp[i1]);
		}
		if (get_matrix) {
			i3 = 0;
			i4 = 0;
			for (i1 = 0; i1 < nd_dof; i1++) {
				for (i2 = 0; i2 < nd_dof; i2++) {
					d_rd_v[i3] += d_rtmp[i4];
					i3++;
					i4++;
				}
				i3 += num_int_dof;
			}
		}
	}

	// Material damping
	d_non_zero = false;
	i1 = 0;
	while (!d_non_zero && i1 < 36) {
		if (pre.dmat[i1].val > 0.0) {
			d_non_zero = true;
		}
		i1++;
	}

	if (d_non_zero) {
		if (cmd.nonlinear_geom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}
		for (i1 = 0; i1 < num_ip; i1++) {
			i2 = i1 * 3;
			vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
			get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
			wt_det_j.set_val(ip_wt[i1]);
			wt_det_j.mult(det_j);
			if (dof_per_nd == 3) {
				ar_to_vec(d_ndx, pre.scr_vec4, 0, 33);
				mat_mul(pre.scr_vec5, pre.glob_disp, pre.scr_vec4, 3, n_dim, 3);
				vec_to_ar(ux, pre.scr_vec5, 0, 9);
				for (i2 = 0; i2 < nd_dof; i2++) {
					get_solid_strain(strain, ux, d_ndx, pre.loc_ori, i2, max_int, cmd.nonlinear_geom);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						//i4 = i3 * nd_dof + i2;
						pre.bmat[i4].set_val(strain[i3]);
						i4 += nd_dof;
					}
				}
			}
			else {
				for (i2 = 0; i2 < nd_dof; i2++) {
					get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, cmd.nonlinear_geom, i2, max_int);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						pre.bmat[i4].set_val(sec_def[i3]);
						i4 += nd_dof;
					}
				}
			}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				bvel[i2].set_val(0.0);
				i4 = 0;
				for (i3 = 0; i3 < nd_dof; i3++) {
					nd = dof_table[i4];
					dof = dof_table[i4 + 1];
					tmp.set_val(pre.bmat[i5]);
					tmp.mult(pre.glob_vel[dof * num_nds + nd]);
					bvel[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			ar_to_vec(bvel, pre.scr_vec4, 0, 9);
			mat_mul(pre.scr_vec5, pre.dmat, pre.scr_vec4, def_dim, def_dim, 1);
			mat_mul(pre.scr_vec4, pre.scr_vec5, pre.bmat, 1, def_dim, nd_dof);
			for (i2 = 0; i2 < nd_dof; i2++) {
				pre.scr_vec4[i2].mult(wt_det_j);
				rvec[i2].add(pre.scr_vec4[i2]);
			}
			if (get_matrix) {
				mat_mul(pre.cbmat, pre.dmat, pre.bmat, def_dim, def_dim, nd_dof);
				i3 = nd_dof * def_dim;
				for (i2 = 0; i2 < i3; i2++) {
					pre.cbmat[i2].mult(wt_det_j);
				}
				i5 = 0;
				for (i2 = 0; i2 < nd_dof; i2++) {
					for (i3 = 0; i3 < nd_dof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							//i5 = i2 * nd_dof + i3;
							//i6 = i4 * nd_dof + i2;
							//i7 = i4 * nd_dof + i3;
							d_rd_v[i5] += pre.bmat[i6].val * pre.cbmat[i7].val;
							i6 += nd_dof;
							i7 += nd_dof;
						}
						i5++;
					}
				}
			}
		}
	}

	return;
}


void Element::get_ru(vector<DiffDoub0>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	int nd_dof;
	int tot_dof;
	DiffDoub0 temp_acc[30];
	vector<DiffDoub0>& rvec = pre.scr_vec1;
	//double d_rdu[1089];
	vector<double>& d_rdu = pre.scr_mat1;
	//double d_rd_t[330];
	vector<double>& d_rd_t = pre.scr_mat2;
	vector<DiffDoub0>& rtmp = pre.scr_vec2;
	//double d_rtmp[1089];
	vector<double>& d_rtmp = pre.scr_mat3;
	double c1;
	double c2;
	DiffDoub0 tmp;
	
	nd_dof = num_nds*dof_per_nd;
	tot_dof = nd_dof + num_int_dof;

	for (i1 = 0; i1 < tot_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix && i1 < nd_dof) {
			i3 = i1 * tot_dof;
			for (i2 = 0; i2 < tot_dof; i2++) {
				d_rdu[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		get_ru_frc_fld(glob_r, globd_rdu, get_matrix, cmd, pre, nd_ar);
		return;
	}

	if (type != 1) {
		get_ruk(rvec, d_rdu, d_rd_t, get_matrix, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
	}

	if (num_int_dof > 0) {
		i2 = nd_dof;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			internal_ru[i1].set_val(rvec[i2]);
			i2++;
		}
		if (get_matrix) {
			i4 = 0;
			for (i1 = 0; i1 < tot_dof; i1++) {
				i3 = i1 * tot_dof + nd_dof;
				for (i2 = nd_dof; i2 < tot_dof; i2++) {
					internal_mat[i4] = d_rdu[i3];
					i3++;
					i4++;
				}
			}
			//condense matrix
			condense_mat(d_rdu,pre.scr_mat4,pre.scr_mat5);
		}
	}
	
	if(cmd.dynamic) {
		if (get_matrix) {
			if (cmd.lump_mass) {
				i2 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					nd = dof_table[i2];
					dof = dof_table[i2 + 1];
					i3 = dof * num_nds + nd;
					temp_acc[i1].set_val(pre.glob_acc[i3]);
					pre.glob_acc[i3].set_val(1.0);
					i2 += 2;
				}
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				c1 = cmd.time_step;
				c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
				i2 = 0;
				i3 = tot_dof + 1;
				for (i1 = 0; i1 < nd_dof; i1++) {
					temp_acc[i1].mult(rtmp[i1]);
					rvec[i1].add(temp_acc[i1]);
					d_rdu[i2] += c1 * rtmp[i1].val;
					i2 += i3;
				}
			}
			else {
				get_rum(rtmp, d_rtmp, true, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
				c1 = cmd.time_step;
				c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					for (i2 = 0; i2 < nd_dof; i2++) {
						d_rdu[i3] += c1 * d_rtmp[i4];
						i3++;
						i4++;
					}
					i3 += num_int_dof;
				}
			}
			if (type != 1) {
				get_rud(rtmp, d_rtmp, true, cmd, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
				c2 = cmd.time_step * cmd.newmark_gamma * c1;
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					for (i2 = 0; i2 < nd_dof; i2++) {
						d_rdu[i3] += c2 * d_rtmp[i4];
						i3++;
						i4++;
					}
					i3 += num_int_dof;
				}
			}
		}
		else {
			if (cmd.lump_mass) {
				i2 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					nd = dof_table[i2];
					dof = dof_table[i2 + 1];
					i3 = dof * num_nds + nd;
					temp_acc[i1].set_val(pre.glob_acc[i3]);
					pre.glob_acc[i3].set_val(1.0);
					i2 += 2;
				}
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					temp_acc[i1].mult(rtmp[i1]);
					rvec[i1].add(temp_acc[i1]);
				}
			}
			else {
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
			}
			if (type != 1) {
				get_rud(rtmp, d_rtmp, false, cmd, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
			}
		}
	}
	
	i4 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i4]];
		dof = dof_table[i4+1];
		glob_ind = nd_ar[nd].dof_index[dof];
		glob_r[glob_ind].add(rvec[i1]);
		if(get_matrix) {
			i3 = i1*tot_dof;
			i5 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd2 = nodes[dof_table[i5]];
				dof2 = dof_table[i5+1];
				glob_ind2 = nd_ar[nd2].dof_index[dof2];
				globd_rdu.add_entry(glob_ind, glob_ind2, d_rdu[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4+= 2;
	}
	
	return;
}

void Element::get_rtk(vector<DiffDoub0>& rvec, vector<double>& d_rd_t, bool get_matrix, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 d_nt[33];
	DiffDoub0 det_j;
	DiffDoub0 tmp;
	DiffDoub0 grad_t[3];
	DiffDoub0 q_vec[3];
	DiffDoub0 rtmp[10];
	DiffDoub0 d_rtmp[100];
	double tmp_s[3];
	DiffDoub0 tmp_t[10];
	DiffDoub0 tmp_tc[9];
	DiffDoub0 tmp_ar[30];

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_t[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		tmp.set_val(ip_wt[i1]);
		det_j.mult(tmp);
		vec_to_ar(tmp_t, pre.glob_temp, 0, 10);
		mat_mul(grad_t, tmp_t, d_ndx, 1, num_nds, 3);
		vec_to_ar(tmp_tc, pre.tcmat, 0, 9);
		mat_mul(q_vec, tmp_tc, grad_t, 3, 3, 1);
		mat_mul(rtmp, d_ndx, q_vec, num_nds, 3, 1);
		for (i2 = 0; i2 < num_nds; i2++) {
			rtmp[i2].mult(det_j);
			rvec[i2].add(rtmp[i2]);
		}
		if (get_matrix) {
			transpose(d_nt, d_ndx, num_nds,3);
			mat_mul(tmp_ar, tmp_tc, d_nt, 3, 3, num_nds);
			mat_mul(d_rtmp, d_ndx, tmp_ar, num_nds, 3, num_nds);
			i3 = num_nds * num_nds;
			for (i2 = 0; i2 < i3; i2++) {
				d_rd_t[i2] += d_rtmp[i2].val*det_j.val;
			}
		}
	}

	return;
}

void Element::get_rtm(vector<DiffDoub0>& rvec, vector<double>& d_rd_tdot, bool get_matrix, bool actual_props, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 tmp;
	DiffDoub0 pt_tdot;
	DiffDoub0 rtmp[10];
	DiffDoub0 d_rtmp[100];
	DiffDoub0 save_cp;
	double tmp_s[3];

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_tdot[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		return;
	}
	else if (type == 1) {
		tmp.set_val(pre.mass_per_el);
		tmp.mult(pre.spec_heat);
		d_rd_tdot[0] = tmp.val;
		tmp.mult(pre.glob_temp[0]);
		rvec[0].set_val(tmp);
		return;
	}

	if (!actual_props) {
		save_cp.set_val(pre.spec_heat);
		pre.spec_heat.set_val(1.0);
	}
	else if (dof_per_nd == 3) {
		pre.spec_heat.mult(pre.mmat[0]);
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		tmp.set_val(ip_wt[i1]);
		det_j.mult(tmp);
		pt_tdot.set_val(0.0);
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(n_vec[i2]);
			tmp.mult(pre.glob_tdot[i2]);
			pt_tdot.add(tmp);
		}
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(n_vec[i2]);
			tmp.mult(pre.spec_heat);
			tmp.mult(pt_tdot);
			tmp.mult(det_j);
			rvec[i2].add(tmp);
		}
		if (get_matrix) {
			mat_mul(d_rtmp,n_vec,n_vec,num_nds,1,num_nds);
			i3 = num_nds * num_nds;
			for (i2 = 0; i2 < i3; i2++) {
				d_rd_tdot[i2] += d_rtmp[i2].val * pre.spec_heat.val * det_j.val;
			}
		}
	}

	if (!actual_props) {
		pre.spec_heat.set_val(save_cp);
	}
	else if (dof_per_nd == 3) {
		pre.spec_heat.dvd(pre.mmat[0]);
	}

	return;
}

void Element::get_rt(vector<DiffDoub0>& glob_r, SparseMat& globd_rd_t, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int glob_ind1;
	int glob_ind2;
	double c1 = 1.0/(cmd.time_step*cmd.newmark_gamma);
	vector<DiffDoub0>& rvec = pre.scr_vec1;
	vector<double>& d_rd_t = pre.scr_mat1;
	vector<DiffDoub0>& rtmp = pre.scr_vec2;
	vector<double>& d_rtmp = pre.scr_mat2;

	if (type == 21) {
		get_rt_frc_fld(glob_r, globd_rd_t, pre.scr_mat1, pre.scr_mat2, get_matrix, cmd, pre, nd_ar);
		return;
	}

	if (type != 1) {
		get_rtk(rvec, d_rd_t, get_matrix, pre);
	}

	if (cmd.dynamic) {
		get_rtm(rtmp, d_rtmp, get_matrix, true, pre);
		i3 = 0;
		for (i1 = 0; i1 < num_nds; i1++) {
			rvec[i1].add(rtmp[i1]);
			if (get_matrix) {
				for (i2 = 0; i2 < num_nds; i2++) {
					d_rd_t[i3] += c1 * d_rtmp[i3];
					i3++;
				}
			}
		}
	}

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		glob_ind1 = nd_ar[nodes[i1]].sorted_rank;
		glob_r[glob_ind1].add(rvec[i1]);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				glob_ind2 = nd_ar[nodes[i2]].sorted_rank;
				globd_rd_t.add_entry(glob_ind1, glob_ind2, d_rd_t[i3]);
				i3++;
			}
		}
	}

	return;
}

void Element::get_ru_frc_fld(vector<DiffDoub0>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd_dof = 6;
	int tot_dof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	DiffDoub0 rvec[6];
	double d_rd_u[36];
	DiffDoub0 d_vec[3];
	DiffDoub0 dist;
	DiffDoub0 dv_vec[3];
	DiffDoub0 d_dvecd_u[18];
	DiffDoub0 d_distd_u[6];
	DiffDoub0 f_n1[3];
	DiffDoub0 df_n1d_u[18];
	DiffDoub0 dto_p;
	DiffDoub0 dto_p1;
	DiffDoub0 dto_p2;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	// potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_nds[i2]);
		tmp.add(pre.glob_disp[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_nds[i2]);
		tmp.sub(pre.glob_disp[i2]);
		d_vec[i1].set_val(tmp);
	}

	dist.set_val(d_vec[0]);
	dist.sqr();
	tmp.set_val(d_vec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.set_val(d_vec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	d_dvecd_u[0].set_val(-1.0);
	d_dvecd_u[3].set_val(1.0);
	d_dvecd_u[7].set_val(-1.0);
	d_dvecd_u[10].set_val(1.0);
	d_dvecd_u[14].set_val(-1.0);
	d_dvecd_u[17].set_val(1.0);

	mat_mul(d_distd_u, d_vec, d_dvecd_u, 1, 3, 6);
	tmp.set_val(1.0);
	tmp.dvd(dist);
	for (i1 = 0; i1 < 6; i1++) {
		d_distd_u[i1].mult(tmp);
	}

	dto_p.set_val(dist);
	i1 = 1;
	while (i1 < pre.frc_fld_exp[0].val) {
		dto_p.mult(dist);
		i1++;
	}
	dto_p1.set_val(dto_p);
	dto_p1.mult(dist);
	dto_p2.set_val(dto_p1);
	dto_p2.mult(dist);

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(d_vec[i1]);
		f_n1[i1].mult(tmp);
	}

	mat_mul(df_n1d_u, d_vec, d_distd_u, 3, 1, 6);
	tmp.set_val(1.0);
	tmp.add(pre.frc_fld_exp[0]);
	tmp.neg();
	tmp.mult(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p2);
	for (i1 = 0; i1 < 18; i1++) {
		df_n1d_u[i1].mult(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 18; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(d_dvecd_u[i1]);
		df_n1d_u[i1].add(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		rvec[i1].set_val(f_n1[i1]);
		rvec[i1].neg();
		rvec[i1 + 3].set_val(f_n1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		d_rd_u[i1] = -df_n1d_u[i1].val;
		d_rd_u[i1 + 18] = df_n1d_u[i1].val;
	}

	// damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_vel[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_vel[i2]);
		dv_vec[i1].set_val(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[1]);
	tmp.dvd(dto_p);

	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(tmp);
		f_n1[i1].mult(dv_vec[i1]);
	}

	/*mat_mul(df_n1d_u, dv_vec, d_distd_u, 3, 1, 6);

	tmp2.set_val(pre.frc_fld_coef[1]);
	tmp2.mult(pre.frc_fld_exp[1]);
	tmp2.neg();
	tmp2.dvd(dto_p1);

	for (i1 = 0; i1 < 18; i1++) {
		df_n1d_u[i1].mult(tmp2);
	}*/

	tmp2.set_val(-cmd.newmark_gamma / (cmd.time_step * (cmd.newmark_beta - cmd.newmark_gamma)));
	tmp.mult(tmp2);

	for (i1 = 0; i1 < 18; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(d_dvecd_u[i1]);
		//df_n1d_u[i1].add(tmp2);
		df_n1d_u[i1].set_val(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		rvec[i1].sub(f_n1[i1]);
		rvec[i1 + 3].add(f_n1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		d_rd_u[i1] -= df_n1d_u[i1].val;
		d_rd_u[i1 + 18] += df_n1d_u[i1].val;
	}

	i4 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i4]];
		dof = dof_table[i4 + 1];
		glob_ind = nd_ar[nd].dof_index[dof];
		glob_r[glob_ind].add(rvec[i1]);
		if (get_matrix) {
			i3 = i1 * tot_dof;
			i5 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd2 = nodes[dof_table[i5]];
				dof2 = dof_table[i5 + 1];
				glob_ind2 = nd_ar[nd2].dof_index[dof2];
				globd_rdu.add_entry(glob_ind, glob_ind2, d_rd_u[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4 += 2;
	}

	return;
}

void Element::get_rt_frc_fld(vector<DiffDoub0>& glob_r, SparseMat& globd_rd_t, vector<double>& d_rd_u, vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd_dof = 6;
	int tot_dof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	DiffDoub0 rvec[2];
	double d_rd_t[4];
	DiffDoub0 d_vec[3];
	DiffDoub0 dist;
	DiffDoub0 dv_vec[3];
	DiffDoub0 d_dvecd_u[18];
	DiffDoub0 d_distd_u[6];
	DiffDoub0 f_n1[3];
	DiffDoub0 df_n1d_u[18];
	DiffDoub0 df_n1d_v[18];
	DiffDoub0 dto_p;
	DiffDoub0 dto_p1;
	DiffDoub0 dto_p2;
	DiffDoub0 t1to3;
	DiffDoub0 t1to4;
	DiffDoub0 t2to3;
	DiffDoub0 t2to4;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 tmp_mat[18];
	DiffDoub0 tmp_mat2[18];

	// potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_nds[i2]);
		tmp.add(pre.glob_disp[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_nds[i2]);
		tmp.sub(pre.glob_disp[i2]);
		d_vec[i1].set_val(tmp);
	}

	dist.set_val(d_vec[0]);
	dist.sqr();
	tmp.set_val(d_vec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.set_val(d_vec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	dto_p.set_val(dist);
	i1 = 1;
	while (i1 < pre.frc_fld_exp[0].val) {
		dto_p.mult(dist);
		i1++;
	}
	dto_p1.set_val(dto_p);
	dto_p1.mult(dist); // ||d||^(p+1)
	dto_p2.set_val(dto_p1);
	dto_p2.mult(dist); // ||d||^(p+2)

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(d_vec[i1]);
		f_n1[i1].mult(tmp);
	}

	// damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_vel[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_vel[i2]);
		dv_vec[i1].set_val(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[1]);
	tmp.dvd(dto_p); // tmp = c_d/||d||^p

	for (i1 = 0; i1 < 3; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(dv_vec[i1]);
		f_n1[i1].add(tmp2);
	}

	// calculate absolute temps
	tmp.set_val(pre.glob_temp[0]);
	tmp.add(pre.ref_temp);
	t1to4.set_val(tmp);
	t1to4.sqr();
	t1to4.sqr();
	t1to3.set_val(t1to4);
	t1to3.dvd(tmp);

	tmp.set_val(pre.glob_temp[1]);
	tmp.add(pre.ref_temp);
	t2to4.set_val(tmp);
	t2to4.sqr();
	t2to4.sqr();
	t2to3.set_val(t2to4);
	t2to3.dvd(tmp);

	// conduction and radiation terms
	tmp.set_val(pre.glob_temp[0]);
	tmp.sub(pre.glob_temp[1]);
	tmp.mult(pre.thrm_fld_coef[0]);
	tmp.dvd(dist);
	rvec[0].set_val(tmp);

	tmp.set_val(t1to4);
	tmp.sub(t2to4);
	tmp.mult(pre.thrm_fld_coef[1]);
	tmp.dvd(dist);
	tmp.dvd(dist);
	rvec[0].add(tmp);

	rvec[1].set_val(rvec[0]);
	rvec[1].neg();

	// work dissipation term

	tmp.set_val(f_n1[0]);
	tmp.mult(dv_vec[0]);
	tmp2.set_val(f_n1[1]);
	tmp2.mult(dv_vec[1]);
	tmp.add(tmp2);
	tmp2.set_val(f_n1[2]);
	tmp2.mult(dv_vec[2]);
	tmp.add(tmp2);
	tmp2.set_val(0.5);
	tmp.mult(tmp2); // tmp = 0.5*dot(fn1,d_v)

	rvec[0].sub(tmp);
	rvec[1].sub(tmp);
	for (i1 = 0; i1 < 2; i1++) {
		nd = nodes[i1];
		glob_ind = nd_ar[nd].sorted_rank;
		glob_r[glob_ind].add(rvec[i1]);
	}

	if (get_matrix) {
		d_dvecd_u[0].set_val(-1.0);
		d_dvecd_u[3].set_val(1.0);
		d_dvecd_u[7].set_val(-1.0);
		d_dvecd_u[10].set_val(1.0);
		d_dvecd_u[14].set_val(-1.0);
		d_dvecd_u[17].set_val(1.0);

		mat_mul(d_distd_u, d_vec, d_dvecd_u, 1, 3, 6);
		tmp.set_val(1.0);
		tmp.dvd(dist);
		for (i1 = 0; i1 < 6; i1++) {
			d_distd_u[i1].mult(tmp);
		}

		mat_mul(df_n1d_u, d_vec, d_distd_u, 3, 1, 6);
		tmp.set_val(1.0);
		tmp.add(pre.frc_fld_exp[0]);
		tmp.neg();
		tmp.mult(pre.frc_fld_coef[0]);
		tmp.dvd(dto_p2);
		for (i1 = 0; i1 < 18; i1++) {
			df_n1d_u[i1].mult(tmp);
		}

		tmp.set_val(pre.frc_fld_coef[0]);
		tmp.dvd(dto_p1);
		for (i1 = 0; i1 < 18; i1++) {
			tmp2.set_val(tmp);
			tmp2.mult(d_dvecd_u[i1]);
			df_n1d_u[i1].add(tmp2);
		}

		mat_mul(tmp_mat, dv_vec, d_distd_u, 3, 1, 6);
		tmp.set_val(pre.frc_fld_coef[1]);
		tmp.mult(pre.frc_fld_exp[1]);
		tmp.dvd(dto_p1);
		tmp.neg();
		for (i1 = 0; i1 < 18; i1++) {
			tmp2.set_val(tmp);
			tmp2.mult(tmp_mat[i1]);
			df_n1d_u[i1].add(tmp2);
		}

		tmp.set_val(pre.frc_fld_coef[1]);
		tmp.dvd(dto_p);
		for (i1 = 0; i1 < 18; i1++) {
			df_n1d_v[i1].set_val(d_dvecd_u[i1]);
			df_n1d_v[i1].mult(tmp);
		}

		// d_rd_t
		tmp.set_val(pre.thrm_fld_coef[0]);
		tmp.dvd(dist);
		tmp2.set_val(pre.thrm_fld_coef[1]);
		tmp2.dvd(dist);
		tmp2.dvd(dist);
		d_rd_t[0] = tmp.val + tmp2.val * t1to3.val;
		d_rd_t[1] = -tmp.val;
		d_rd_t[2] = -tmp.val;
		d_rd_t[3] = tmp.val + tmp2.val * t2to3.val;

		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[nd].sorted_rank;
			for (i2 = 0; i2 < 2; i2++) {
				nd2 = nodes[i2];
				glob_ind2 = nd_ar[nd2].sorted_rank;
				globd_rd_t.add_entry(glob_ind, glob_ind2, d_rd_t[i3]);
				i3++;
			}
		}
        
		// d_rd_u
		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] = 0.0;
		}

		// conduction term
		tmp.set_val(pre.glob_temp[0]);
		tmp.sub(pre.glob_temp[1]);
		tmp.mult(pre.thrm_fld_coef[0]);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.neg();
		tmp_mat[0].set_val(tmp);
		tmp_mat[1].set_val(tmp);
		tmp_mat[1].neg();

		mat_mul(tmp_mat2, tmp_mat, d_distd_u, 2, 1, 6);
		
		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] += tmp_mat2[i1].val;
		}

		// radiation term
		tmp.set_val(t1to4);
		tmp.sub(t2to4);
		tmp.mult(pre.thrm_fld_coef[1]);
		tmp2.set_val(2.0);
		tmp.mult(tmp2);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.neg();
		tmp_mat[0].set_val(tmp);
		tmp_mat[1].set_val(tmp);
		tmp_mat[1].neg();

		mat_mul(tmp_mat2, tmp_mat, d_distd_u, 2, 1, 6);

		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] += tmp_mat2[i1].val;
		}

		// work dissipation term
		mat_mul(tmp_mat, dv_vec, df_n1d_u, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_u[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

		// d_rd_v

		for (i1 = 0; i1 < 12; i1++) {
			d_rd_v[i1] = 0.0;
		}

		mat_mul(tmp_mat, dv_vec, df_n1d_v, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_v[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

		mat_mul(tmp_mat, f_n1, d_dvecd_u, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_v[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

	}

	return;
}

void Element::get_app_load(vector<DiffDoub0>& app_ld, Load& ld_pt, bool n_lgeom, DiffDoub0StressPrereq& pre, vector<Section>& sec_ar, vector<Face>& fc_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int glob_ind;
	int nd;
	int dof;
	int nd_dof = num_nds*dof_per_nd;
	int num_lay = sec_ar[sect_ptr].layers.size();
	string ld_type = ld_pt.type;
	vector<double>& d_rd_a = pre.scr_mat1;
	int fc_num_nds;
	bool nd_in_face;
	vector<DiffDoub0>& el_app_ld = pre.scr_vec1;
	DiffDoub0 tot_nd_f[6];
	DiffDoub0 inp_mag;
	DiffDoub0 vec_mag;
	DiffDoub0 el_vol;
	DiffDoub0 sc_fact;
	DiffDoub0 el_cent[3];
	DiffDoub0 cent_to_el[3];
	DiffDoub0 ax_to_el[3];
	DiffDoub0 ang_vel2;
	DiffDoub0 nn_inv;
	DiffDoub0 fc_area;
	DiffDoub0 fc_norm[3];
	DiffDoub0 trac[3];
	DiffDoub0 dp;
	DiffDoub0 tmp;

	for (i1 = 0; i1 < nd_dof; i1++) {
		pre.glob_acc[i1].set_val(0.0);
	}

	if (ld_type == "bodyForce") {
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			pre.glob_acc[i3].set_val(ld_pt.load[dof]);
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, false, n_lgeom, pre, nd_ar, dv_ar);
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			dof = dof_table[i2 + 1];
			tot_nd_f[dof].add(el_app_ld[i1]);
			i2 += 2;
		}
		for (i1 = 0; i1 < dof_per_nd; i1++) {
			tmp.set_val(ld_pt.load[i1]);
			tmp.sqr();
			inp_mag.add(tmp);
			tmp.set_val(tot_nd_f[i1]);
			tmp.sqr();
			vec_mag.add(tmp);
		}
		inp_mag.sqt();
		vec_mag.sqt();
		if (type == 41 || type == 3) {
			el_vol.set_val(0.0);
			for (i1 = 0; i1 < num_lay; i1++) {
				get_volume(tmp, pre, i1, sec_ar, dv_ar);
				el_vol.add(tmp);
			}
		}
		else {
			get_volume(el_vol, pre, 0, sec_ar, dv_ar);
		}
		sc_fact.set_val(inp_mag);
		sc_fact.mult(el_vol);
		sc_fact.dvd(vec_mag);
		for (i1 = 0; i1 < nd_dof; i1++) {
			el_app_ld[i1].mult(sc_fact);
		}
	}
	else if (ld_type == "gravitational") {
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			if (dof < 3) {
				pre.glob_acc[i3].set_val(ld_pt.load[dof]);
			}
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, true, n_lgeom, pre, nd_ar, dv_ar);
	}
	else if (ld_type == "centrifugal") {
		ang_vel2.set_val(ld_pt.angular_vel);
		ang_vel2.sqr();
		i3 = 0;
		nn_inv.set_val(1.0 / num_nds);
		dp.set_val(0.0);
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < num_nds; i2++) {
				el_cent[i1].add(pre.glob_nds[i3]);
				i3++;
			}
			el_cent[i1].mult(nn_inv);
			cent_to_el[i1].set_val(el_cent[i1]);
			tmp.set_val(ld_pt.center[i1]);
			cent_to_el[i1].sub(tmp);
			tmp.set_val(ld_pt.axis[i1]);
			tmp.mult(cent_to_el[i1]);
			dp.add(tmp);
		}
		for (i1 = 0; i1 < 3; i1++) {
			ax_to_el[i1].set_val(cent_to_el[i1]);
			tmp.set_val(ld_pt.axis[i1]);
			tmp.mult(dp);
			ax_to_el[i1].sub(tmp);
			ax_to_el[i1].mult(ang_vel2);
		}
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			if (dof < 3) {
				pre.glob_acc[i3].set_val(ax_to_el[dof]);
			}
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, true, n_lgeom, pre, nd_ar, dv_ar);
	}
	else if (ld_type == "surfaceTraction" || ld_type == "surfacePressure") {
		for (auto& fi : faces) {
			Face& this_fc = fc_ar[fi];
			if (this_fc.on_surf) {
				this_fc.get_area_normal(fc_area, fc_norm, nd_ar, dv_ar);
				dp.set_val(0.0);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.set_val(ld_pt.normal_dir[i1]);
					tmp.mult(fc_norm[i1]);
					dp.add(tmp);
				}
				tmp.set_val(r_pio180* ld_pt.norm_tol);
				tmp.cs();
				if (dp.val > tmp.val) {
					if (ld_type == "surfacePressure") {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].set_val(ld_pt.load[0]);
							trac[i1].mult(fc_norm[i1]);
							trac[i1].neg();
						}
					}
					else {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].set_val(ld_pt.load[i1]);
						}
					}
					fc_num_nds = this_fc.num_nds;
					for (i1 = 0; i1 < fc_num_nds; i1++) {
						i4 = this_fc.loc_nodes[i1];
						for (i3 = 0; i3 < 3; i3++) {
							//i4 = i3 * num_nds + fc_loc_nd[i1];
							pre.glob_acc[i4].set_val(trac[i3]);
							i4 += num_nds;
						}
					}
					get_rum(el_app_ld, d_rd_a, false, false, n_lgeom, pre, nd_ar, dv_ar);
					i2 = 0;
					for (i1 = 0; i1 < nd_dof; i1++) {
						nd = dof_table[i2];
						dof = dof_table[i2 + 1];
						nd_in_face = false;
						for (i3 = 0; i3 < fc_num_nds; i3++) {
							if (this_fc.loc_nodes[i3] == nd) {
								nd_in_face = true;
							}
						}
						if (!nd_in_face) {
							el_app_ld[i1].set_val(0.0);
						}
						tot_nd_f[dof].add(el_app_ld[i1]);
						i2 += 2;
					}
					inp_mag.set_val(0.0);
					vec_mag.set_val(0.0);
					for (i1 = 0; i1 < 3; i1++) {
						tmp.set_val(tot_nd_f[i1]);
						tmp.sqr();
						vec_mag.add(tmp);
						tmp.set_val(trac[i1]);
						tmp.sqr();
						inp_mag.add(tmp);
					}
					inp_mag.sqt();
					vec_mag.sqt();
					tmp.set_val(fc_area);
					tmp.mult(inp_mag);
					tmp.dvd(vec_mag);
					for (i1 = 0; i1 < nd_dof; i1++) {
						el_app_ld[i1].mult(tmp);
					}
				}
			}
		}
	}

	i2 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i2]];
		dof = dof_table[i2 + 1];
		glob_ind = nd_ar[nd].dof_index[dof];
		app_ld[glob_ind].add(el_app_ld[i1]);
		i2 += 2;
	}

	return;
}

void Element::get_app_therm_load(vector<DiffDoub0>& app_ld, Load& ld_pt, DiffDoub0StressPrereq& pre, vector<Section>& sec_ar, vector<Face>& fc_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int glob_ind;
	int num_lay = sec_ar[sect_ptr].layers.size();
	string ld_type = ld_pt.type;
	vector<DiffDoub0>& el_app_ld = pre.scr_vec1;
	vector<double>& d_rd_t = pre.scr_mat1;
	int fc_num_nds;
	bool nd_in_face;
	DiffDoub0 fc_area;
	DiffDoub0 fc_norm[3];
	DiffDoub0 tot_hg;
	DiffDoub0 el_vol;
	DiffDoub0 dp;
	DiffDoub0 tmp;

	for (i1 = 0; i1 < num_nds; i1++) {
		pre.glob_tdot[i1].set_val(0.0);
	}

	if (ld_type == "bodyHeatGen") {
		for (i1 = 0; i1 < num_nds; i1++) {
			pre.glob_tdot[i1].set_val(ld_pt.load[0]);
		}
		get_rtm(el_app_ld, d_rd_t,false,false,pre);
		tot_hg.set_val(0.0);
		for (i1 = 0; i1 < num_nds; i1++) {
			tot_hg.add(el_app_ld[i1]);
		}
		if (num_lay > 0) {
			el_vol.set_val(0.0);
			for (i1 = 0; i1 < num_lay; i1++) {
				get_volume(tmp, pre, i1, sec_ar, dv_ar);
				el_vol.add(tmp);
			}
		}
		else {
			get_volume(el_vol, pre, 0, sec_ar, dv_ar);
		}
		tmp.set_val(ld_pt.load[0]);
		tmp.mult(el_vol);
		tmp.dvd(tot_hg);
		for (i1 = 0; i1 < num_nds; i1++) {
			el_app_ld[i1].mult(tmp);
		}
	}
	else if (ld_type == "surfaceFlux") {
		for (auto& fi : faces) {
			Face& this_fc = fc_ar[fi];
			if (this_fc.on_surf) {
				this_fc.get_area_normal(fc_area, fc_norm, nd_ar, dv_ar);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.set_val(ld_pt.normal_dir[i1]);
					tmp.mult(fc_norm[i1]);
					dp.add(tmp);
				}
				tmp.set_val(r_pio180 * ld_pt.norm_tol);
				tmp.cs();
				if (dp.val > tmp.val) {
					fc_num_nds = this_fc.num_nds;
					for (i1 = 0; i1 < fc_num_nds; i1++) {
						i2 = this_fc.loc_nodes[i1];
						pre.glob_tdot[i2].set_val(ld_pt.load[0]);
					}
					get_rtm(el_app_ld, d_rd_t, false, false, pre);
					tot_hg.set_val(0.0);
					for (i1 = 0; i1 < num_nds; i1++) {
						nd_in_face = false;
						for (i2 = 0; i2 < fc_num_nds; i2++) {
							if (this_fc.loc_nodes[i2] == i1) {
								nd_in_face = true;
							}
						}
						if (!nd_in_face) {
							el_app_ld[i1].set_val(0.0);
						}
						tot_hg.add(el_app_ld[i1]);
					}
					tmp.set_val(ld_pt.load[0]);
					tmp.mult(fc_area);
					tmp.dvd(tot_hg);
					for (i1 = 0; i1 < num_nds; i1++) {
						el_app_ld[i1].mult(tmp);
					}
				}
			}
		}
	}

	for (i1 = 0; i1 < num_nds; i1++) {
		glob_ind = nd_ar[nodes[i1]].sorted_rank;
		app_ld[glob_ind].add(el_app_ld[i1]);
	}

	return;
}

//end dup
 
//skip 
 
//diff_doub1 versions: 
//dup1

void Element::get_ruk(vector<DiffDoub1>& rvec, vector<double>& d_rdu, vector<double>& d_rd_t, bool get_matrix, bool n_lgeom, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int tot_dof;
	
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 d_jwt;
	
	DiffDoub1 strain[6];
	DiffDoub1 stress[6];
	DiffDoub1 thrm_stn[6];
	DiffDoub1 ip_temp;
	DiffDoub1 cte[6];
	DiffDoub1 cten[90];
	DiffDoub1 ux[9];
	DiffDoub1 sec_def[9];
	DiffDoub1 sec_fc_mom[9];
	
	DiffDoub1 tmp;
	DiffDoub1 tmp_c[81];
	DiffDoub1 tmp_gd[60];
	DiffDoub1 tmp_te[6];

	vec_to_ar(tmp_c, pre.cmat, 0, 81);
	vec_to_ar(tmp_gd, pre.glob_disp, 0, 60);
	vec_to_ar(tmp_te, pre.therm_exp, 0, 6);
	
	tot_dof = num_nds*dof_per_nd + num_int_dof;

	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < tot_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < tot_dof; i2++) {
				d_rdu[i3] = 0.0;
				i3++;
			}
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_t[i4] = 0.0;
				i4++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (dof_per_nd == 6) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}
	}
	
	for (i1 = 0; i1 < num_ip; i1++) {
		get_ip_data(n_vec,d_ndx,det_j,pre.loc_nds,&int_pts[3*i1]);
		d_jwt.set_val(det_j);
		tmp.set_val(ip_wt[i1]);
		d_jwt.mult(tmp);
		ip_temp.set_val(0.0);
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(pre.glob_temp[i2]);
			tmp.mult(n_vec[i2]);
			ip_temp.add(tmp);
		}
		if(dof_per_nd == 3) {
			mat_mul(ux,tmp_gd,d_ndx,3,n_dim,3);
			get_solid_strain(strain,ux,d_ndx,pre.loc_ori,max_int,max_int,n_lgeom);
			for (i2 = 0; i2 < 6; i2++) {
				thrm_stn[i2].set_val(pre.therm_exp[i2]);
				thrm_stn[i2].mult(ip_temp);
				strain[i2].sub(thrm_stn[i2]);
				strain[i2].sub(pre.einit[i2]);
			}
			mat_mul(stress,tmp_c,strain,6,6,1);
			if (get_matrix) {
				mat_mul(cte, tmp_c, tmp_te, 6, 6, 1);
				mat_mul(cten, cte, n_vec, 6, 1, num_nds);
			}
			for (i2 = 0; i2 < tot_dof; i2++) {
				get_solid_strain(strain,ux,d_ndx,pre.loc_ori,i2,max_int,n_lgeom);
				i4 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					tmp.set_val(stress[i3]);
					tmp.mult(strain[i3]);
					tmp.mult(d_jwt);
					rvec[i2].add(tmp);
					if(get_matrix) {
						pre.bmat[i4].set_val(strain[i3]);
						i4+= tot_dof;
					}
				}
				if(n_lgeom && get_matrix) {
					i5 = (tot_dof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < tot_dof; i3++) {
						get_solid_strain(strain,ux,d_ndx,pre.loc_ori,i2,i3,n_lgeom);
						for (i4 = 0; i4 < 6; i4++) {
							d_rdu[i5]+= stress[i4].val*strain[i4].val*d_jwt.val;
						}
						d_rdu[i6] = d_rdu[i5];
						i5++;
						i6+= tot_dof;
					}
				}
			}
		} else {
			get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,max_int,max_int);
			mat_mul(sec_fc_mom,tmp_c,sec_def,def_dim,def_dim,1);
			for (i2 = 0; i2 < 6; i2++) {
				tmp.set_val(pre.therm_exp[i2]);
				tmp.mult(ip_temp);
				tmp.add(pre.einit[i2]);
				sec_fc_mom[i2].sub(tmp);
			}
			if (get_matrix) {
				mat_mul(cten, tmp_te, n_vec, 6, 1, num_nds);
			}
			for (i2 = 0; i2 < tot_dof; i2++) {
				get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,i2,max_int);
				i4 = i2;
				for (i3 = 0; i3 < def_dim; i3++) {
					tmp.set_val(sec_fc_mom[i3]);
					tmp.mult(sec_def[i3]);
					tmp.mult(d_jwt);
					rvec[i2].add(tmp);
					if(get_matrix) {
						pre.bmat[i4].set_val(sec_def[i3]);
						i4+= tot_dof;
					}
				}
				if(n_lgeom && get_matrix) {
					i5 = (tot_dof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < tot_dof; i3++) {
						get_section_def(sec_def,pre.glob_disp,pre.inst_ori,pre.loc_ori,pre.glob_nds,d_ndx,n_vec,n_lgeom,i2,i3);
						for (i4 = 0; i4 < def_dim; i4++) {
							d_rdu[i5]+= sec_fc_mom[i4].val*sec_def[i4].val*d_jwt.val;
						}
						d_rdu[i6] = d_rdu[i5];
						i5++;
						i6+= tot_dof;
					}					
				}
			}
		}
		if(get_matrix) {
			mat_mul(pre.cbmat,pre.cmat,pre.bmat,def_dim,def_dim,tot_dof);
			i3 = def_dim*tot_dof;
			for (i2 = 0; i2 < i3; i2++) {
				pre.cbmat[i2].mult(d_jwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < tot_dof; i2++) {
				for (i3 = 0; i3 < tot_dof; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < def_dim; i4++) {
						d_rdu[i5]+= pre.bmat[i6].val*pre.cbmat[i7].val;
						i6+= tot_dof;
						i7+= tot_dof;
					}
					i5++;
				}
			}
			for (i2 = 0; i2 < 6 * num_nds; i2++) {
				cten[i2].mult(d_jwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < tot_dof; i2++) {
				for (i3 = 0; i3 < num_nds; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < 6; i4++) {
						d_rd_t[i5] -= pre.bmat[i6].val * cten[i7].val;
						i6 += tot_dof;
						i7 += num_nds;
					}
					i5++;
				}
			}
		}		
	}
	
	return;
}

void Element::get_rum(vector<DiffDoub1>& rvec, vector<double>& d_rd_a, bool get_matrix, bool actual_props, bool n_lgeom, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int nd1;
	int dof1;
	int nd2;
	int dof2;
	int nd_dof = num_nds * dof_per_nd;

	DiffDoub1 inst_disp[60];

	double tmp_s[3];
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 d_jwt;

	DiffDoub1 tmp;
	DiffDoub1 tmp61[6];
	DiffDoub1 tmp62[6];
	DiffDoub1 det_jwt;

	DiffDoub1 rtmp[30];

	DiffDoub1 save_m[36];


	i3 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < nd_dof; i2++) {
				d_rd_a[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1) {
		i2 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			rvec[i1].set_val(pre.glob_acc[i1]);
			rvec[i1].mult(pre.mass_per_el);
			if (get_matrix) {
				d_rd_a[i2] = pre.mass_per_el.val;
				i2 += 4;
			}
		}
		return;
	}
	else if (type == 21) {
		return;
	}

	if (dof_per_nd == 6) {
		if (n_lgeom) {
		    get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 1);
		}
		if (!actual_props) {
			for (i1 = 0; i1 < 36; i1++) {
				save_m[i1].set_val(pre.mmat[i1]);
				pre.mmat[i1].set_val(0.0);
			}
			for (i1 = 0; i1 < 36; i1 += 7) {
				pre.mmat[i1].set_val(1.0);
			}
		}
	}
	else {
		if (!actual_props) {
			save_m[0].set_val(pre.mmat[0]);
			pre.mmat[0].set_val(1.0);
		}
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		det_jwt.set_val(ip_wt[i1]);
		det_jwt.mult(det_j);
		if (dof_per_nd == 6) {
			// build matrix [d_u^i/d_u^g] * {n}
			i7 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd1 = dof_table[i7];
				dof1 = dof_table[i7 + 1];
				get_inst_disp(inst_disp, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, n_lgeom, i2, max_int);
				i6 = 0;
				i5 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					pre.bmat[i5].set_val(0.0);
					for (i4 = 0; i4 < num_nds; i4++) {
						tmp.set_val(n_vec[i4]);
						tmp.mult(inst_disp[i6]);
						pre.bmat[i5].add(tmp);
						i6++;
					}
					i5 += nd_dof;
					i6 += (n_dim - num_nds);
				}
				i7 += 2;
			}
			// tmp61 = bmat*{acc}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				i4 = 0;
				tmp61[i2].set_val(0.0);
				for (i3 = 0; i3 < nd_dof; i3++) {
					nd1 = dof_table[i4];
					dof1 = dof_table[i4 + 1];
					i6 = dof1 * num_nds + nd1;
					tmp.set_val(pre.bmat[i5]);
					tmp.mult(pre.glob_acc[i6]);
					tmp61[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			// pre.scr_vec5 = [m][b]{a}, originally tmp62
			ar_to_vec(tmp61, pre.scr_vec4, 0, 6);
			mat_mul(pre.scr_vec5, pre.mmat, pre.scr_vec4, 6, 6, 1);
			mat_mul(pre.scr_vec4, pre.scr_vec5, pre.bmat, 1, 6, nd_dof);
			// update rvec
			for (i2 = 0; i2 < nd_dof; i2++) {
				tmp.set_val(pre.scr_vec4[i2]);
				tmp.mult(det_jwt);
				rvec[i2].add(tmp);
			}
			if (get_matrix) {
				mat_mul(pre.cbmat, pre.mmat, pre.bmat, 6, 6, nd_dof);
				i4 = 6 * nd_dof;
				for (i2 = 0; i2 < i4; i2++) {
					pre.cbmat[i2].mult(det_jwt);
				}
				i5 = 0;
				for (i2 = 0; i2 < nd_dof; i2++) {
					i7 = 0;
					for (i3 = 0; i3 < nd_dof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							d_rd_a[i5] += pre.cbmat[i6].val * pre.bmat[i7].val;
							i6 += nd_dof;
							i7 += nd_dof;
						}
						i5++;
					}
				}
			}
		}
		else {
			ar_to_vec(n_vec, pre.scr_vec4, 0, 11);
			mat_mul(pre.bmat, pre.glob_acc, pre.scr_vec4, 3, num_nds, 1);
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd1 = dof_table[2 * i2];
				dof1 = dof_table[2 * i2 + 1];
				tmp.set_val(n_vec[nd1]);
				tmp.mult(pre.mmat[0]);
				tmp.mult(pre.bmat[dof1]);
				tmp.mult(det_jwt);
				rvec[i2].add(tmp);
				if (get_matrix) {
					for (i3 = 0; i3 < nd_dof; i3++) {
						nd2 = dof_table[2 * i3];
						dof2 = dof_table[2 * i3 + 1];
						if (dof2 == dof1) {
							tmp.set_val(n_vec[nd1]);
							tmp.mult(n_vec[nd2]);
							tmp.mult(pre.mmat[0]);
							tmp.mult(det_jwt);
							i4 = i2 * nd_dof + i3;
							d_rd_a[i4] += tmp.val;
						}
					}
				}
			}
		}
	}

	if (dof_per_nd == 6) {
		if (!actual_props) {
			for (i1 = 0; i1 < 36; i1++) {
				pre.mmat[i1].set_val(save_m[i1]);
			}
		}
	}
	else {
		if (!actual_props) {
			pre.mmat[0].set_val(save_m[0]);
		}
	}

	return;
}

void Element::get_rud(vector<DiffDoub1>& rvec, vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd_dof = num_nds * dof_per_nd;
	DiffDoub1 tmp;
	vector<DiffDoub1>& rtmp = pre.scr_vec3;
	//double d_rtmp[1089];
	vector<double>& d_rtmp = pre.scr_mat4;
	//double d_rd_t[330];
	vector<double>& d_rd_t = pre.scr_mat5;
	bool d_non_zero;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 wt_det_j;
	DiffDoub1 strain[6];
	DiffDoub1 sec_def[9];
	DiffDoub1 ux[9];
	DiffDoub1 bvel[9];
	DiffDoub1 dbvel[9];

	double tmp_s[3];

	i3 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < nd_dof; i2++) {
				d_rd_v[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (cmd.ray_damp_cm > 0.0) {
		tmp.set_val(cmd.ray_damp_cm);
		for (i1 = 0; i1 < 36; i1++) {
			pre.mmat[i1].mult(tmp);
		}
		for (i1 = 0; i1 < nd_dof; i1++) {
			pre.glob_acc[i1].set_val(pre.glob_vel[i1]);
		}
		get_rum(rtmp, d_rtmp, get_matrix, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
		for (i1 = 0; i1 < nd_dof; i1++) {
			rvec[i1].add(rtmp[i1]);
		}
		if (get_matrix) {
			i3 = 0;
			for (i1 = 0; i1 < nd_dof; i1++) {
				for (i2 = 0; i2 < nd_dof; i2++) {
					d_rd_v[i3] += d_rtmp[i3];
					i3++;
				}
			}
		}
	}
	if (cmd.ray_damp_ck > 0.0) {
		tmp.set_val(cmd.ray_damp_ck);
		for (i1 = 0; i1 < 81; i1++) {
			pre.cmat[i1].mult(tmp);
		}
		i3 = 0;
		i4 = 0;
		for (i1 = 0; i1 < dof_per_nd; i1++) {
			for (i2 = 0; i2 < n_dim; i2++) {
				if (i2 >= num_nds) {
					pre.glob_disp[i3].set_val(0.0);
				}
				else {
					pre.glob_disp[i3].set_val(pre.glob_vel[i4]);
					i4++;
				}
				i3++;
			}
		}
		get_ruk(rtmp, d_rtmp, d_rd_t, get_matrix, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
		for (i1 = 0; i1 < nd_dof; i1++) {
			rvec[i1].add(rtmp[i1]);
		}
		if (get_matrix) {
			i3 = 0;
			i4 = 0;
			for (i1 = 0; i1 < nd_dof; i1++) {
				for (i2 = 0; i2 < nd_dof; i2++) {
					d_rd_v[i3] += d_rtmp[i4];
					i3++;
					i4++;
				}
				i3 += num_int_dof;
			}
		}
	}

	// Material damping
	d_non_zero = false;
	i1 = 0;
	while (!d_non_zero && i1 < 36) {
		if (pre.dmat[i1].val > 0.0) {
			d_non_zero = true;
		}
		i1++;
	}

	if (d_non_zero) {
		if (cmd.nonlinear_geom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}
		for (i1 = 0; i1 < num_ip; i1++) {
			i2 = i1 * 3;
			vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
			get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
			wt_det_j.set_val(ip_wt[i1]);
			wt_det_j.mult(det_j);
			if (dof_per_nd == 3) {
				ar_to_vec(d_ndx, pre.scr_vec4, 0, 33);
				mat_mul(pre.scr_vec5, pre.glob_disp, pre.scr_vec4, 3, n_dim, 3);
				vec_to_ar(ux, pre.scr_vec5, 0, 9);
				for (i2 = 0; i2 < nd_dof; i2++) {
					get_solid_strain(strain, ux, d_ndx, pre.loc_ori, i2, max_int, cmd.nonlinear_geom);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						//i4 = i3 * nd_dof + i2;
						pre.bmat[i4].set_val(strain[i3]);
						i4 += nd_dof;
					}
				}
			}
			else {
				for (i2 = 0; i2 < nd_dof; i2++) {
					get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, cmd.nonlinear_geom, i2, max_int);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						pre.bmat[i4].set_val(sec_def[i3]);
						i4 += nd_dof;
					}
				}
			}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				bvel[i2].set_val(0.0);
				i4 = 0;
				for (i3 = 0; i3 < nd_dof; i3++) {
					nd = dof_table[i4];
					dof = dof_table[i4 + 1];
					tmp.set_val(pre.bmat[i5]);
					tmp.mult(pre.glob_vel[dof * num_nds + nd]);
					bvel[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			ar_to_vec(bvel, pre.scr_vec4, 0, 9);
			mat_mul(pre.scr_vec5, pre.dmat, pre.scr_vec4, def_dim, def_dim, 1);
			mat_mul(pre.scr_vec4, pre.scr_vec5, pre.bmat, 1, def_dim, nd_dof);
			for (i2 = 0; i2 < nd_dof; i2++) {
				pre.scr_vec4[i2].mult(wt_det_j);
				rvec[i2].add(pre.scr_vec4[i2]);
			}
			if (get_matrix) {
				mat_mul(pre.cbmat, pre.dmat, pre.bmat, def_dim, def_dim, nd_dof);
				i3 = nd_dof * def_dim;
				for (i2 = 0; i2 < i3; i2++) {
					pre.cbmat[i2].mult(wt_det_j);
				}
				i5 = 0;
				for (i2 = 0; i2 < nd_dof; i2++) {
					for (i3 = 0; i3 < nd_dof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							//i5 = i2 * nd_dof + i3;
							//i6 = i4 * nd_dof + i2;
							//i7 = i4 * nd_dof + i3;
							d_rd_v[i5] += pre.bmat[i6].val * pre.cbmat[i7].val;
							i6 += nd_dof;
							i7 += nd_dof;
						}
						i5++;
					}
				}
			}
		}
	}

	return;
}


void Element::get_ru(vector<DiffDoub1>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	int nd_dof;
	int tot_dof;
	DiffDoub1 temp_acc[30];
	vector<DiffDoub1>& rvec = pre.scr_vec1;
	//double d_rdu[1089];
	vector<double>& d_rdu = pre.scr_mat1;
	//double d_rd_t[330];
	vector<double>& d_rd_t = pre.scr_mat2;
	vector<DiffDoub1>& rtmp = pre.scr_vec2;
	//double d_rtmp[1089];
	vector<double>& d_rtmp = pre.scr_mat3;
	double c1;
	double c2;
	DiffDoub1 tmp;
	
	nd_dof = num_nds*dof_per_nd;
	tot_dof = nd_dof + num_int_dof;

	for (i1 = 0; i1 < tot_dof; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix && i1 < nd_dof) {
			i3 = i1 * tot_dof;
			for (i2 = 0; i2 < tot_dof; i2++) {
				d_rdu[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		get_ru_frc_fld(glob_r, globd_rdu, get_matrix, cmd, pre, nd_ar);
		return;
	}

	if (type != 1) {
		get_ruk(rvec, d_rdu, d_rd_t, get_matrix, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
	}

	if (num_int_dof > 0) {
		i2 = nd_dof;
		for (i1 = 0; i1 < num_int_dof; i1++) {
			internal_ru[i1].set_val(rvec[i2]);
			i2++;
		}
		if (get_matrix) {
			i4 = 0;
			for (i1 = 0; i1 < tot_dof; i1++) {
				i3 = i1 * tot_dof + nd_dof;
				for (i2 = nd_dof; i2 < tot_dof; i2++) {
					internal_mat[i4] = d_rdu[i3];
					i3++;
					i4++;
				}
			}
			//condense matrix
			condense_mat(d_rdu,pre.scr_mat4,pre.scr_mat5);
		}
	}
	
	if(cmd.dynamic) {
		if (get_matrix) {
			if (cmd.lump_mass) {
				i2 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					nd = dof_table[i2];
					dof = dof_table[i2 + 1];
					i3 = dof * num_nds + nd;
					temp_acc[i1].set_val(pre.glob_acc[i3]);
					pre.glob_acc[i3].set_val(1.0);
					i2 += 2;
				}
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				c1 = cmd.time_step;
				c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
				i2 = 0;
				i3 = tot_dof + 1;
				for (i1 = 0; i1 < nd_dof; i1++) {
					temp_acc[i1].mult(rtmp[i1]);
					rvec[i1].add(temp_acc[i1]);
					d_rdu[i2] += c1 * rtmp[i1].val;
					i2 += i3;
				}
			}
			else {
				get_rum(rtmp, d_rtmp, true, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
				c1 = cmd.time_step;
				c1 = -1.0 / (c1 * c1 * (cmd.newmark_beta - cmd.newmark_gamma));
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					for (i2 = 0; i2 < nd_dof; i2++) {
						d_rdu[i3] += c1 * d_rtmp[i4];
						i3++;
						i4++;
					}
					i3 += num_int_dof;
				}
			}
			if (type != 1) {
				get_rud(rtmp, d_rtmp, true, cmd, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
				c2 = cmd.time_step * cmd.newmark_gamma * c1;
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					for (i2 = 0; i2 < nd_dof; i2++) {
						d_rdu[i3] += c2 * d_rtmp[i4];
						i3++;
						i4++;
					}
					i3 += num_int_dof;
				}
			}
		}
		else {
			if (cmd.lump_mass) {
				i2 = 0;
				for (i1 = 0; i1 < nd_dof; i1++) {
					nd = dof_table[i2];
					dof = dof_table[i2 + 1];
					i3 = dof * num_nds + nd;
					temp_acc[i1].set_val(pre.glob_acc[i3]);
					pre.glob_acc[i3].set_val(1.0);
					i2 += 2;
				}
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					temp_acc[i1].mult(rtmp[i1]);
					rvec[i1].add(temp_acc[i1]);
				}
			}
			else {
				get_rum(rtmp, d_rtmp, false, true, cmd.nonlinear_geom, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
			}
			if (type != 1) {
				get_rud(rtmp, d_rtmp, false, cmd, pre, nd_ar, dv_ar);
				for (i1 = 0; i1 < nd_dof; i1++) {
					rvec[i1].add(rtmp[i1]);
				}
			}
		}
	}
	
	i4 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i4]];
		dof = dof_table[i4+1];
		glob_ind = nd_ar[nd].dof_index[dof];
		glob_r[glob_ind].add(rvec[i1]);
		if(get_matrix) {
			i3 = i1*tot_dof;
			i5 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd2 = nodes[dof_table[i5]];
				dof2 = dof_table[i5+1];
				glob_ind2 = nd_ar[nd2].dof_index[dof2];
				globd_rdu.add_entry(glob_ind, glob_ind2, d_rdu[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4+= 2;
	}
	
	return;
}

void Element::get_rtk(vector<DiffDoub1>& rvec, vector<double>& d_rd_t, bool get_matrix, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 d_nt[33];
	DiffDoub1 det_j;
	DiffDoub1 tmp;
	DiffDoub1 grad_t[3];
	DiffDoub1 q_vec[3];
	DiffDoub1 rtmp[10];
	DiffDoub1 d_rtmp[100];
	double tmp_s[3];
	DiffDoub1 tmp_t[10];
	DiffDoub1 tmp_tc[9];
	DiffDoub1 tmp_ar[30];

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_t[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		tmp.set_val(ip_wt[i1]);
		det_j.mult(tmp);
		vec_to_ar(tmp_t, pre.glob_temp, 0, 10);
		mat_mul(grad_t, tmp_t, d_ndx, 1, num_nds, 3);
		vec_to_ar(tmp_tc, pre.tcmat, 0, 9);
		mat_mul(q_vec, tmp_tc, grad_t, 3, 3, 1);
		mat_mul(rtmp, d_ndx, q_vec, num_nds, 3, 1);
		for (i2 = 0; i2 < num_nds; i2++) {
			rtmp[i2].mult(det_j);
			rvec[i2].add(rtmp[i2]);
		}
		if (get_matrix) {
			transpose(d_nt, d_ndx, num_nds,3);
			mat_mul(tmp_ar, tmp_tc, d_nt, 3, 3, num_nds);
			mat_mul(d_rtmp, d_ndx, tmp_ar, num_nds, 3, num_nds);
			i3 = num_nds * num_nds;
			for (i2 = 0; i2 < i3; i2++) {
				d_rd_t[i2] += d_rtmp[i2].val*det_j.val;
			}
		}
	}

	return;
}

void Element::get_rtm(vector<DiffDoub1>& rvec, vector<double>& d_rd_tdot, bool get_matrix, bool actual_props, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 tmp;
	DiffDoub1 pt_tdot;
	DiffDoub1 rtmp[10];
	DiffDoub1 d_rtmp[100];
	DiffDoub1 save_cp;
	double tmp_s[3];

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		rvec[i1].set_val(0.0);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				d_rd_tdot[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		return;
	}
	else if (type == 1) {
		tmp.set_val(pre.mass_per_el);
		tmp.mult(pre.spec_heat);
		d_rd_tdot[0] = tmp.val;
		tmp.mult(pre.glob_temp[0]);
		rvec[0].set_val(tmp);
		return;
	}

	if (!actual_props) {
		save_cp.set_val(pre.spec_heat);
		pre.spec_heat.set_val(1.0);
	}
	else if (dof_per_nd == 3) {
		pre.spec_heat.mult(pre.mmat[0]);
	}

	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(tmp_s, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, tmp_s);
		tmp.set_val(ip_wt[i1]);
		det_j.mult(tmp);
		pt_tdot.set_val(0.0);
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(n_vec[i2]);
			tmp.mult(pre.glob_tdot[i2]);
			pt_tdot.add(tmp);
		}
		for (i2 = 0; i2 < num_nds; i2++) {
			tmp.set_val(n_vec[i2]);
			tmp.mult(pre.spec_heat);
			tmp.mult(pt_tdot);
			tmp.mult(det_j);
			rvec[i2].add(tmp);
		}
		if (get_matrix) {
			mat_mul(d_rtmp,n_vec,n_vec,num_nds,1,num_nds);
			i3 = num_nds * num_nds;
			for (i2 = 0; i2 < i3; i2++) {
				d_rd_tdot[i2] += d_rtmp[i2].val * pre.spec_heat.val * det_j.val;
			}
		}
	}

	if (!actual_props) {
		pre.spec_heat.set_val(save_cp);
	}
	else if (dof_per_nd == 3) {
		pre.spec_heat.dvd(pre.mmat[0]);
	}

	return;
}

void Element::get_rt(vector<DiffDoub1>& glob_r, SparseMat& globd_rd_t, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int glob_ind1;
	int glob_ind2;
	double c1 = 1.0/(cmd.time_step*cmd.newmark_gamma);
	vector<DiffDoub1>& rvec = pre.scr_vec1;
	vector<double>& d_rd_t = pre.scr_mat1;
	vector<DiffDoub1>& rtmp = pre.scr_vec2;
	vector<double>& d_rtmp = pre.scr_mat2;

	if (type == 21) {
		get_rt_frc_fld(glob_r, globd_rd_t, pre.scr_mat1, pre.scr_mat2, get_matrix, cmd, pre, nd_ar);
		return;
	}

	if (type != 1) {
		get_rtk(rvec, d_rd_t, get_matrix, pre);
	}

	if (cmd.dynamic) {
		get_rtm(rtmp, d_rtmp, get_matrix, true, pre);
		i3 = 0;
		for (i1 = 0; i1 < num_nds; i1++) {
			rvec[i1].add(rtmp[i1]);
			if (get_matrix) {
				for (i2 = 0; i2 < num_nds; i2++) {
					d_rd_t[i3] += c1 * d_rtmp[i3];
					i3++;
				}
			}
		}
	}

	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		glob_ind1 = nd_ar[nodes[i1]].sorted_rank;
		glob_r[glob_ind1].add(rvec[i1]);
		if (get_matrix) {
			for (i2 = 0; i2 < num_nds; i2++) {
				glob_ind2 = nd_ar[nodes[i2]].sorted_rank;
				globd_rd_t.add_entry(glob_ind1, glob_ind2, d_rd_t[i3]);
				i3++;
			}
		}
	}

	return;
}

void Element::get_ru_frc_fld(vector<DiffDoub1>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd_dof = 6;
	int tot_dof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	DiffDoub1 rvec[6];
	double d_rd_u[36];
	DiffDoub1 d_vec[3];
	DiffDoub1 dist;
	DiffDoub1 dv_vec[3];
	DiffDoub1 d_dvecd_u[18];
	DiffDoub1 d_distd_u[6];
	DiffDoub1 f_n1[3];
	DiffDoub1 df_n1d_u[18];
	DiffDoub1 dto_p;
	DiffDoub1 dto_p1;
	DiffDoub1 dto_p2;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	// potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_nds[i2]);
		tmp.add(pre.glob_disp[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_nds[i2]);
		tmp.sub(pre.glob_disp[i2]);
		d_vec[i1].set_val(tmp);
	}

	dist.set_val(d_vec[0]);
	dist.sqr();
	tmp.set_val(d_vec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.set_val(d_vec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	d_dvecd_u[0].set_val(-1.0);
	d_dvecd_u[3].set_val(1.0);
	d_dvecd_u[7].set_val(-1.0);
	d_dvecd_u[10].set_val(1.0);
	d_dvecd_u[14].set_val(-1.0);
	d_dvecd_u[17].set_val(1.0);

	mat_mul(d_distd_u, d_vec, d_dvecd_u, 1, 3, 6);
	tmp.set_val(1.0);
	tmp.dvd(dist);
	for (i1 = 0; i1 < 6; i1++) {
		d_distd_u[i1].mult(tmp);
	}

	dto_p.set_val(dist);
	i1 = 1;
	while (i1 < pre.frc_fld_exp[0].val) {
		dto_p.mult(dist);
		i1++;
	}
	dto_p1.set_val(dto_p);
	dto_p1.mult(dist);
	dto_p2.set_val(dto_p1);
	dto_p2.mult(dist);

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(d_vec[i1]);
		f_n1[i1].mult(tmp);
	}

	mat_mul(df_n1d_u, d_vec, d_distd_u, 3, 1, 6);
	tmp.set_val(1.0);
	tmp.add(pre.frc_fld_exp[0]);
	tmp.neg();
	tmp.mult(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p2);
	for (i1 = 0; i1 < 18; i1++) {
		df_n1d_u[i1].mult(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 18; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(d_dvecd_u[i1]);
		df_n1d_u[i1].add(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		rvec[i1].set_val(f_n1[i1]);
		rvec[i1].neg();
		rvec[i1 + 3].set_val(f_n1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		d_rd_u[i1] = -df_n1d_u[i1].val;
		d_rd_u[i1 + 18] = df_n1d_u[i1].val;
	}

	// damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_vel[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_vel[i2]);
		dv_vec[i1].set_val(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[1]);
	tmp.dvd(dto_p);

	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(tmp);
		f_n1[i1].mult(dv_vec[i1]);
	}

	/*mat_mul(df_n1d_u, dv_vec, d_distd_u, 3, 1, 6);

	tmp2.set_val(pre.frc_fld_coef[1]);
	tmp2.mult(pre.frc_fld_exp[1]);
	tmp2.neg();
	tmp2.dvd(dto_p1);

	for (i1 = 0; i1 < 18; i1++) {
		df_n1d_u[i1].mult(tmp2);
	}*/

	tmp2.set_val(-cmd.newmark_gamma / (cmd.time_step * (cmd.newmark_beta - cmd.newmark_gamma)));
	tmp.mult(tmp2);

	for (i1 = 0; i1 < 18; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(d_dvecd_u[i1]);
		//df_n1d_u[i1].add(tmp2);
		df_n1d_u[i1].set_val(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		rvec[i1].sub(f_n1[i1]);
		rvec[i1 + 3].add(f_n1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		d_rd_u[i1] -= df_n1d_u[i1].val;
		d_rd_u[i1 + 18] += df_n1d_u[i1].val;
	}

	i4 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i4]];
		dof = dof_table[i4 + 1];
		glob_ind = nd_ar[nd].dof_index[dof];
		glob_r[glob_ind].add(rvec[i1]);
		if (get_matrix) {
			i3 = i1 * tot_dof;
			i5 = 0;
			for (i2 = 0; i2 < nd_dof; i2++) {
				nd2 = nodes[dof_table[i5]];
				dof2 = dof_table[i5 + 1];
				glob_ind2 = nd_ar[nd2].dof_index[dof2];
				globd_rdu.add_entry(glob_ind, glob_ind2, d_rd_u[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4 += 2;
	}

	return;
}

void Element::get_rt_frc_fld(vector<DiffDoub1>& glob_r, SparseMat& globd_rd_t, vector<double>& d_rd_u, vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd_dof = 6;
	int tot_dof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int glob_ind;
	int glob_ind2;
	DiffDoub1 rvec[2];
	double d_rd_t[4];
	DiffDoub1 d_vec[3];
	DiffDoub1 dist;
	DiffDoub1 dv_vec[3];
	DiffDoub1 d_dvecd_u[18];
	DiffDoub1 d_distd_u[6];
	DiffDoub1 f_n1[3];
	DiffDoub1 df_n1d_u[18];
	DiffDoub1 df_n1d_v[18];
	DiffDoub1 dto_p;
	DiffDoub1 dto_p1;
	DiffDoub1 dto_p2;
	DiffDoub1 t1to3;
	DiffDoub1 t1to4;
	DiffDoub1 t2to3;
	DiffDoub1 t2to4;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 tmp_mat[18];
	DiffDoub1 tmp_mat2[18];

	// potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_nds[i2]);
		tmp.add(pre.glob_disp[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_nds[i2]);
		tmp.sub(pre.glob_disp[i2]);
		d_vec[i1].set_val(tmp);
	}

	dist.set_val(d_vec[0]);
	dist.sqr();
	tmp.set_val(d_vec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.set_val(d_vec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	dto_p.set_val(dist);
	i1 = 1;
	while (i1 < pre.frc_fld_exp[0].val) {
		dto_p.mult(dist);
		i1++;
	}
	dto_p1.set_val(dto_p);
	dto_p1.mult(dist); // ||d||^(p+1)
	dto_p2.set_val(dto_p1);
	dto_p2.mult(dist); // ||d||^(p+2)

	tmp.set_val(pre.frc_fld_coef[0]);
	tmp.dvd(dto_p1);
	for (i1 = 0; i1 < 3; i1++) {
		f_n1[i1].set_val(d_vec[i1]);
		f_n1[i1].mult(tmp);
	}

	// damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.set_val(pre.glob_vel[i2]);
		i2 -= 1;
		tmp.sub(pre.glob_vel[i2]);
		dv_vec[i1].set_val(tmp);
	}

	tmp.set_val(pre.frc_fld_coef[1]);
	tmp.dvd(dto_p); // tmp = c_d/||d||^p

	for (i1 = 0; i1 < 3; i1++) {
		tmp2.set_val(tmp);
		tmp2.mult(dv_vec[i1]);
		f_n1[i1].add(tmp2);
	}

	// calculate absolute temps
	tmp.set_val(pre.glob_temp[0]);
	tmp.add(pre.ref_temp);
	t1to4.set_val(tmp);
	t1to4.sqr();
	t1to4.sqr();
	t1to3.set_val(t1to4);
	t1to3.dvd(tmp);

	tmp.set_val(pre.glob_temp[1]);
	tmp.add(pre.ref_temp);
	t2to4.set_val(tmp);
	t2to4.sqr();
	t2to4.sqr();
	t2to3.set_val(t2to4);
	t2to3.dvd(tmp);

	// conduction and radiation terms
	tmp.set_val(pre.glob_temp[0]);
	tmp.sub(pre.glob_temp[1]);
	tmp.mult(pre.thrm_fld_coef[0]);
	tmp.dvd(dist);
	rvec[0].set_val(tmp);

	tmp.set_val(t1to4);
	tmp.sub(t2to4);
	tmp.mult(pre.thrm_fld_coef[1]);
	tmp.dvd(dist);
	tmp.dvd(dist);
	rvec[0].add(tmp);

	rvec[1].set_val(rvec[0]);
	rvec[1].neg();

	// work dissipation term

	tmp.set_val(f_n1[0]);
	tmp.mult(dv_vec[0]);
	tmp2.set_val(f_n1[1]);
	tmp2.mult(dv_vec[1]);
	tmp.add(tmp2);
	tmp2.set_val(f_n1[2]);
	tmp2.mult(dv_vec[2]);
	tmp.add(tmp2);
	tmp2.set_val(0.5);
	tmp.mult(tmp2); // tmp = 0.5*dot(fn1,d_v)

	rvec[0].sub(tmp);
	rvec[1].sub(tmp);
	for (i1 = 0; i1 < 2; i1++) {
		nd = nodes[i1];
		glob_ind = nd_ar[nd].sorted_rank;
		glob_r[glob_ind].add(rvec[i1]);
	}

	if (get_matrix) {
		d_dvecd_u[0].set_val(-1.0);
		d_dvecd_u[3].set_val(1.0);
		d_dvecd_u[7].set_val(-1.0);
		d_dvecd_u[10].set_val(1.0);
		d_dvecd_u[14].set_val(-1.0);
		d_dvecd_u[17].set_val(1.0);

		mat_mul(d_distd_u, d_vec, d_dvecd_u, 1, 3, 6);
		tmp.set_val(1.0);
		tmp.dvd(dist);
		for (i1 = 0; i1 < 6; i1++) {
			d_distd_u[i1].mult(tmp);
		}

		mat_mul(df_n1d_u, d_vec, d_distd_u, 3, 1, 6);
		tmp.set_val(1.0);
		tmp.add(pre.frc_fld_exp[0]);
		tmp.neg();
		tmp.mult(pre.frc_fld_coef[0]);
		tmp.dvd(dto_p2);
		for (i1 = 0; i1 < 18; i1++) {
			df_n1d_u[i1].mult(tmp);
		}

		tmp.set_val(pre.frc_fld_coef[0]);
		tmp.dvd(dto_p1);
		for (i1 = 0; i1 < 18; i1++) {
			tmp2.set_val(tmp);
			tmp2.mult(d_dvecd_u[i1]);
			df_n1d_u[i1].add(tmp2);
		}

		mat_mul(tmp_mat, dv_vec, d_distd_u, 3, 1, 6);
		tmp.set_val(pre.frc_fld_coef[1]);
		tmp.mult(pre.frc_fld_exp[1]);
		tmp.dvd(dto_p1);
		tmp.neg();
		for (i1 = 0; i1 < 18; i1++) {
			tmp2.set_val(tmp);
			tmp2.mult(tmp_mat[i1]);
			df_n1d_u[i1].add(tmp2);
		}

		tmp.set_val(pre.frc_fld_coef[1]);
		tmp.dvd(dto_p);
		for (i1 = 0; i1 < 18; i1++) {
			df_n1d_v[i1].set_val(d_dvecd_u[i1]);
			df_n1d_v[i1].mult(tmp);
		}

		// d_rd_t
		tmp.set_val(pre.thrm_fld_coef[0]);
		tmp.dvd(dist);
		tmp2.set_val(pre.thrm_fld_coef[1]);
		tmp2.dvd(dist);
		tmp2.dvd(dist);
		d_rd_t[0] = tmp.val + tmp2.val * t1to3.val;
		d_rd_t[1] = -tmp.val;
		d_rd_t[2] = -tmp.val;
		d_rd_t[3] = tmp.val + tmp2.val * t2to3.val;

		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[nd].sorted_rank;
			for (i2 = 0; i2 < 2; i2++) {
				nd2 = nodes[i2];
				glob_ind2 = nd_ar[nd2].sorted_rank;
				globd_rd_t.add_entry(glob_ind, glob_ind2, d_rd_t[i3]);
				i3++;
			}
		}
        
		// d_rd_u
		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] = 0.0;
		}

		// conduction term
		tmp.set_val(pre.glob_temp[0]);
		tmp.sub(pre.glob_temp[1]);
		tmp.mult(pre.thrm_fld_coef[0]);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.neg();
		tmp_mat[0].set_val(tmp);
		tmp_mat[1].set_val(tmp);
		tmp_mat[1].neg();

		mat_mul(tmp_mat2, tmp_mat, d_distd_u, 2, 1, 6);
		
		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] += tmp_mat2[i1].val;
		}

		// radiation term
		tmp.set_val(t1to4);
		tmp.sub(t2to4);
		tmp.mult(pre.thrm_fld_coef[1]);
		tmp2.set_val(2.0);
		tmp.mult(tmp2);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.dvd(dist);
		tmp.neg();
		tmp_mat[0].set_val(tmp);
		tmp_mat[1].set_val(tmp);
		tmp_mat[1].neg();

		mat_mul(tmp_mat2, tmp_mat, d_distd_u, 2, 1, 6);

		for (i1 = 0; i1 < 12; i1++) {
			d_rd_u[i1] += tmp_mat2[i1].val;
		}

		// work dissipation term
		mat_mul(tmp_mat, dv_vec, df_n1d_u, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_u[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

		// d_rd_v

		for (i1 = 0; i1 < 12; i1++) {
			d_rd_v[i1] = 0.0;
		}

		mat_mul(tmp_mat, dv_vec, df_n1d_v, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_v[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

		mat_mul(tmp_mat, f_n1, d_dvecd_u, 1, 3, 6);
		i3 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				d_rd_v[i3] -= 0.5 * tmp_mat[i2].val;
				i3++;
			}
		}

	}

	return;
}

void Element::get_app_load(vector<DiffDoub1>& app_ld, Load& ld_pt, bool n_lgeom, DiffDoub1StressPrereq& pre, vector<Section>& sec_ar, vector<Face>& fc_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int glob_ind;
	int nd;
	int dof;
	int nd_dof = num_nds*dof_per_nd;
	int num_lay = sec_ar[sect_ptr].layers.size();
	string ld_type = ld_pt.type;
	vector<double>& d_rd_a = pre.scr_mat1;
	int fc_num_nds;
	bool nd_in_face;
	vector<DiffDoub1>& el_app_ld = pre.scr_vec1;
	DiffDoub1 tot_nd_f[6];
	DiffDoub1 inp_mag;
	DiffDoub1 vec_mag;
	DiffDoub1 el_vol;
	DiffDoub1 sc_fact;
	DiffDoub1 el_cent[3];
	DiffDoub1 cent_to_el[3];
	DiffDoub1 ax_to_el[3];
	DiffDoub1 ang_vel2;
	DiffDoub1 nn_inv;
	DiffDoub1 fc_area;
	DiffDoub1 fc_norm[3];
	DiffDoub1 trac[3];
	DiffDoub1 dp;
	DiffDoub1 tmp;

	for (i1 = 0; i1 < nd_dof; i1++) {
		pre.glob_acc[i1].set_val(0.0);
	}

	if (ld_type == "bodyForce") {
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			pre.glob_acc[i3].set_val(ld_pt.load[dof]);
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, false, n_lgeom, pre, nd_ar, dv_ar);
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			dof = dof_table[i2 + 1];
			tot_nd_f[dof].add(el_app_ld[i1]);
			i2 += 2;
		}
		for (i1 = 0; i1 < dof_per_nd; i1++) {
			tmp.set_val(ld_pt.load[i1]);
			tmp.sqr();
			inp_mag.add(tmp);
			tmp.set_val(tot_nd_f[i1]);
			tmp.sqr();
			vec_mag.add(tmp);
		}
		inp_mag.sqt();
		vec_mag.sqt();
		if (type == 41 || type == 3) {
			el_vol.set_val(0.0);
			for (i1 = 0; i1 < num_lay; i1++) {
				get_volume(tmp, pre, i1, sec_ar, dv_ar);
				el_vol.add(tmp);
			}
		}
		else {
			get_volume(el_vol, pre, 0, sec_ar, dv_ar);
		}
		sc_fact.set_val(inp_mag);
		sc_fact.mult(el_vol);
		sc_fact.dvd(vec_mag);
		for (i1 = 0; i1 < nd_dof; i1++) {
			el_app_ld[i1].mult(sc_fact);
		}
	}
	else if (ld_type == "gravitational") {
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			if (dof < 3) {
				pre.glob_acc[i3].set_val(ld_pt.load[dof]);
			}
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, true, n_lgeom, pre, nd_ar, dv_ar);
	}
	else if (ld_type == "centrifugal") {
		ang_vel2.set_val(ld_pt.angular_vel);
		ang_vel2.sqr();
		i3 = 0;
		nn_inv.set_val(1.0 / num_nds);
		dp.set_val(0.0);
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < num_nds; i2++) {
				el_cent[i1].add(pre.glob_nds[i3]);
				i3++;
			}
			el_cent[i1].mult(nn_inv);
			cent_to_el[i1].set_val(el_cent[i1]);
			tmp.set_val(ld_pt.center[i1]);
			cent_to_el[i1].sub(tmp);
			tmp.set_val(ld_pt.axis[i1]);
			tmp.mult(cent_to_el[i1]);
			dp.add(tmp);
		}
		for (i1 = 0; i1 < 3; i1++) {
			ax_to_el[i1].set_val(cent_to_el[i1]);
			tmp.set_val(ld_pt.axis[i1]);
			tmp.mult(dp);
			ax_to_el[i1].sub(tmp);
			ax_to_el[i1].mult(ang_vel2);
		}
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = dof_table[i2];
			dof = dof_table[i2 + 1];
			i3 = dof * num_nds + nd;
			if (dof < 3) {
				pre.glob_acc[i3].set_val(ax_to_el[dof]);
			}
			i2 += 2;
		}
		get_rum(el_app_ld, d_rd_a, false, true, n_lgeom, pre, nd_ar, dv_ar);
	}
	else if (ld_type == "surfaceTraction" || ld_type == "surfacePressure") {
		for (auto& fi : faces) {
			Face& this_fc = fc_ar[fi];
			if (this_fc.on_surf) {
				this_fc.get_area_normal(fc_area, fc_norm, nd_ar, dv_ar);
				dp.set_val(0.0);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.set_val(ld_pt.normal_dir[i1]);
					tmp.mult(fc_norm[i1]);
					dp.add(tmp);
				}
				tmp.set_val(r_pio180* ld_pt.norm_tol);
				tmp.cs();
				if (dp.val > tmp.val) {
					if (ld_type == "surfacePressure") {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].set_val(ld_pt.load[0]);
							trac[i1].mult(fc_norm[i1]);
							trac[i1].neg();
						}
					}
					else {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].set_val(ld_pt.load[i1]);
						}
					}
					fc_num_nds = this_fc.num_nds;
					for (i1 = 0; i1 < fc_num_nds; i1++) {
						i4 = this_fc.loc_nodes[i1];
						for (i3 = 0; i3 < 3; i3++) {
							//i4 = i3 * num_nds + fc_loc_nd[i1];
							pre.glob_acc[i4].set_val(trac[i3]);
							i4 += num_nds;
						}
					}
					get_rum(el_app_ld, d_rd_a, false, false, n_lgeom, pre, nd_ar, dv_ar);
					i2 = 0;
					for (i1 = 0; i1 < nd_dof; i1++) {
						nd = dof_table[i2];
						dof = dof_table[i2 + 1];
						nd_in_face = false;
						for (i3 = 0; i3 < fc_num_nds; i3++) {
							if (this_fc.loc_nodes[i3] == nd) {
								nd_in_face = true;
							}
						}
						if (!nd_in_face) {
							el_app_ld[i1].set_val(0.0);
						}
						tot_nd_f[dof].add(el_app_ld[i1]);
						i2 += 2;
					}
					inp_mag.set_val(0.0);
					vec_mag.set_val(0.0);
					for (i1 = 0; i1 < 3; i1++) {
						tmp.set_val(tot_nd_f[i1]);
						tmp.sqr();
						vec_mag.add(tmp);
						tmp.set_val(trac[i1]);
						tmp.sqr();
						inp_mag.add(tmp);
					}
					inp_mag.sqt();
					vec_mag.sqt();
					tmp.set_val(fc_area);
					tmp.mult(inp_mag);
					tmp.dvd(vec_mag);
					for (i1 = 0; i1 < nd_dof; i1++) {
						el_app_ld[i1].mult(tmp);
					}
				}
			}
		}
	}

	i2 = 0;
	for (i1 = 0; i1 < nd_dof; i1++) {
		nd = nodes[dof_table[i2]];
		dof = dof_table[i2 + 1];
		glob_ind = nd_ar[nd].dof_index[dof];
		app_ld[glob_ind].add(el_app_ld[i1]);
		i2 += 2;
	}

	return;
}

void Element::get_app_therm_load(vector<DiffDoub1>& app_ld, Load& ld_pt, DiffDoub1StressPrereq& pre, vector<Section>& sec_ar, vector<Face>& fc_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int glob_ind;
	int num_lay = sec_ar[sect_ptr].layers.size();
	string ld_type = ld_pt.type;
	vector<DiffDoub1>& el_app_ld = pre.scr_vec1;
	vector<double>& d_rd_t = pre.scr_mat1;
	int fc_num_nds;
	bool nd_in_face;
	DiffDoub1 fc_area;
	DiffDoub1 fc_norm[3];
	DiffDoub1 tot_hg;
	DiffDoub1 el_vol;
	DiffDoub1 dp;
	DiffDoub1 tmp;

	for (i1 = 0; i1 < num_nds; i1++) {
		pre.glob_tdot[i1].set_val(0.0);
	}

	if (ld_type == "bodyHeatGen") {
		for (i1 = 0; i1 < num_nds; i1++) {
			pre.glob_tdot[i1].set_val(ld_pt.load[0]);
		}
		get_rtm(el_app_ld, d_rd_t,false,false,pre);
		tot_hg.set_val(0.0);
		for (i1 = 0; i1 < num_nds; i1++) {
			tot_hg.add(el_app_ld[i1]);
		}
		if (num_lay > 0) {
			el_vol.set_val(0.0);
			for (i1 = 0; i1 < num_lay; i1++) {
				get_volume(tmp, pre, i1, sec_ar, dv_ar);
				el_vol.add(tmp);
			}
		}
		else {
			get_volume(el_vol, pre, 0, sec_ar, dv_ar);
		}
		tmp.set_val(ld_pt.load[0]);
		tmp.mult(el_vol);
		tmp.dvd(tot_hg);
		for (i1 = 0; i1 < num_nds; i1++) {
			el_app_ld[i1].mult(tmp);
		}
	}
	else if (ld_type == "surfaceFlux") {
		for (auto& fi : faces) {
			Face& this_fc = fc_ar[fi];
			if (this_fc.on_surf) {
				this_fc.get_area_normal(fc_area, fc_norm, nd_ar, dv_ar);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.set_val(ld_pt.normal_dir[i1]);
					tmp.mult(fc_norm[i1]);
					dp.add(tmp);
				}
				tmp.set_val(r_pio180 * ld_pt.norm_tol);
				tmp.cs();
				if (dp.val > tmp.val) {
					fc_num_nds = this_fc.num_nds;
					for (i1 = 0; i1 < fc_num_nds; i1++) {
						i2 = this_fc.loc_nodes[i1];
						pre.glob_tdot[i2].set_val(ld_pt.load[0]);
					}
					get_rtm(el_app_ld, d_rd_t, false, false, pre);
					tot_hg.set_val(0.0);
					for (i1 = 0; i1 < num_nds; i1++) {
						nd_in_face = false;
						for (i2 = 0; i2 < fc_num_nds; i2++) {
							if (this_fc.loc_nodes[i2] == i1) {
								nd_in_face = true;
							}
						}
						if (!nd_in_face) {
							el_app_ld[i1].set_val(0.0);
						}
						tot_hg.add(el_app_ld[i1]);
					}
					tmp.set_val(ld_pt.load[0]);
					tmp.mult(fc_area);
					tmp.dvd(tot_hg);
					for (i1 = 0; i1 < num_nds; i1++) {
						el_app_ld[i1].mult(tmp);
					}
				}
			}
		}
	}

	for (i1 = 0; i1 < num_nds; i1++) {
		glob_ind = nd_ar[nodes[i1]].sorted_rank;
		app_ld[glob_ind].add(el_app_ld[i1]);
	}

	return;
}

//end dup
 
//end skip 
 
 
