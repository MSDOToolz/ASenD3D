#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"

using namespace std;

const double r_2o3 = 0.666666666666666667;
const double r_1o3 = 0.333333333333333333;
const double r_1o6 = 0.166666666666666667;
const double r_pi = 3.14159265358979324;
const double r_pio180 = 0.0174532925199432958;
const double r_1ort3 = 0.577350269189625765;
const int max_int = 2000000000;

//dup1
void Element::get_nd_disp(vector<DiffDoub0>& glob_disp, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_disp[i3].set_val(this_nd.displacement[i2]);
			i3+= n_dim;
		}
	}
	
	if(num_int_dof > 0) {
		i2 = 2*num_nds*dof_per_nd;
		for (i3 = 0; i3 < num_int_dof; i3++) {
			i4 = dof_table[i2];
			i5 = dof_table[i2+1];
			i6 = n_dim*i5 + i4;
			glob_disp[i6].set_val(internal_disp[i3]);
			i2+= 2;
		}
	}
	
	return;
}

void Element::get_nd_vel(vector<DiffDoub0>& glob_vel, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_vel[i3].set_val(this_nd.velocity[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_acc(vector<DiffDoub0>& glob_acc, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_acc[i3].set_val(this_nd.acceleration[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_fl_vel(vector<DiffDoub0>& fl_vel, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			fl_vel[i3].set_val(this_nd.fl_vel[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_fl_vdot(vector<DiffDoub0>& fl_vdot, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			fl_vdot[i3].set_val(this_nd.fl_vel_dot[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_temp(vector<DiffDoub0>& glob_temp, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		glob_temp[i1].set_val(this_nd.temperature);
	}
	return;
}

void Element::get_nd_tdot(vector<DiffDoub0>& glob_tdot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		glob_tdot[i1].set_val(this_nd.temp_change_rate);
	}
	return;
}

void Element::get_nd_fl_den(vector<DiffDoub0>& fl_den, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		fl_den[i1].set_val(this_nd.fl_den);
	}

	return;
}

void Element::get_nd_fl_den_dot(vector<DiffDoub0>& fl_den_dot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		fl_den_dot[i1].set_val(this_nd.fl_den_dot);
	}

	return;
}

void Element::get_nd_turb_e(vector<DiffDoub0>& turb_e, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		turb_e[i1].set_val(this_nd.turb_e);
	}

	return;
}

void Element::get_nd_turb_edot(vector<DiffDoub0>& turb_edot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		turb_edot[i1].set_val(this_nd.turb_edot);
	}

	return;
}

void Element::eval_n(DiffDoub0 n_vec[], DiffDoub0 d_nds[], double spt[]) {
	if(type == 4 || type == 400) {
		n_vec[0].set_val(1.0-spt[0]-spt[1]-spt[2]);
		d_nds[0].set_val(-1.0);
		d_nds[1].set_val(-1.0);
		d_nds[2].set_val(-1.0);
        
		n_vec[1].set_val(spt[0]);
		d_nds[3].set_val(1.0);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
        
	    n_vec[2].set_val(spt[1]);
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(1.0);
		d_nds[8].set_val(0.0);
		
        n_vec[3].set_val(spt[2]);
		d_nds[9].set_val(0.0);
		d_nds[10].set_val(0.0);
		d_nds[11].set_val(1.0);
	} else if(type == 6 || type == 600) {
		n_vec[0].set_val(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
		d_nds[0].set_val(-0.5*(1.0-spt[2]));
	    d_nds[1].set_val(-0.5*(1.0-spt[2]));
		d_nds[2].set_val(-0.5*(1.0-spt[0]-spt[1]));
		
	    n_vec[1].set_val(0.5*spt[0]*(1.0-spt[2]));
		d_nds[3].set_val(0.5*(1.0-spt[2]));
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(-0.5*spt[0]);
		
	    n_vec[2].set_val(0.5*spt[1]*(1.0-spt[2]));
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(0.5*(1.0-spt[2]));
		d_nds[8].set_val(-0.5*spt[1]);
		
		n_vec[3].set_val(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
		d_nds[9].set_val(-0.5*(1.0+spt[2]));
		d_nds[10].set_val(-0.5*(1.0+spt[2]));
		d_nds[11].set_val(0.5*(1.0-spt[0]-spt[1]));
		
	    n_vec[4].set_val(0.5*spt[0]*(1.0+spt[2]));
		d_nds[12].set_val(0.5*(1.0+spt[2]));
		d_nds[13].set_val(0.0);
		d_nds[14].set_val(0.5*spt[0]);
		
	    n_vec[5].set_val(0.5*spt[1]*(1.0+spt[2]));
		d_nds[15].set_val(0.0);
		d_nds[16].set_val(0.5*(1.0+spt[2]));
		d_nds[17].set_val(0.5*spt[1]);
	} else if(type == 8 || type == 81 || type == 800) {
		n_vec[0].set_val(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
		d_nds[0].set_val(-0.125*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[1].set_val(-0.125*(1.0-spt[0])*(1.0-spt[2]));
		d_nds[2].set_val(-0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        n_vec[1].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[3].set_val(0.125*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[4].set_val(-0.125*(1.0+spt[0])*(1.0-spt[2]));
		d_nds[5].set_val(-0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        n_vec[2].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[6].set_val(0.125*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[7].set_val(0.125*(1.0+spt[0])*(1.0-spt[2]));
		d_nds[8].set_val(-0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        n_vec[3].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[9].set_val(-0.125*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[10].set_val(0.125*(1.0-spt[0])*(1.0-spt[2]));
		d_nds[11].set_val(-0.125*(1.0-spt[0])*(1.0+spt[1]));
		
        n_vec[4].set_val(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[12].set_val(-0.125*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[13].set_val(-0.125*(1.0-spt[0])*(1.0+spt[2]));
		d_nds[14].set_val(0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        n_vec[5].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[15].set_val(0.125*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[16].set_val(-0.125*(1.0+spt[0])*(1.0+spt[2]));
		d_nds[17].set_val(0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        n_vec[6].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[18].set_val(0.125*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[19].set_val(0.125*(1.0+spt[0])*(1.0+spt[2]));
		d_nds[20].set_val(0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        n_vec[7].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[21].set_val(-0.125*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[22].set_val(0.125*(1.0-spt[0])*(1.0+spt[2]));
		d_nds[23].set_val(0.125*(1.0-spt[0])*(1.0+spt[1]));
		
		if(type == 81) {
			n_vec[8].set_val(1.0-spt[0]*spt[0]);
			d_nds[24].set_val(-2.0*spt[0]);
		    d_nds[25].set_val(0.0);
		    d_nds[26].set_val(0.0);
			
            n_vec[9].set_val(1.0-spt[1]*spt[1]);
			d_nds[27].set_val(0.0);
		    d_nds[28].set_val(-2.0*spt[1]);
		    d_nds[29].set_val(0.0);
			
            n_vec[10].set_val(1.0-spt[2]*spt[2]);
			d_nds[30].set_val(0.0);
		    d_nds[31].set_val(0.0);
		    d_nds[32].set_val(-2.0*spt[2]);
		}
	}
	else if (type == 10 || type == 1000) {
		double p1 = 1.0 - spt[0] - spt[1] - spt[2];
		double p2 = p1 - 0.5;
		n_vec[0].set_val(2.0 * p1 * p2);
		d_nds[0].set_val(2.0 * (-p2 - p1));
		d_nds[1].set_val(2.0 * (-p2 - p1));
		d_nds[2].set_val(2.0 * (-p2 - p1));

		n_vec[1].set_val(-2.0 * spt[0] * (0.5 - spt[0]));
		d_nds[3].set_val(-2.0 * (0.5 - spt[0] - spt[0]));
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
		n_vec[2].set_val(-2.0 * spt[1] * (0.5 - spt[1]));
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(-2.0 * (0.5 - spt[1] - spt[1]));
		d_nds[8].set_val(0.0);
		
		n_vec[3].set_val(-2.0 * spt[2] * (0.5 - spt[2]));
		d_nds[9].set_val(0.0);
		d_nds[10].set_val(0.0);
		d_nds[11].set_val(-2.0 * (0.5 - spt[2] - spt[2]));
		
		n_vec[4].set_val(4.0 * spt[0] * p1);
		d_nds[12].set_val(4.0 * (p1 - spt[0]));
		d_nds[13].set_val(-4.0 * spt[0]);
		d_nds[14].set_val(-4.0 * spt[0]);
		
		n_vec[5].set_val(4.0 * spt[0] * spt[1]);
		d_nds[15].set_val(4.0 * spt[1]);
		d_nds[16].set_val(4.0 * spt[0]);
		d_nds[17].set_val(0.0);
		
		n_vec[6].set_val(4.0 * spt[1] * p1);
		d_nds[18].set_val(-4.0 * spt[1]);
		d_nds[19].set_val(4.0 * (p1 - spt[1]));
		d_nds[20].set_val(-4.0 * spt[1]);
		
		n_vec[7].set_val(4.0 * spt[2] * p1);
		d_nds[21].set_val(-4.0 * spt[2]);
		d_nds[22].set_val(-4.0 * spt[2]);
		d_nds[23].set_val(4.0 * (p1 - spt[2]));
		
		n_vec[8].set_val(4.0 * spt[0] * spt[2]);
		d_nds[24].set_val(4.0 * spt[2]);
		d_nds[25].set_val(0.0);
		d_nds[26].set_val(4.0 * spt[0]);
		
		n_vec[9].set_val(4.0 * spt[1] * spt[2]);
		d_nds[27].set_val(0.0);
		d_nds[28].set_val(4.0 * spt[2]);
		d_nds[29].set_val(4.0 * spt[1]);
	}
	else if (type == 3) {
		n_vec[0].set_val(1.0-spt[0]-spt[1]);
		d_nds[0].set_val(-1.0);
		d_nds[1].set_val(-1.0);
		d_nds[2].set_val(0.0);
		
		n_vec[1].set_val(spt[0]);
		d_nds[3].set_val(1.0);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
		n_vec[2].set_val(spt[1]);
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(1.0);
		d_nds[8].set_val(0.0);
		
		n_vec[3].set_val(spt[0]*(1.0-spt[0]-spt[1]));
		d_nds[9].set_val(1.0-spt[0]-spt[1] - spt[0]);
		d_nds[10].set_val(-spt[0]);
		d_nds[11].set_val(0.0);
		
		n_vec[4].set_val(spt[0]*spt[1]);
		d_nds[12].set_val(spt[1]);
		d_nds[13].set_val(spt[0]);
		d_nds[14].set_val(0.0);
		
		n_vec[5].set_val(spt[1]*(1.0-spt[0]-spt[1]));
		d_nds[15].set_val(-spt[1]);
		d_nds[16].set_val(1.0-spt[0]-spt[1]-spt[1]);
		d_nds[17].set_val(0.0);
		
	} else if(type == 41) {
		n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));
		d_nds[0].set_val(-0.25*(1.0-spt[1]));
		d_nds[1].set_val(-0.25*(1.0-spt[0]));
		d_nds[2].set_val(0.0);
		
        n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
		d_nds[3].set_val(0.25*(1.0-spt[1]));
		d_nds[4].set_val(-0.25*(1.0+spt[0]));
		d_nds[5].set_val(0.0);
		
        n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
		d_nds[6].set_val(0.25*(1.0+spt[1]));
		d_nds[7].set_val(0.25*(1.0+spt[0]));
		d_nds[8].set_val(0.0);
		
        n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
		d_nds[9].set_val(-0.25*(1.0+spt[1]));
		d_nds[10].set_val(0.25*(1.0-spt[0]));
		d_nds[11].set_val(0.0);
		
		n_vec[4].set_val(1.0-spt[0]*spt[0]);
		d_nds[12].set_val(-2.0*spt[0]);
		d_nds[13].set_val(0.0);
		d_nds[14].set_val(0.0);
		
        n_vec[5].set_val(1.0-spt[1]*spt[1]);
		d_nds[15].set_val(0.0);
		d_nds[16].set_val(-2.0*spt[1]);
		d_nds[17].set_val(0.0);
		
	    n_vec[6].set_val((1.0-spt[0]*spt[0])*(1.0-spt[1]));
		d_nds[18].set_val(-2.0*spt[0]*(1.0-spt[1]));
		d_nds[19].set_val(-1.0+spt[0]*spt[0]);
		d_nds[20].set_val(0.0);
		
	    n_vec[7].set_val((1.0-spt[1]*spt[1])*(1.0+spt[0]));
		d_nds[21].set_val(1.0-spt[1]*spt[1]);
		d_nds[22].set_val(-2.0*spt[1]*(1.0+spt[0]));
		d_nds[23].set_val(0.0);
		
	    n_vec[8].set_val((1.0-spt[0]*spt[0])*(1.0+spt[1]));
		d_nds[24].set_val(-2.0*spt[0]*(1.0+spt[1]));
		d_nds[25].set_val(1.0-spt[0]*spt[0]);
		d_nds[26].set_val(0.0);
		
	    n_vec[9].set_val((1.0-spt[1]*spt[1])*(1.0-spt[0]));
		d_nds[27].set_val(-1.0+spt[1]*spt[1]);
		d_nds[28].set_val(-2.0*spt[1]*(1.0-spt[0]));
		d_nds[29].set_val(0.0);
		
	} else if(type == 2) {
		n_vec[0].set_val(0.5*(1.0 - spt[0]));
		d_nds[0].set_val(-0.5);
		d_nds[1].set_val(0.0);
		d_nds[2].set_val(0.0);
		
        n_vec[1].set_val(0.5*(1.0 + spt[0]));
		d_nds[3].set_val(0.5);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
        n_vec[2].set_val(1.0 - spt[0]*spt[0]);
		d_nds[0].set_val(-2.0*spt[0]);
		d_nds[1].set_val(0.0);
		d_nds[2].set_val(0.0);
	}
	return;
}

void Element::get_ip_data(DiffDoub0 n_vec[], DiffDoub0 d_ndx[], DiffDoub0& det_j, vector<DiffDoub0>& loc_nds, double spt[]) {
	int i1;
	int i2;
	DiffDoub0 n_cent[11];
	DiffDoub0 d_nds[33];
	DiffDoub0 d_nds_cent[33];
	DiffDoub0 j_mat[9];
	DiffDoub0 tmp_nds[30];
	DiffDoub0 j_cent[9];
	DiffDoub0 det_cent;
	DiffDoub0 j_inv[9];
	DiffDoub0 j_inv_cent[9];
	double s_cent[3] = { 0.0,0.0,0.0 };
	DiffDoub0 x_vec[3];
	DiffDoub0 b_vec[3];
	DiffDoub0 z_dir;
	DiffDoub0 tmp;
	
	eval_n(n_vec, d_nds, spt);
	vec_to_ar(tmp_nds, loc_nds, 0, 30);
	mat_mul(j_mat,tmp_nds,d_nds,3,num_nds,3);

	if(type == 41 || type == 3) {
		z_dir.set_val(j_mat[0]);
		z_dir.mult(j_mat[4]);
		tmp.set_val(j_mat[3]);
		tmp.mult(j_mat[1]);
		z_dir.sub(tmp);
		if(z_dir.val > 0.0) {
			j_mat[8].set_val(1.0);
		} else {
			j_mat[8].set_val(-1.0);
		}
	} else if(type == 2) {
		j_mat[4].set_val(1.0);
		if(j_mat[0].val > 0.0) {
			j_mat[8].set_val(1.0);
		} else {
			j_mat[8].set_val(-1.0);
		}
	}
	
	get_det_inv(det_j,j_inv,j_mat,3,0,x_vec,b_vec);
	
	// mat_mul(d_nds,d_nds,j_inv,n_dim,3,3);
	mat_mul(d_ndx,d_nds,j_inv,num_nds,3,3);

	if (n_dim > num_nds) {
		eval_n(n_cent, d_nds_cent, s_cent);
		mat_mul(j_cent, tmp_nds, d_nds_cent, 3, num_nds, 3);

		if (type == 41 || type == 3) {
			z_dir.set_val(j_cent[0]);
			z_dir.mult(j_cent[4]);
			tmp.set_val(j_cent[3]);
			tmp.mult(j_cent[1]);
			z_dir.sub(tmp);
			if (z_dir.val > 0.0) {
				j_cent[8].set_val(1.0);
			}
			else {
				j_cent[8].set_val(-1.0);
			}
		}
		else if (type == 2) {
			j_cent[4].set_val(1.0);
			if (j_cent[0].val > 0.0) {
				j_cent[8].set_val(1.0);
			}
			else {
				j_cent[8].set_val(-1.0);
			}
		}

		get_det_inv(det_cent, j_inv_cent, j_cent, 3, 0, x_vec, b_vec);

		i1 = 3 * num_nds;
		i2 = n_dim - num_nds;
		mat_mul(&d_ndx[i1], &d_nds[i1], j_inv_cent, i2, 3, 3);
	}
	
	return;
}

void Element::get_inst_ori(vector<DiffDoub0>& inst_ori_mat, vector<DiffDoub0>& loc_ori, vector<DiffDoub0>& glob_disp, int stat) {
	// stat = 1: nonlinear geometry, Element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub0 version, 1st order for DiffDoub1 version
	DiffDoub0 rot[3];
	DiffDoub0 nnds;
	DiffDoub0 one;
	DiffDoub0 tmp_ori[9];
	DiffDoub0 tmp_inst[9];
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int st_index;
	int i_ori_size = (num_nds+1)*144;
	for (i1 = 0; i1 < i_ori_size; i1++) {
		inst_ori_mat[i1].set_val(0.0);
	}
	bool is_diff = loc_ori[0].diff_type();
	
	nnds.set_val(0.0);
	one.set_val(1.0);
	rot[0].set_val(0.0);
	rot[1].set_val(0.0);
	rot[2].set_val(0.0);
	i2 = 3 * n_dim;
	for (i1 = 0; i1 < num_nds; i1++) {
		rot[0].add(glob_disp[i2]);
		rot[1].add(glob_disp[i2 + n_dim]);
		rot[2].add(glob_disp[i2 + 2 * n_dim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);

	vec_to_ar(tmp_ori, loc_ori, 0, 9);
	
	if(stat == 1) {
		d_orid_thet(tmp_inst, tmp_ori, rot, 0, 0);
		ar_to_vec(tmp_inst, inst_ori_mat, 0, 9);
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			st_index = 144 * i1;
			d_orid_thet(tmp_inst, tmp_ori, rot, 0, 0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			for (i2 = 1; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else if(is_diff) {
		for (i2 = 0; i2 < 4; i2++) {
			st_index = 36*i2;
			d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			i3 = st_index;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
				i3++;
				i4++;
			}
		}
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			for (i2 = 0; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else {
		for (i2 = 0; i2 < 4; i2++) {
			st_index = 36*i2;
			d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			i3 = st_index;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				st_index = 36*i2 + 9*i6;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,i6);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 36*i6 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			for (i2 = 0; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					st_index = 144*i1 + 36*i2 + 9*i6;
					d_orid_thet(tmp_inst,tmp_ori,rot,i2,i6);
					ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
					i3 = st_index;
					i4 = 144*i1 + 36*i6 + 9*i2;
					for (i5 = 0; i5 < 9; i5++) {
						inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
						i3++;
						i4++;
					}
				}
			}
		}
	}

	return;
}

void Element::get_inst_disp(DiffDoub0 inst_disp[], vector<DiffDoub0>& glob_disp, vector<DiffDoub0>& inst_ori_mat, vector<DiffDoub0>& loc_ori, vector<DiffDoub0>& x_glob, bool n_lgeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int nd_oind;
	int dof;
	int dof_oind;
	int nd2;
	int dof2;
	int dof2_oind;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 nn_inv;
	DiffDoub0 nn_inv2;
	
	i2 = 6*n_dim;
	for (i1 = 0; i1 < i2; i1++) {
		inst_disp[i1].set_val(0.0);
	}
	
	if (!n_lgeom) {
		if (dv1 == max_int && dv2 == max_int) {
			i7 = 3 * n_dim;
			for (i1 = 0; i1 < 3; i1++) {
				i4 = i1 * n_dim;
				for (i2 = 0; i2 < num_nds; i2++) {
					i5 = i1 * 3;
					i6 = i2;
					for (i3 = 0; i3 < 3; i3++) {
						//i4 = i1 * n_dim + i2;
						//i5 = i1 * 3 + i3;
						//i6 = i3 * n_dim + i2;
						tmp.set_val(loc_ori[i5]);
						tmp.mult(glob_disp[i6]);
						inst_disp[i4].add(tmp);
						tmp.set_val(loc_ori[i5]);
						tmp.mult(glob_disp[i6 + i7]);
						inst_disp[i4 + i7].add(tmp);
						i5++;
						i6 += n_dim;
					}
					i4++;
				}
			}
			i2 = 2 * num_nds * dof_per_nd;
			for (i1 = 0; i1 < num_int_dof; i1++) {
				nd = dof_table[i2];
				dof = dof_table[i2 + 1];
				i3 = dof * n_dim + nd;
				inst_disp[i3].set_val(glob_disp[i3]);
				i2 += 2;
			}
		}
		else if ((dv1 + dv2) >= max_int) {
			dv1 = dv1 + dv2 - max_int;
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < num_nds) {
					i2 = nd;
					i3 = dof;
					for (i1 = 0; i1 < 3; i1++) {
						//i2 = i1 * n_dim + nd;
						//i3 = i1 * 3 + dof;
						inst_disp[i2].set_val(loc_ori[i3]);
						i2 += n_dim;
						i3 += 3;
					}
				}
				else {
					i1 = dof * n_dim + nd;
					inst_disp[i1].set_val(1.0);
				}
			}
			else {
				i2 = 3 * n_dim + nd;
				i3 = dof - 3;
				for (i1 = 0; i1 < 3; i1++) {
					//i2 = (i1 + 3) * n_dim + nd;
					//i3 = i1 * 3 + (dof - 3);
					inst_disp[i2].set_val(loc_ori[i3]);
					i2 += n_dim;
					i3 += 3;
				}
			}
		}
	}
	else {
		if (dv1 == max_int && dv2 == max_int) {
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < num_nds; i2++) {
					i4 = i1 * n_dim + i2;
					i5 = i2;
					i6 = i2;
					i7 = i1 * 3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.set_val(glob_disp[i5]);
						tmp.add(x_glob[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp2.set_val(loc_ori[i7]);
						tmp2.mult(x_glob[i6]);
						tmp.sub(tmp2);
						inst_disp[i4].add(tmp);
						i5 += n_dim;
						i6 += num_nds;
						i7++;
					}
				}
			}

			i2 = num_nds * dof_per_nd;
			i3 = i2 + num_int_dof;
			for (i1 = i2; i1 < i3; i1++) {
				i4 = dof_table[2 * i1];
				i5 = dof_table[2 * i1 + 1];
				i6 = i5 * n_dim + i4;
				inst_disp[i6].set_val(glob_disp[i6]);
			}

			for (i1 = 0; i1 < num_nds; i1++) {
				i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
				i4 = 4 * n_dim + i1;
				i5 = 5 * n_dim + i1;
				nd_oind = 144 * (i1 + 1);
				for (i2 = 0; i2 < 3; i2++) {
					i6 = nd_oind + 3 + i2;
					i7 = 6 + i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i3].add(tmp);
					i6 = nd_oind + 6 + i2;
					i7 = i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i4].add(tmp);
					i6 = nd_oind + i2;
					i7 = 3 + i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i5].add(tmp);
				}
			}
		}
		else if ((dv1 + dv2) >= max_int) {
			dv1 = dv1 + dv2 - max_int;
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < num_nds) {
					i1 = nd; // index in inst_disp
					i2 = dof; //index in inst_ori_mat
					inst_disp[i1].set_val(inst_ori_mat[i2]);
					i1 = n_dim + nd;
					i2 = 3 + dof;
					inst_disp[i1].set_val(inst_ori_mat[i2]);
					i1 = 2 * n_dim + nd;
					i2 = 6 + dof;
					inst_disp[i1].set_val(inst_ori_mat[i2]);
				}
				else {
					i1 = dof * n_dim + nd;
					inst_disp[i1].set_val(1.0);
				}
			}
			else { // dof is rotation
				nn_inv.set_val(1.0 / num_nds);
				dof_oind = 36 * (dof - 2);
				for (i1 = 0; i1 < 3; i1++) {
					for (i2 = 0; i2 < num_nds; i2++) {
						i4 = i1 * n_dim + i2;
						i5 = i2;
						i6 = i2;
						i7 = dof_oind + i1 * 3;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.set_val(glob_disp[i5]);
							tmp.add(x_glob[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i5 += n_dim;
							i6 += num_nds;
							i7++;
						}
					}
				}

				for (i1 = 0; i1 < num_nds; i1++) {
					i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
					i4 = 4 * n_dim + i1;
					i5 = 5 * n_dim + i1;
					nd_oind = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = nd_oind + 3 + i2;
						i7 = dof_oind + 6 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i3].add(tmp);
						i6 = nd_oind + 6 + i2;
						i7 = dof_oind + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i4].add(tmp);
						i6 = nd_oind + i2;
						i7 = dof_oind + 3 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i5].add(tmp);
						if (i1 == nd) {
							i6 = nd_oind + dof_oind + 3 + i2;
							i7 = 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + 6 + i2;
							i7 = i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + i2;
							i7 = 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i5].add(tmp);
						}
					}
				}
			}
		}
		else {
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			nd2 = dof_table[2 * dv2];
			dof2 = dof_table[2 * dv2 + 1];
			nn_inv.set_val(1.0 / num_nds);
			nn_inv2.set_val(nn_inv);
			nn_inv2.sqr();
			if (dof > 2 && dof2 > 2) {
				for (i1 = 0; i1 < 3; i1++) {
					dof_oind = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
					for (i2 = 0; i2 < num_nds; i2++) {
						i4 = n_dim * i1 + i2;
						i5 = i2;
						i6 = i2;
						i7 = dof_oind;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.set_val(glob_disp[i5]);
							tmp.add(x_glob[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv2);
							inst_disp[i4].add(tmp);
							i5 += n_dim;
							i6 += num_nds;
							i7++;
						}
					}
				}

				dof_oind = 36 * (dof - 2);
				dof2_oind = 9 * (dof2 - 2);
				for (i1 = 0; i1 < num_nds; i1++) {
					i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
					i4 = 4 * n_dim + i1;
					i5 = 5 * n_dim + i1;
					nd_oind = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = nd_oind + 3 + i2;
						i7 = dof_oind + dof2_oind + 6 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i3].add(tmp);
						i6 = nd_oind + 6 + i2;
						i7 = dof_oind + dof2_oind + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i4].add(tmp);
						i6 = nd_oind + i2;
						i7 = dof_oind + dof2_oind + 3 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i5].add(tmp);
						if (i1 == nd) {
							i6 = nd_oind + dof_oind + 3 + i2;
							i7 = dof2_oind + 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + 6 + i2;
							i7 = dof2_oind + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + i2;
							i7 = dof2_oind + 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i5].add(tmp);
						}
						if (i1 == nd2) {
							i6 = nd_oind + dof2_oind + 3 + i2;
							i7 = dof_oind + 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof2_oind + 6 + i2;
							i7 = dof_oind + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof2_oind + i2;
							i7 = dof_oind + 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i5].add(tmp);
						}
						if (i1 == nd && i1 == nd2) {
							i6 = nd_oind + dof_oind + dof2_oind + 3 + i2;
							i7 = 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + dof2_oind + 6 + i2;
							i7 = i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + dof2_oind + i2;
							i7 = 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i5].add(tmp);
						}
					}
				}
			}
			else if (dof < 3 && dof2 < 3) {
				return;
			}
			else {
				if (dof > 2) {
					i1 = dof;
					dof = dof2;
					dof2 = i1;
					i1 = nd;
					nd = nd2;
					nd2 = i1;
				}
				i1 = 36 * (dof2 - 2) + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = nd;
				inst_disp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 3 + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = n_dim + nd;
				inst_disp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 6 + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = 2 * n_dim + nd;
				inst_disp[i2].add(tmp);
			}
		}
	}

	
	return;
}

void Element::get_stress_prereq(DiffDoub0StressPrereq& pre, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int num_lay;
	DiffDoub0 offset;
	get_nd_crds(pre.glob_nds, nd_ar, dv_ar);
	get_loc_ori(pre.loc_ori, sec_ar, dv_ar);
	get_nd_disp(pre.glob_disp, nd_ar);
	get_nd_vel(pre.glob_vel, nd_ar);
	get_nd_acc(pre.glob_acc, nd_ar);
	get_nd_temp(pre.glob_temp, nd_ar);
	get_nd_tdot(pre.glob_tdot, nd_ar);
	if (dof_per_nd == 6) {
		correct_orient(pre.loc_ori, pre.glob_nds);
		if (type != 2) {
			get_layer_thk_z(pre.layer_thk, pre.layer_z, offset, sec_ar, dv_ar);
			get_layer_angle(pre.layer_ang, sec_ar, dv_ar);
			get_layer_q(pre.layer_q, sec_ar, mat_ar, dv_ar);
			get_layer_d(pre.layer_d, sec_ar, mat_ar, dv_ar);
			get_layer_th_exp(pre.layer_te, sec_ar, mat_ar, dv_ar);
			get_layer_einit(pre.layer_e0, sec_ar, dv_ar);
			get_layer_den(pre.layer_den, sec_ar, mat_ar, dv_ar);
			get_layer_cond(pre.layer_tc, sec_ar, mat_ar, dv_ar);
			get_layer_spec_heat(pre.layer_sh, sec_ar, mat_ar, dv_ar);
			get_abd(pre.cmat, pre.layer_thk, pre.layer_z, pre.layer_q, pre.layer_ang, sec_ar);
			get_shell_damp(pre.dmat, pre.layer_thk, pre.layer_z, pre.layer_d, pre.layer_ang, sec_ar);
			get_shell_exp_load(pre.therm_exp, pre.einit, pre.layer_thk, pre.layer_z, pre.layer_q, pre.layer_te, pre.layer_e0, pre.layer_ang, sec_ar);
			get_shell_mass(pre.mmat, pre.layer_thk, pre.layer_z, pre.layer_den, sec_ar, dv_ar);
			get_shell_cond(pre.tcmat, pre.layer_thk, pre.layer_ang, pre.layer_tc, sec_ar, dv_ar);
			get_shell_spec_heat(pre.spec_heat, pre.layer_thk, pre.layer_sh, pre.layer_den, sec_ar);
		}
		else {
			get_beam_stiff(pre.cmat, sec_ar, mat_ar, dv_ar);
			get_beam_damp(pre.dmat, sec_ar, mat_ar, dv_ar);
			get_beam_exp_load(pre.therm_exp, pre.einit, sec_ar, mat_ar, dv_ar);
			get_beam_mass(pre.mmat, sec_ar, mat_ar, dv_ar);
			get_beam_cond(pre.tcmat, sec_ar, mat_ar, dv_ar);
			get_beam_spec_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
		}
	}
	else if (type == 21) {
		get_frc_fld_const(pre.frc_fld_coef, pre.frc_fld_exp, sec_ar, dv_ar);
		get_thrm_fld_const(pre.thrm_fld_coef, pre.ref_temp, sec_ar, dv_ar);
	}
	else if (type == 1) {
		get_mass_per_el(pre.mass_per_el, sec_ar, dv_ar);
		get_specific_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
	}
	else {
		get_solid_stiff(pre.cmat, sec_ar, mat_ar, dv_ar);
		get_solid_damp(pre.dmat, sec_ar, mat_ar, dv_ar);
		get_thermal_exp(pre.therm_exp, pre.einit, sec_ar, mat_ar, dv_ar);
		get_density(pre.mmat[0], 0, sec_ar, mat_ar, dv_ar);
		get_conductivity(pre.tcmat, sec_ar, mat_ar, dv_ar);
		get_specific_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
	}
	mat_mul(pre.loc_nds, pre.loc_ori, pre.glob_nds, 3, 3, num_nds);


	return;
}

void Element::get_fluid_prereq(DiffDoub0FlPrereq& pre, vector<Section>& sec_ar, vector<Fluid>& fl_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	double sec_prop;
	
	get_nd_crds(pre.glob_nds, nd_ar, dv_ar);
	get_nd_disp(pre.glob_disp, nd_ar);
	i2 = 3 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		pre.glob_nds[i1].add(pre.glob_disp[i1]);
	}
	get_nd_vel(pre.glob_vel, nd_ar);
	get_nd_fl_den(pre.fl_den, nd_ar);
	get_nd_fl_vel(pre.fl_vel, nd_ar);
	get_nd_temp(pre.fl_temp, nd_ar);
	get_nd_turb_e(pre.fl_turb_e, nd_ar);
	get_nd_fl_den_dot(pre.fl_den_dot, nd_ar);
	get_nd_fl_vdot(pre.fl_vel_dot, nd_ar);
	get_nd_tdot(pre.fl_tdot, nd_ar);
	get_nd_turb_edot(pre.fl_turb_edot, nd_ar);

	Section& this_sec = sec_ar[sect_ptr];
	Fluid& this_fl = fl_ar[this_sec.fl_ptr];

	sec_prop = this_sec.ref_turb_e;
	pre.ref_turb_e.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "refTurbE", dv_ar);

	sec_prop = this_sec.grad_vturb_coef;
	pre.grad_vturb_coef.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "gradVTurbCoef", dv_ar);

	sec_prop = this_sec.diss_turb_coef;
	pre.diss_turb_coef.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "dissTurbCoef", dv_ar);

	sec_prop = this_fl.viscosity;
	pre.ref_visc.set_val(sec_prop);
	get_gen_prop(pre.ref_visc, "viscosity", dv_ar);

	sec_prop = this_sec.den_vis_coef;
	pre.den_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.den_vis_coef, "denVisCoef", dv_ar);

	sec_prop = this_sec.temp_vis_coef;
	pre.temp_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.temp_vis_coef, "tempVisCoef", dv_ar);

	sec_prop = this_sec.turb_vis_coef;
	pre.turb_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.turb_vis_coef, "turbVisCoef", dv_ar);

	sec_prop = this_sec.ref_enth;
	pre.ref_enth.set_val(sec_prop);
	get_gen_prop(pre.ref_enth, "refEnth", dv_ar);

	sec_prop = this_sec.enth_coef;
	pre.den_enth_coef.set_val(sec_prop);
	get_gen_prop(pre.den_enth_coef, "denEnthCoef", dv_ar);

	sec_prop = this_sec.enth_exp;
	pre.den_enth_exp.set_val(sec_prop);
	get_gen_prop(pre.den_enth_exp, "denEnthExp", dv_ar);

	sec_prop = this_sec.pres_coef;
	pre.den_pres_coef.set_val(sec_prop);
	get_gen_prop(pre.den_pres_coef, "denPresCoef", dv_ar);

	sec_prop = this_sec.pres_exp;
	pre.den_pres_exp.set_val(sec_prop);
	get_gen_prop(pre.den_pres_exp, "denPresExp", dv_ar);

	sec_prop = this_sec.ref_den;
	pre.ref_den.set_val(sec_prop);
	get_gen_prop(pre.ref_den, "refDen", dv_ar);

	sec_prop = this_sec.ref_temp;
	pre.ref_temp.set_val(sec_prop);
	get_gen_prop(pre.ref_temp, "refTemp", dv_ar);

	sec_prop = this_fl.therm_cond;
	pre.therm_cond.set_val(sec_prop);
	get_gen_prop(pre.therm_cond, "thermCond", dv_ar);

	sec_prop = this_fl.spec_heat;
	pre.spec_heat.set_val(sec_prop);
	get_gen_prop(pre.spec_heat, "specHeat", dv_ar);

	sec_prop = this_fl.ideal_gas;
	pre.i_gconst.set_val(sec_prop);
	get_gen_prop(pre.i_gconst, "iGConst", dv_ar);

	return;
}

void Element::get_volume(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 thk;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;
    double s_tmp[3];

	if (type == 2) {
		thk.set_val(sec_ar[sect_ptr].area); //rem: update to factor in dvars;
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value(dv_val);
				dv_val.mult(tmp);
				thk.add(dv_val);
			}
		}
	} else if (type == 3 || type == 41) {
		thk.set_val(pre.layer_thk[layer]);
	} else {
		thk.set_val(1.0);
	}

	vol.set_val(0.0);
	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(s_tmp, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, s_tmp);
		tmp.set_val(ip_wt[i1]);
		tmp.mult(det_j);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::get_section_def(DiffDoub0 sec_def[], vector<DiffDoub0>& glob_disp,  vector<DiffDoub0>& inst_ori_mat, vector<DiffDoub0>& loc_ori, vector<DiffDoub0>& x_glob, DiffDoub0 d_ndx[], DiffDoub0 n_vec[], bool n_lgeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub0 inst_disp[60];
	DiffDoub0 ux[9];
	DiffDoub0 rx[9];
	DiffDoub0 rot[3];
	DiffDoub0 tmp;

	if (dof_per_nd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			sec_def[i1].set_val(0.0);
		}
		return;
	}
	
	get_inst_disp(inst_disp, glob_disp, inst_ori_mat, loc_ori, x_glob, n_lgeom, dv1, dv2);
	
	if(dv1 == max_int && dv2 == max_int) {
		mat_mul(ux,inst_disp,d_ndx,3,n_dim,3);
		i1 = 3*n_dim;
		mat_mul(rx,&inst_disp[i1],d_ndx,3,n_dim,3);
		mat_mul(rot,&inst_disp[i1],n_vec,3,n_dim,1);
	} else if((dv1 + dv2) >= max_int) {
		dv1 = dv1 + dv2 - max_int;
		nd = dof_table[2*dv1];
		dof = dof_table[2*dv1+1];
		if(dof < 3) {
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].set_val(inst_disp[i4]);
					ux[i3].mult(d_ndx[i5]);
					rx[i3].set_val(0.0);
					i3++;
					i5++;
				}
				rot[i1].set_val(0.0);
				i4+= n_dim;
			}
		} else {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = n_dim*i1;
					i6 = i2;
					i7 = n_dim*(i1+3);
					ux[i4].set_val(0.0);
					rx[i4].set_val(0.0);
					for (i3 = 0; i3 < num_nds; i3++) {
						tmp.set_val(inst_disp[i5]);
						tmp.mult(d_ndx[i6]);
						ux[i4].add(tmp);
						tmp.set_val(inst_disp[i7]);
						tmp.mult(d_ndx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].set_val(0.0);
				i3 = n_dim*(i1+3);
				for (i2 = 0; i2 < num_nds; i2++) {
					tmp.set_val(inst_disp[i3]);
					tmp.mult(n_vec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		}
	} else {
		nd = dof_table[2*dv1];
		dof = dof_table[2*dv1+1];
		nd2 = dof_table[2*dv2];
		dof2 = dof_table[2*dv2+1];
		if(dof < 3 && dof2 < 3) {
			for (i1 = 0; i1 < 9; i1++) {
				sec_def[i1].set_val(0.0);
			}
			return;
		} else if(dof > 2 && dof2 > 2) {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = n_dim*i1;
					i6 = i2;
					i7 = n_dim*(i1+3);
					ux[i4].set_val(0.0);
					for (i3 = 0; i3 < num_nds; i3++) {
						tmp.set_val(inst_disp[i5]);
						tmp.mult(d_ndx[i6]);
						ux[i4].add(tmp);
						tmp.set_val(inst_disp[i7]);
						tmp.mult(d_ndx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].set_val(0.0);
				i3 = n_dim*(i1+3);
				for (i2 = 0; i2 < num_nds; i2++) {
					tmp.set_val(inst_disp[i3]);
					tmp.mult(n_vec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		} else {
			if(dof > dof2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].set_val(inst_disp[i4]);
					ux[i3].mult(d_ndx[i5]);
					rx[i3].set_val(0.0);
					i3++;
					i5++;
				}
				rot[i1].set_val(0.0);
				i4+= n_dim;
			}
		}
	}
	
	if(type == 2) {
		sec_def[0].set_val(ux[0]);
		sec_def[1].set_val(ux[3]);
		sec_def[1].sub(rot[2]);
		sec_def[2].set_val(ux[6]);
		sec_def[2].add(rot[1]);
		sec_def[3].set_val(rx[0]);
		sec_def[4].set_val(rx[3]);
		sec_def[5].set_val(rx[6]);
	} else {
		sec_def[0].set_val(ux[0]);
		sec_def[1].set_val(ux[4]);
		sec_def[2].set_val(ux[1]);
		sec_def[2].add(ux[3]);
		sec_def[3].set_val(rx[3]);
		sec_def[4].set_val(rx[1]);
		sec_def[4].neg();
		sec_def[5].set_val(rx[4]);
		sec_def[5].sub(rx[0]);
		sec_def[6].set_val(ux[6]);
		sec_def[6].add(rot[1]);
		sec_def[7].set_val(ux[7]);
		sec_def[7].sub(rot[0]);
		sec_def[8].set_val(2.0);
		sec_def[8].mult(rot[2]);
		sec_def[8].sub(ux[3]);
		sec_def[8].add(ux[1]);
	}
	
	return;
}

void Element::get_solid_strain(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 d_ndx[], vector<DiffDoub0>& loc_ori, int dv1, int dv2, bool n_lgeom) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub0 ux_l[9];
	DiffDoub0 strn_mat[9];
	DiffDoub0 tmp_ori[9];
	DiffDoub0 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strn_mat[i1].set_val(0.0);
	}
	if(dv1 == max_int && dv2 == max_int) {
		vec_to_ar(tmp_ori, loc_ori, 0, 9);
		mat_mul(ux_l,tmp_ori,ux,3,3,3);
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			i5 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				strn_mat[i4].add(ux_l[i4]);
				strn_mat[i4].add(ux_l[i5]);
				if(n_lgeom) {
					i6 = i1;
					i7 = i2;
				    for (i3 = 0; i3 < 3; i3++) {
						tmp.set_val(ux[i6]);
						tmp.mult(ux[i7]);
						strn_mat[i4].add(tmp);
						i6+= 3;
						i7+= 3;
				    }
				}
				i4++;
				i5+= 3;
			}
		}
	} else if((dv1 + dv2) >= max_int) {
		if(dv1 < max_int) {
			nd = dof_table[2*dv1];
			dof = dof_table[2*dv1+1];
		} else {
			nd = dof_table[2*dv2];
			dof = dof_table[2*dv2+1];
		}
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				i5 = 3*i1 + dof;
				i6 = 3*nd + i2;
				tmp.set_val(loc_ori[i5]);
				tmp.mult(d_ndx[i6]);
				strn_mat[i4].add(tmp);
				i5 = 3*i2 + dof;
				i6 = 3*nd + i1;
				tmp.set_val(loc_ori[i5]);
				tmp.mult(d_ndx[i6]);
				strn_mat[i4].add(tmp);
				if(n_lgeom) {
					i5 = 3*nd + i1;
					i6 = 3*dof + i2;
					tmp.set_val(d_ndx[i5]);
					tmp.mult(ux[i6]);
					strn_mat[i4].add(tmp);
					i5 = 3*nd + i2;
					i6 = 3*dof + i1;
					tmp.set_val(d_ndx[i5]);
					tmp.mult(ux[i6]);
					strn_mat[i4].add(tmp);
				}
				i4++;
			}
		}
	} else {
		if(n_lgeom) {
			nd = dof_table[2*dv1];
			dof = dof_table[2*dv1+1];
			nd2 = dof_table[2*dv2];
			dof2 = dof_table[2*dv2+1];
			if(dof == dof2) {
				for (i1 = 0; i1 < 3; i1++) {
					i4 = 4*i1;
					for (i2 = i1; i2 < 3; i2++) {
						i5 = 3*nd + i1;
						i6 = 3*nd2 + i2;
						tmp.set_val(d_ndx[i5]);
						tmp.mult(d_ndx[i6]);
						strn_mat[i4].add(tmp);
						i5 = 3*nd2 + i1;
						i6 = 3*nd + i2;
						tmp.set_val(d_ndx[i5]);
						tmp.mult(d_ndx[i6]);
						strn_mat[i4].add(tmp);
						i4++;
					}
				}
			}
		}
	}
	
	tmp.set_val(0.5);
	strain[0].set_val(strn_mat[0]);
	strain[0].mult(tmp);
	strain[1].set_val(strn_mat[4]);
	strain[1].mult(tmp);
	strain[2].set_val(strn_mat[8]);
	strain[2].mult(tmp);
	strain[3].set_val(strn_mat[1]);
	strain[4].set_val(strn_mat[2]);
	strain[5].set_val(strn_mat[5]);
	
	return;
}

void Element::get_stress_strain(DiffDoub0 stress[], DiffDoub0 strain[], double spt[], int layer, bool n_lgeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 ux[9];
	DiffDoub0 sec_def[9];
	DiffDoub0 sect_strn[3];
	DiffDoub0 adj_stn[6];
	DiffDoub0 tmp;
	DiffDoub0 ip_temp;
	DiffDoub0 tmp_ar[60];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);

	ip_temp.set_val(0.0);
	for (i1 = 0; i1 < num_nds; i1++) {
		tmp.set_val(pre.glob_temp[i1]);
		tmp.mult(n_vec[i1]);
		ip_temp.add(tmp);
	}

	if (type == 41 || type == 3) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}

		get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, max_int, max_int);
		
		sect_strn[0].set_val(sec_def[0]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[3]);
		sect_strn[0].add(tmp);
		
		sect_strn[1].set_val(sec_def[1]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[4]);
		sect_strn[1].sub(tmp);

		sect_strn[2].set_val(sec_def[2]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[5]);
		sect_strn[2].add(tmp);

		tmp.set_val(pre.layer_ang[layer]);
		tmp.neg();
		transform_strain(strain, sect_strn, tmp);
		i2 = 3 * layer;
		for (i1 = 0; i1 < 3; i1++) {
			tmp.set_val(pre.layer_te[i2]);
			tmp.mult(ip_temp);
			adj_stn[i1].set_val(strain[i1]);
			adj_stn[i1].sub(tmp);
			adj_stn[i1].sub(pre.layer_e0[i2]);
			i2++;
		}
		mat_mul(stress, &pre.layer_q[9 * layer], adj_stn, 3, 3, 1);

		strain[3].set_val(strain[2]);
		strain[2].set_val(0.0);
		strain[4].set_val(0.0);
		strain[5].set_val(0.0);

		stress[3].set_val(stress[2]);
		stress[2].set_val(0.0);
		stress[4].set_val(0.0);
		stress[5].set_val(0.0);

	} else if(type != 2) {
		vec_to_ar(tmp_ar, pre.glob_disp, 0, 60);
		mat_mul(ux, tmp_ar, d_ndx, 3, n_dim, 3);
		get_solid_strain(strain, ux, d_ndx, pre.loc_ori, max_int, max_int, n_lgeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.set_val(pre.therm_exp[i1]);
			tmp.mult(ip_temp);
			adj_stn[i1].set_val(strain[i1]);
			adj_stn[i1].sub(tmp);
			adj_stn[i1].sub(pre.einit[i1]);
		}
		vec_to_ar(tmp_ar, pre.cmat, 0, 36);
		mat_mul(stress, tmp_ar, adj_stn, 6, 6, 1);
	}
	return;
}

void Element::d_stress_straind_u(vector<DiffDoub0>& dsd_u, vector<DiffDoub0>& ded_u, vector<DiffDoub0>& dsd_t, double spt[], int layer, bool n_lgeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 ux[9];
	DiffDoub0 sec_def[9];
	DiffDoub0 sect_strn[3];
	DiffDoub0 d_strain[6];
	DiffDoub0 d_stress[6];
	DiffDoub0 cte[6];
	DiffDoub0 cten[60];
	DiffDoub0 tmp;
	DiffDoub0 tmp_ar[60];
	DiffDoub0 tmp_ar2[6];

	int tot_dof = num_nds * dof_per_nd + num_int_dof;
	i2 = 6 * tot_dof;
	for (i1 = 0; i1 < i2; i1++) {
		dsd_u[i1].set_val(0.0);
		ded_u[i1].set_val(0.0);
	}
	i2 = 6 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		dsd_t[i1].set_val(0.0);
	}

	if (type == 41 || type == 3) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}

		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
		for (i1 = 0; i1 < tot_dof; i1++) {
			get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, i1, max_int);

			sect_strn[0].set_val(sec_def[0]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[3]);
			sect_strn[0].add(tmp);

			sect_strn[1].set_val(sec_def[1]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[4]);
			sect_strn[1].sub(tmp);

			sect_strn[2].set_val(sec_def[2]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[5]);
			sect_strn[2].add(tmp);

			tmp.set_val(pre.layer_ang[layer]);
			tmp.neg();
			transform_strain(d_strain, sect_strn, tmp);
			mat_mul(d_stress, &pre.layer_q[9 * layer], d_strain, 3, 3, 1);

			ded_u[i1].set_val(d_strain[0]);
			ded_u[i1 + tot_dof].set_val(d_strain[1]);
			ded_u[i1 + 3 * tot_dof].set_val(d_strain[2]);

			dsd_u[i1].set_val(d_stress[0]);
			dsd_u[tot_dof + i1].set_val(d_stress[1]);
			dsd_u[3 * tot_dof + i1].set_val(d_stress[2]);
		}
		mat_mul(cte, &pre.layer_q[9 * layer], &pre.layer_te[3 * layer], 3, 3, 1);
		cte[0].neg();
		cte[1].neg();
		cte[2].neg();
		mat_mul(cten, cte, n_vec, 3, 1, num_nds);
		for (i1 = 0; i1 < num_nds; i1++) {
			dsd_t[i1].set_val(cten[i1]);
			dsd_t[i1 + num_nds].set_val(cten[i1 + num_nds]);
			dsd_t[i1 + 3 * num_nds].set_val(cten[i1 + 2 * num_nds]);
		}
	}
	else if (type != 2) {
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
		vec_to_ar(tmp_ar, pre.glob_disp, 0, 60);
		mat_mul(ux, tmp_ar, d_ndx, 3, n_dim, 3);
		vec_to_ar(tmp_ar, pre.cmat, 0, 36);
		for (i1 = 0; i1 < tot_dof; i1++) {
			get_solid_strain(d_strain, ux, d_ndx, pre.loc_ori, i1, max_int, n_lgeom);
			mat_mul(d_stress, tmp_ar, d_strain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				ded_u[i3].set_val(d_strain[i2]);
				dsd_u[i3].set_val(d_stress[i2]);
				i3 += tot_dof;
			}
		}
		vec_to_ar(tmp_ar2, pre.therm_exp, 0, 6);
		mat_mul(cte, tmp_ar, tmp_ar2, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			cte[i1].neg();
		}
		mat_mul(tmp_ar, cte, n_vec, 6, 1, num_nds);
		ar_to_vec(tmp_ar, dsd_t, 0, 60);
	}

	return;
}

void Element::get_def_frc_mom(DiffDoub0 def[], DiffDoub0 frc_mom[], double spt[], bool n_lgeom, DiffDoub0StressPrereq& pre) {
	int i1;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 pt_temp;
	DiffDoub0 tmp;
	DiffDoub0 tmp_ar[36];

	if (dof_per_nd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			def[i1].set_val(0.0);
		}
		return;
	}

	if (n_lgeom) {
		get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
	}

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	get_section_def(def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, max_int, max_int);

	pt_temp.set_val(0.0);
	for (i1 = 0; i1 < num_nds; i1++) {
		tmp.set_val(pre.glob_temp[i1]);
		tmp.mult(n_vec[i1]);
		pt_temp.add(tmp);
	}

	vec_to_ar(tmp_ar, pre.cmat, 0, 36);
	mat_mul(frc_mom, tmp_ar, def, def_dim, def_dim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frc_mom[i1].sub(pre.einit[i1]);
		tmp.set_val(pt_temp);
		tmp.mult(pre.therm_exp[i1]);
		frc_mom[i1].sub(tmp);
	}

	return;
}

void Element::d_def_frc_momd_u(vector<DiffDoub0>& d_defd_u, vector<DiffDoub0>& d_frc_momd_u, vector<DiffDoub0>& d_frc_momd_t, double spt[], bool n_lgeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	int tot_dof;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 pt_temp;
	DiffDoub0 tmp;
	DiffDoub0 def[9];
	DiffDoub0 tmp_ar[60];
	DiffDoub0 tmp_ar2[6];

	if (n_lgeom) {
		get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
	}


	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	tot_dof = num_nds * dof_per_nd + num_int_dof;
	for (i1 = 0; i1 < tot_dof; i1++) {
		get_section_def(def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, i1, max_int);
		i2 = i1;
		for (i3 = 0; i3 < def_dim; i3++) {
			d_defd_u[i2].set_val(def[i3]);
			i2 += tot_dof;
		}
	}

	mat_mul(d_frc_momd_u, pre.cmat, d_defd_u, def_dim, def_dim, tot_dof);

	vec_to_ar(tmp_ar2, pre.therm_exp, 0, 6);
	mat_mul(tmp_ar, tmp_ar2, n_vec, 6, 1, num_nds);
	ar_to_vec(tmp_ar, d_frc_momd_t, 0, 60);

	return;
}

void Element::get_flux_tgrad(DiffDoub0 flux[], DiffDoub0 t_grad[], double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 tmp_ar[10];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	vec_to_ar(tmp_ar, pre.glob_temp, 0, 10);
	mat_mul(t_grad, tmp_ar, d_ndx, 1, num_nds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vec_to_ar(tmp_ar, pre.layer_tc, i1, i1 + 9);
		mat_mul(flux, tmp_ar, t_grad, 3, 3, 1);
	}
	else {
		vec_to_ar(tmp_ar, pre.tcmat, 0, 9);
		mat_mul(flux, tmp_ar, t_grad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::d_flux_tgradd_t(vector<DiffDoub0>& d_fd_t, vector<DiffDoub0>& d_tg, double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub0 n_vec[11];
	DiffDoub0 d_ndx[33];
	DiffDoub0 det_j;
	DiffDoub0 tmp_ar1[33];
	DiffDoub0 tmp_ar2[9];
	DiffDoub0 tmp_ar3[30];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	transpose(tmp_ar1, d_ndx, num_nds, 3);
	ar_to_vec(tmp_ar1, d_tg, 0, 33);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vec_to_ar(tmp_ar2, pre.layer_tc, i1, i1 + 9);
		mat_mul(tmp_ar3, tmp_ar2, tmp_ar1, 3, 3, num_nds);
		ar_to_vec(tmp_ar3, d_fd_t, 0, 30);
	}
	else {
		mat_mul(d_fd_t, pre.tcmat, d_tg, 3, 3, num_nds);
	}
	i2 = 3 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		d_fd_t[i1].neg();
	}

	return;
}

void Element::put_vec_to_glob_mat(SparseMat& q_mat, vector<DiffDoub0>& el_qvec, bool for_therm, int mat_row, vector<Node>& nd_ar) {
	int i1;
	int nd_dof = num_nds*dof_per_nd;
	int tot_dof = nd_dof + num_int_dof;
	int nd;
	int dof;
	int glob_ind;

	if (for_therm) {
		for (i1 = 0; i1 < num_nds; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[i1].sorted_rank;
			q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < tot_dof; i1++) {
			if (i1 < nd_dof) {
				nd = nodes[dof_table[2 * i1]];
				dof = dof_table[2 * i1 + 1];
				glob_ind = nd_ar[nd].dof_index[dof];
				q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
			}
			else {
				dof = i1 - nd_dof;
				glob_ind = int_dof_index + dof;
				q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
			}
		}
	}

	return;
}

//end dup
 
//skip 
 
//diff_doub1 versions: 
//dup1
void Element::get_nd_disp(vector<DiffDoub1>& glob_disp, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_disp[i3].set_val(this_nd.displacement[i2]);
			i3+= n_dim;
		}
	}
	
	if(num_int_dof > 0) {
		i2 = 2*num_nds*dof_per_nd;
		for (i3 = 0; i3 < num_int_dof; i3++) {
			i4 = dof_table[i2];
			i5 = dof_table[i2+1];
			i6 = n_dim*i5 + i4;
			glob_disp[i6].set_val(internal_disp[i3]);
			i2+= 2;
		}
	}
	
	return;
}

void Element::get_nd_vel(vector<DiffDoub1>& glob_vel, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_vel[i3].set_val(this_nd.velocity[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_acc(vector<DiffDoub1>& glob_acc, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			glob_acc[i3].set_val(this_nd.acceleration[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_fl_vel(vector<DiffDoub1>& fl_vel, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			fl_vel[i3].set_val(this_nd.fl_vel[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_fl_vdot(vector<DiffDoub1>& fl_vdot, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			fl_vdot[i3].set_val(this_nd.fl_vel_dot[i2]);
			i3 += num_nds;
		}
	}

	return;
}

void Element::get_nd_temp(vector<DiffDoub1>& glob_temp, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		glob_temp[i1].set_val(this_nd.temperature);
	}
	return;
}

void Element::get_nd_tdot(vector<DiffDoub1>& glob_tdot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		glob_tdot[i1].set_val(this_nd.temp_change_rate);
	}
	return;
}

void Element::get_nd_fl_den(vector<DiffDoub1>& fl_den, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		fl_den[i1].set_val(this_nd.fl_den);
	}

	return;
}

void Element::get_nd_fl_den_dot(vector<DiffDoub1>& fl_den_dot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		fl_den_dot[i1].set_val(this_nd.fl_den_dot);
	}

	return;
}

void Element::get_nd_turb_e(vector<DiffDoub1>& turb_e, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		turb_e[i1].set_val(this_nd.turb_e);
	}

	return;
}

void Element::get_nd_turb_edot(vector<DiffDoub1>& turb_edot, vector<Node>& nd_ar) {
	int i1;

	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		turb_edot[i1].set_val(this_nd.turb_edot);
	}

	return;
}

void Element::eval_n(DiffDoub1 n_vec[], DiffDoub1 d_nds[], double spt[]) {
	if(type == 4 || type == 400) {
		n_vec[0].set_val(1.0-spt[0]-spt[1]-spt[2]);
		d_nds[0].set_val(-1.0);
		d_nds[1].set_val(-1.0);
		d_nds[2].set_val(-1.0);
        
		n_vec[1].set_val(spt[0]);
		d_nds[3].set_val(1.0);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
        
	    n_vec[2].set_val(spt[1]);
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(1.0);
		d_nds[8].set_val(0.0);
		
        n_vec[3].set_val(spt[2]);
		d_nds[9].set_val(0.0);
		d_nds[10].set_val(0.0);
		d_nds[11].set_val(1.0);
	} else if(type == 6 || type == 600) {
		n_vec[0].set_val(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
		d_nds[0].set_val(-0.5*(1.0-spt[2]));
	    d_nds[1].set_val(-0.5*(1.0-spt[2]));
		d_nds[2].set_val(-0.5*(1.0-spt[0]-spt[1]));
		
	    n_vec[1].set_val(0.5*spt[0]*(1.0-spt[2]));
		d_nds[3].set_val(0.5*(1.0-spt[2]));
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(-0.5*spt[0]);
		
	    n_vec[2].set_val(0.5*spt[1]*(1.0-spt[2]));
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(0.5*(1.0-spt[2]));
		d_nds[8].set_val(-0.5*spt[1]);
		
		n_vec[3].set_val(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
		d_nds[9].set_val(-0.5*(1.0+spt[2]));
		d_nds[10].set_val(-0.5*(1.0+spt[2]));
		d_nds[11].set_val(0.5*(1.0-spt[0]-spt[1]));
		
	    n_vec[4].set_val(0.5*spt[0]*(1.0+spt[2]));
		d_nds[12].set_val(0.5*(1.0+spt[2]));
		d_nds[13].set_val(0.0);
		d_nds[14].set_val(0.5*spt[0]);
		
	    n_vec[5].set_val(0.5*spt[1]*(1.0+spt[2]));
		d_nds[15].set_val(0.0);
		d_nds[16].set_val(0.5*(1.0+spt[2]));
		d_nds[17].set_val(0.5*spt[1]);
	} else if(type == 8 || type == 81 || type == 800) {
		n_vec[0].set_val(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
		d_nds[0].set_val(-0.125*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[1].set_val(-0.125*(1.0-spt[0])*(1.0-spt[2]));
		d_nds[2].set_val(-0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        n_vec[1].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[3].set_val(0.125*(1.0-spt[1])*(1.0-spt[2]));
		d_nds[4].set_val(-0.125*(1.0+spt[0])*(1.0-spt[2]));
		d_nds[5].set_val(-0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        n_vec[2].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[6].set_val(0.125*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[7].set_val(0.125*(1.0+spt[0])*(1.0-spt[2]));
		d_nds[8].set_val(-0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        n_vec[3].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[9].set_val(-0.125*(1.0+spt[1])*(1.0-spt[2]));
		d_nds[10].set_val(0.125*(1.0-spt[0])*(1.0-spt[2]));
		d_nds[11].set_val(-0.125*(1.0-spt[0])*(1.0+spt[1]));
		
        n_vec[4].set_val(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[12].set_val(-0.125*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[13].set_val(-0.125*(1.0-spt[0])*(1.0+spt[2]));
		d_nds[14].set_val(0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        n_vec[5].set_val(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[15].set_val(0.125*(1.0-spt[1])*(1.0+spt[2]));
		d_nds[16].set_val(-0.125*(1.0+spt[0])*(1.0+spt[2]));
		d_nds[17].set_val(0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        n_vec[6].set_val(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[18].set_val(0.125*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[19].set_val(0.125*(1.0+spt[0])*(1.0+spt[2]));
		d_nds[20].set_val(0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        n_vec[7].set_val(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[21].set_val(-0.125*(1.0+spt[1])*(1.0+spt[2]));
		d_nds[22].set_val(0.125*(1.0-spt[0])*(1.0+spt[2]));
		d_nds[23].set_val(0.125*(1.0-spt[0])*(1.0+spt[1]));
		
		if(type == 81) {
			n_vec[8].set_val(1.0-spt[0]*spt[0]);
			d_nds[24].set_val(-2.0*spt[0]);
		    d_nds[25].set_val(0.0);
		    d_nds[26].set_val(0.0);
			
            n_vec[9].set_val(1.0-spt[1]*spt[1]);
			d_nds[27].set_val(0.0);
		    d_nds[28].set_val(-2.0*spt[1]);
		    d_nds[29].set_val(0.0);
			
            n_vec[10].set_val(1.0-spt[2]*spt[2]);
			d_nds[30].set_val(0.0);
		    d_nds[31].set_val(0.0);
		    d_nds[32].set_val(-2.0*spt[2]);
		}
	}
	else if (type == 10 || type == 1000) {
		double p1 = 1.0 - spt[0] - spt[1] - spt[2];
		double p2 = p1 - 0.5;
		n_vec[0].set_val(2.0 * p1 * p2);
		d_nds[0].set_val(2.0 * (-p2 - p1));
		d_nds[1].set_val(2.0 * (-p2 - p1));
		d_nds[2].set_val(2.0 * (-p2 - p1));

		n_vec[1].set_val(-2.0 * spt[0] * (0.5 - spt[0]));
		d_nds[3].set_val(-2.0 * (0.5 - spt[0] - spt[0]));
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
		n_vec[2].set_val(-2.0 * spt[1] * (0.5 - spt[1]));
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(-2.0 * (0.5 - spt[1] - spt[1]));
		d_nds[8].set_val(0.0);
		
		n_vec[3].set_val(-2.0 * spt[2] * (0.5 - spt[2]));
		d_nds[9].set_val(0.0);
		d_nds[10].set_val(0.0);
		d_nds[11].set_val(-2.0 * (0.5 - spt[2] - spt[2]));
		
		n_vec[4].set_val(4.0 * spt[0] * p1);
		d_nds[12].set_val(4.0 * (p1 - spt[0]));
		d_nds[13].set_val(-4.0 * spt[0]);
		d_nds[14].set_val(-4.0 * spt[0]);
		
		n_vec[5].set_val(4.0 * spt[0] * spt[1]);
		d_nds[15].set_val(4.0 * spt[1]);
		d_nds[16].set_val(4.0 * spt[0]);
		d_nds[17].set_val(0.0);
		
		n_vec[6].set_val(4.0 * spt[1] * p1);
		d_nds[18].set_val(-4.0 * spt[1]);
		d_nds[19].set_val(4.0 * (p1 - spt[1]));
		d_nds[20].set_val(-4.0 * spt[1]);
		
		n_vec[7].set_val(4.0 * spt[2] * p1);
		d_nds[21].set_val(-4.0 * spt[2]);
		d_nds[22].set_val(-4.0 * spt[2]);
		d_nds[23].set_val(4.0 * (p1 - spt[2]));
		
		n_vec[8].set_val(4.0 * spt[0] * spt[2]);
		d_nds[24].set_val(4.0 * spt[2]);
		d_nds[25].set_val(0.0);
		d_nds[26].set_val(4.0 * spt[0]);
		
		n_vec[9].set_val(4.0 * spt[1] * spt[2]);
		d_nds[27].set_val(0.0);
		d_nds[28].set_val(4.0 * spt[2]);
		d_nds[29].set_val(4.0 * spt[1]);
	}
	else if (type == 3) {
		n_vec[0].set_val(1.0-spt[0]-spt[1]);
		d_nds[0].set_val(-1.0);
		d_nds[1].set_val(-1.0);
		d_nds[2].set_val(0.0);
		
		n_vec[1].set_val(spt[0]);
		d_nds[3].set_val(1.0);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
		n_vec[2].set_val(spt[1]);
		d_nds[6].set_val(0.0);
		d_nds[7].set_val(1.0);
		d_nds[8].set_val(0.0);
		
		n_vec[3].set_val(spt[0]*(1.0-spt[0]-spt[1]));
		d_nds[9].set_val(1.0-spt[0]-spt[1] - spt[0]);
		d_nds[10].set_val(-spt[0]);
		d_nds[11].set_val(0.0);
		
		n_vec[4].set_val(spt[0]*spt[1]);
		d_nds[12].set_val(spt[1]);
		d_nds[13].set_val(spt[0]);
		d_nds[14].set_val(0.0);
		
		n_vec[5].set_val(spt[1]*(1.0-spt[0]-spt[1]));
		d_nds[15].set_val(-spt[1]);
		d_nds[16].set_val(1.0-spt[0]-spt[1]-spt[1]);
		d_nds[17].set_val(0.0);
		
	} else if(type == 41) {
		n_vec[0].set_val(0.25*(1.0-spt[0])*(1.0-spt[1]));
		d_nds[0].set_val(-0.25*(1.0-spt[1]));
		d_nds[1].set_val(-0.25*(1.0-spt[0]));
		d_nds[2].set_val(0.0);
		
        n_vec[1].set_val(0.25*(1.0+spt[0])*(1.0-spt[1]));
		d_nds[3].set_val(0.25*(1.0-spt[1]));
		d_nds[4].set_val(-0.25*(1.0+spt[0]));
		d_nds[5].set_val(0.0);
		
        n_vec[2].set_val(0.25*(1.0+spt[0])*(1.0+spt[1]));
		d_nds[6].set_val(0.25*(1.0+spt[1]));
		d_nds[7].set_val(0.25*(1.0+spt[0]));
		d_nds[8].set_val(0.0);
		
        n_vec[3].set_val(0.25*(1.0-spt[0])*(1.0+spt[1]));
		d_nds[9].set_val(-0.25*(1.0+spt[1]));
		d_nds[10].set_val(0.25*(1.0-spt[0]));
		d_nds[11].set_val(0.0);
		
		n_vec[4].set_val(1.0-spt[0]*spt[0]);
		d_nds[12].set_val(-2.0*spt[0]);
		d_nds[13].set_val(0.0);
		d_nds[14].set_val(0.0);
		
        n_vec[5].set_val(1.0-spt[1]*spt[1]);
		d_nds[15].set_val(0.0);
		d_nds[16].set_val(-2.0*spt[1]);
		d_nds[17].set_val(0.0);
		
	    n_vec[6].set_val((1.0-spt[0]*spt[0])*(1.0-spt[1]));
		d_nds[18].set_val(-2.0*spt[0]*(1.0-spt[1]));
		d_nds[19].set_val(-1.0+spt[0]*spt[0]);
		d_nds[20].set_val(0.0);
		
	    n_vec[7].set_val((1.0-spt[1]*spt[1])*(1.0+spt[0]));
		d_nds[21].set_val(1.0-spt[1]*spt[1]);
		d_nds[22].set_val(-2.0*spt[1]*(1.0+spt[0]));
		d_nds[23].set_val(0.0);
		
	    n_vec[8].set_val((1.0-spt[0]*spt[0])*(1.0+spt[1]));
		d_nds[24].set_val(-2.0*spt[0]*(1.0+spt[1]));
		d_nds[25].set_val(1.0-spt[0]*spt[0]);
		d_nds[26].set_val(0.0);
		
	    n_vec[9].set_val((1.0-spt[1]*spt[1])*(1.0-spt[0]));
		d_nds[27].set_val(-1.0+spt[1]*spt[1]);
		d_nds[28].set_val(-2.0*spt[1]*(1.0-spt[0]));
		d_nds[29].set_val(0.0);
		
	} else if(type == 2) {
		n_vec[0].set_val(0.5*(1.0 - spt[0]));
		d_nds[0].set_val(-0.5);
		d_nds[1].set_val(0.0);
		d_nds[2].set_val(0.0);
		
        n_vec[1].set_val(0.5*(1.0 + spt[0]));
		d_nds[3].set_val(0.5);
		d_nds[4].set_val(0.0);
		d_nds[5].set_val(0.0);
		
        n_vec[2].set_val(1.0 - spt[0]*spt[0]);
		d_nds[0].set_val(-2.0*spt[0]);
		d_nds[1].set_val(0.0);
		d_nds[2].set_val(0.0);
	}
	return;
}

void Element::get_ip_data(DiffDoub1 n_vec[], DiffDoub1 d_ndx[], DiffDoub1& det_j, vector<DiffDoub1>& loc_nds, double spt[]) {
	int i1;
	int i2;
	DiffDoub1 n_cent[11];
	DiffDoub1 d_nds[33];
	DiffDoub1 d_nds_cent[33];
	DiffDoub1 j_mat[9];
	DiffDoub1 tmp_nds[30];
	DiffDoub1 j_cent[9];
	DiffDoub1 det_cent;
	DiffDoub1 j_inv[9];
	DiffDoub1 j_inv_cent[9];
	double s_cent[3] = { 0.0,0.0,0.0 };
	DiffDoub1 x_vec[3];
	DiffDoub1 b_vec[3];
	DiffDoub1 z_dir;
	DiffDoub1 tmp;
	
	eval_n(n_vec, d_nds, spt);
	vec_to_ar(tmp_nds, loc_nds, 0, 30);
	mat_mul(j_mat,tmp_nds,d_nds,3,num_nds,3);

	if(type == 41 || type == 3) {
		z_dir.set_val(j_mat[0]);
		z_dir.mult(j_mat[4]);
		tmp.set_val(j_mat[3]);
		tmp.mult(j_mat[1]);
		z_dir.sub(tmp);
		if(z_dir.val > 0.0) {
			j_mat[8].set_val(1.0);
		} else {
			j_mat[8].set_val(-1.0);
		}
	} else if(type == 2) {
		j_mat[4].set_val(1.0);
		if(j_mat[0].val > 0.0) {
			j_mat[8].set_val(1.0);
		} else {
			j_mat[8].set_val(-1.0);
		}
	}
	
	get_det_inv(det_j,j_inv,j_mat,3,0,x_vec,b_vec);
	
	// mat_mul(d_nds,d_nds,j_inv,n_dim,3,3);
	mat_mul(d_ndx,d_nds,j_inv,num_nds,3,3);

	if (n_dim > num_nds) {
		eval_n(n_cent, d_nds_cent, s_cent);
		mat_mul(j_cent, tmp_nds, d_nds_cent, 3, num_nds, 3);

		if (type == 41 || type == 3) {
			z_dir.set_val(j_cent[0]);
			z_dir.mult(j_cent[4]);
			tmp.set_val(j_cent[3]);
			tmp.mult(j_cent[1]);
			z_dir.sub(tmp);
			if (z_dir.val > 0.0) {
				j_cent[8].set_val(1.0);
			}
			else {
				j_cent[8].set_val(-1.0);
			}
		}
		else if (type == 2) {
			j_cent[4].set_val(1.0);
			if (j_cent[0].val > 0.0) {
				j_cent[8].set_val(1.0);
			}
			else {
				j_cent[8].set_val(-1.0);
			}
		}

		get_det_inv(det_cent, j_inv_cent, j_cent, 3, 0, x_vec, b_vec);

		i1 = 3 * num_nds;
		i2 = n_dim - num_nds;
		mat_mul(&d_ndx[i1], &d_nds[i1], j_inv_cent, i2, 3, 3);
	}
	
	return;
}

void Element::get_inst_ori(vector<DiffDoub1>& inst_ori_mat, vector<DiffDoub1>& loc_ori, vector<DiffDoub1>& glob_disp, int stat) {
	// stat = 1: nonlinear geometry, Element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub1 version, 1st order for DiffDoub1 version
	DiffDoub1 rot[3];
	DiffDoub1 nnds;
	DiffDoub1 one;
	DiffDoub1 tmp_ori[9];
	DiffDoub1 tmp_inst[9];
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int st_index;
	int i_ori_size = (num_nds+1)*144;
	for (i1 = 0; i1 < i_ori_size; i1++) {
		inst_ori_mat[i1].set_val(0.0);
	}
	bool is_diff = loc_ori[0].diff_type();
	
	nnds.set_val(0.0);
	one.set_val(1.0);
	rot[0].set_val(0.0);
	rot[1].set_val(0.0);
	rot[2].set_val(0.0);
	i2 = 3 * n_dim;
	for (i1 = 0; i1 < num_nds; i1++) {
		rot[0].add(glob_disp[i2]);
		rot[1].add(glob_disp[i2 + n_dim]);
		rot[2].add(glob_disp[i2 + 2 * n_dim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);

	vec_to_ar(tmp_ori, loc_ori, 0, 9);
	
	if(stat == 1) {
		d_orid_thet(tmp_inst, tmp_ori, rot, 0, 0);
		ar_to_vec(tmp_inst, inst_ori_mat, 0, 9);
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			st_index = 144 * i1;
			d_orid_thet(tmp_inst, tmp_ori, rot, 0, 0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			for (i2 = 1; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else if(is_diff) {
		for (i2 = 0; i2 < 4; i2++) {
			st_index = 36*i2;
			d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			i3 = st_index;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
				i3++;
				i4++;
			}
		}
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			for (i2 = 0; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else {
		for (i2 = 0; i2 < 4; i2++) {
			st_index = 36*i2;
			d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
			ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
			i3 = st_index;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				st_index = 36*i2 + 9*i6;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,i6);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 36*i6 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
			}
		}
		for (i1 = 1; i1 <= num_nds; i1++) {
			i2 = 3*n_dim + i1 - 1;
			rot[0].set_val(glob_disp[i2]);
			rot[1].set_val(glob_disp[i2+n_dim]);
			rot[2].set_val(glob_disp[i2+2*n_dim]);
			for (i2 = 0; i2 < 4; i2++) {
				st_index = 144*i1 + 36*i2;
				d_orid_thet(tmp_inst,tmp_ori,rot,i2,0);
				ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
				i3 = st_index;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					st_index = 144*i1 + 36*i2 + 9*i6;
					d_orid_thet(tmp_inst,tmp_ori,rot,i2,i6);
					ar_to_vec(tmp_inst, inst_ori_mat, st_index, st_index + 9);
					i3 = st_index;
					i4 = 144*i1 + 36*i6 + 9*i2;
					for (i5 = 0; i5 < 9; i5++) {
						inst_ori_mat[i4].set_val(inst_ori_mat[i3]);
						i3++;
						i4++;
					}
				}
			}
		}
	}

	return;
}

void Element::get_inst_disp(DiffDoub1 inst_disp[], vector<DiffDoub1>& glob_disp, vector<DiffDoub1>& inst_ori_mat, vector<DiffDoub1>& loc_ori, vector<DiffDoub1>& x_glob, bool n_lgeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int nd_oind;
	int dof;
	int dof_oind;
	int nd2;
	int dof2;
	int dof2_oind;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 nn_inv;
	DiffDoub1 nn_inv2;
	
	i2 = 6*n_dim;
	for (i1 = 0; i1 < i2; i1++) {
		inst_disp[i1].set_val(0.0);
	}
	
	if (!n_lgeom) {
		if (dv1 == max_int && dv2 == max_int) {
			i7 = 3 * n_dim;
			for (i1 = 0; i1 < 3; i1++) {
				i4 = i1 * n_dim;
				for (i2 = 0; i2 < num_nds; i2++) {
					i5 = i1 * 3;
					i6 = i2;
					for (i3 = 0; i3 < 3; i3++) {
						//i4 = i1 * n_dim + i2;
						//i5 = i1 * 3 + i3;
						//i6 = i3 * n_dim + i2;
						tmp.set_val(loc_ori[i5]);
						tmp.mult(glob_disp[i6]);
						inst_disp[i4].add(tmp);
						tmp.set_val(loc_ori[i5]);
						tmp.mult(glob_disp[i6 + i7]);
						inst_disp[i4 + i7].add(tmp);
						i5++;
						i6 += n_dim;
					}
					i4++;
				}
			}
			i2 = 2 * num_nds * dof_per_nd;
			for (i1 = 0; i1 < num_int_dof; i1++) {
				nd = dof_table[i2];
				dof = dof_table[i2 + 1];
				i3 = dof * n_dim + nd;
				inst_disp[i3].set_val(glob_disp[i3]);
				i2 += 2;
			}
		}
		else if ((dv1 + dv2) >= max_int) {
			dv1 = dv1 + dv2 - max_int;
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < num_nds) {
					i2 = nd;
					i3 = dof;
					for (i1 = 0; i1 < 3; i1++) {
						//i2 = i1 * n_dim + nd;
						//i3 = i1 * 3 + dof;
						inst_disp[i2].set_val(loc_ori[i3]);
						i2 += n_dim;
						i3 += 3;
					}
				}
				else {
					i1 = dof * n_dim + nd;
					inst_disp[i1].set_val(1.0);
				}
			}
			else {
				i2 = 3 * n_dim + nd;
				i3 = dof - 3;
				for (i1 = 0; i1 < 3; i1++) {
					//i2 = (i1 + 3) * n_dim + nd;
					//i3 = i1 * 3 + (dof - 3);
					inst_disp[i2].set_val(loc_ori[i3]);
					i2 += n_dim;
					i3 += 3;
				}
			}
		}
	}
	else {
		if (dv1 == max_int && dv2 == max_int) {
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < num_nds; i2++) {
					i4 = i1 * n_dim + i2;
					i5 = i2;
					i6 = i2;
					i7 = i1 * 3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.set_val(glob_disp[i5]);
						tmp.add(x_glob[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp2.set_val(loc_ori[i7]);
						tmp2.mult(x_glob[i6]);
						tmp.sub(tmp2);
						inst_disp[i4].add(tmp);
						i5 += n_dim;
						i6 += num_nds;
						i7++;
					}
				}
			}

			i2 = num_nds * dof_per_nd;
			i3 = i2 + num_int_dof;
			for (i1 = i2; i1 < i3; i1++) {
				i4 = dof_table[2 * i1];
				i5 = dof_table[2 * i1 + 1];
				i6 = i5 * n_dim + i4;
				inst_disp[i6].set_val(glob_disp[i6]);
			}

			for (i1 = 0; i1 < num_nds; i1++) {
				i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
				i4 = 4 * n_dim + i1;
				i5 = 5 * n_dim + i1;
				nd_oind = 144 * (i1 + 1);
				for (i2 = 0; i2 < 3; i2++) {
					i6 = nd_oind + 3 + i2;
					i7 = 6 + i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i3].add(tmp);
					i6 = nd_oind + 6 + i2;
					i7 = i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i4].add(tmp);
					i6 = nd_oind + i2;
					i7 = 3 + i2;
					tmp.set_val(inst_ori_mat[i6]);
					tmp.mult(inst_ori_mat[i7]);
					inst_disp[i5].add(tmp);
				}
			}
		}
		else if ((dv1 + dv2) >= max_int) {
			dv1 = dv1 + dv2 - max_int;
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < num_nds) {
					i1 = nd; // index in inst_disp
					i2 = dof; //index in inst_ori_mat
					inst_disp[i1].set_val(inst_ori_mat[i2]);
					i1 = n_dim + nd;
					i2 = 3 + dof;
					inst_disp[i1].set_val(inst_ori_mat[i2]);
					i1 = 2 * n_dim + nd;
					i2 = 6 + dof;
					inst_disp[i1].set_val(inst_ori_mat[i2]);
				}
				else {
					i1 = dof * n_dim + nd;
					inst_disp[i1].set_val(1.0);
				}
			}
			else { // dof is rotation
				nn_inv.set_val(1.0 / num_nds);
				dof_oind = 36 * (dof - 2);
				for (i1 = 0; i1 < 3; i1++) {
					for (i2 = 0; i2 < num_nds; i2++) {
						i4 = i1 * n_dim + i2;
						i5 = i2;
						i6 = i2;
						i7 = dof_oind + i1 * 3;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.set_val(glob_disp[i5]);
							tmp.add(x_glob[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i5 += n_dim;
							i6 += num_nds;
							i7++;
						}
					}
				}

				for (i1 = 0; i1 < num_nds; i1++) {
					i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
					i4 = 4 * n_dim + i1;
					i5 = 5 * n_dim + i1;
					nd_oind = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = nd_oind + 3 + i2;
						i7 = dof_oind + 6 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i3].add(tmp);
						i6 = nd_oind + 6 + i2;
						i7 = dof_oind + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i4].add(tmp);
						i6 = nd_oind + i2;
						i7 = dof_oind + 3 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv);
						inst_disp[i5].add(tmp);
						if (i1 == nd) {
							i6 = nd_oind + dof_oind + 3 + i2;
							i7 = 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + 6 + i2;
							i7 = i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + i2;
							i7 = 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i5].add(tmp);
						}
					}
				}
			}
		}
		else {
			nd = dof_table[2 * dv1];
			dof = dof_table[2 * dv1 + 1];
			nd2 = dof_table[2 * dv2];
			dof2 = dof_table[2 * dv2 + 1];
			nn_inv.set_val(1.0 / num_nds);
			nn_inv2.set_val(nn_inv);
			nn_inv2.sqr();
			if (dof > 2 && dof2 > 2) {
				for (i1 = 0; i1 < 3; i1++) {
					dof_oind = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
					for (i2 = 0; i2 < num_nds; i2++) {
						i4 = n_dim * i1 + i2;
						i5 = i2;
						i6 = i2;
						i7 = dof_oind;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.set_val(glob_disp[i5]);
							tmp.add(x_glob[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv2);
							inst_disp[i4].add(tmp);
							i5 += n_dim;
							i6 += num_nds;
							i7++;
						}
					}
				}

				dof_oind = 36 * (dof - 2);
				dof2_oind = 9 * (dof2 - 2);
				for (i1 = 0; i1 < num_nds; i1++) {
					i3 = 3 * n_dim + i1; // indexes of thetax, y and z for Node i1
					i4 = 4 * n_dim + i1;
					i5 = 5 * n_dim + i1;
					nd_oind = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = nd_oind + 3 + i2;
						i7 = dof_oind + dof2_oind + 6 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i3].add(tmp);
						i6 = nd_oind + 6 + i2;
						i7 = dof_oind + dof2_oind + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i4].add(tmp);
						i6 = nd_oind + i2;
						i7 = dof_oind + dof2_oind + 3 + i2;
						tmp.set_val(inst_ori_mat[i6]);
						tmp.mult(inst_ori_mat[i7]);
						tmp.mult(nn_inv2);
						inst_disp[i5].add(tmp);
						if (i1 == nd) {
							i6 = nd_oind + dof_oind + 3 + i2;
							i7 = dof2_oind + 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + 6 + i2;
							i7 = dof2_oind + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + i2;
							i7 = dof2_oind + 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i5].add(tmp);
						}
						if (i1 == nd2) {
							i6 = nd_oind + dof2_oind + 3 + i2;
							i7 = dof_oind + 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof2_oind + 6 + i2;
							i7 = dof_oind + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof2_oind + i2;
							i7 = dof_oind + 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							tmp.mult(nn_inv);
							inst_disp[i5].add(tmp);
						}
						if (i1 == nd && i1 == nd2) {
							i6 = nd_oind + dof_oind + dof2_oind + 3 + i2;
							i7 = 6 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i3].add(tmp);
							i6 = nd_oind + dof_oind + dof2_oind + 6 + i2;
							i7 = i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i4].add(tmp);
							i6 = nd_oind + dof_oind + dof2_oind + i2;
							i7 = 3 + i2;
							tmp.set_val(inst_ori_mat[i6]);
							tmp.mult(inst_ori_mat[i7]);
							inst_disp[i5].add(tmp);
						}
					}
				}
			}
			else if (dof < 3 && dof2 < 3) {
				return;
			}
			else {
				if (dof > 2) {
					i1 = dof;
					dof = dof2;
					dof2 = i1;
					i1 = nd;
					nd = nd2;
					nd2 = i1;
				}
				i1 = 36 * (dof2 - 2) + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = nd;
				inst_disp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 3 + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = n_dim + nd;
				inst_disp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 6 + dof;
				tmp.set_val(inst_ori_mat[i1]);
				tmp.mult(nn_inv);
				i2 = 2 * n_dim + nd;
				inst_disp[i2].add(tmp);
			}
		}
	}

	
	return;
}

void Element::get_stress_prereq(DiffDoub1StressPrereq& pre, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int num_lay;
	DiffDoub1 offset;
	get_nd_crds(pre.glob_nds, nd_ar, dv_ar);
	get_loc_ori(pre.loc_ori, sec_ar, dv_ar);
	get_nd_disp(pre.glob_disp, nd_ar);
	get_nd_vel(pre.glob_vel, nd_ar);
	get_nd_acc(pre.glob_acc, nd_ar);
	get_nd_temp(pre.glob_temp, nd_ar);
	get_nd_tdot(pre.glob_tdot, nd_ar);
	if (dof_per_nd == 6) {
		correct_orient(pre.loc_ori, pre.glob_nds);
		if (type != 2) {
			get_layer_thk_z(pre.layer_thk, pre.layer_z, offset, sec_ar, dv_ar);
			get_layer_angle(pre.layer_ang, sec_ar, dv_ar);
			get_layer_q(pre.layer_q, sec_ar, mat_ar, dv_ar);
			get_layer_d(pre.layer_d, sec_ar, mat_ar, dv_ar);
			get_layer_th_exp(pre.layer_te, sec_ar, mat_ar, dv_ar);
			get_layer_einit(pre.layer_e0, sec_ar, dv_ar);
			get_layer_den(pre.layer_den, sec_ar, mat_ar, dv_ar);
			get_layer_cond(pre.layer_tc, sec_ar, mat_ar, dv_ar);
			get_layer_spec_heat(pre.layer_sh, sec_ar, mat_ar, dv_ar);
			get_abd(pre.cmat, pre.layer_thk, pre.layer_z, pre.layer_q, pre.layer_ang, sec_ar);
			get_shell_damp(pre.dmat, pre.layer_thk, pre.layer_z, pre.layer_d, pre.layer_ang, sec_ar);
			get_shell_exp_load(pre.therm_exp, pre.einit, pre.layer_thk, pre.layer_z, pre.layer_q, pre.layer_te, pre.layer_e0, pre.layer_ang, sec_ar);
			get_shell_mass(pre.mmat, pre.layer_thk, pre.layer_z, pre.layer_den, sec_ar, dv_ar);
			get_shell_cond(pre.tcmat, pre.layer_thk, pre.layer_ang, pre.layer_tc, sec_ar, dv_ar);
			get_shell_spec_heat(pre.spec_heat, pre.layer_thk, pre.layer_sh, pre.layer_den, sec_ar);
		}
		else {
			get_beam_stiff(pre.cmat, sec_ar, mat_ar, dv_ar);
			get_beam_damp(pre.dmat, sec_ar, mat_ar, dv_ar);
			get_beam_exp_load(pre.therm_exp, pre.einit, sec_ar, mat_ar, dv_ar);
			get_beam_mass(pre.mmat, sec_ar, mat_ar, dv_ar);
			get_beam_cond(pre.tcmat, sec_ar, mat_ar, dv_ar);
			get_beam_spec_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
		}
	}
	else if (type == 21) {
		get_frc_fld_const(pre.frc_fld_coef, pre.frc_fld_exp, sec_ar, dv_ar);
		get_thrm_fld_const(pre.thrm_fld_coef, pre.ref_temp, sec_ar, dv_ar);
	}
	else if (type == 1) {
		get_mass_per_el(pre.mass_per_el, sec_ar, dv_ar);
		get_specific_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
	}
	else {
		get_solid_stiff(pre.cmat, sec_ar, mat_ar, dv_ar);
		get_solid_damp(pre.dmat, sec_ar, mat_ar, dv_ar);
		get_thermal_exp(pre.therm_exp, pre.einit, sec_ar, mat_ar, dv_ar);
		get_density(pre.mmat[0], 0, sec_ar, mat_ar, dv_ar);
		get_conductivity(pre.tcmat, sec_ar, mat_ar, dv_ar);
		get_specific_heat(pre.spec_heat, sec_ar, mat_ar, dv_ar);
	}
	mat_mul(pre.loc_nds, pre.loc_ori, pre.glob_nds, 3, 3, num_nds);


	return;
}

void Element::get_fluid_prereq(DiffDoub1FlPrereq& pre, vector<Section>& sec_ar, vector<Fluid>& fl_ar, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	double sec_prop;
	
	get_nd_crds(pre.glob_nds, nd_ar, dv_ar);
	get_nd_disp(pre.glob_disp, nd_ar);
	i2 = 3 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		pre.glob_nds[i1].add(pre.glob_disp[i1]);
	}
	get_nd_vel(pre.glob_vel, nd_ar);
	get_nd_fl_den(pre.fl_den, nd_ar);
	get_nd_fl_vel(pre.fl_vel, nd_ar);
	get_nd_temp(pre.fl_temp, nd_ar);
	get_nd_turb_e(pre.fl_turb_e, nd_ar);
	get_nd_fl_den_dot(pre.fl_den_dot, nd_ar);
	get_nd_fl_vdot(pre.fl_vel_dot, nd_ar);
	get_nd_tdot(pre.fl_tdot, nd_ar);
	get_nd_turb_edot(pre.fl_turb_edot, nd_ar);

	Section& this_sec = sec_ar[sect_ptr];
	Fluid& this_fl = fl_ar[this_sec.fl_ptr];

	sec_prop = this_sec.ref_turb_e;
	pre.ref_turb_e.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "refTurbE", dv_ar);

	sec_prop = this_sec.grad_vturb_coef;
	pre.grad_vturb_coef.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "gradVTurbCoef", dv_ar);

	sec_prop = this_sec.diss_turb_coef;
	pre.diss_turb_coef.set_val(sec_prop);
	get_gen_prop(pre.ref_turb_e, "dissTurbCoef", dv_ar);

	sec_prop = this_fl.viscosity;
	pre.ref_visc.set_val(sec_prop);
	get_gen_prop(pre.ref_visc, "viscosity", dv_ar);

	sec_prop = this_sec.den_vis_coef;
	pre.den_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.den_vis_coef, "denVisCoef", dv_ar);

	sec_prop = this_sec.temp_vis_coef;
	pre.temp_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.temp_vis_coef, "tempVisCoef", dv_ar);

	sec_prop = this_sec.turb_vis_coef;
	pre.turb_vis_coef.set_val(sec_prop);
	get_gen_prop(pre.turb_vis_coef, "turbVisCoef", dv_ar);

	sec_prop = this_sec.ref_enth;
	pre.ref_enth.set_val(sec_prop);
	get_gen_prop(pre.ref_enth, "refEnth", dv_ar);

	sec_prop = this_sec.enth_coef;
	pre.den_enth_coef.set_val(sec_prop);
	get_gen_prop(pre.den_enth_coef, "denEnthCoef", dv_ar);

	sec_prop = this_sec.enth_exp;
	pre.den_enth_exp.set_val(sec_prop);
	get_gen_prop(pre.den_enth_exp, "denEnthExp", dv_ar);

	sec_prop = this_sec.pres_coef;
	pre.den_pres_coef.set_val(sec_prop);
	get_gen_prop(pre.den_pres_coef, "denPresCoef", dv_ar);

	sec_prop = this_sec.pres_exp;
	pre.den_pres_exp.set_val(sec_prop);
	get_gen_prop(pre.den_pres_exp, "denPresExp", dv_ar);

	sec_prop = this_sec.ref_den;
	pre.ref_den.set_val(sec_prop);
	get_gen_prop(pre.ref_den, "refDen", dv_ar);

	sec_prop = this_sec.ref_temp;
	pre.ref_temp.set_val(sec_prop);
	get_gen_prop(pre.ref_temp, "refTemp", dv_ar);

	sec_prop = this_fl.therm_cond;
	pre.therm_cond.set_val(sec_prop);
	get_gen_prop(pre.therm_cond, "thermCond", dv_ar);

	sec_prop = this_fl.spec_heat;
	pre.spec_heat.set_val(sec_prop);
	get_gen_prop(pre.spec_heat, "specHeat", dv_ar);

	sec_prop = this_fl.ideal_gas;
	pre.i_gconst.set_val(sec_prop);
	get_gen_prop(pre.i_gconst, "iGConst", dv_ar);

	return;
}

void Element::get_volume(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 thk;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;
    double s_tmp[3];

	if (type == 2) {
		thk.set_val(sec_ar[sect_ptr].area); //rem: update to factor in dvars;
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value(dv_val);
				dv_val.mult(tmp);
				thk.add(dv_val);
			}
		}
	} else if (type == 3 || type == 41) {
		thk.set_val(pre.layer_thk[layer]);
	} else {
		thk.set_val(1.0);
	}

	vol.set_val(0.0);
	for (i1 = 0; i1 < num_ip; i1++) {
		i2 = 3 * i1;
		vec_to_ar(s_tmp, int_pts, i2, i2 + 3);
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, s_tmp);
		tmp.set_val(ip_wt[i1]);
		tmp.mult(det_j);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::get_section_def(DiffDoub1 sec_def[], vector<DiffDoub1>& glob_disp,  vector<DiffDoub1>& inst_ori_mat, vector<DiffDoub1>& loc_ori, vector<DiffDoub1>& x_glob, DiffDoub1 d_ndx[], DiffDoub1 n_vec[], bool n_lgeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub1 inst_disp[60];
	DiffDoub1 ux[9];
	DiffDoub1 rx[9];
	DiffDoub1 rot[3];
	DiffDoub1 tmp;

	if (dof_per_nd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			sec_def[i1].set_val(0.0);
		}
		return;
	}
	
	get_inst_disp(inst_disp, glob_disp, inst_ori_mat, loc_ori, x_glob, n_lgeom, dv1, dv2);
	
	if(dv1 == max_int && dv2 == max_int) {
		mat_mul(ux,inst_disp,d_ndx,3,n_dim,3);
		i1 = 3*n_dim;
		mat_mul(rx,&inst_disp[i1],d_ndx,3,n_dim,3);
		mat_mul(rot,&inst_disp[i1],n_vec,3,n_dim,1);
	} else if((dv1 + dv2) >= max_int) {
		dv1 = dv1 + dv2 - max_int;
		nd = dof_table[2*dv1];
		dof = dof_table[2*dv1+1];
		if(dof < 3) {
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].set_val(inst_disp[i4]);
					ux[i3].mult(d_ndx[i5]);
					rx[i3].set_val(0.0);
					i3++;
					i5++;
				}
				rot[i1].set_val(0.0);
				i4+= n_dim;
			}
		} else {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = n_dim*i1;
					i6 = i2;
					i7 = n_dim*(i1+3);
					ux[i4].set_val(0.0);
					rx[i4].set_val(0.0);
					for (i3 = 0; i3 < num_nds; i3++) {
						tmp.set_val(inst_disp[i5]);
						tmp.mult(d_ndx[i6]);
						ux[i4].add(tmp);
						tmp.set_val(inst_disp[i7]);
						tmp.mult(d_ndx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].set_val(0.0);
				i3 = n_dim*(i1+3);
				for (i2 = 0; i2 < num_nds; i2++) {
					tmp.set_val(inst_disp[i3]);
					tmp.mult(n_vec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		}
	} else {
		nd = dof_table[2*dv1];
		dof = dof_table[2*dv1+1];
		nd2 = dof_table[2*dv2];
		dof2 = dof_table[2*dv2+1];
		if(dof < 3 && dof2 < 3) {
			for (i1 = 0; i1 < 9; i1++) {
				sec_def[i1].set_val(0.0);
			}
			return;
		} else if(dof > 2 && dof2 > 2) {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = n_dim*i1;
					i6 = i2;
					i7 = n_dim*(i1+3);
					ux[i4].set_val(0.0);
					for (i3 = 0; i3 < num_nds; i3++) {
						tmp.set_val(inst_disp[i5]);
						tmp.mult(d_ndx[i6]);
						ux[i4].add(tmp);
						tmp.set_val(inst_disp[i7]);
						tmp.mult(d_ndx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].set_val(0.0);
				i3 = n_dim*(i1+3);
				for (i2 = 0; i2 < num_nds; i2++) {
					tmp.set_val(inst_disp[i3]);
					tmp.mult(n_vec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		} else {
			if(dof > dof2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].set_val(inst_disp[i4]);
					ux[i3].mult(d_ndx[i5]);
					rx[i3].set_val(0.0);
					i3++;
					i5++;
				}
				rot[i1].set_val(0.0);
				i4+= n_dim;
			}
		}
	}
	
	if(type == 2) {
		sec_def[0].set_val(ux[0]);
		sec_def[1].set_val(ux[3]);
		sec_def[1].sub(rot[2]);
		sec_def[2].set_val(ux[6]);
		sec_def[2].add(rot[1]);
		sec_def[3].set_val(rx[0]);
		sec_def[4].set_val(rx[3]);
		sec_def[5].set_val(rx[6]);
	} else {
		sec_def[0].set_val(ux[0]);
		sec_def[1].set_val(ux[4]);
		sec_def[2].set_val(ux[1]);
		sec_def[2].add(ux[3]);
		sec_def[3].set_val(rx[3]);
		sec_def[4].set_val(rx[1]);
		sec_def[4].neg();
		sec_def[5].set_val(rx[4]);
		sec_def[5].sub(rx[0]);
		sec_def[6].set_val(ux[6]);
		sec_def[6].add(rot[1]);
		sec_def[7].set_val(ux[7]);
		sec_def[7].sub(rot[0]);
		sec_def[8].set_val(2.0);
		sec_def[8].mult(rot[2]);
		sec_def[8].sub(ux[3]);
		sec_def[8].add(ux[1]);
	}
	
	return;
}

void Element::get_solid_strain(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 d_ndx[], vector<DiffDoub1>& loc_ori, int dv1, int dv2, bool n_lgeom) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub1 ux_l[9];
	DiffDoub1 strn_mat[9];
	DiffDoub1 tmp_ori[9];
	DiffDoub1 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strn_mat[i1].set_val(0.0);
	}
	if(dv1 == max_int && dv2 == max_int) {
		vec_to_ar(tmp_ori, loc_ori, 0, 9);
		mat_mul(ux_l,tmp_ori,ux,3,3,3);
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			i5 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				strn_mat[i4].add(ux_l[i4]);
				strn_mat[i4].add(ux_l[i5]);
				if(n_lgeom) {
					i6 = i1;
					i7 = i2;
				    for (i3 = 0; i3 < 3; i3++) {
						tmp.set_val(ux[i6]);
						tmp.mult(ux[i7]);
						strn_mat[i4].add(tmp);
						i6+= 3;
						i7+= 3;
				    }
				}
				i4++;
				i5+= 3;
			}
		}
	} else if((dv1 + dv2) >= max_int) {
		if(dv1 < max_int) {
			nd = dof_table[2*dv1];
			dof = dof_table[2*dv1+1];
		} else {
			nd = dof_table[2*dv2];
			dof = dof_table[2*dv2+1];
		}
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				i5 = 3*i1 + dof;
				i6 = 3*nd + i2;
				tmp.set_val(loc_ori[i5]);
				tmp.mult(d_ndx[i6]);
				strn_mat[i4].add(tmp);
				i5 = 3*i2 + dof;
				i6 = 3*nd + i1;
				tmp.set_val(loc_ori[i5]);
				tmp.mult(d_ndx[i6]);
				strn_mat[i4].add(tmp);
				if(n_lgeom) {
					i5 = 3*nd + i1;
					i6 = 3*dof + i2;
					tmp.set_val(d_ndx[i5]);
					tmp.mult(ux[i6]);
					strn_mat[i4].add(tmp);
					i5 = 3*nd + i2;
					i6 = 3*dof + i1;
					tmp.set_val(d_ndx[i5]);
					tmp.mult(ux[i6]);
					strn_mat[i4].add(tmp);
				}
				i4++;
			}
		}
	} else {
		if(n_lgeom) {
			nd = dof_table[2*dv1];
			dof = dof_table[2*dv1+1];
			nd2 = dof_table[2*dv2];
			dof2 = dof_table[2*dv2+1];
			if(dof == dof2) {
				for (i1 = 0; i1 < 3; i1++) {
					i4 = 4*i1;
					for (i2 = i1; i2 < 3; i2++) {
						i5 = 3*nd + i1;
						i6 = 3*nd2 + i2;
						tmp.set_val(d_ndx[i5]);
						tmp.mult(d_ndx[i6]);
						strn_mat[i4].add(tmp);
						i5 = 3*nd2 + i1;
						i6 = 3*nd + i2;
						tmp.set_val(d_ndx[i5]);
						tmp.mult(d_ndx[i6]);
						strn_mat[i4].add(tmp);
						i4++;
					}
				}
			}
		}
	}
	
	tmp.set_val(0.5);
	strain[0].set_val(strn_mat[0]);
	strain[0].mult(tmp);
	strain[1].set_val(strn_mat[4]);
	strain[1].mult(tmp);
	strain[2].set_val(strn_mat[8]);
	strain[2].mult(tmp);
	strain[3].set_val(strn_mat[1]);
	strain[4].set_val(strn_mat[2]);
	strain[5].set_val(strn_mat[5]);
	
	return;
}

void Element::get_stress_strain(DiffDoub1 stress[], DiffDoub1 strain[], double spt[], int layer, bool n_lgeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 ux[9];
	DiffDoub1 sec_def[9];
	DiffDoub1 sect_strn[3];
	DiffDoub1 adj_stn[6];
	DiffDoub1 tmp;
	DiffDoub1 ip_temp;
	DiffDoub1 tmp_ar[60];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);

	ip_temp.set_val(0.0);
	for (i1 = 0; i1 < num_nds; i1++) {
		tmp.set_val(pre.glob_temp[i1]);
		tmp.mult(n_vec[i1]);
		ip_temp.add(tmp);
	}

	if (type == 41 || type == 3) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}

		get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, max_int, max_int);
		
		sect_strn[0].set_val(sec_def[0]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[3]);
		sect_strn[0].add(tmp);
		
		sect_strn[1].set_val(sec_def[1]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[4]);
		sect_strn[1].sub(tmp);

		sect_strn[2].set_val(sec_def[2]);
		tmp.set_val(pre.layer_z[layer]);
		tmp.mult(sec_def[5]);
		sect_strn[2].add(tmp);

		tmp.set_val(pre.layer_ang[layer]);
		tmp.neg();
		transform_strain(strain, sect_strn, tmp);
		i2 = 3 * layer;
		for (i1 = 0; i1 < 3; i1++) {
			tmp.set_val(pre.layer_te[i2]);
			tmp.mult(ip_temp);
			adj_stn[i1].set_val(strain[i1]);
			adj_stn[i1].sub(tmp);
			adj_stn[i1].sub(pre.layer_e0[i2]);
			i2++;
		}
		mat_mul(stress, &pre.layer_q[9 * layer], adj_stn, 3, 3, 1);

		strain[3].set_val(strain[2]);
		strain[2].set_val(0.0);
		strain[4].set_val(0.0);
		strain[5].set_val(0.0);

		stress[3].set_val(stress[2]);
		stress[2].set_val(0.0);
		stress[4].set_val(0.0);
		stress[5].set_val(0.0);

	} else if(type != 2) {
		vec_to_ar(tmp_ar, pre.glob_disp, 0, 60);
		mat_mul(ux, tmp_ar, d_ndx, 3, n_dim, 3);
		get_solid_strain(strain, ux, d_ndx, pre.loc_ori, max_int, max_int, n_lgeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.set_val(pre.therm_exp[i1]);
			tmp.mult(ip_temp);
			adj_stn[i1].set_val(strain[i1]);
			adj_stn[i1].sub(tmp);
			adj_stn[i1].sub(pre.einit[i1]);
		}
		vec_to_ar(tmp_ar, pre.cmat, 0, 36);
		mat_mul(stress, tmp_ar, adj_stn, 6, 6, 1);
	}
	return;
}

void Element::d_stress_straind_u(vector<DiffDoub1>& dsd_u, vector<DiffDoub1>& ded_u, vector<DiffDoub1>& dsd_t, double spt[], int layer, bool n_lgeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 ux[9];
	DiffDoub1 sec_def[9];
	DiffDoub1 sect_strn[3];
	DiffDoub1 d_strain[6];
	DiffDoub1 d_stress[6];
	DiffDoub1 cte[6];
	DiffDoub1 cten[60];
	DiffDoub1 tmp;
	DiffDoub1 tmp_ar[60];
	DiffDoub1 tmp_ar2[6];

	int tot_dof = num_nds * dof_per_nd + num_int_dof;
	i2 = 6 * tot_dof;
	for (i1 = 0; i1 < i2; i1++) {
		dsd_u[i1].set_val(0.0);
		ded_u[i1].set_val(0.0);
	}
	i2 = 6 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		dsd_t[i1].set_val(0.0);
	}

	if (type == 41 || type == 3) {
		if (n_lgeom) {
			get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
		}

		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
		for (i1 = 0; i1 < tot_dof; i1++) {
			get_section_def(sec_def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, i1, max_int);

			sect_strn[0].set_val(sec_def[0]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[3]);
			sect_strn[0].add(tmp);

			sect_strn[1].set_val(sec_def[1]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[4]);
			sect_strn[1].sub(tmp);

			sect_strn[2].set_val(sec_def[2]);
			tmp.set_val(pre.layer_z[layer]);
			tmp.mult(sec_def[5]);
			sect_strn[2].add(tmp);

			tmp.set_val(pre.layer_ang[layer]);
			tmp.neg();
			transform_strain(d_strain, sect_strn, tmp);
			mat_mul(d_stress, &pre.layer_q[9 * layer], d_strain, 3, 3, 1);

			ded_u[i1].set_val(d_strain[0]);
			ded_u[i1 + tot_dof].set_val(d_strain[1]);
			ded_u[i1 + 3 * tot_dof].set_val(d_strain[2]);

			dsd_u[i1].set_val(d_stress[0]);
			dsd_u[tot_dof + i1].set_val(d_stress[1]);
			dsd_u[3 * tot_dof + i1].set_val(d_stress[2]);
		}
		mat_mul(cte, &pre.layer_q[9 * layer], &pre.layer_te[3 * layer], 3, 3, 1);
		cte[0].neg();
		cte[1].neg();
		cte[2].neg();
		mat_mul(cten, cte, n_vec, 3, 1, num_nds);
		for (i1 = 0; i1 < num_nds; i1++) {
			dsd_t[i1].set_val(cten[i1]);
			dsd_t[i1 + num_nds].set_val(cten[i1 + num_nds]);
			dsd_t[i1 + 3 * num_nds].set_val(cten[i1 + 2 * num_nds]);
		}
	}
	else if (type != 2) {
		get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
		vec_to_ar(tmp_ar, pre.glob_disp, 0, 60);
		mat_mul(ux, tmp_ar, d_ndx, 3, n_dim, 3);
		vec_to_ar(tmp_ar, pre.cmat, 0, 36);
		for (i1 = 0; i1 < tot_dof; i1++) {
			get_solid_strain(d_strain, ux, d_ndx, pre.loc_ori, i1, max_int, n_lgeom);
			mat_mul(d_stress, tmp_ar, d_strain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				ded_u[i3].set_val(d_strain[i2]);
				dsd_u[i3].set_val(d_stress[i2]);
				i3 += tot_dof;
			}
		}
		vec_to_ar(tmp_ar2, pre.therm_exp, 0, 6);
		mat_mul(cte, tmp_ar, tmp_ar2, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			cte[i1].neg();
		}
		mat_mul(tmp_ar, cte, n_vec, 6, 1, num_nds);
		ar_to_vec(tmp_ar, dsd_t, 0, 60);
	}

	return;
}

void Element::get_def_frc_mom(DiffDoub1 def[], DiffDoub1 frc_mom[], double spt[], bool n_lgeom, DiffDoub1StressPrereq& pre) {
	int i1;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 pt_temp;
	DiffDoub1 tmp;
	DiffDoub1 tmp_ar[36];

	if (dof_per_nd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			def[i1].set_val(0.0);
		}
		return;
	}

	if (n_lgeom) {
		get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
	}

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	get_section_def(def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, max_int, max_int);

	pt_temp.set_val(0.0);
	for (i1 = 0; i1 < num_nds; i1++) {
		tmp.set_val(pre.glob_temp[i1]);
		tmp.mult(n_vec[i1]);
		pt_temp.add(tmp);
	}

	vec_to_ar(tmp_ar, pre.cmat, 0, 36);
	mat_mul(frc_mom, tmp_ar, def, def_dim, def_dim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frc_mom[i1].sub(pre.einit[i1]);
		tmp.set_val(pt_temp);
		tmp.mult(pre.therm_exp[i1]);
		frc_mom[i1].sub(tmp);
	}

	return;
}

void Element::d_def_frc_momd_u(vector<DiffDoub1>& d_defd_u, vector<DiffDoub1>& d_frc_momd_u, vector<DiffDoub1>& d_frc_momd_t, double spt[], bool n_lgeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	int tot_dof;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 pt_temp;
	DiffDoub1 tmp;
	DiffDoub1 def[9];
	DiffDoub1 tmp_ar[60];
	DiffDoub1 tmp_ar2[6];

	if (n_lgeom) {
		get_inst_ori(pre.inst_ori, pre.loc_ori, pre.glob_disp, 2);
	}


	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	tot_dof = num_nds * dof_per_nd + num_int_dof;
	for (i1 = 0; i1 < tot_dof; i1++) {
		get_section_def(def, pre.glob_disp, pre.inst_ori, pre.loc_ori, pre.glob_nds, d_ndx, n_vec, n_lgeom, i1, max_int);
		i2 = i1;
		for (i3 = 0; i3 < def_dim; i3++) {
			d_defd_u[i2].set_val(def[i3]);
			i2 += tot_dof;
		}
	}

	mat_mul(d_frc_momd_u, pre.cmat, d_defd_u, def_dim, def_dim, tot_dof);

	vec_to_ar(tmp_ar2, pre.therm_exp, 0, 6);
	mat_mul(tmp_ar, tmp_ar2, n_vec, 6, 1, num_nds);
	ar_to_vec(tmp_ar, d_frc_momd_t, 0, 60);

	return;
}

void Element::get_flux_tgrad(DiffDoub1 flux[], DiffDoub1 t_grad[], double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 tmp_ar[10];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	vec_to_ar(tmp_ar, pre.glob_temp, 0, 10);
	mat_mul(t_grad, tmp_ar, d_ndx, 1, num_nds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vec_to_ar(tmp_ar, pre.layer_tc, i1, i1 + 9);
		mat_mul(flux, tmp_ar, t_grad, 3, 3, 1);
	}
	else {
		vec_to_ar(tmp_ar, pre.tcmat, 0, 9);
		mat_mul(flux, tmp_ar, t_grad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::d_flux_tgradd_t(vector<DiffDoub1>& d_fd_t, vector<DiffDoub1>& d_tg, double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub1 n_vec[11];
	DiffDoub1 d_ndx[33];
	DiffDoub1 det_j;
	DiffDoub1 tmp_ar1[33];
	DiffDoub1 tmp_ar2[9];
	DiffDoub1 tmp_ar3[30];

	get_ip_data(n_vec, d_ndx, det_j, pre.loc_nds, spt);
	transpose(tmp_ar1, d_ndx, num_nds, 3);
	ar_to_vec(tmp_ar1, d_tg, 0, 33);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vec_to_ar(tmp_ar2, pre.layer_tc, i1, i1 + 9);
		mat_mul(tmp_ar3, tmp_ar2, tmp_ar1, 3, 3, num_nds);
		ar_to_vec(tmp_ar3, d_fd_t, 0, 30);
	}
	else {
		mat_mul(d_fd_t, pre.tcmat, d_tg, 3, 3, num_nds);
	}
	i2 = 3 * num_nds;
	for (i1 = 0; i1 < i2; i1++) {
		d_fd_t[i1].neg();
	}

	return;
}

void Element::put_vec_to_glob_mat(SparseMat& q_mat, vector<DiffDoub1>& el_qvec, bool for_therm, int mat_row, vector<Node>& nd_ar) {
	int i1;
	int nd_dof = num_nds*dof_per_nd;
	int tot_dof = nd_dof + num_int_dof;
	int nd;
	int dof;
	int glob_ind;

	if (for_therm) {
		for (i1 = 0; i1 < num_nds; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[i1].sorted_rank;
			q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < tot_dof; i1++) {
			if (i1 < nd_dof) {
				nd = nodes[dof_table[2 * i1]];
				dof = dof_table[2 * i1 + 1];
				glob_ind = nd_ar[nd].dof_index[dof];
				q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
			}
			else {
				dof = i1 - nd_dof;
				glob_ind = int_dof_index + dof;
				q_mat.add_entry(mat_row, glob_ind, el_qvec[i1].val);
			}
		}
	}

	return;
}

//end dup
 
//end skip 
 
 
void Element::get_el_vec(vector<double>& el_vec, vector<double>& glob_vec, bool for_therm, bool intnl, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int nd;
	int dof;
	int glob_ind;
	int nd_dof;

	if (for_therm) {
		for (i1 = 0; i1 < num_nds; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[nd].sorted_rank;
			el_vec[i1] = glob_vec[glob_ind];
		}
	}
	else {
		nd_dof = num_nds * dof_per_nd;
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = nodes[dof_table[i2]];
			dof = dof_table[i2 + 1];
			glob_ind = nd_ar[nd].dof_index[dof];
			el_vec[i1] = glob_vec[glob_ind];
			i2 += 2;
		}
		if (intnl) {
			for (i1 = 0; i1 < num_int_dof; i1++) {
				el_vec[i1 + nd_dof] = glob_vec[i1 + int_dof_index];
			}
		}
	}

	return;
}

void Element::add_to_glob_vec(vector<double>& el_vec, vector<double>& glob_vec, bool for_therm, bool intnl, vector<Node>& nd_ar) {
	int i1;
	int i2;
	int nd;
	int dof;
	int glob_ind;
	int nd_dof;

	if (for_therm) {
		for (i1 = 0; i1 < num_nds; i1++) {
			nd = nodes[i1];
			glob_ind = nd_ar[nd].sorted_rank;
			glob_vec[glob_ind] += el_vec[i1];
		}
	}
	else {
		nd_dof = num_nds * dof_per_nd;
		i2 = 0;
		for (i1 = 0; i1 < nd_dof; i1++) {
			nd = nodes[dof_table[i2]];
			dof = dof_table[i2 + 1];
			glob_ind = nd_ar[nd].dof_index[dof];
			glob_vec[glob_ind] += el_vec[i1];
			i2 += 2;
		}
		if (intnl) {
			for (i1 = 0; i1 < num_int_dof; i1++) {
				glob_vec[i1 + int_dof_index] += el_vec[i1 + nd_dof];
			}
		}
	}

	return;
}
 
 
// 

