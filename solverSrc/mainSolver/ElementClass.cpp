#include <vector>
#include "ElementClass.h"
#include "constants.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "FaceClass.h"
#include "matrixFunctions.h"

using namespace std;

// stress prerequisite classes

//dup1
DiffDoub0StressPrereq::DiffDoub0StressPrereq() {
	glob_nds = vector<DiffDoub0>(30);
	loc_nds = vector<DiffDoub0>(30);
	loc_ori = vector<DiffDoub0>(9);
	inst_ori = vector<DiffDoub0>(720);
	glob_disp = vector<DiffDoub0>(60);
	glob_vel = vector<DiffDoub0>(30);
	glob_acc = vector<DiffDoub0>(30);
	glob_temp = vector<DiffDoub0>(10);
	glob_tdot = vector<DiffDoub0>(10);
	cmat = vector<DiffDoub0>(81);
	mmat = vector<DiffDoub0>(36);
	dmat = vector<DiffDoub0>(81);
	therm_exp = vector<DiffDoub0>(6);
	einit = vector<DiffDoub0>(6);
	tcmat = vector<DiffDoub0>(9);
	bmat = vector<DiffDoub0>(288);
	cbmat = vector<DiffDoub0>(288);
	frc_fld_coef = vector<DiffDoub0>(2);
	frc_fld_exp = vector<DiffDoub0>(2);
	thrm_fld_coef = vector<DiffDoub0>(2);
	scr_mat1 = vector<double>(3600);
	scr_mat2 = vector<double>(3600);
	scr_mat3 = vector<double>(3600);
	scr_mat4 = vector<double>(3600);
	scr_mat5 = vector<double>(3600);
	scr_vec1 = vector<DiffDoub0>(60);
	scr_vec2 = vector<DiffDoub0>(60);
	scr_vec3 = vector<DiffDoub0>(60);
	scr_vec4 = vector<DiffDoub0>(60);
	scr_vec5 = vector<DiffDoub0>(60);
	current_lay_len = 0;
	return;
}

void DiffDoub0StressPrereq::allocate_layers_dfd0(int num_layers) {
	if (num_layers != 0) {
		layer_z = vector<DiffDoub0>(num_layers);
		layer_thk = vector<DiffDoub0>(num_layers);
		layer_ang = vector<DiffDoub0>(num_layers);
		layer_q = vector<DiffDoub0>(9 * num_layers);
		layer_d = vector<DiffDoub0>(9 * num_layers);
		layer_te = vector<DiffDoub0>(3 * num_layers);
		layer_e0 = vector<DiffDoub0>(3 * num_layers);
		layer_den = vector<DiffDoub0>(num_layers);
		layer_tc = vector<DiffDoub0>(9 * num_layers);
		layer_sh = vector<DiffDoub0>(num_layers);
	}
	current_lay_len = num_layers;
	return;
}

DiffDoub0FlPrereq::DiffDoub0FlPrereq() {
	glob_nds = vector<DiffDoub0>(30);
	glob_disp = vector<DiffDoub0>(30);
	glob_vel = vector<DiffDoub0>(30);
	fl_den = vector<DiffDoub0>(10);
	fl_vel = vector<DiffDoub0>(30);
	fl_temp = vector<DiffDoub0>(10);
	fl_turb_e = vector<DiffDoub0>(10);
	fl_den_dot = vector<DiffDoub0>(10);
	fl_vel_dot = vector<DiffDoub0>(30);
	fl_tdot = vector<DiffDoub0>(10);
	fl_turb_edot = vector<DiffDoub0>(10);
	scratch = vector<double>(3600);

	return;
}

//end dup 
 
//skip 
 
//DiffDoub1 versions: 
//dup1
DiffDoub1StressPrereq::DiffDoub1StressPrereq() {
	glob_nds = vector<DiffDoub1>(30);
	loc_nds = vector<DiffDoub1>(30);
	loc_ori = vector<DiffDoub1>(9);
	inst_ori = vector<DiffDoub1>(720);
	glob_disp = vector<DiffDoub1>(60);
	glob_vel = vector<DiffDoub1>(30);
	glob_acc = vector<DiffDoub1>(30);
	glob_temp = vector<DiffDoub1>(10);
	glob_tdot = vector<DiffDoub1>(10);
	cmat = vector<DiffDoub1>(81);
	mmat = vector<DiffDoub1>(36);
	dmat = vector<DiffDoub1>(81);
	therm_exp = vector<DiffDoub1>(6);
	einit = vector<DiffDoub1>(6);
	tcmat = vector<DiffDoub1>(9);
	bmat = vector<DiffDoub1>(288);
	cbmat = vector<DiffDoub1>(288);
	frc_fld_coef = vector<DiffDoub1>(2);
	frc_fld_exp = vector<DiffDoub1>(2);
	thrm_fld_coef = vector<DiffDoub1>(2);
	scr_mat1 = vector<double>(3600);
	scr_mat2 = vector<double>(3600);
	scr_mat3 = vector<double>(3600);
	scr_mat4 = vector<double>(3600);
	scr_mat5 = vector<double>(3600);
	scr_vec1 = vector<DiffDoub1>(60);
	scr_vec2 = vector<DiffDoub1>(60);
	scr_vec3 = vector<DiffDoub1>(60);
	scr_vec4 = vector<DiffDoub1>(60);
	scr_vec5 = vector<DiffDoub1>(60);
	current_lay_len = 0;
	return;
}

void DiffDoub1StressPrereq::allocate_layers_dfd1(int num_layers) {
	if (num_layers != 0) {
		layer_z = vector<DiffDoub1>(num_layers);
		layer_thk = vector<DiffDoub1>(num_layers);
		layer_ang = vector<DiffDoub1>(num_layers);
		layer_q = vector<DiffDoub1>(9 * num_layers);
		layer_d = vector<DiffDoub1>(9 * num_layers);
		layer_te = vector<DiffDoub1>(3 * num_layers);
		layer_e0 = vector<DiffDoub1>(3 * num_layers);
		layer_den = vector<DiffDoub1>(num_layers);
		layer_tc = vector<DiffDoub1>(9 * num_layers);
		layer_sh = vector<DiffDoub1>(num_layers);
	}
	current_lay_len = num_layers;
	return;
}

DiffDoub1FlPrereq::DiffDoub1FlPrereq() {
	glob_nds = vector<DiffDoub1>(30);
	glob_disp = vector<DiffDoub1>(30);
	glob_vel = vector<DiffDoub1>(30);
	fl_den = vector<DiffDoub1>(10);
	fl_vel = vector<DiffDoub1>(30);
	fl_temp = vector<DiffDoub1>(10);
	fl_turb_e = vector<DiffDoub1>(10);
	fl_den_dot = vector<DiffDoub1>(10);
	fl_vel_dot = vector<DiffDoub1>(30);
	fl_tdot = vector<DiffDoub1>(10);
	fl_turb_edot = vector<DiffDoub1>(10);
	scratch = vector<double>(3600);

	return;
}

//end dup 
 
//end skip 
 
 
 
Element::Element() {
	return;
}

void Element::initialize_type(int new_type) {
	type = new_type;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;

	s_cent[0] = 0.0;
	s_cent[1] = 0.0;
	s_cent[2] = 0.0;
	if(type == 4 || type == 400) {
		num_nds = 4;
		dof_per_nd = 3;
		num_int_dof = 0;
		n_dim = 4;
		def_dim = 6;
		num_ip = 1;
		num_faces = 4;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(i1 * 2); //int[i1][2];
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = 0.25;
		int_pts[1] = 0.25;
		int_pts[2] = 0.25;
		ip_wt[0] = r_1o6;
		nd_spts = vector<double>(12);
		nd_spts[0] = 0.0;
		nd_spts[1] = 0.0;
		nd_spts[2] = 0.0;
		nd_spts[3] = 1.0;
		nd_spts[4] = 0.0;
		nd_spts[5] = 0.0;
		nd_spts[6] = 0.0;
		nd_spts[7] = 1.0;
		nd_spts[8] = 0.0;
		nd_spts[9] = 0.0;
		nd_spts[10] = 0.0;
		nd_spts[11] = 1.0;
	} else if(type == 6 || type == 600) {
		num_nds = 6;
		dof_per_nd = 3;
		num_int_dof = 0;
		n_dim = 6;
		def_dim = 6;
		num_ip = 2;
		num_faces = 5;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(i1 * 2);// int[i1][2];
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = r_1o3;
		int_pts[1] = r_1o3;
		int_pts[2] = -r_1ort3;
		int_pts[3] = r_1o3;
		int_pts[4] = r_1o3;
		int_pts[5] = r_1ort3;
		ip_wt[0] = 0.5;
		ip_wt[1] = 0.5;
		s_cent[0] = r_1o3;
		s_cent[1] = r_1o3;
		s_cent[2] = 0.0;
		nd_spts = vector<double>(18);
		nd_spts[0] = 0.0;
		nd_spts[1] = 0.0;
		nd_spts[2] = -1.0;
		nd_spts[3] = 1.0;
		nd_spts[4] = 0.0;
		nd_spts[5] = -1.0;
		nd_spts[6] = 0.0;
		nd_spts[7] = 1.0;
		nd_spts[8] = -1.0;
		nd_spts[9] = 0.0;
		nd_spts[10] = 0.0;
		nd_spts[11] = 1.0;
		nd_spts[12] = 1.0;
		nd_spts[13] = 0.0;
		nd_spts[14] = 1.0;
		nd_spts[15] = 0.0;
		nd_spts[16] = 1.0;
		nd_spts[17] = 1.0;
	} else if(type == 8 || type == 81 || type == 800) {
		num_nds = 8;
		dof_per_nd = 3;
		if(type == 8 || type == 800) {
		    num_int_dof = 0;
			n_dim = 8;
		} else {
			num_int_dof = 9;
			n_dim = 11;
		}
		def_dim = 6;
		num_ip = 8;
		num_faces = 6;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(i1 * 2);//int[i1][2];
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		double s_val[2] = {-r_1ort3,r_1ort3};
		i4 = 0;
		i5 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 2; i2++) {
				for (i3 = 0; i3 < 2; i3++) {
					int_pts[i5] = s_val[i3];
					i5++;
					int_pts[i5] = s_val[i2];
					i5++;
					int_pts[i5] = s_val[i1];
					i5++;
					ip_wt[i4] = 1.0;
					i4++;
				}
			}
		}
		if(type == 81) {
			i3 = 24;
			for (i1 = 8; i1 < 11; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					dof_table[2*i3] = i1;
					dof_table[2*i3+1] = i2;
					i3++;
				}
			}
		}
		nd_spts = vector<double>(24);
		nd_spts[0] = -1.0;
		nd_spts[1] = -1.0;
		nd_spts[2] = -1.0;
		nd_spts[3] = 1.0;
		nd_spts[4] = -1.0;
		nd_spts[5] = -1.0;
		nd_spts[6] = 1.0;
		nd_spts[7] = 1.0;
		nd_spts[8] = -1.0;
		nd_spts[9] = -1.0;
		nd_spts[10] = 1.0;
		nd_spts[11] = -1.0;
		nd_spts[12] = -1.0;
		nd_spts[13] = -1.0;
		nd_spts[14] = 1.0;
		nd_spts[15] = 1.0;
		nd_spts[16] = -1.0;
		nd_spts[17] = 1.0;
		nd_spts[18] = 1.0;
		nd_spts[19] = 1.0;
		nd_spts[20] = 1.0;
		nd_spts[21] = -1.0;
		nd_spts[22] = 1.0;
		nd_spts[23] = 1.0;
	}
	else if (type == 10 || type == 1000) {
		num_nds = 10;
		dof_per_nd = 3;
		num_int_dof = 0;
		n_dim = 10;
		def_dim = 6;
		num_ip = 4;
		num_faces = 4;
		nodes = vector<int>(num_nds);
		i1 = num_nds * dof_per_nd + num_int_dof;
		dof_table = vector<int>(2 * i1);
		int_pts = vector<double>(num_ip * 3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = r_tet1;
		int_pts[1] = r_tet1;
		int_pts[2] = r_tet1;
		int_pts[3] = r_tet2;
		int_pts[4] = r_tet1;
		int_pts[5] = r_tet1;
		int_pts[6] = r_tet1;
		int_pts[7] = r_tet2;
		int_pts[8] = r_tet1;
		int_pts[9] = r_tet1;
		int_pts[10] = r_tet1;
		int_pts[11] = r_tet2;
		ip_wt[0] = r_1o24;
		ip_wt[1] = r_1o24;
		ip_wt[2] = r_1o24;
		ip_wt[3] = r_1o24;
		s_cent[0] = 0.25;
		s_cent[1] = 0.25;
		s_cent[2] = 0.25;
		nd_spts = vector<double>(30);
		nd_spts[0] = 0.0;
		nd_spts[1] = 0.0;
		nd_spts[2] = 0.0;
		nd_spts[3] = 1.0;
		nd_spts[4] = 0.0;
		nd_spts[5] = 0.0;
		nd_spts[6] = 0.0;
		nd_spts[7] = 1.0;
		nd_spts[8] = 0.0;
		nd_spts[9] = 0.0;
		nd_spts[10] = 0.0;
		nd_spts[11] = 1.0;
		nd_spts[12] = 0.5;
		nd_spts[13] = 0.0;
		nd_spts[14] = 0.0;
		nd_spts[15] = 0.5;
		nd_spts[16] = 0.5;
		nd_spts[17] = 0.0;
		nd_spts[18] = 0.0;
		nd_spts[19] = 0.5;
		nd_spts[20] = 0.0;
		nd_spts[21] = 0.0;
		nd_spts[22] = 0.0;
		nd_spts[23] = 0.5;
		nd_spts[24] = 0.5;
		nd_spts[25] = 0.0;
		nd_spts[26] = 0.5;
		nd_spts[27] = 0.0;
		nd_spts[28] = 0.5;
		nd_spts[29] = 0.5;
	}
	else if (type == 3) {
		num_nds = 3;
		dof_per_nd = 6;
		num_int_dof = 3;
		n_dim = 6;
		def_dim = 9;
		num_ip = 3;
		num_faces = 2;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(2*i1);
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = r_1o6;
		int_pts[1] = r_1o6;
		int_pts[2] = 0.0;
		int_pts[3] = r_2o3;
		int_pts[4] = r_1o6;
		int_pts[5] = 0.0;
		int_pts[6] = r_1o6;
		int_pts[7] = r_2o3;
		int_pts[8] = 0.0;
		ip_wt[0] = r_1o6;
		ip_wt[1] = r_1o6;
		ip_wt[2] = r_1o6;
		dof_table[36] = 3;
		dof_table[37] = 2;
		dof_table[38] = 4;
		dof_table[39] = 2;
		dof_table[40] = 5;
		dof_table[41] = 2;
		s_cent[0] = r_1o3;
		s_cent[1] = r_1o3;
		s_cent[2] = r_1o3;
	} else if(type == 41) {
		num_nds = 4;
		dof_per_nd = 6;
		num_int_dof = 8;
		n_dim = 10;
		def_dim = 9;
		num_ip = 4;
		num_faces = 2;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(2*i1);
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = -r_1ort3;
		int_pts[1] = -r_1ort3;
		int_pts[2] = 0.0;
		int_pts[3] = r_1ort3;
		int_pts[4] = -r_1ort3;
		int_pts[5] = 0.0;
		int_pts[6] = -r_1ort3;
		int_pts[7] = r_1ort3;
		int_pts[8] = 0.0;
		int_pts[9] = r_1ort3;
		int_pts[10] = r_1ort3;
		int_pts[11] = 0.0;
		ip_wt[0] = 1.0;
		ip_wt[1] = 1.0;
		ip_wt[2] = 1.0;
		ip_wt[3] = 1.0;
		dof_table[48] = 4;
		dof_table[49] = 0;
		dof_table[50] = 5;
		dof_table[51] = 0;
		dof_table[52] = 4;
		dof_table[53] = 1;
		dof_table[54] = 5;
		dof_table[55] = 1;
		dof_table[56] = 6;
		dof_table[57] = 2;
		dof_table[58] = 7;
		dof_table[59] = 2;
		dof_table[60] = 8;
		dof_table[61] = 2;
		dof_table[62] = 9;
		dof_table[63] = 2;
	} else if(type == 2) {
		num_nds = 2;
		dof_per_nd = 6;
		num_int_dof = 2;
		n_dim = 3;
		def_dim = 6;
		num_ip = 2;
		num_faces = 0;
		nodes = vector<int>(num_nds);
		i1 = num_nds*dof_per_nd + num_int_dof;
		dof_table = vector<int>(i1*2);
		int_pts = vector<double>(num_ip*3);
		ip_wt = vector<double>(num_ip);
		int_pts[0] = -r_1ort3;
		int_pts[1] = 0.0;
		int_pts[2] = 0.0;
		int_pts[3] = r_1ort3;
		int_pts[4] = 0.0;
		int_pts[5] = 0.0;
		ip_wt[0] = 1.0;
		ip_wt[1] = 1.0;
		dof_table[24] = 2;
		dof_table[25] = 1;
		dof_table[26] = 2;
		dof_table[27] = 2;
	}
	else if (type == 21) { //force field
		num_nds = 2;
		dof_per_nd = 3;
		num_int_dof = 0;
		n_dim = 2;
		def_dim = 0;
		num_ip = 0;
		num_faces = 0;
		nodes = vector<int>(num_nds);
		dof_table = vector<int>(12);
		int_pts = vector<double>(1);
		ip_wt = vector<double>(1);
	}
	else if (type == 1) { // mass
		num_nds = 1;
		dof_per_nd = 3;
		num_int_dof = 0;
		n_dim = 1;
		def_dim = 0;
		num_ip = 0;
		num_faces = 0;
		nodes = vector<int>(num_nds);
		dof_table = vector<int>(6);
		int_pts = vector<double>(1);
		int_pts = vector<double>(1);
	}
	
	i3 = 0;
	for (i1 = 0; i1 < num_nds; i1++) {
		for (i2 = 0; i2 < dof_per_nd; i2++) {
			dof_table[i3] = i1;
			dof_table[i3+1] = i2;
			i3+= 2;
		}
	}
	
	if(num_int_dof != 0) {
		internal_disp = vector<double>(num_int_dof);
		int_prev_disp = vector<double>(num_int_dof);
		internald_ldu = vector<double>(num_int_dof);
		internal_adj = vector<double>(num_int_dof);
		internal_ru = vector<DiffDoub1>(num_int_dof);
		i1 = (num_nds*dof_per_nd + num_int_dof)*num_int_dof;
		internal_mat = vector<double>(i1);
	}

	int_dof_index = 0;
	
	sect_ptr = max_int;
	
	return;
}

void Element::set_nodes(int new_nds[]) {
	int i1;
	for (i1 = 0; i1 < num_nds; i1++) {
		nodes[i1] = new_nds[i1];
	}
	return;
}

void Element::initialize_faces(vector<Face>& glob_fc_lst, int& fi) {
	//fi = the current number of faces that have been written into glob_fc_lst
	Face* new_fc = nullptr;
	if(type == 4 || type == 400) {
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,1,nodes[1]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,1,nodes[1]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,3,nodes[3]);
		new_fc->set_node(2,2,nodes[2]);
		faces.push_back(fi);
		fi++;
	} else if(type == 6 || type == 600) {
        glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,1,nodes[1]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,3,nodes[3]);
		new_fc->set_node(1,4,nodes[4]);
		new_fc->set_node(2,5,nodes[5]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,4,nodes[4]);
		new_fc->set_node(3,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,1,nodes[1]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,5,nodes[5]);
		new_fc->set_node(3,4,nodes[4]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,3,nodes[3]);
		new_fc->set_node(2,5,nodes[5]);
		new_fc->set_node(3,2,nodes[2]);
		faces.push_back(fi);
		fi++;
	} else if(type == 8 || type == 81 || type == 800) {
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,3,nodes[3]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,1,nodes[1]);
		new_fc->set_node(3,0,nodes[0]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,4,nodes[4]);
		new_fc->set_node(1,5,nodes[5]);
		new_fc->set_node(2,6,nodes[6]);
		new_fc->set_node(3,7,nodes[7]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,5,nodes[5]);
		new_fc->set_node(3,4,nodes[4]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,1,nodes[1]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,6,nodes[6]);
		new_fc->set_node(3,5,nodes[5]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,2,nodes[2]);
		new_fc->set_node(1,3,nodes[3]);
		new_fc->set_node(2,7,nodes[7]);
		new_fc->set_node(3,6,nodes[6]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,3,nodes[3]);
		new_fc->set_node(1,0,nodes[0]);
		new_fc->set_node(2,4,nodes[4]);
		new_fc->set_node(3,7,nodes[7]);
		faces.push_back(fi);
		fi++;
	}
	else if (type == 10 || type == 1000) {
		glob_fc_lst[fi].num_nds = 6;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0, 0, nodes[0]);
		new_fc->set_node(1, 2, nodes[2]);
		new_fc->set_node(2, 1, nodes[1]);
		new_fc->set_node(3, 6, nodes[6]);
		new_fc->set_node(4, 5, nodes[5]);
		new_fc->set_node(5, 4, nodes[4]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 6;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0, 0, nodes[0]);
		new_fc->set_node(1, 1, nodes[1]);
		new_fc->set_node(2, 3, nodes[3]);
		new_fc->set_node(3, 4, nodes[4]);
		new_fc->set_node(4, 8, nodes[8]);
		new_fc->set_node(5, 7, nodes[7]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 6;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0, 1, nodes[1]);
		new_fc->set_node(1, 2, nodes[2]);
		new_fc->set_node(2, 3, nodes[3]);
		new_fc->set_node(3, 5, nodes[5]);
		new_fc->set_node(4, 9, nodes[9]);
		new_fc->set_node(5, 8, nodes[8]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 6;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0, 0, nodes[0]);
		new_fc->set_node(1, 3, nodes[3]);
		new_fc->set_node(2, 2, nodes[2]);
		new_fc->set_node(3, 7, nodes[7]);
		new_fc->set_node(4, 9, nodes[9]);
		new_fc->set_node(5, 6, nodes[6]);
		faces.push_back(fi);
		fi++;
	}
	else if (type == 3) {
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,2,nodes[2]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 3;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,2,nodes[2]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,0,nodes[0]);
		faces.push_back(fi);
		fi++;
	} else if(type == 41) {
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,0,nodes[0]);
		new_fc->set_node(1,1,nodes[1]);
		new_fc->set_node(2,2,nodes[2]);
		new_fc->set_node(3,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		glob_fc_lst[fi].num_nds = 4;
		new_fc = &glob_fc_lst[fi];
		new_fc->set_node(0,3,nodes[3]);
		new_fc->set_node(1,2,nodes[2]);
		new_fc->set_node(2,1,nodes[1]);
		new_fc->set_node(3,0,nodes[0]);
		faces.push_back(fi);
		fi++;
	} 
	
	return;
}

void Element::set_int_disp(double new_disp[]) {
	int i1;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		internal_disp[i1] = new_disp[i1];
	}
	return;
}

void Element::set_int_prev_disp(double new_disp[]) {
	int i1;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		int_prev_disp[i1] = new_disp[i1];
	}
	return;
}

void Element::advance_int_disp() {
	int i1;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		int_prev_disp[i1] = internal_disp[i1];
	}
	return;
}

void Element::backstep_int_disp() {
	int i1;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		internal_disp[i1] = int_prev_disp[i1];
	}
	return;
}

void Element::set_intd_ld_u(vector<double>& globd_ld_u) {
	int i1;
	int i2 = int_dof_index;
	for (i1 = 0; i1 < num_int_dof; i1++) {
		internald_ldu[i1] = globd_ld_u[i2];
		i2++;
	}
	return;
}

int Element::get_num_layers(vector<Section>& sec_lst) {
	return sec_lst[sect_ptr].layers.size();
}


void Element::add_design_variable(int d_index, double coef) {
	int i1 = design_vars.size();
	IDCapsule dv;
	dv.int_dat = d_index;
	dv.doub_dat = coef;
	design_vars.push_back(dv);
	return;
}

void Element::add_comp_dvar(int d_index) {
	for (auto& dv : comp_dvars) {
		if (dv == d_index) {
			return;
		}
	}
	comp_dvars.push_back(d_index);
	return;
}


