#include "MeshFace.h"
#include "constants.h"
#include "MeshNode.h"
#include "MeshElement.h"
#include "utilities.h"
#include <iostream>
#include <cmath>

using namespace std;

MeshFace::MeshFace() {
	nodes[0] = max_int;
	nodes[1] = max_int;
	nodes[2] = max_int;
	elements[0] = max_int;
	elements[1] = max_int;
	return;
}

void MeshFace::copy_data(MeshFace& f_in) {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		nodes[i1] = f_in.nodes[i1];
		norm_dir[i1] = f_in.norm_dir[i1];
	}
	elements[0] = f_in.elements[0];
	elements[1] = f_in.elements[1];
	proj_dist = f_in.proj_dist;
	return;
}

void MeshFace::init_norm_dir(vector<MeshNode>& nd_ar) {
	double v1[3];
	double v2[3];
	double cp[3];
	double mag;
	double area;
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		v1[i1] = n2[i1] - n1[i1];
		v2[i1] = n3[i1] - n1[i1];
	}
	cross_prod(cp, v1, v2);
	mag = sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]);
	area = 0.5 * mag;
	mag = 1.0 / mag;
	norm_dir[0] = mag * cp[0];
	norm_dir[1] = mag * cp[1];
	norm_dir[2] = mag * cp[2];
	proj_dist = 1.2408064788027997 * sqrt(area); // (2^1.5/3^0.75)*sqrt(area)
	return;
}

void MeshFace::norm_dir_from_el_cent(double cent[], vector<MeshNode>& nd_ar) {
	double fc_cent[3];
	double d_vec[3];
	get_centroid(fc_cent,nd_ar);
	d_vec[0] = fc_cent[0] - cent[0];
	d_vec[1] = fc_cent[1] - cent[1];
	d_vec[2] = fc_cent[2] - cent[2];
	double dp = d_vec[0] * norm_dir[0] + d_vec[1] * norm_dir[1] + d_vec[2] * norm_dir[2];
	if (dp < 0.0) {
		norm_dir[0] *= -1.0;
		norm_dir[1] *= -1.0;
		norm_dir[2] *= -1.0;
	}
	return;
}

void MeshFace::get_centroid(double cent[], vector<MeshNode>& nd_ar) {
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	cent[0] = 0.333333333333333 * (n1[0] + n2[0] + n3[0]);
	cent[1] = 0.333333333333333 * (n1[1] + n2[1] + n3[1]);
	cent[2] = 0.333333333333333 * (n1[2] + n2[2] + n3[2]);
	return;
}

double MeshFace::get_longest_edge_len(vector<MeshNode>& nd_ar) {
	double d_vec[3];
	double dist;
	double max_dist = 0.0;
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	d_vec[0] = n2[0] - n1[0];
	d_vec[1] = n2[1] - n1[1];
	d_vec[2] = n2[2] - n1[2];
	dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
	if (dist > max_dist) {
		max_dist = dist;
	}
	d_vec[0] = n3[0] - n2[0];
	d_vec[1] = n3[1] - n2[1];
	d_vec[2] = n3[2] - n2[2];
	dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
	if (dist > max_dist) {
		max_dist = dist;
	}
	d_vec[0] = n1[0] - n3[0];
	d_vec[1] = n1[1] - n3[1];
	d_vec[2] = n1[2] - n3[2];
	dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
	if (dist > max_dist) {
		max_dist = dist;
	}
	return max_dist;
}

bool MeshFace::get_intersection(double out_param[], double pt[], double vec[], vector<MeshNode>& nd_ar) {
	int i1;
	double mat[9];
	double mat_mag;
	double rhs[3];
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	mat[0] = vec[0];
	mat[3] = vec[1];
	mat[6] = vec[2];
	mat[1] = n1[0] - n2[0];
	mat[4] = n1[1] - n2[1];
	mat[7] = n1[2] - n2[2];
	mat[2] = n1[0] - n3[0];
	mat[5] = n1[1] - n3[1];
	mat[8] = n1[2] - n3[2];
	mat_mag = 0.0;
	for (i1 = 0; i1 < 9; i1++) {
		mat_mag += (mat[i1] * mat[i1]);
	}
	mat_mag = sqrt(mat_mag);
	q_rfactor(mat, 3, 0, 2, 0, 2, 0);
	double det = mat[0] * mat[4] * mat[8];
	if (abs(det) < 1.0e-6*mat_mag*mat_mag*mat_mag) {
		return false;
	}
	else {
		rhs[0] = n1[0] - pt[0];
		rhs[1] = n1[1] - pt[1];
		rhs[2] = n1[2] - pt[2];
		solveq_rx_eqb(out_param, mat, rhs, 3, 0, 2, 0, 2, 0);
		if (out_param[1] > 0.0000000001 && out_param[2] > 0.0000000001) {
			if ((out_param[1] + out_param[2]) < 0.9999999999) {
				return true;
			}
		}
		return false;
	}
}

bool MeshFace::edges_intersect(int fc, double dist_tol, vector<MeshNode>& nd_ar, vector<MeshFace>& fc_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double mat[6];
	double rhs[3];
	double soln[2];
	double det;
	double d_vec[3];
	double dist;
	double v1[3];
	double v2[3];
	double mat_mag;
	double* n11;
	double* n21;
	double* n12;
	double* n22;
	int* fc_nodes = &fc_ar[fc].nodes[0];
	
	for (i1 = 0; i1 < 3; i1++) {
		n11 = &nd_ar[nodes[i1]].coord[0];
		i3 = i1 + 1;
		if (i3 > 2) {
			i3 = 0;
		}
		n21 = &nd_ar[nodes[i3]].coord[0];
		for (i2 = 0; i2 < 3; i2++) {
			n12 = &nd_ar[fc_nodes[i2]].coord[0];
			i4 = i2 + 1;
			if (i4 > 2) {
				i4 = 0;
			}
			n22 = &nd_ar[fc_nodes[i4]].coord[0];
			v1[0] = n21[0] - n11[0];
			v1[1] = n21[1] - n11[1];
			v1[2] = n21[2] - n11[2];
			v2[0] = n22[0] - n12[0];
			v2[1] = n22[1] - n12[1];
			v2[2] = n22[2] - n12[2];
			mat[0] = v1[0];
			mat[1] = -v2[0];
			mat[2] = v1[1];
			mat[3] = -v2[1];
			mat[4] = v1[2];
			mat[5] = -v2[2];
			rhs[0] = n12[0] - n11[0];
			rhs[1] = n12[1] - n11[1];
			rhs[2] = n12[2] - n11[2];
			mat_mag = 0.0;
			for (i5 = 0; i5 < 6; i5++) {
				mat_mag += (mat[i5] * mat[i5]);
			}
			mat_mag = sqrt(mat_mag);
			q_rfactor(mat, 2, 0, 2, 0, 1, 0);
			det = mat[0] * mat[3];
			if (abs(det) > 1.0e-6*mat_mag*mat_mag) {
				solveq_rx_eqb(soln, mat, rhs, 2, 0, 2, 0, 1, 0);
				if (soln[0] > 0.0000000001 && soln[0] < 0.9999999999 && soln[1] > 0.0000000001 && soln[1] < 0.9999999999) {
					d_vec[0] = n11[0] + soln[0] * v1[0] - n12[0] - soln[1] * v2[0];
					d_vec[1] = n11[1] + soln[0] * v1[1] - n12[1] - soln[1] * v2[1];
					d_vec[2] = n11[2] + soln[0] * v1[2] - n12[2] - soln[1] * v2[2];
					dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
					if (dist < dist_tol) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

int MeshFace::get_shared_nodes(int nd_pts[], bool shared[], int fc, vector<MeshNode>& nd_ar, vector<MeshFace>& fc_ar) {
	int i1;
	int i2;
	int* fc_nds = &fc_ar[fc].nodes[0];
	nd_pts[0] = nodes[0];
	nd_pts[1] = nodes[1];
	nd_pts[2] = nodes[2];
	shared[0] = false;
	shared[1] = false;
	shared[2] = false;
	int num_shared = 0;
	for (i1 = 0; i1 < 3; i1++) {
		for (i2 = 0; i2 < 3; i2++) {
			if (nd_pts[i1] == fc_nds[i2]) {
				shared[i1] = true;
				num_shared++;
			}
		}
	}
	return num_shared;
}

void MeshFace::print_info(vector<MeshNode>& nd_ar) {
	int i1;
	double* crd;
	cout << "nodes:" << endl;
	for (i1 = 0; i1 < 3; i1++) {
		crd = &nd_ar[nodes[i1]].coord[0];
		cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
	}
	return;
}