#include "MeshElement.h"
#include "MeshNode.h"
#include "utilities.h"
#include <vector>

using namespace std;

MeshElement::MeshElement() {
	return;
}

void MeshElement::get_centroid(double cent[], vector<MeshNode>& nd_ar) {
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	double* n4 = &nd_ar[nodes[3]].coord[0];
	cent[0] = 0.25 * (n1[0] + n2[0] + n3[0] + n4[0]);
	cent[1] = 0.25 * (n1[1] + n2[1] + n3[1] + n4[1]);
	cent[2] = 0.25 * (n1[2] + n2[2] + n3[2] + n4[2]);
	return;
}

double MeshElement::get_volume(vector<MeshNode>& nd_ar) {
	double mat[9];
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	double* n4 = &nd_ar[nodes[3]].coord[0];
	mat[0] = n2[0] - n1[0];
	mat[3] = n2[1] - n1[1];
	mat[6] = n2[2] - n1[2];
	mat[1] = n3[0] - n1[0];
	mat[4] = n3[1] - n1[1];
	mat[7] = n3[2] - n1[2];
	mat[2] = n4[0] - n1[0];
	mat[5] = n4[1] - n1[1];
	mat[8] = n4[2] - n1[2];
	q_rfactor(mat, 3, 0, 2, 0, 2, 0);
	double vol = mat[0] * mat[4] * mat[8];
	return vol;
}

bool MeshElement::point_in(double pt[], vector<MeshNode>& nd_ar) {
	double mat[9];
	double rhs[3];
	double soln[3];
	double sol_sum;
	double* n1 = &nd_ar[nodes[0]].coord[0];
	double* n2 = &nd_ar[nodes[1]].coord[0];
	double* n3 = &nd_ar[nodes[2]].coord[0];
	double* n4 = &nd_ar[nodes[3]].coord[0];
	mat[0] = n2[0] - n1[0];
	mat[3] = n2[1] - n1[1];
	mat[6] = n2[2] - n1[2];
	mat[1] = n3[0] - n1[0];
	mat[4] = n3[1] - n1[1];
	mat[7] = n3[2] - n1[2];
	mat[2] = n4[0] - n1[0];
	mat[5] = n4[1] - n1[1];
	mat[8] = n4[2] - n1[2];
	rhs[0] = pt[0] - n1[0];
	rhs[1] = pt[1] - n1[1];
	rhs[2] = pt[2] - n1[2];
	q_rfactor(mat, 3, 0, 2, 0, 2, 0);
	solveq_rx_eqb(soln, mat, rhs, 3, 0, 2, 0, 2, 0);
	if (soln[0] > 0.0000000001 && soln[1] > 0.0000000001 && soln[2] > 0.0000000001) {
		sol_sum = soln[0] + soln[1] + soln[2];
		if (sol_sum < 0.9999999999) {
			return true;
		}
	}
	return false;
}
