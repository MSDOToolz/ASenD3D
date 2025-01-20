#include "ConstraintClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

const int max_int = 2000000000;

ConstraintTerm::ConstraintTerm() {
	node_set = "";
	ns_ptr = max_int;
	int dof = max_int;
	coef = 0.0;
	return;
}


Constraint::Constraint() {
	terms.clear();
	rhs = 0.0;
	scale_fact = 0.0;
	return;
}

void Constraint::build_mat(vector<Node>& nd_ar, vector<Set>& set_ar) {
	int set_len = 1;
	int seti_len;
	int nd_index;
	int row;
	int col;
	int dof;
	double coef;
	string set_nm;
	bool set_found;
	int this_set;

	for (auto& tm : terms) {
		this_set = tm.ns_ptr;
		seti_len = set_ar[this_set].labels.size();
		if (seti_len > set_len) {
			set_len = seti_len;
		}
	}

	mat.set_dim(set_len);
	for (auto& tm : terms) {
		coef = tm.coef;
		dof = tm.dof;
		this_set = tm.ns_ptr;
		list<int>& set_labs = set_ar[this_set].labels;
		if (set_labs.size() == 1) {
			nd_index = *set_labs.begin();
			if (type == "displacement") {
				col = nd_ar[nd_index].dof_index[dof - 1];
			}
			else {
				col = nd_ar[nd_index].sorted_rank;
			}
			for (row = 0; row < set_len; row++) {
				mat.add_entry(row, col, coef);
			}
		}
		else {
			row = 0;
			for (auto& nd : set_labs) {
				if (type == "displacement") {
					col = nd_ar[nd].dof_index[dof - 1];
				}
				else {
					col = nd_ar[nd].sorted_rank;
				}
				mat.add_entry(row, col, coef);
				row++;
			}
		}
	}
	return;
}

void Constraint::full_vec_multiply(vector<double>& prod, vector<double>& vec, vector<double>& tmp_v) {
	// compute scale_fact*[mat]^t*([mat]*vec + q_vec)
	int i1;
	for (auto& tv : tmp_v) {
		tv = 0.0;
	}
	mat.vector_multiply(tmp_v, vec, false);
	for (auto& tv : tmp_v) {
		tv *= scale_fact;
	}
	mat.vector_multiply(prod, tmp_v, true);
	return;
}

void Constraint::get_load(vector<double>& c_ld, vector<double>& u_vec, vector<double>& q_vec, int res_dim) {
	int i1;
	int dim = mat.dim;
	for (i1 = 0; i1 < dim; i1++) {
		q_vec[i1] = -rhs;
	}
	mat.vector_multiply(q_vec, u_vec, false);
	for (i1 = 0; i1 < dim; i1++) {
		q_vec[i1] *= -scale_fact;
	}
	mat.vector_multiply(c_ld, q_vec, true);

	return;
}

void Constraint::write_to_file(ofstream& out_file) {
	out_file << "type: " << type << "\n";
	out_file << "scale factor: " << scale_fact << "\n";
	out_file << "rhs: " << rhs << "\n";
	out_file << "matrix: \n";
	mat.write_to_file(out_file);
	return;
}

ConstraintList::ConstraintList() {
	const_vec.clear();
	return;
}

void ConstraintList::set_scale_fact(double new_sf) {
	for (auto& cnst : const_vec) {
		cnst.scale_fact = new_sf;
	}
	return;
}

void ConstraintList::get_total_vec_mult(vector<double>& prod, vector<double>& vec, vector<double>& tmp_v) {
	for (auto& cnst : const_vec) {
		cnst.full_vec_multiply(prod, vec, tmp_v);
	}
	return;
}

void ConstraintList::get_total_load(vector<double>& c_ld, vector<double>& u_vec, vector<double>& q_vec, int res_dim) {
	for (auto& cnst : const_vec) {
		cnst.get_load(c_ld,u_vec,q_vec,res_dim);
	}
	return;
}

void ConstraintList::write_all_to_file(string file_name) {
	ofstream out_file;
	out_file.open(file_name);
	for (auto& cnst : const_vec) {
		cnst.write_to_file(out_file);
		out_file << "##------------------------------------------------------\n";
	}
	out_file.close();
	return;
}