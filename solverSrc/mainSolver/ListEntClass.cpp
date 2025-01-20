#include <cmath>
#include <cstddef>
#include <string>
#include <fstream>
#include <vector>
#include "ListEntClass.h"

using namespace std;

IDCapsule::IDCapsule() {
	int_dat = 0;
	doub_dat = 0.0;
}

MatrixEnt::MatrixEnt() {
	row = 0;
	col = 0;
	value = 0.0;
}

MatrixRow::MatrixRow() {
	row_vec.clear();
}

void MatrixRow::add_entry(int row, int col, double val) {
	auto iter = row_vec.begin();
	bool inserted = false;
	while (iter != row_vec.end() && !inserted) {
		if (iter->col == col) {
			iter->value += val;
			inserted = true;
		}
		++iter;
	}
	if (!inserted) {
		MatrixEnt new_ent;
		new_ent.row = row;
		new_ent.col = col;
		new_ent.value = val;
		row_vec.push_back(new_ent);
	}
	return;
}

SparseMat::SparseMat() {
	dim = 0;
	matrix.clear();
	return;
}

void SparseMat::set_dim(int new_dim) {
	int i1;
	dim = new_dim;
	matrix = vector<MatrixRow>(new_dim);
	return;
}

void SparseMat::zero_all() {
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.row_vec) {
			i2.value = 0.0;
		}
	}
	
	return;
}

void SparseMat::add_entry(int row, int col, double val) {
	matrix[row].add_entry(row, col, val);
	return;
}

void SparseMat::add_matrix(SparseMat& inp_mat) {
	for (auto& i1 : inp_mat.matrix) {
		for (auto& i2 : i1.row_vec) {
			add_entry(i2.row, i2.col, i2.value);
		}
	}
	return;
}

void SparseMat::vector_multiply(vector<double>& prod, vector<double>& inp_vec, bool transpose) {
	if(transpose) {
		for (auto& i1 : matrix) {
			for (auto& i2 : i1.row_vec) {
				prod[i2.col] += i2.value * inp_vec[i2.row];
			}
		}
	} else {
		for (auto& i1 : matrix) {
			for (auto& i2 : i1.row_vec) {
				prod[i2.row] += i2.value * inp_vec[i2.col];
			}
		}
	}
	return;
}

double SparseMat::get_max_abs_val() {
	double this_val;
	double max_val = 0.0;
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.row_vec) {
			this_val = abs(i2.value);
			if (this_val > max_val) {
				max_val = this_val;
			}
		}
	}

	return max_val;
}

void SparseMat::write_to_file(ofstream& out_file) {
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.row_vec) {
			out_file << "    - [" << i2.row << ", " << i2.col << ", " << i2.value << "]\n";
		}
	}
	return;
}
