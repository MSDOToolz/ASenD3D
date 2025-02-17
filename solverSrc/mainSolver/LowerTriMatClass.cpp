#include "LowerTriMatClass.h"
#include "ListEntClass.h"
#include "ConstraintClass.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

LowerTriMat::LowerTriMat() {
	mat.clear();
	range.clear();
	min_col.clear();
	z_vec.clear();
	dim = 0;
	size = 0;
	max_bandwidth = 0;
	allocated = false;
	return;
}

void LowerTriMat::set_dim(int new_dim) {
	dim = new_dim;
	range = vector<int>(new_dim+1);
	min_col = vector<int>(new_dim);
	z_vec = vector<double>(new_dim);
	return;
}

void LowerTriMat::allocate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list, int block_dim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int row;
	int curr_block;
	int blk_mc;
	int const_dim;
	
	if(!allocated) {
		set_dim(sp_mat.dim);
	}
	
	for (auto& mr : sp_mat.matrix) {
		row = mr.row_vec.begin()->row;
		range[row] = 1;
		curr_block = row / block_dim;
		blk_mc = block_dim * curr_block;
		for (auto& me : mr.row_vec) {
			i2 = me.col;
			if(i2 <= row && i2 >= blk_mc) {
				i3 = row - i2 + 1;
				if(i3 > range[row]) {
					range[row] = i3;
				}
			}
		}
	}
	
	for (auto& cnst : c_list.const_vec) {
		SparseMat& this_mat = cnst.mat;
		for (auto& mr : this_mat.matrix) {
			for (auto& me : mr.row_vec) {
				i2 = me.col;
				curr_block = i2 / block_dim;
				blk_mc = block_dim * curr_block;
				for (auto& me2 : mr.row_vec) {
					i3 = me2.col;
					if(i3 <= i2 && i3 >= blk_mc) {
						i4 = i2 - i3 + 1;
						if(i4 > range[i2]) {
							range[i2] = i4;
						}
					}
				}
			}
		}
	}

	size = 0;
	max_bandwidth = 0;
	for (i1 = 0; i1 < dim; i1++) {
		size+= range[i1];
		min_col[i1] = i1 - range[i1] + 1;
		if(range[i1] > max_bandwidth) {
			max_bandwidth = range[i1];
		}
	}
	range[dim] = size;
	for (i1 = (dim-1); i1 >= 0; i1--) {
		range[i1] = range[i1+1] - range[i1];
	}
	
	mat = vector<double>(size);
	
	allocated = true;

	cout << "maxBandwidth: " << max_bandwidth << endl;
	
	return;
}

bool LowerTriMat::is_allocated() {
	return allocated;
}

void LowerTriMat::populate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list) {
	int i1;
	int i2;
	int i3;
	int i4;
	int const_dim;
	double const_sf;
	
	for (auto& me : mat) {
		me = 0.0;
	}
	
	i1 = 0;
	for (auto& mr : sp_mat.matrix) {
		for (auto& me : mr.row_vec) {
			i2 = me.col;
			if(i2 <= i1 && i2 >= min_col[i1]) {
				i3 = range[i1] + (i2 - min_col[i1]);
				mat[i3]+= me.value;
			}
		}
		i1++;
	}
	
	for (auto& cnst : c_list.const_vec) {
		SparseMat& this_mat = cnst.mat;
		const_sf = cnst.scale_fact;
		for (auto& mr : this_mat.matrix) {
			for (auto& me : mr.row_vec) {
				i2 = me.col;
				for (auto& me2 : mr.row_vec) {
					i3 = me2.col;
					if(i3 <= i2 && i3 >= min_col[i2]) {
						i4 = range[i2] + (i3 - min_col[i2]);
						mat[i4] += const_sf * me.value * me2.value;
					}
				}
			}
		}
	}		
    	
	return;
}

void LowerTriMat::ldl_factor() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int st_col;
	double sum;
	
	vector<double> ld_vec(dim);
	
	for (i1 = 0; i1 < dim; i1++) { // i1 = column in l
		//form ld_vec
		i3 = range[i1];
		for (i2 = min_col[i1]; i2 < i1; i2++) {
			i4 = range[i2+1] - 1;
			ld_vec[i2] = mat[i3]*mat[i4];
			i3++;
		}
		//get d term for column i1
		st_col = min_col[i1];				
		i2 = range[i1];
		sum = 0.0;
		for (i3 = st_col; i3 < i1; i3++) {
			sum+= ld_vec[i3]*mat[i2];
			i2++;
		}
		i2 = range[i1+1] - 1;
		mat[i2] = mat[i2] - sum;
		//get l terms for column i1
		i2 = i1 + 1;
		while(i2 < (i1 + max_bandwidth) && i2 < dim) { // i2 = row in l
		    if(min_col[i2] <= i1) {
				st_col = min_col[i2];
				if(min_col[i1] > st_col) {
					st_col = min_col[i1];
				}
				i4 = range[i2] + (st_col - min_col[i2]);
				sum = 0.0;
				for (i5 = st_col; i5 < i1; i5++) {
					sum+= ld_vec[i5]*mat[i4];
					i4++;
				}
				i3 = range[i2] + (i1 - min_col[i2]);
				i4 = range[i1+1] - 1;
				mat[i3] = (mat[i3] - sum)/mat[i4];
			}
			i2++;
		}
	}
	
	return;
}

void LowerTriMat::ldl_solve(vector<double>& soln_vec, vector<double>& rhs) {
	int i1;
	int i2;
	int i3;
	int stop_row;
	double sum;
	
	z_vec[0] = rhs[0];
	for (i1 = 1; i1 < dim; i1++) {
		i2 = range[i1];
		sum = 0.0;
		for (i3 = min_col[i1]; i3 < i1; i3++) {
			sum+= mat[i2]*z_vec[i3];
			i2++;
		}
		z_vec[i1] = rhs[i1] - sum;
	}
	
	for (i1 = dim-1; i1 >= 0; i1--) {
		sum = 0.0;
		stop_row = i1 + max_bandwidth + 1;
		if(stop_row > dim) {
			stop_row = dim;
		}
		for (i2 = (i1+1); i2 < stop_row; i2++) {
			i3 = range[i2] + (i1 - min_col[i2]);
			if(i3 >= range[i2] && i3 < range[i2+1]) {
				sum+= mat[i3]*soln_vec[i2];
			}
		}
		i2 = range[i1+1] - 1;
		soln_vec[i1] = (z_vec[i1]/mat[i2]) - sum;
	}
	return;
}

bool LowerTriMat::pos_def() {
	int i1;
	double d_val;
	if (!allocated) {
		return false;
	}
	for (i1 = 0; i1 < dim; i1++) {
		d_val = mat[range[i1 + 1] - 1];
		if (d_val < 0.0) {
			return false;
		}
	}
	return true;
}

void LowerTriMat::write_to_file(ofstream& out_file) {
	int i1;
	int i2;
	int col;
	
	for (i1 = 0; i1 < dim; i1++) {
		col = min_col[i1];
		for (i2 = range[i1]; i2 < range[i1 + 1]; i2++) {
			out_file << "    - [" << i1 << ", " << col << ", " << mat[i2] << "]\n";
			col++;
		}
	}

	return;
}