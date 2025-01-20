#include <iostream>
#include <iomanip>
#include <vector>
#include "LUMatClass.h"
#include "ListEntClass.h"
#include "ConstraintClass.h"

using namespace std;

LUMat::LUMat() {
	l_mat.clear();
	l_range.clear();
	l_min_col.clear();
	l_size = 0;
	u_mat.clear();
	u_range.clear();
	u_min_row.clear();
	u_size = 0;
	z_vec.clear();
	dim = 0;
	max_bandwidth = 0;
	allocated = false;
	return;
}

void LUMat::set_dim(int new_dim) {
	dim = new_dim;
	l_range = vector<int>(new_dim+1);
	l_min_col = vector<int>(new_dim);
	u_range = vector<int>(new_dim+1);
	u_min_row = vector<int>(new_dim);
	z_vec = vector<double>(new_dim);
	return;
}

void LUMat::allocate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list, int block_dim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int curr_block;
	int const_dim;
	int blk_max_col;
	int blk_min_col;
	
	if(!allocated) {
		set_dim(sp_mat.dim);
	}

	for (i1 = 0; i1 < dim; i1++) {
		l_range[i1] = 0;
		u_range[i1] = 1;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		curr_block = i1/block_dim;
		blk_min_col = curr_block*block_dim;
		blk_max_col = blk_min_col + block_dim;
		MatrixRow& mr = sp_mat.matrix[i1];
		for (auto& me : mr.row_vec) {
			i2 = me.col;
			if(i2 >= blk_min_col && i2 < blk_max_col) {
				if(i2 < i1) { // in l
				    i3 = i1 - i2;
					if(i3 > l_range[i1]) {
						l_range[i1] = i3;
					}
				}
				else {
					i3 = i2 - i1 + 1;
					if(i3 > u_range[i2]) {
						u_range[i2] = i3;
					}
				}
			}
		}
	}
	
	for (auto& cnst : c_list.const_vec) {
		SparseMat& this_mat = cnst.mat;
		for (auto& mr : this_mat.matrix) {
			for (auto& me : mr.row_vec) {
				i2 = me.col;
				curr_block = i2/block_dim;
				blk_min_col = curr_block*block_dim;
				blk_max_col = blk_min_col + block_dim;
				for (auto& me2 : mr.row_vec) {
					i3 = me2.col;
					if(i3 >= blk_min_col && i3 < blk_max_col) {
						if(i3 < i2) { // in lu
						    i4 = i2 - i3;
							if(i4 > l_range[i2]) {
								l_range[i2] = i4;
							}
						}
						else {
							i4 = i3 - i2 + 1;
							if(i4 > u_range[i3]) {
								u_range[i3] = i4;
							}
						}
					}
				}
			}
		}
	}
	
	l_size = 0;
	u_size = 0;
	max_bandwidth = 0;
	for (i1 = 0; i1 < dim; i1++) {
		l_size += l_range[i1];
		u_size += u_range[i1];
		l_min_col[i1] = i1 - l_range[i1];
		if(l_range[i1] > max_bandwidth) {
			max_bandwidth = l_range[i1];
		}
		u_min_row[i1] = i1 - u_range[i1] + 1;
		if(u_range[i1] > max_bandwidth) {
			max_bandwidth = u_range[i1];
		}
	}
	l_range[dim] = l_size;
	u_range[dim] = u_size;
	for (i1 = (dim-1); i1 >= 0; i1--) {
		l_range[i1] = l_range[i1+1] - l_range[i1];
		u_range[i1] = u_range[i1+1] - u_range[i1];
	}
	
	l_mat = vector<double>(l_size);
	u_mat = vector<double>(u_size);
	
	allocated = true;
	
	cout << "maxBandwidth: " << max_bandwidth << endl;
	
	return;
}

void LUMat::populate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list) {
	int i1;
	int i2;
	int i3;
	int i4;
	int const_dim;
	double const_sf;
	
	for (auto& lm : l_mat) {
		lm = 0.0;
	}
	
	for (auto& um : u_mat) {
		um = 0.0;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		MatrixRow& mr = sp_mat.matrix[i1];
		for (auto& me : mr.row_vec) {
			i2 = me.col;
			if(i2 >= l_min_col[i1] && i2 < i1) { // in l
				i3 = l_range[i1] + i2 - l_min_col[i1];
			    l_mat[i3] += me.value;
			}
			else if(i1 >= u_min_row[i2] && i1 <= i2) {
				i3 = u_range[i2] + i1 - u_min_row[i2];
				u_mat[i3] += me.value;
			}
		}
	}
	
	for (auto& cnst : c_list.const_vec) {
		SparseMat& this_mat = cnst.mat;
		const_sf = cnst.scale_fact;
		for (auto& mr : this_mat.matrix) {
			for (auto& me : mr.row_vec) {
				i2 = me.col;
				for (auto& me2 : mr.row_vec) {
					i3 = me2.col;
					if(i3 >= l_min_col[i2] && i3 < i2) { // in lu
						i4 = l_range[i2] + i3 - l_min_col[i2];
						l_mat[i4] += const_sf*me.value*me2.value;
					}
					else if(i2 >= u_min_row[i3] && i2 <= i3) {
						i4 = u_range[i3] + i2 - u_min_row[i3];
						u_mat[i4] += const_sf*me.value*me2.value;
					}
				}
			}
		}
	}
	
	return;
}

void LUMat::lu_factor() {
	int i1;
	int i2;
	int i3;
	int this_ind;
	int l_ind;
	int u_ind;
	int u_piv;
	int max_col;
	int max_row;
	double tmp;
	
	for (i1 = 0; i1 < dim; i1++) {
		max_col = i1 + max_bandwidth;
		if(max_col > dim) {
			max_col = dim;
		}
		for (i2 = i1; i2 < max_col; i2++) {
			if(u_min_row[i2] <= i1) {
				this_ind = u_range[i2] + i1 - u_min_row[i2];
				i3 = u_min_row[i2];
				if(l_min_col[i1] > i3) {
					i3 = l_min_col[i1];
				}
				l_ind = l_range[i1] + i3 - l_min_col[i1];
				u_ind = u_range[i2] + i3 - u_min_row[i2];
				tmp = u_mat[this_ind];
				while(l_ind < l_range[i1+1]) {
					tmp -= l_mat[l_ind]*u_mat[u_ind];
					l_ind++;
					u_ind++;
				}
				u_mat[this_ind] = tmp;
			}
		}
		max_row = max_col;
		for (i2 = (i1+1); i2 < max_row; i2++) {
			if(l_min_col[i2] <= i1) {
				this_ind = l_range[i2] + i1 - l_min_col[i2];
				i3 = l_min_col[i2];
				if(u_min_row[i1] > i3) {
					i3 = u_min_row[i1];
				}
				l_ind = l_range[i2] + i3 - l_min_col[i2];
				u_ind = u_range[i1] + i3 - u_min_row[i1];
				tmp = l_mat[this_ind];
				u_piv = u_range[i1+1] - 1;
				while(u_ind < u_piv) {
					tmp -= l_mat[l_ind]*u_mat[u_ind];
					l_ind++;
					u_ind++;
				}
				l_mat[this_ind] = tmp/u_mat[u_piv];
			}
		}
	}
	
	return;
}

void LUMat::lu_solve(vector<double>& soln_vec, vector<double>& rhs, bool transpose) {
	int i1;
	int i2;
	int i3;
	int max_col;
	int max_row;
	int u_piv;
	double tmp;
	
	if(!transpose) {
		for (i1 = 0; i1 < dim; i1++) {
			tmp = rhs[i1];
			i2 = l_min_col[i1];
			for (i3 = l_range[i1]; i3 < l_range[i1+1]; i3++) {
				tmp -= l_mat[i3]*z_vec[i2];
				i2++;
			}
			z_vec[i1] = tmp;
		}
		for (i1 = (dim-1); i1 >= 0; i1--) {
			tmp = z_vec[i1];
			max_col = i1 + max_bandwidth;
			if(max_col > dim) {
				max_col = dim;
			}
			for (i2 = (i1+1); i2 < max_col; i2++) {
				if(u_min_row[i2] <= i1) {
					i3 = u_range[i2] + i1 - u_min_row[i2];
					tmp -= u_mat[i3]*soln_vec[i2];
				}
			}
			u_piv = u_range[i1+1] - 1;
			soln_vec[i1] = tmp/u_mat[u_piv];
		}
	}
	else {
		for (i1 = 0; i1 < dim; i1++) {
			tmp = rhs[i1];
			i2 = u_min_row[i1];
			for (i3 = u_range[i1]; i3 < (u_range[i1+1]-1); i3++) {
				tmp -= u_mat[i3]*z_vec[i2];
				i2++;
			}
			u_piv = u_range[i1+1] - 1;
			z_vec[i1] = tmp/u_mat[u_piv];
		}
		for (i1 = (dim-1); i1 >= 0; i1--) {
			tmp = z_vec[i1];
			max_row = i1 + max_bandwidth;
			if(max_row > dim) {
				max_row = dim;
			}
			for (i2 = (i1+1); i2 < max_row; i2++) {
				if(l_min_col[i2] <= i1) {
 					i3 = l_range[i2] + i1 - l_min_col[i2];
                    tmp -= l_mat[i3]*soln_vec[i2];
				}
			}
			soln_vec[i1] = tmp;
		}
	}
	
	return;
}