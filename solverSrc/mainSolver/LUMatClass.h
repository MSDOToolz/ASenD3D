#ifndef lumat
#define lumat
#include "ListEntClass.h"
#include "ConstraintClass.h"
#include <vector>

class LUMat {
public:
    std::vector<double> l_mat;
	std::vector<int> l_range;
	std::vector<int> l_min_col;
	int l_size;
	std::vector<double> u_mat;
	std::vector<int> u_range;
	std::vector<int> u_min_row;
	int u_size;
	std::vector<double> z_vec;
	int dim;
	int max_bandwidth;
	bool allocated;
	
    LUMat();
	
	void set_dim(int new_dim);
	
	void allocate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list, int block_dim);
	
	void populate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list);
	
	void lu_factor();
	
	void lu_solve(std::vector<double>& soln_vec, std::vector<double>& rhs, bool transpose);
};

#endif