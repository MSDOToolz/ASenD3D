#ifndef lowertrimat
#define lowertrimat
#include <fstream>
#include <vector>
#include "ListEntClass.h"
#include "ConstraintClass.h"

class LowerTriMat {
	public:
	    std::vector<double> mat;
		std::vector<int> range;
		std::vector<int> min_col;
		std::vector<double> z_vec;
		int dim;
		int size;
		int max_bandwidth;
		bool allocated;
	
        LowerTriMat();
		
		void set_dim(int new_dim);
		
		void allocate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list, int block_dim);
		
		bool is_allocated();
		
		void populate_from_sparse_mat(SparseMat& sp_mat, ConstraintList& c_list);

		//void copy_to_full_mat(double mat_fp[]);
		
		void ldl_factor();
		
		void ldl_solve(std::vector<double>& soln_vec, std::vector<double>& rhs);

		bool pos_def();

		void write_to_file(std::ofstream& out_file);
};

#endif