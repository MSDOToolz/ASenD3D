#ifndef constraints
#define constraints
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <fstream>

class ConstraintTerm {
	public:
	    std::string node_set;
		int ns_ptr;
		int dof;
		double coef;
		
	    ConstraintTerm();
		
};

class Constraint {
	public:
		std::string type;
		std::list<ConstraintTerm> terms;
		double rhs;
		SparseMat mat;
		double scale_fact;
		
	    Constraint();
		
		void build_mat(std::vector<Node>& nd_ar, std::vector<Set>& set_ar);

		void full_vec_multiply(std::vector<double>& prod, std::vector<double>& vec, std::vector<double>& tmp_v);

		void get_load(std::vector<double>& c_ld, std::vector<double>& u_vec, std::vector<double>& q_vec,int res_dim);

		void write_to_file(std::ofstream& out_file);
};

class ConstraintList {
	public:
		std::vector<Constraint> const_vec;
	
	    ConstraintList();

		void set_scale_fact(double new_sf);

		void get_total_vec_mult(std::vector<double>& prod, std::vector<double>& vec, std::vector<double>& tmp_v);

		void get_total_load(std::vector<double>& c_ld, std::vector<double>& u_vec, std::vector<double>& q_vec, int res_dim);

		void write_all_to_file(std::string file_name);
};

#endif