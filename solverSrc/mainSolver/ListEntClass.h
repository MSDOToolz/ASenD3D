#ifndef listent
#define listent
#include <string>
#include <fstream>
#include <vector>
#include <list>

class IDCapsule {
public:
	int int_dat;
	double doub_dat;

	IDCapsule();
};

class MatrixEnt {
    public:
	    int row;
		int col;
		double value;
		
	    MatrixEnt();
};

class MatrixRow {
public:
	std::list<MatrixEnt> row_vec;

	MatrixRow();

	void add_entry(int row, int col, double val);
};

class SparseMat {
	public:
	    int dim;
	    std::vector<MatrixRow> matrix;
		
	    SparseMat();
		
		void set_dim(int new_dim);
		
		void zero_all();
		
		void add_entry(int row, int col, double val);

		void add_matrix(SparseMat& inp_mat);
		
		void vector_multiply(std::vector<double>& prod, std::vector<double>& inp_vec, bool transpose);

		double get_max_abs_val();

		void write_to_file(std::ofstream& out_file);
};
#endif