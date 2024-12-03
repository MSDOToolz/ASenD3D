#ifndef LISTENT
#define LISTENT
#include <string>
#include <fstream>
#include <vector>
#include <list>

class IDCapsule {
public:
	int intDat;
	double doubDat;

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
	std::list<MatrixEnt> rowVec;

	MatrixRow();

	void addEntry(int row, int col, double val);
};

class SparseMat {
	public:
	    int dim;
	    std::vector<MatrixRow> matrix;
		
	    SparseMat();
		
		void setDim(int newDim);
		
		void zeroAll();
		
		void addEntry(int row, int col, double val);

		void addMatrix(SparseMat& inpMat);
		
		void vectorMultiply(std::vector<double>& prod, std::vector<double>& inpVec, bool transpose);

		double getMaxAbsVal();

		void writeToFile(std::ofstream& outFile);
};
#endif