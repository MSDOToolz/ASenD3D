#ifndef LUMAT
#define LUMAT
#include "ListEntClass.h"
#include "ConstraintClass.h"

class LUMat {
private:
    double* lMat;
	int* lRange;
	int* lMinCol;
	int lSize;
	double* uMat;
	int* uRange;
	int* uMinRow;
	int uSize;
	double* zVec;
	int dim;
	int maxBandwidth;
	bool allocated;
	
public:
    LUMat();
	
	void setDim(int newDim);
	
	void allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim);
	
	bool isAllocated();
	
	void populateFromSparseMat(SparseMat& spMat, ConstraintList& cList);
	
	void luFactor();
	
	void luSolve(double solnVec[], double rhs[], bool transpose);
	
	~LUMat();
};

#endif