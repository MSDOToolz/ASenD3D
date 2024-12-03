#ifndef LUMAT
#define LUMAT
#include "ListEntClass.h"
#include "ConstraintClass.h"
#include <vector>

class LUMat {
public:
    std::vector<double> lMat;
	std::vector<int> lRange;
	std::vector<int> lMinCol;
	int lSize;
	std::vector<double> uMat;
	std::vector<int> uRange;
	std::vector<int> uMinRow;
	int uSize;
	std::vector<double> zVec;
	int dim;
	int maxBandwidth;
	bool allocated;
	
    LUMat();
	
	void setDim(int newDim);
	
	void allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim);
	
	void populateFromSparseMat(SparseMat& spMat, ConstraintList& cList);
	
	void luFactor();
	
	void luSolve(std::vector<double>& solnVec, std::vector<double>& rhs, bool transpose);
};

#endif