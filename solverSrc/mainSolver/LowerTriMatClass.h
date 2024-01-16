#ifndef LOWERTRIMAT
#define LOWERTRIMAT
#include <fstream>
#include "ListEntClass.h"
#include "ConstraintClass.h"

class LowerTriMat {
	private:
	    double *mat;
		int *range;
		int *minCol;
		double *zVec;
		int dim;
		int size;
		int maxBandwidth;
		bool allocated;
	
	public:
        LowerTriMat();
		
		void setDim(int newDim);
		
		void allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim);
		
		bool isAllocated();
		
		void populateFromSparseMat(SparseMat& spMat, ConstraintList& cList);

		void copyToFullMat(double matFP[]);
		
		void ldlFactor();
		
		void ldlSolve(double solnVec[], double rhs[]);

		void writeToFile(std::ofstream& outFile);
		
		void destroy();
};

#endif