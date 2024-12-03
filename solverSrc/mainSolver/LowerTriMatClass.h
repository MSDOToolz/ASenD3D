#ifndef LOWERTRIMAT
#define LOWERTRIMAT
#include <fstream>
#include <vector>
#include "ListEntClass.h"
#include "ConstraintClass.h"

class LowerTriMat {
	public:
	    std::vector<double> mat;
		std::vector<int> range;
		std::vector<int> minCol;
		std::vector<double> zVec;
		int dim;
		int size;
		int maxBandwidth;
		bool allocated;
	
        LowerTriMat();
		
		void setDim(int newDim);
		
		void allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim);
		
		bool isAllocated();
		
		void populateFromSparseMat(SparseMat& spMat, ConstraintList& cList);

		//void copyToFullMat(double matFP[]);
		
		void ldlFactor();
		
		void ldlSolve(std::vector<double>& solnVec, std::vector<double>& rhs);

		bool posDef();

		void writeToFile(std::ofstream& outFile);
};

#endif