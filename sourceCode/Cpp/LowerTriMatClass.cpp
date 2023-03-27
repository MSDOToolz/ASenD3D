#include "LowerTriMatClass.h"
#include "ListEntClass.h"
#include "ConstraintClass.h"
#include <iostream>

using namespace std;

LowerTriMat::LowerTriMat() {
	mat = NULL;
	range = NULL;
	minCol = NULL;
	zVec = NULL;
	allocated = false;
	return;
}

void LowerTriMat::setDim(int newDim) {
	dim = newDim;
	range = new int[newDim+1];
	minCol = new int[newDim];
	zVec = new double[newDim];
	return;
}


void LowerTriMat::allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int maxBw) {
	int i1;
	int i2;
	int i3;
	int i4;
	int constDim;
	
	MatrixEnt *mEnt1;
	MatrixEnt *mEnt2;
	Constraint *thisConst;
	
	for (i1 = 0; i1 < dim; i1++) {
		range[i1] = 1;
		mEnt1 = spMat.getFirstEnt(i1);
		while(mEnt1) {
			i2 = mEnt1->col;
			if(i2 <= i1) {
				i3 = i1 - i2 + 1;
				if(i3 > range[i1] && i3 <= maxBw) {
					range[i1] = i3;
				}
			}
			mEnt1 = mEnt1->nextEnt;
		}
	}
	
	thisConst = cList.getFirst();
	while(thisConst) {
		constDim = thisConst->getMatDim();
		for (i1 = 0; i1 < constDim; i1++) {
			mEnt1 = thisConst->getMatFirst(i1);
			while(mEnt1) {
				i2 = mEnt1->col;
				mEnt2 = thisConst->getMatFirst(i1);
				while(mEnt2) {
					i3 = mEnt2->col;
					if(i3 <= i2) {
						i4 = i2 - i3 + 1;
						if(i4 > range[i2] && i4 <= maxBw) {
							range[i2] = i4;
						}
					}
					mEnt2 = mEnt2->nextEnt;
				}
				mEnt1 = mEnt1->nextEnt;
			}
		}
		thisConst = thisConst->getNext();
	}

	size = 0;
	maxBandwidth = 0;
	for (i1 = 0; i1 < dim; i1++) {
		size+= range[i1];
		minCol[i1] = i1 - range[i1] + 1;
		if(range[i1] > maxBandwidth) {
			maxBandwidth = range[i1];
		}
	}
	range[dim] = size;
	for (i1 = (dim-1); i1 >= 0; i1--) {
		range[i1] = range[i1+1] - range[i1];
	}
	
	mat = new double[size];
	
	allocated = true;
	
	return;
}

bool LowerTriMat::isAllocated() {
	return allocated;
}

void LowerTriMat::populateFromSparseMat(SparseMat& spMat, ConstraintList& cList) {
	int i1;
	int i2;
	int i3;
	int i4;
	int constDim;
	
	MatrixEnt *mEnt1;
	MatrixEnt *mEnt2;
	Constraint *thisConst;
	
	for (i1 = 0; i1 < size; i1++) {
		mat[i1] = 0.0;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		mEnt1 = spMat.getFirstEnt(i1);
		while(mEnt1) {
			i2 = mEnt1->col;
			if(i2 <= i1 && i2 >= minCol[i1]) {
				i3 = range[i1] + (i2 - minCol[i1]);
				mat[i3]+= mEnt1->value;
			}
			mEnt1 = mEnt1->nextEnt;
		}
	}
	
	thisConst = cList.getFirst();
	while(thisConst) {
		constDim = thisConst->getMatDim();
		for (i1 = 0; i1 < constDim; i1++) {
			mEnt1 = thisConst->getMatFirst(i1);
			while(mEnt1) {
				i2 = mEnt1->col;
				mEnt2 = thisConst->getMatFirst(i1);
				while(mEnt2) {
					i3 = mEnt2->col;
					if(i3 <= i2 && i3 >= minCol[i2]) {
						i4 = range[i2] + (i3 - minCol[i2]);
						mat[i4]+= mEnt1->value*mEnt2->value;
					}
					mEnt2 = mEnt2->nextEnt;
				}
				mEnt1 = mEnt1->nextEnt;
			}
		}
		thisConst = thisConst->getNext();
	}		
    	
	return;
}

void LowerTriMat::ldlFactor() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int stCol;
	double sum;
	
	double *ldVec = new double[dim];
	
	for (i1 = 0; i1 < dim; i1++) { // i1 = column in L
		//form ldVec
		i3 = range[i1];
		for (i2 = minCol[i1]; i2 < i1; i2++) {
			i4 = range[i2+1] - 1;
			ldVec[i2] = mat[i3]*mat[i4];
			i3++;
		}
		//get d term for column i1
		stCol = minCol[i1];				
		i2 = range[i1];
		sum = 0.0;
		for (i3 = stCol; i3 < i1; i3++) {
			sum+= ldVec[i3]*mat[i2];
			i2++;
		}
		i2 = range[i1+1] - 1;
		mat[i2] = mat[i2] - sum;
		//get L terms for column i1
		i2 = i1 + 1;
		while(i2 < i1 + maxBandwidth && i2 < dim) { // i2 = row in L
		    if(minCol[i2] <= i1) {
				stCol = minCol[i2];
				if(minCol[i1] > stCol) {
					stCol = minCol[i1];
				}
				i4 = range[i2] + (stCol - minCol[i2]);
				sum = 0.0;
				for (i5 = stCol; i5 < i1; i5++) {
					sum+= ldVec[i5]*mat[i4];
					i4++;
				}
				i3 = range[i2] + (i1 - minCol[i2]);
				i4 = range[i1+1] - 1;
				mat[i3] = (mat[i3] - sum)/mat[i4];
			}
			i2++;
		}
	}
	
	delete[] ldVec;
	return;
}

void LowerTriMat::ldlSolve(double solnVec[], double rhs[]) {
	int i1;
	int i2;
	int i3;
	int stopRow;
	double sum;
	
	zVec[0] = rhs[0];
	for (i1 = 1; i1 < dim; i1++) {
		i2 = range[i1];
		sum = 0.0;
		for (i3 = minCol[i1]; i3 < i1; i3++) {
			sum+= mat[i2]*zVec[i3];
			i2++;
		}
		zVec[i1] = rhs[i1] - sum;
	}
	
	for (i1 = dim-1; i1 >= 0; i1--) {
		sum = 0.0;
		stopRow = i1 + maxBandwidth + 1;
		if(stopRow > dim) {
			stopRow = dim;
		}
		for (i2 = (i1+1); i2 < stopRow; i2++) {
			i3 = range[i2] + (i1 - minCol[i2]);
			if(i3 >= range[i2] && i3 < range[i2+1]) {
				sum+= mat[i3]*solnVec[i2];
			}
		}
		i2 = range[i1+1] - 1;
		solnVec[i1] = (zVec[i1]/mat[i2]) - sum;
	}
	return;
}

void LowerTriMat::destroy() {
	delete[] mat;
	delete[] range;
	delete[] minCol;
	delete[] zVec;
	return;
}