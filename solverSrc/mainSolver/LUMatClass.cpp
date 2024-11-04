#include "LUMatClass.h"
#include "ListEntClass.h"
#include "ConstraintClass.h"

using namespace std;

LUMat::LUMat() {
	lMat = nullptr;
	lRange = nullptr;
	lMinCol = nullptr;
	lSize = 0;
	uMat = nullptr;
	uRange = nullptr;
	uMinRow = nullptr;
	uSize = 0;
	zVec = nullptr;
	dim = 0;
	maxBandwidth = 0;
	allocated = false;
	return;
}

void LUMat::setDim(int newDim) {
	dim = newDim;
	lRange = new int[newDim+1];
	lMinCol = new int[newDim];
	uRange = new int[newDim+1];
	uMinRow = new int[newDim];
	zVec = new double[newDim];
	return;
}

void LUMat::allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int currBlock;
	int blkMaxCol;
	int blkMinCol;
	
	MatrixEnt* mEnt1;
	MatrixEnt* mEnt2;
	Constraint* thisConst;
	
	if(!lMat) {
		setDim(spMat.getDim());
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		lRange[i1] = 0;
		uRange[i1] = 1;
		currBlock = i1/blockDim;
		blkMinCol = currBlock*blockDim;
		blkMaxCol = blkMinCol + blockDim;
		mEnt1 = spMat.getFirstEnt(i1);
		while(mEnt1) {
			i2 = mEnt1->col;
			if(i2 >= blkMinCol && i2 < blkMaxCol) {
				if(i2 < i1) { // in L
				    i3 = i1 - i2;
					if(i3 > lRange[i1]) {
						lRange[i1] = i3;
					}
				}
				else {
					i3 = i2 - i1 + 1;
					if(i3 > uRange[i2]) {
						uRange[i2] = i3;
					}
				}
			}
			mEnt1 = mEnt1->nextEnt;
		}
	}
	
	thisConst = cList.getFirst();
	while(thisConst) {
		constDim = thisConst->getMatDim();
		for (i1 = 0; i1 < constDim; i1++) {
			mEnt1 = thisConst-> getMatFirst(i1);
			while(mEnt1) {
				i2 = mEnt1->col;
				currBlock = i2/blockDim;
				blkMinCol = currBlock*blockDim;
				blkMaxCol = blkMinCol + blockDim;
				mEnt2 = thisConst->getMatFirst(i1);
				while(mEnt2) {
					i3 = mEnt2->col;
					if(i3 >= blkMinCol && i3 < blkMaxCol) {
						if(i3 < i2) { // in LU
						    i4 = i2 - i3;
							if(i4 > lRange[i2]) {
								lRange[i2] = i4;
							}
						}
						else {
							i4 = i3 - i2 + 1;
							if(i4 > uRange[i3]) {
								uRange[i3] = i4;
							}
						}
					}
					mEnt2 = mEnt2->nextEnt;
				}
				mEnt1 = mEnt1->nextEnt;
			}
		}
		thisConst = thisConst->getNext();
	}
	
	lSize = 0;
	uSize = 0;
	maxBandwidth = 0;
	for (i1 = 0; i1 < dim; i1++) {
		lSize += lRange[i1];
		uSize += uRange[i1];
		lMinCol[i1] = i1 - lRange[i1];
		if(lRange[i1] > maxBandwidth) {
			maxBandwidth = lRange[i1];
		}
		uMinRow[i1] = i1 - uRange[i1] + 1;
		if(uRange[i1] > maxBandwidth) {
			maxBandwidth = uRange[i1];
		}
	}
	lRange[dim] = lSize;
	uRange[dim] = uSize;
	for (i1 = (dim-1); i1 >= 0; i1--) {
		lRange[i1] = lRange[i1+1] - lRange[i1];
		uRange[i1] = uRange[i1+1] - uRange[i1];
	}
	
	lMat = new double[lSize];
	uMat = new double[uSize];
	
	allocated = true;
	
	cout << "maxBandwidth: " << maxBandwidth << endl;
	
	return;
}

bool LUMat::isAllocated() {
	return allocated;
}

void LUMat::populateFromSparseMat(SparseMat& spMat, ConstraintList& cList) {
	int i1;
	int i2;
	int i3;
	int i4;
	int constDim;
	double constSF;
	
	MatrixEnt* mEnt1;
	MatrixEnt* mEnt2;
	Constraint* thisConst;
	
	for (i1 = 0; i1 < lSize; i1++) {
		lMat[i1] = 0.0;
	}
	
	for (i1 = 0; i1 < uSize; i1++) {
		uMat[i1] = 0.0;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		mEnt1 = spMat.getFirstEnt(i1);
		while(mEnt1) {
			i2 = mEnt1->col;
			if(i2 >= lMinCol[i1] && i2 < i1) { // in L
				i3 = lRange[i1] + i2 - lMinCol[i1];
			    lMat[i3] += mEnt1->value;
			}
			else if(i1 >= uMinRow[i2] && i1 <= i2) {
				i3 = uRange[i2] + i1 - uMinRow[i2];
				uMat[i3] += mEnt1->value;
			}
			mEnt1 = mEnt1->nextEnt;
		}
	}
	
	thisConst = cList.getFirst();
	while(thisConst) {
		constDim = thisConst->getMatDim();
		constSF = thisConst->getScaleFact();
		for (i1 = 0; i1 < constDim; i1++) {
			mEnt1 = thisConst-> getMatFirst(i1);
			while(mEnt1) {
				i2 = mEnt1->col;
				mEnt2 = thisConst->getMatFirst(i1);
				while(mEnt2) {
					i3 = mEnt2->col;
					if(i3 >= lMinCol[i2] && i3 < i2) { // in LU
						i4 = lRange[i2] + i3 - lMinCol[i2];
						lMat[i4] += constSF*mEnt1->value*mEnt2->value;
					}
					else if(i2 >= uMinRow[i3] && i2 <= i3) {
						i4 = uRange[i3] + i2 - uMinRow[i3];
						uMat[i4] += constSF*mEnt1->value*mEnt2->value;
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

void LUMat::luFactor() {
	int i1;
	int i2;
	int i3;
	int thisInd;
	int lInd;
	int uInd;
	int uPiv;
	int maxCol;
	int maxRow;
	double tmp;
	
	for (i1 = 0; i1 < dim; i1++) {
		maxCol = i1 + maxBandwidth;
		if(maxCol > dim) {
			maxCol = dim;
		}
		for (i2 = i1; i2 < maxCol; i2++) {
			if(uMinRow[i2] <= i1) {
				thisInd = uRange[i2] + i1 - uMinRow[i2];
				i3 = uMinRow[i2];
				if(lMinCol[i1] > i3) {
					i3 = lMinCol[i1];
				}
				lInd = lRange[i1] + i3 - lMinCol[i1];
				uInd = uRange[i2] + i3 - uMinRow[i2];
				tmp = uMat[thisInd];
				while(lInd < lRange[i1+1]) {
					tmp -= lMat[lInd]*uMat[uInd];
					lInd++;
					uInd++;
				}
				uMat[thisInd] = tmp;
			}
		}
		maxRow = maxCol;
		for (i2 = (i1+1); i2 < maxRow; i2++) {
			if(lMinCol[i2] <= i1) {
				thisInd = lRange[i2] + i1 - lMinCol[i2];
				i3 = lMinCol[i2];
				if(uMinRow[i1] > i3) {
					i3 = uMinRow[i1];
				}
				lInd = lRange[i2] + i3 - lMinCol[i2];
				uInd = uRange[i1] + i3 - uMinRow[i1];
				tmp = lMat[thisInd];
				uPiv = uRange[i1+1] - 1;
				while(uInd < uPiv) {
					tmp -= lMat[lInd]*uMat[uInd];
					lInd++;
					uInd++;
				}
				lMat[thisInd] = tmp/uMat[uPiv];
			}
		}
	}
	
	return;
}

void LUMat::luSolve(double solnVec[], double rhs[], bool transpose) {
	int i1;
	int i2;
	int i3;
	int maxCol;
	int maxRow;
	int uPiv;
	double tmp;
	
	if(!transpose) {
		for (i1 = 0; i1 < dim; i1++) {
			tmp = rhs[i1];
			i2 = lMinCol[i1];
			for (i3 = lRange[i1]; i3 < lRange[i1+1]; i3++) {
				tmp -= lMat[i3]*zVec[i2];
				i2++;
			}
			zVec[i1] = tmp;
		}
		for (i1 = (dim-1); i1 >= 0; i1--) {
			tmp = zVec[i1];
			maxCol = i1 + maxBandwidth;
			if(maxCol > dim) {
				maxCol = dim;
			}
			for (i2 = (i1+1); i2 < maxCol; i2++) {
				if(uMinRow[i2] <= i1) {
					i3 = uRange[i2] + i1 - uMinRow[i2];
					tmp -= uMat[i3]*solnVec[i2];
				}
			}
			uPiv = uRange[i1+1] - 1;
			solnVec[i1] = tmp/uMat[uPiv];
		}
	}
	else {
		for (i1 = 0; i1 < dim; i1++) {
			tmp = rhs[i1];
			i2 = uMinRow[i1];
			for (i3 = uRange[i1]; i3 < (uRange[i1+1]-1); i3++) {
				tmp -= uMat[i3]*zVec[i2];
				i2++;
			}
			uPiv = uRange[i1+1] - 1;
			zVec[i1] = tmp/uMat[uPiv];
		}
		for (i1 = (dim-1); i1 >= 0; i1--) {
			tmp = zVec[i1];
			maxRow = i1 + maxBandwidth;
			if(maxRow > dim) {
				maxRow = dim;
			}
			for (i2 = (i1+1); i2 < maxRow; i2++) {
				if(lMinCol[i2] <= i1) {
 l					i3 = lRange[i2] + i1 - lMinCol[i2];
                    tmp -= lMat[i3]*solnVec[i2];
				}
			}
			solnVec[i1] = tmp;
		}
	}
	
	return;
}

LUMat::~LUMat() {
	if (dim > 0) {
		delete[] lRange;
		delete[] lMinCol;
		delete[] uRange;
		delete[] uMinRow;
		delete[] zVec;
	}
	if(lMat) {
		delete[] lMat;
		delete[] uMat;
	}
}