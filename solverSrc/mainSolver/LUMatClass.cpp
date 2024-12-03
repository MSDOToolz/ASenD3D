#include <iostream>
#include <iomanip>
#include <vector>
#include "LUMatClass.h"
#include "ListEntClass.h"
#include "ConstraintClass.h"

using namespace std;

LUMat::LUMat() {
	lMat.clear();
	lRange.clear();
	lMinCol.clear();
	lSize = 0;
	uMat.clear();
	uRange.clear();
	uMinRow.clear();
	uSize = 0;
	zVec.clear();
	dim = 0;
	maxBandwidth = 0;
	allocated = false;
	return;
}

void LUMat::setDim(int newDim) {
	dim = newDim;
	lRange = vector<int>(newDim+1);
	lMinCol = vector<int>(newDim);
	uRange = vector<int>(newDim+1);
	uMinRow = vector<int>(newDim);
	zVec = vector<double>(newDim);
	return;
}

void LUMat::allocateFromSparseMat(SparseMat& spMat, ConstraintList& cList, int blockDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int currBlock;
	int constDim;
	int blkMaxCol;
	int blkMinCol;
	
	if(!allocated) {
		setDim(spMat.dim);
	}

	for (i1 = 0; i1 < dim; i1++) {
		lRange[i1] = 0;
		uRange[i1] = 1;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		currBlock = i1/blockDim;
		blkMinCol = currBlock*blockDim;
		blkMaxCol = blkMinCol + blockDim;
		MatrixRow& mr = spMat.matrix[i1];
		for (auto& me : mr.rowVec) {
			i2 = me.col;
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
		}
	}
	
	for (auto& cnst : cList.constVec) {
		SparseMat& thisMat = cnst.mat;
		for (auto& mr : thisMat.matrix) {
			for (auto& me : mr.rowVec) {
				i2 = me.col;
				currBlock = i2/blockDim;
				blkMinCol = currBlock*blockDim;
				blkMaxCol = blkMinCol + blockDim;
				for (auto& me2 : mr.rowVec) {
					i3 = me2.col;
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
				}
			}
		}
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
	
	lMat = vector<double>(lSize);
	uMat = vector<double>(uSize);
	
	allocated = true;
	
	cout << "maxBandwidth: " << maxBandwidth << endl;
	
	return;
}

void LUMat::populateFromSparseMat(SparseMat& spMat, ConstraintList& cList) {
	int i1;
	int i2;
	int i3;
	int i4;
	int constDim;
	double constSF;
	
	for (auto& lm : lMat) {
		lm = 0.0;
	}
	
	for (auto& um : uMat) {
		um = 0.0;
	}
	
	for (i1 = 0; i1 < dim; i1++) {
		MatrixRow& mr = spMat.matrix[i1];
		for (auto& me : mr.rowVec) {
			i2 = me.col;
			if(i2 >= lMinCol[i1] && i2 < i1) { // in L
				i3 = lRange[i1] + i2 - lMinCol[i1];
			    lMat[i3] += me.value;
			}
			else if(i1 >= uMinRow[i2] && i1 <= i2) {
				i3 = uRange[i2] + i1 - uMinRow[i2];
				uMat[i3] += me.value;
			}
		}
	}
	
	for (auto& cnst : cList.constVec) {
		SparseMat& thisMat = cnst.mat;
		constSF = cnst.scaleFact;
		for (auto& mr : thisMat.matrix) {
			for (auto& me : mr.rowVec) {
				i2 = me.col;
				for (auto& me2 : mr.rowVec) {
					i3 = me2.col;
					if(i3 >= lMinCol[i2] && i3 < i2) { // in LU
						i4 = lRange[i2] + i3 - lMinCol[i2];
						lMat[i4] += constSF*me.value*me2.value;
					}
					else if(i2 >= uMinRow[i3] && i2 <= i3) {
						i4 = uRange[i3] + i2 - uMinRow[i3];
						uMat[i4] += constSF*me.value*me2.value;
					}
				}
			}
		}
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

void LUMat::luSolve(vector<double>& solnVec, vector<double>& rhs, bool transpose) {
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
 					i3 = lRange[i2] + i1 - lMinCol[i2];
                    tmp -= lMat[i3]*solnVec[i2];
				}
			}
			solnVec[i1] = tmp;
		}
	}
	
	return;
}