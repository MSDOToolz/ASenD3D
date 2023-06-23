#include <iostream>
#include "matrixFunctions.h"
#include "ListEntClass.h"
#include "LowerTriMatClass.h"

using namespace std;

int main() {
	double mat1[9] = {2.0,-1.0,-0.5,-1.0,3.0,-1.0,-0.5,-1.0,4.0};
	qRFactor(mat1,3,0,2,0,2,0);
	
	cout << "Factored Matrix" << endl;
	
	cout << mat1[0] << ", " << mat1[1] << ", " << mat1[2] << endl;
	cout << mat1[3] << ", " << mat1[4] << ", " << mat1[5] << endl;
	cout << mat1[6] << ", " << mat1[7] << ", " << mat1[8] << endl;
	
	double xVec[3];
	double bVec[3] = {1.0,0.0,0.0};
	
	solveqRxEqb(xVec,mat1,bVec,3,0,2,0,2,0);
	
	cout << "xVec solution: " << xVec[0] << ", " << xVec[1] << ", " << xVec[2] << endl;
	cout << "Check: 0.6769, 0.2769, 0.1538" << endl;
	
 	mat1[0] = 2.0;
	mat1[1] = -1.0;
	mat1[2] = -0.5;
	mat1[3] = -1.0;
	mat1[4] = 3.0;
	mat1[5] = -1.0;
	mat1[6] = -0.5;
	mat1[7] = -1.0;
	mat1[8] = 4.0;
	
	double eVecs[9];
	double eVals[3];
	
	eigenSolve(eVals, eVecs, mat1, 3, 0);
	
	cout << "Eigenvalues: " << endl;
	
	cout << eVals[0] << ", " << eVals[1] << ", " << eVals[2] << endl;
	
	cout << "Check: 1.0596, 3.3175, 4.6229" << endl;
	
	cout << "Eigenvectors: " << endl;
	
	cout << eVecs[0] << ", " << eVecs[1] << ", " << eVecs[2] << endl;
	cout << eVecs[3] << ", " << eVecs[4] << ", " << eVecs[5] << endl;
	cout << eVecs[6] << ", " << eVecs[7] << ", " << eVecs[8] << endl;
	
	cout << "Check:" << endl;
	
	cout << "0.7648, 0.6425, 0.0484" << endl;
	cout << "0.5591, -0.6244, -0.5454" << endl;
	cout << "0.3202, -0.4442, 0.8368" << endl; 
	
	Doub mat2[9];
	mat2[0].setVal(2.0);
	mat2[1].setVal(-1.0);
	mat2[2].setVal(-0.5);
	mat2[3].setVal(-1.0);
	mat2[4].setVal(3.0);
	mat2[5].setVal(-1.0);
	mat2[6].setVal(-0.5);
	mat2[7].setVal(-1.0);
	mat2[8].setVal(4.0);
	
	Doub det;
	Doub invMat2[9];
	Doub dbX[3];
	Doub dbB[3];
	getDetInv(det,invMat2,mat2,3,0,dbX,dbB);
	
	cout << "Inverse Matrix: " << endl;
	
	cout << invMat2[0].val << ", " << invMat2[1].val << ", " << invMat2[2].val << endl;
	cout << invMat2[3].val << ", " << invMat2[4].val << ", " << invMat2[5].val << endl;
	cout << invMat2[6].val << ", " << invMat2[7].val << ", " << invMat2[8].val << endl;
	
	cout << "Check: " << endl;
	cout << "0.6769, 0.2769, 0.1538" << endl;
	cout << "0.2769, 0.4769, 0.1538" << endl;
	cout << "0.1538, 0.1538, 0.3077" << endl;
	
	mat2[0].setVal(2.0);
	mat2[1].setVal(-1.0);
	mat2[2].setVal(-0.5);
	mat2[3].setVal(-1.0);
	mat2[4].setVal(3.0);
	mat2[5].setVal(-1.0);
	mat2[6].setVal(-0.5);
	mat2[7].setVal(-1.0);
	mat2[8].setVal(4.0);
	
	Doub prod[9];
	
	matMul(prod,mat2,invMat2,3,3,3);
	
	cout << "Matrix Product: " << endl;
	
	cout << prod[0].val << ", " << prod[1].val << ", " << prod[2].val << endl;
	cout << prod[3].val << ", " << prod[4].val << ", " << prod[5].val << endl;
	cout << prod[6].val << ", " << prod[7].val << ", " << prod[8].val << endl;
	
// Test with sparse and lower triangular matrices	
	SparseMat spMat;
	spMat.setDim(3);
	spMat.addEntry(0,0,2.0);
	spMat.addEntry(1,0,-1.0);
	spMat.addEntry(1,1,3.0);
	spMat.addEntry(2,0,-0.5);
	spMat.addEntry(2,1,-1.0);
	spMat.addEntry(2,2,4.0);
	
	ConstraintList cList;
	
	LowerTriMat ltMat;
	
	ltMat.setDim(3);
	ltMat.allocateFromSparseMat(spMat,cList,3);
	ltMat.populateFromSparseMat(spMat,cList);
	
    ltMat.ldlFactor();
	bVec[0] = 1.0;
	bVec[1] = 0.0;
	bVec[2] = 0.0;
	ltMat.ldlSolve(xVec,bVec);
	
	cout << "xVec solution: " << xVec[0] << ", " << xVec[1] << ", " << xVec[2] << endl;
	cout << "Check: 0.6769, 0.2769, 0.1538" << endl;
	
	spMat.destroy();
	cList.destroy();
	ltMat.destroy();
	
	return 0;
}
