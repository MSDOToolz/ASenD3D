#include "matrixFunctions.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LowerTriMatClass.h"
#include "ConstraintClass.h"

using namespace std;

const double pi_2 = 1.57079632679489662;
const double tol = 1.0e-12;
const double magtol = 1.0e-8;
const double maxMag = 1.0e+12;
const double minMag = 1.0e-12;

double getDist(double p1[], double p2[]) {
	double dist;
	double vec[3];
	vec[0] = p1[0] - p2[0];
	vec[1] = p1[1] - p2[1];
	vec[2] = p1[2] - p2[2];
	dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	return dist;
}

void crossProd(double prod[], double v1[], double v2[]) {
	prod[0] = v1[1]*v2[2] - v1[2]*v2[1];
	prod[1] = v1[2]*v2[0] - v1[0]*v2[2];
	prod[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return;
}

void qRFactor(double mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Min;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	double theta;
	double sth;
	double cth;
	double p1;
	double p2;

	for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k11 = (i2Min-1)*colDim + i1;
			k12 = i2*colDim + i1;
			if(abs(mat[k11]) < tol) {
				mat[k11] = tol;
			}
			theta = atan(mat[k12]/mat[k11]);
			sth = sin(theta);
			cth = cos(theta);
			i3Min = i1;
			if(triDiag == 2) {
				i3Max = i1 + 2;
				if(i3Max > endCol) {
					i3Max = endCol;
				}
			} else {
				i3Max = endCol;
			}
			for (i3 = i3Min; i3 <= i3Max; i3++) {
				k22 = (i2Min-1)*colDim + i3;
				k23 = i2*colDim + i3;
				p1 = cth*mat[k22] + sth*mat[k23];
				p2 = -sth*mat[k22] + cth*mat[k23];
				mat[k22] = p1;
				mat[k23] = p2;
			}
			mat[k12] = theta;
		}
	}
	return;
}

void solveqRxEqb(double xVec[], double mat[], double bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int k11;
	int k12;
	double theta;
	double sth;
	double cth;
	double p1;
	double p2;
	double rowSum;

    for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		i3 = i2Min - 1;
		for (i2 = i2Min; i2 <=i2Max; i2++) {
			k12 = i2*colDim + i1;
			theta = mat[k12];
			sth = sin(theta);
			cth = cos(theta);
			p1 = cth*bVec[i3] + sth*bVec[i2];
			p2 = -sth*bVec[i3] + cth*bVec[i2];
		    bVec[i3] = p1;
			bVec[i2] = p2;
		}
		xVec[i1] = 0.0;
	}
	
	for (i1 = endCol; i1 >= stCol; i1--) {
		i3 = stRow + (i1 - stCol);
		i2Min = i1 + 1;
		if(triDiag == 2) {
			i2Max = i1 + 2;
			if(i2Max > endCol) {
				i2Max = endCol;
			}
		} else {
			i2Max = endCol;
		}
		rowSum = 0.0;
		k11 = i3*colDim + i2Min;
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			rowSum+= mat[k11]*xVec[i2];
			k11++;
		}
		k11 = i3*colDim + i1;
		xVec[i1] = (bVec[i3] - rowSum)/mat[k11];
	}
	
	return;
}

void gMResSparse(double soln[], SparseMat& mat, ConstraintList& cnst, LowerTriMat& pcMat, double rhs[], double convTol, int maxIt, int restart) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int itCt;
	double resNrm;
	double tmp;
	int dim = mat.getDim();

	double* resVec = new double[dim];
	double* tmpV = new double[dim];
	double* tmpV2 = new double[dim];
	i1 = dim * (restart + 1);
	double* hMat = new double[i1];
	i1 = restart * (restart + 1);
	double* phiMat = new double[i1];

	for (i1 = 0; i1 < dim; i1++) {
		soln[i1] = 0.0;
		tmpV[i1] = -rhs[i1];
	}

	pcMat.ldlSolve(resVec, tmpV);
	resNrm = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		resNrm += resVec[i1] * resVec[i1];
	}
	resNrm = sqrt(resNrm);

	itCt = 0;
	while (resNrm > convTol && itCt < maxIt) {
		tmp = 1.0 / resNrm;
		for (i1 = 0; i1 < dim; i1++) {
			hMat[i1] = tmp * resVec[i1];
		}
		i2 = restart * (restart + 1);
		for (i1 = 0; i1 < i2; i1++) {
			phiMat[i1] = 0.0;
		}
		//Generate basis vectors
		for (i1 = 1; i1 < (restart + 1); i1++) {
			//Multiply by previous vector
			for (i2 = 0; i2 < dim; i2++) {
				tmpV[i2] = 0.0;
			}
			i2 = dim * (i1 - 1);
			mat.vectorMultiply(tmpV, &hMat[i2], false);
			cnst.getTotalVecMult(tmpV, &hMat[i2], tmpV2);
			pcMat.ldlSolve(tmpV2, tmpV);
			//Orthogonalize with all previous vectors
			for (i2 = 0; i2 < i1; i2++) {
				i4 = dim * i2;
				tmp = 0.0;
				for (i3 = 0; i3 < dim; i3++) {
					tmp += tmpV2[i3] * hMat[i4];
					i4++;
				}
				i5 = i2 * restart + (i1 - 1);
				phiMat[i5] = tmp;
				i4 = dim * i2;
				for (i3 = 0; i3 < dim; i3++) {
					tmpV2[i3] -= tmp * hMat[i4];
					i4++;
				}
			}
			//Calculate magnitude and set new unit basis vector
			tmp = 0.0;
			for (i3 = 0; i3 < dim; i3++) {
				tmp += tmpV2[i3] * tmpV2[i3];
			}
			tmp = sqrt(tmp);
			i5 = i1 * restart + (i1 - 1);
			phiMat[i5] = tmp;
			tmp = 1.0 / tmp;
			i4 = dim * i1;
			for (i2 = 0; i2 < dim; i2++) {
				hMat[i4] = tmp * tmpV2[i2];
				i4++;
			}
		}
		//Find the least squares solution
		i3 = 0;
		for (i1 = 0; i1 < restart; i1++) {
			tmpV[i1] = 0.0;
			for (i2 = 0; i2 < dim; i2++) {
				tmpV[i1] -= hMat[i3] * resVec[i2];
				i3++;
			}
		}
		qRFactor(phiMat, restart, 0, (restart-1), 0, (restart-1), 1);
		solveqRxEqb(tmpV2, phiMat, tmpV, restart, 0, (restart-1), 0, (restart-1), 1);
		//Update the solution vector
		i3 = 0;
		for (i1 = 0; i1 < restart; i1++) {
			for (i2 = 0; i2 < dim; i2++) {
				soln[i2] += hMat[i3] * tmpV2[i1];
				i3++;
			}
		}
		//Update residual vector
		for (i1 = 0; i1 < dim; i1++) {
			tmpV[i1] = -rhs[i1];
		}
		mat.vectorMultiply(tmpV, soln, false);
		cnst.getTotalVecMult(tmpV, soln, tmpV2);
		pcMat.ldlSolve(resVec, tmpV);
		resNrm = 0.0;
		for (i1 = 0; i1 < dim; i1++) {
			resNrm += resVec[i1] * resVec[i1];
		}
		resNrm = sqrt(resNrm);
		itCt += restart;
	}

	if (resNrm > convTol) {
		cout << "Warning: GMRES solver did not converge to the requested tolerance of " << convTol << endl;
		cout << "Residual norm after " << itCt << " iterations: " << resNrm << endl;
	}

	delete[] resVec;
	delete[] tmpV;
	delete[] tmpV2;
	delete[] hMat;
	delete[] phiMat;

	return;
}

void symFactor(double mat[], double qMat[], int matDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double theta;
	double sth;
	double cth;
	double a1;
	double a2;
	
	i3 = 0;
	for (i1 = 0; i1 < matDim; i1++) {
		for (i2 = 0; i2 < matDim; i2++) {
			if(i2 == i1) {
				qMat[i3] = 1.0;
			} else {
				qMat[i3] = 0.0;
			}
			i3++;
		}
	}
	
	for (i1 = 0; i1 < matDim; i1++) {
		for (i2 = i1+2; i2 < matDim; i2++) {
			i3 = (i1+1)*matDim + i1;
			if(abs(mat[i3]) < tol) {
				theta = pi_2;
			} else {
				i4 = i2*matDim + i1;
				theta = atan(mat[i4]/mat[i3]);
			}
			sth = sin(theta);
			cth = cos(theta);
			i4 = (i1+1)*matDim + i1;
			i5 = i2*matDim + i1;
			for (i3 = i1; i3 < matDim; i3++) {
				a1 = cth*mat[i4] + sth*mat[i5];
				a2 = -sth*mat[i4] + cth*mat[i5];
				mat[i4] = a1;
				mat[i5] = a2;
				i4++;
				i5++;
 			}
			i4 = i1 + 1;
			i5 = i2;
			for (i3 = 0; i3 < matDim; i3++) {
				a1 = cth*mat[i4] + sth*mat[i5];
				a2 = -sth*mat[i4] + cth*mat[i5];
				mat[i4] = a1;
				mat[i5] = a2;
				a1 = cth*qMat[i4] + sth*qMat[i5];
				a2 = -sth*qMat[i4] + cth*qMat[i5];
				qMat[i4] = a1;
				qMat[i5] = a2;
				i4+= matDim;
				i5+= matDim;
			}
		}
	}		
	
	return;
}

void getCharFun(DiffDoub1& cFun, DiffDoub1 mat[], int matDim, DoubList& eVals, double lam, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Min;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;

	for (i1 = 0; i1 <= (matDim - 1); i1++) {
		i2Min = i1 + 1;
		if(triDiag == 0) {
			i2Max = matDim - 1;
		} else {
			i2Max = i1 + 1;
			if(i2Max > (matDim - 1)) {
				i2Max = matDim - 1;
			}
		}
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k11 = (i2Min-1)*matDim + i1;
			k12 = i2*matDim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.setVal(mat[k12]);
			tmp.setVal(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.setVal(theta);
			sth.sn();
			cth.setVal(theta);
			cth.cs();
			i3Min = i1;
			if(triDiag == 2) {
				i3Max = i1 + 2;
				if(i3Max > (matDim - 1)) {
					i3Max = matDim - 1;
				}
			} else {
				i3Max = matDim - 1;
			}
			for (i3 = i3Min; i3 <= i3Max; i3++) {
				k22 = (i2Min-1)*matDim + i3;
				k23 = i2*matDim + i3;
				p1.setVal(cth);
				p1.mult(mat[k22]);
				tmp.setVal(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.setVal(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.setVal(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].setVal(p1);
				mat[k23].setVal(p2);
			}
			mat[k12].setVal(theta);
		}
	}
	
	cFun.setVal(mat[0]);
	i2 = matDim + 1;
	for (i1 = 1; i1 < matDim; i1++) {
		cFun.mult(mat[i2]);
		while(abs(cFun.val) > maxMag) {
			cFun.val*= minMag;
			cFun.dval*= minMag;
		}
		while(abs(cFun.val) < minMag) {
			cFun.val*= maxMag;
			cFun.dval*= maxMag;
		}
		i2+= (matDim + 1);
	}
	
	DoubListEnt *thisEval = eVals.getFirst();
	double term;
	while(thisEval) {
		term = thisEval->value - lam;
		tmp.setVal(term,-1.0);
		cFun.dvd(tmp);
		while(abs(cFun.val) > maxMag) {
			cFun.val*= minMag;
			cFun.dval*= minMag;
		}
		while(abs(cFun.val) < minMag) {
			cFun.val*= maxMag;
			cFun.dval*= maxMag;
		}
		thisEval = thisEval->next;
	}
	
	return;
}

void getEvals(double eVals[], double mat[], int matDim, double lamInit, double convTol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int matSize = matDim*matDim;
	DoubList eValList;
	DiffDoub1 *matCopy = new DiffDoub1[matSize];
	DiffDoub1 cFun;
	double lam = lamInit;
	double dLam;
	double dTrm;
	double backStep = 1.0e+6*convTol;
	int numFound = 0;
	int maxIt = 100*matDim;
	int it = 0;
	while(numFound < matDim && it < maxIt) {
		i3 = 0;
		for (i1 = 0; i1 < matDim; i1++) {
			for (i2 = 0; i2 < matDim; i2++) {
				if(i1 == i2) {
					dTrm = mat[i3] - lam;
					matCopy[i3].setVal(dTrm,-1.0);
				} else {
					matCopy[i3].setVal(mat[i3]);
				}
				i3++;
			}
		}
		getCharFun(cFun,matCopy,matDim,eValList,lam,triDiag);
		dLam = -cFun.val/cFun.dval;
		lam+= dLam;
		if(abs(dLam) < convTol) {
			eValList.addEntry(lam);
			numFound++;
			lam-= backStep;
		}
		it++;
	}
	DoubListEnt *thisEnt = eValList.getFirst();
	i1 = 0;
	while(thisEnt) {
		eVals[i1] = thisEnt->value;
		thisEnt = thisEnt->next;
		i1++;
	}
	delete[] matCopy;
	return;
}

double getLowEval(double mat[], int matDim) {
	double *vec = new double[matDim];
	double *prevVec = new double[matDim];
	double high;
	double tmp;
	double mag;
	double dp;
	int i1;
	int i2;
	int i3;
	int it;
	int maxIt;
	
	mag = 0.0;
	for (i1 = 0; i1 < matDim; i1++) {
		tmp = sin(1.0*i1);
		prevVec[i1] = tmp;
		mag+= tmp*tmp;
	}
	mag = 1.0/sqrt(mag);
	for (i1 = 0; i1 < matDim; i1++) {
		prevVec[i1]*= mag;
	}
	
	bool conv = false;
	maxIt = 100*matDim;
	it = 0;
	while(!conv && it < maxIt) {
		i3 = 0;
		mag = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			vec[i1] = 0.0;
			for (i2 = 0; i2 < matDim; i2++) {
				vec[i1]+= mat[i3]*prevVec[i2];
				i3++;
			}
			mag+= vec[i1]*vec[i1];
		}
		mag = sqrt(mag);
		dp = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			dp+= vec[i1]*prevVec[i1];
		}
		if(abs(dp) > 0.99999999*mag) {
			conv = true;
		} else {
			mag = 1.0/mag;
			for (i1 = 0; i1 < matDim; i1++) {
				prevVec[i1] = mag*vec[i1];
			}
		}
		it++;
	}
	
	if(dp < 0.0) {
		delete[] vec;
	    delete[] prevVec;
		return dp;
	}
	
	high = dp;
    
    i3 = 0;	
	for (i1 = 0; i1 < matDim; i1++) {
		mat[i3]-= high;
		i3+= (matDim + 1);
	}
	
	mag = 0.0;
	for (i1 = 0; i1 < matDim; i1++) {
		tmp = sin(1.0*i1);
		prevVec[i1] = tmp;
		mag+= tmp*tmp;
	}
	mag = 1.0/sqrt(mag);
	for (i1 = 0; i1 < matDim; i1++) {
		prevVec[i1]*= mag;
	}
	
	conv = false;
	it = 0;
	while(!conv && it < maxIt) {
		i3 = 0;
		mag = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			vec[i1] = 0.0;
			for (i2 = 0; i2 < matDim; i2++) {
				vec[i1]+= mat[i3]*prevVec[i2];
				i3++;
			}
			mag+= vec[i1]*vec[i1];
		}
		mag = sqrt(mag);
		dp = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			dp+= vec[i1]*prevVec[i1];
		}
		if(abs(dp) > 0.99999999*mag) {
			conv = true;
		} else {
			mag = 1.0/mag;
			for (i1 = 0; i1 < matDim; i1++) {
				prevVec[i1] = mag*vec[i1];
			}
		}
		it++;
	}
	
	i3 = 0;	
	for (i1 = 0; i1 < matDim; i1++) {
		mat[i3]+= high;
		i3+= (matDim + 1);
	}
	
	dp+= high;
	delete[] vec;
	delete[] prevVec;
	return dp;
}

void eigenSolve(double eVals[], double eVecs[], double mat[], int matDim, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double tmp;
	int matSize = matDim*matDim;
	int newTd = triDiag;
	double* qMat = nullptr;
	if(triDiag == 0) {
		qMat = new double[matSize];
		symFactor(mat,qMat,matDim);
		newTd = 1;
	}
	
    double low = getLowEval(mat,matDim);
	
	double matMag = 0.0;
	for (i1 = 0; i1 < matSize; i1++) {
		tmp = abs(mat[i1]);
		if(tmp > matMag) {
			matMag = tmp;
		}
	}
	
	low-= 0.001*matMag;
	double convTol = 1.0e-12*matMag;
	getEvals(eVals,mat,matDim,low,convTol,newTd);
	
	//DoubListEnt *thisVal = eVals.getFirst();
	double *matCopy = new double[matSize];
	double *vec = new double[matDim];
	double *prevVec = new double[matDim];
	double *bVec = new double[matDim];
	double shift;
	bool conv;
	double mag;
	double dp;
	int it;
	int maxIt = 100*matDim;
	for (i5 = 0; i5 < matDim; i5++) {
		shift = eVals[i5] + 1.0e-8*matMag;
		i3 = 0;
		for (i1 = 0; i1 < matDim; i1++) {
			for (i2 = 0; i2 < matDim; i2++) {
				if (i1 == i2) {
					matCopy[i3] = mat[i3] - shift;
				}
				else {
					matCopy[i3] = mat[i3];
				}
				i3++;
			}
		}
		qRFactor(matCopy, matDim, 0, (matDim - 1), 0, (matDim - 1), newTd);
		mag = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			tmp = sin(1.0 * i1 * (i5 + 1));
			prevVec[i1] = tmp;
			mag += tmp * tmp;
		}
		mag = 1.0 / sqrt(mag);
		for (i1 = 0; i1 < matDim; i1++) {
			prevVec[i1] = mag * prevVec[i1];
		}
		conv = false;
		it = 0;
		while (!conv && it < maxIt) {
			for (i1 = 0; i1 < matDim; i1++) {
				bVec[i1] = prevVec[i1];
			}
			solveqRxEqb(vec, matCopy, bVec, matDim, 0, (matDim - 1), 0, (matDim - 1), newTd);
			mag = 0.0;
			dp = 0.0;
			for (i1 = 0; i1 < matDim; i1++) {
				mag += vec[i1] * vec[i1];
				dp += vec[i1] * prevVec[i1];
			}
			mag = sqrt(mag);
			if (abs(dp) > 0.999999 * mag) {
				conv = true;
				if (triDiag == 0) {
					i3 = 0;
					//i4 = i5;
					i4 = i5 * matDim;
					for (i1 = 0; i1 < matDim; i1++) {
						//i4 = i5*matDim + i1
						eVecs[i4] = 0.0;
						for (i2 = 0; i2 < matDim; i2++) {
							eVecs[i4] += qMat[i3] * prevVec[i2];
							i3++;
						}
						//i4 += matDim;
						i4++;
					}
				}
				else {
					//i4 = i5;
					i4 = i5 * matDim;
					for (i1 = 0; i1 < matDim; i1++) {
						//i4 = i5*matDim + i1;
						eVecs[i4] = prevVec[i1];
						//i4 += matDim;
						i4++;
					}
				}
			}
			else {
				mag = 1.0 / mag;
				for (i1 = 0; i1 < matDim; i1++) {
					prevVec[i1] = mag * vec[i1];
				}
			}
			it++;
		}
	}

	if (triDiag == 0) {
		delete[] qMat;
	}

	delete[] matCopy;
	delete[] vec;
	delete[] prevVec;
	delete[] bVec;
	return;
}

void symEigenSolve(double eVals[], double eVecs[], double mat[], int matDim, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i3Min;
	int i3Max;
	int i4;
	int loopCt;
	double tmp;
	double mag;
	double dp;

	double* tempV1 = new double[matDim];
	double* tempV2 = new double[matDim];
	double* eVtmp = new double[matDim * matDim];
	int* srtOrder = new int[matDim];

	for (i1 = 0; i1 < matDim; i1++) {
		for (i2 = 0; i2 < matDim; i2++) {
			tmp = sin(1.0 * i2 * (i1+1));
			tempV1[i2] = tmp;
		}
		for (i2 = 0; i2 < i1; i2++) {
			dp = 0.0;
			i4 = i2 * matDim;
			for (i3 = 0; i3 < matDim; i3++) {
				dp += tempV1[i3] * eVtmp[i4];
				i4++;
			}
			i4 = i2 * matDim;
			for (i3 = 0; i3 < matDim; i3++) {
				tempV1[i3] -= dp*eVtmp[i4];
				i4++;
			}
		}
		mag = 0.0;
		for (i2 = 0; i2 < matDim; i2++) {
			mag += tempV1[i2] * tempV1[i2];
		}
		mag = 1.0 / sqrt(mag);
		for (i2 = 0; i2 < matDim; i2++) {
			tempV1[i2] *= mag;
		}
		dp = 0.0;
		loopCt = 0;
		while (abs(dp) < 0.99999999*mag && loopCt < 10000) {
			// Multiply by mat
			for (i2 = 0; i2 < matDim; i2++) {
				tempV2[i2] = 0.0;
				if (triDiag == 1) {
					i3Min = i2 - 1;
					if (i3Min < 0) {
						i3Min = 0;
					}
					i3Max = i2 + 2;
					if (i3Max > matDim) {
						i3Max = matDim;
					}
				}
				else {
					i3Min = 0;
					i3Max = matDim;
				}
				i4 = i2 * matDim + i3Min;
				for (i3 = i3Min; i3 < i3Max; i3++) {
					//i4 = i2 * matDim + i3;
					tempV2[i2] += mat[i4] * tempV1[i3];
					i4++;
				}
			}
			// orthogonalize with previous vectors
			for (i2 = 0; i2 < i1; i2++) {
				dp = 0.0;
				i4 = i2 * matDim;
				for (i3 = 0; i3 < matDim; i3++) {
					dp += tempV2[i3] * eVtmp[i4];
					i4++;
				}
				i4 = i2 * matDim;
				for (i3 = 0; i3 < matDim; i3++) {
					tempV2[i3] -= dp * eVtmp[i4];
					i4++;
				}
			}
			// get mag, dp
			mag = 0.0;
			dp = 0.0;
			for (i2 = 0; i2 < matDim; i2++) {
				mag += tempV2[i2] * tempV2[i2];
				dp += tempV1[i2] * tempV2[i2];
			}
			mag = sqrt(mag);
			// update tempV1
			tmp = 1.0 / mag;
			for (i2 = 0; i2 < matDim; i2++) {
				tempV1[i2] = tmp * tempV2[i2];
			}
			loopCt++;
		}
		eVals[i1] = dp;
		i3 = i1 * matDim;
		for (i2 = 0; i2 < matDim; i2++) {
			eVtmp[i3] = tempV1[i2];
			i3++;
		}
	}

	// sort pairs

	for (i1 = 0; i1 < matDim; i1++) {
		srtOrder[i1] = i1;
	}

	for (i1 = 0; i1 < matDim; i1++) {
		for (i2 = 0; i2 < (matDim - 1); i2++) {
			if (eVals[i2 + 1] < eVals[i2]) {
				tmp = eVals[i2];
				eVals[i2] = eVals[i2 + 1];
				eVals[i2 + 1] = tmp;
				i3 = srtOrder[i2];
				srtOrder[i2] = srtOrder[i2 + 1];
				srtOrder[i2 + 1] = i3;
			}
		}
	}

	for (i1 = 0; i1 < matDim; i1++) {
		i2 = matDim * i1;
		i3 = matDim * srtOrder[i1];
		for (i4 = 0; i4 < matDim; i4++) {
			eVecs[i2] = eVtmp[i3];
			i2++;
			i3++;
		}
	}

	delete[] tempV1;
	delete[] tempV2;
	delete[] eVtmp;
	delete[] srtOrder;

	return;
}

void eigenFull(double eVals[], double eVecs[], int numPairs, LowerTriMat& mat, double massMat[], int matDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double tmp;
	double mag;

	i3 = matDim * matDim;
	double* matFP = new double[i3];
	double* qMat = new double[i3];
	double* eVectmp = new double[i3];
	double* eValtmp = new double[matDim];
	double* vtmp = new double[matDim];

	for (i1 = 0; i1 < i3; i1++) {
		matFP[i1] = 0.0;
	}

	mat.copyToFullMat(matFP);

	for (i1 = 0; i1 < matDim; i1++) {
		massMat[i1] = 1.0 / sqrt(massMat[i1]);
		i3 = i1 * matDim;
		i4 = i1;
		for (i2 = 0; i2 < matDim; i2++) {
			matFP[i3] *= massMat[i1];
			matFP[i4] *= massMat[i1];
			i3++;
			i4 += matDim;
		}
	}

	//symFactor(matFP, qMat, matDim);

	symEigenSolve(eValtmp, eVectmp, matFP, matDim, 1);

	for (i1 = 0; i1 < numPairs; i1++) {
		i4 = 0;
		i5 = i1 * matDim;
		for (i2 = 0; i2 < matDim; i2++) {
			/*vtmp[i2] = 0.0;
			i5 = i1 * matDim;
			for (i3 = 0; i3 < matDim; i3++) {
				vtmp[i2] += qMat[i4] * eVectmp[i5];
				i4++;
				i5++;
			}*/
			vtmp[i2] = eVectmp[i5];
			i5++;
		}
		mag = 0.0;
		for (i2 = 0; i2 < matDim; i2++) {
			tmp = vtmp[i2] * massMat[i2];
			vtmp[i2] = tmp;
			mag += tmp * tmp;
		}
		mag = 1.0/sqrt(mag);
		i3 = i1 * matDim;
		for (i2 = 0; i2 < matDim; i2++) {
			eVecs[i3] = mag * vtmp[i2];
			i3++;
		}
		eVals[i1] = eValtmp[i1];
	}

	for (i1 = 0; i1 < matDim; i1++) {
		tmp = 1.0 / massMat[i1];
		massMat[i1] = tmp * tmp;
	}

	delete[] matFP;
	delete[] qMat;
	delete[] eVectmp;
	delete[] eValtmp;
	delete[] vtmp;

	return;
}

void eigenSparseDirect(double eVals[], double eVecs[], int numPairs, LowerTriMat& mat, double massMat[], int matDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	int hCols = 5*numPairs;
	i1 = hCols*matDim;
	double *hMat = new double[i1];
	i1 = hCols*(hCols - 1);
	double *coefMat = new double[i1];
	
	for (i2 = 0; i2 < i1; i2++) {
		coefMat[i2] = 0.0;
	}
	
	double mag = 0.0;
	double tmp;
	for (i2 = 0; i2 < matDim; i2++) {
		massMat[i2] = sqrt(massMat[i2]);
	}

	for (i2 = matDim; i2 < (2 * matDim); i2++) {
		hMat[i2] = sin(1.0 * i2);
	}

	mat.ldlSolve(&hMat[0], &hMat[matDim]);

	mag = 0.0;
	for (i2 = 0; i2 < matDim; i2++) {
		mag += hMat[i2] * hMat[i2];
	}
	mag = 1.0 / sqrt(mag);

	for (i2 = 0; i2 < matDim; i2++) {
		hMat[i2] *= mag;
	}
	
	double *tVec1 = new double[matDim];
	double *tVec2 = new double[matDim];
	double dp;
	double stVec;
	for (i1 = 1; i1 < hCols; i1++) {
		i3 = matDim*(i1-1);
		for (i2 = 0; i2 < matDim; i2++) {
			tVec1[i2] = massMat[i2]*hMat[i3];
			i3++;
		}
		mat.ldlSolve(tVec2,tVec1);
		for (i2 = 0; i2 < matDim; i2++) {
			tVec1[i2] = tVec2[i2]*massMat[i2];
		}
		if(i1 == 1) {
			stVec = 0;
		} else {
			stVec = i1 - 2;
		}
		for (i2 = stVec; i2 < i1; i2++) {
			dp = 0.0;
			i3 = matDim*i2;
			for (i4 = 0; i4 < matDim; i4++) {
				dp+= tVec1[i4]*hMat[i3];
				i3++;
			}
			i3 = matDim*i2;
			for (i4 = 0; i4 < matDim; i4++) {
				tVec1[i4]-= dp*hMat[i3];
				i3++;
			}
			i3 = (hCols-1)*i2 + (i1-1);
			coefMat[i3] = dp;
		}
		mag = 0.0;
		for (i2 = 0; i2 < matDim; i2++) {
			mag+= tVec1[i2]*tVec1[i2];
		}
		mag = sqrt(mag);
		i3 = hCols*i1 - 1; // (hCols-1)*i1 + (i1 - 1)
		coefMat[i3] = mag;
		mag = 1.0/mag;
		i3 = matDim*i1;
		for (i2 = 0; i2 < matDim; i2++) {
			hMat[i3] = mag*tVec1[i2];
			i3++;
		}
	}
	
	i1 = hCols - 1;
	double *coefVals = new double[i1];
	double *coefVecs = new double[i1*i1];
	//eigenSolve(coefVals, coefVecs, coefMat, i1, 2);
	symEigenSolve(coefVals, coefVecs, coefMat, i1, 1);

	for (i3 = 0; i3 < matDim; i3++) {
		massMat[i3] = 1.0/massMat[i3];
	}
	
	i1 = 0;
	i2 = hCols - 2;
	while (i1 < numPairs && i2 >= 0) {
		for (i3 = 0; i3 < matDim; i3++) {
			tVec1[i3] = 0.0;
		}
		i5 = 0;
		for (i3 = 0; i3 < (hCols-1); i3++) {
			//i6 = (hCols-1)*i3 + i2;
			i6 = i2 * (hCols - 1) + i3;
			for (i4 = 0; i4 < matDim; i4++) {
				tVec1[i4] += hMat[i5] * coefVecs[i6];
				i5++;
			}
		}
		mag = 0.0;
		for (i3 = 0; i3 < matDim; i3++) {
			tmp = tVec1[i3]*massMat[i3];
			tVec2[i3] = tmp;
			mag+= tmp*tmp;
		}
		mag = 1.0/sqrt(mag);
		for (i3 = 0; i3 < matDim; i3++) {
			tVec2[i3] *= mag;
		}
		if (i1 > 0) {
			tmp = 0.0;
			i4 = matDim * (i1 - 1);
			for (i3 = 0; i3 < matDim; i3++) {
				tmp += eVecs[i4] * tVec2[i3];
				i4++;
			}
			if (abs(tmp) < 0.9) {
				for (i3 = 0; i3 < matDim; i3++) {
					eVecs[i4] = tVec2[i3];
					i4++;
				}
				eVals[i1] = 1.0 / coefVals[i2];
				i1++;
				i2--;
			}
			else {
				i2--;
			}
		}
		else {
			i4 = 0;
			for (i3 = 0; i3 < matDim; i3++) {
				//i4 = i1 * matDim + i3;
				eVecs[i4] = tVec2[i3];
				i4++;
			}
			eVals[i1] = 1.0 / coefVals[i2];
			i1++;
			i2--;
		}
	}
	
	for (i1 = 0; i1 < matDim; i1++) {
		tmp = 1.0/massMat[i1];
		massMat[i1] = tmp*tmp;
	}
	
	delete[] hMat;
	delete[] coefMat;
	delete[] tVec1;
	delete[] tVec2;
	delete[] coefVals;
	delete[] coefVecs;
	
	return;
}

//dup1
void qRFactor(DiffDoub0 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Min;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;

	for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k11 = (i2Min-1)*colDim + i1;
			k12 = i2*colDim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.setVal(mat[k12]);
			tmp.setVal(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.setVal(theta);
			sth.sn();
			cth.setVal(theta);
			cth.cs();
			i3Min = i1;
			if(triDiag == 2) {
				i3Max = i1 + 2;
				if(i3Max > endCol) {
					i3Max = endCol;
				}
			} else {
				i3Max = endCol;
			}
			for (i3 = i3Min; i3 <= i3Max; i3++) {
				k22 = (i2Min-1)*colDim + i3;
				k23 = i2*colDim + i3;
				p1.setVal(cth);
				p1.mult(mat[k22]);
				tmp.setVal(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.setVal(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.setVal(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].setVal(p1);
				mat[k23].setVal(p2);
			}
			mat[k12].setVal(theta);
		}
	}
	return;
}

void solveqRxEqb(DiffDoub0 xVec[], DiffDoub0 mat[], DiffDoub0 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int k11;
	int k12;
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;
	DiffDoub0 rowSum;

    for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		i3 = i2Min - 1;
		for (i2 = i2Min; i2 <=i2Max; i2++) {
			k12 = i2*colDim + i1;
			theta.setVal(mat[k12]);
			sth.setVal(theta);
			sth.sn();
			cth.setVal(theta);
			cth.cs();
			p1.setVal(cth);
			p1.mult(bVec[i3]);
			tmp.setVal(sth);
			tmp.mult(bVec[i2]);
			p1.add(tmp);
			p2.setVal(sth);
			p2.neg();
			p2.mult(bVec[i3]);
			tmp.setVal(cth);
			tmp.mult(bVec[i2]);
			p2.add(tmp);
		    bVec[i3].setVal(p1);
			bVec[i2].setVal(p2);
		}
		xVec[i1].setVal(0.0);
	}
	
	
	for (i1 = endCol; i1 >= stCol; i1--) {
		i3 = stRow + (i1 - stCol);
		i2Min = i1 + 1;
		if(triDiag == 2) {
			i2Max = i1 + 2;
			if(i2Max > endCol) {
				i2Max = endCol;
			}
		} else {
			i2Max = endCol;
		}
		rowSum.setVal(0.0);
		k11 = i3*colDim + i2Min;
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			tmp.setVal(mat[k11]);
			tmp.mult(xVec[i2]);
			rowSum.add(tmp);
			k11++;
		}
		tmp.setVal(bVec[i3]);
		tmp.sub(rowSum);
		k11 = i3*colDim + i1;
		tmp.dvd(mat[k11]);
		xVec[i1].setVal(tmp);
	}
	
	return;
}

void getDetInv(DiffDoub0& det, DiffDoub0 inv[], DiffDoub0 mat[], int colDim, int triDiag, DiffDoub0 xVec[], DiffDoub0 bVec[]) {
	qRFactor(mat, colDim, 0, (colDim-1), 0, (colDim-1), triDiag);
	int i1;
    int i2;
	int i3;
	det.setVal(1.0);
	for (i1 = 0; i1 < colDim; i1++) {
		for (i2 = 0; i2 < colDim; i2++) {
			if(i1 == i2) {
				bVec[i2].setVal(1.0);
			} else {
			    bVec[i2].setVal(0.0);
			}
		}
		solveqRxEqb(xVec,mat,bVec,colDim,0,(colDim-1),0,(colDim-1),triDiag);
		for (i2 = 0; i2 < colDim; i2++) {
			i3 = i2*colDim + i1;
			inv[i3].setVal(xVec[i2]);
		}
		i3 = i1*colDim + i1;
		det.mult(mat[i3]);
	}
	return;
}

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
void qRFactor(DiffDoub1 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Min;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;

	for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k11 = (i2Min-1)*colDim + i1;
			k12 = i2*colDim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.setVal(mat[k12]);
			tmp.setVal(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.setVal(theta);
			sth.sn();
			cth.setVal(theta);
			cth.cs();
			i3Min = i1;
			if(triDiag == 2) {
				i3Max = i1 + 2;
				if(i3Max > endCol) {
					i3Max = endCol;
				}
			} else {
				i3Max = endCol;
			}
			for (i3 = i3Min; i3 <= i3Max; i3++) {
				k22 = (i2Min-1)*colDim + i3;
				k23 = i2*colDim + i3;
				p1.setVal(cth);
				p1.mult(mat[k22]);
				tmp.setVal(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.setVal(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.setVal(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].setVal(p1);
				mat[k23].setVal(p2);
			}
			mat[k12].setVal(theta);
		}
	}
	return;
}

void solveqRxEqb(DiffDoub1 xVec[], DiffDoub1 mat[], DiffDoub1 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int k11;
	int k12;
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;
	DiffDoub1 rowSum;

    for (i1 = stCol; i1 <= endCol; i1++) {
		i2Min = stRow + (i1 - stCol) + 1;
		if(triDiag == 0) {
			i2Max = endRow;
		} else {
			i2Max = stRow + (i1 - stCol) + 1;
			if(i2Max > endRow) {
				i2Max = endRow;
			}
		}
		i3 = i2Min - 1;
		for (i2 = i2Min; i2 <=i2Max; i2++) {
			k12 = i2*colDim + i1;
			theta.setVal(mat[k12]);
			sth.setVal(theta);
			sth.sn();
			cth.setVal(theta);
			cth.cs();
			p1.setVal(cth);
			p1.mult(bVec[i3]);
			tmp.setVal(sth);
			tmp.mult(bVec[i2]);
			p1.add(tmp);
			p2.setVal(sth);
			p2.neg();
			p2.mult(bVec[i3]);
			tmp.setVal(cth);
			tmp.mult(bVec[i2]);
			p2.add(tmp);
		    bVec[i3].setVal(p1);
			bVec[i2].setVal(p2);
		}
		xVec[i1].setVal(0.0);
	}
	
	
	for (i1 = endCol; i1 >= stCol; i1--) {
		i3 = stRow + (i1 - stCol);
		i2Min = i1 + 1;
		if(triDiag == 2) {
			i2Max = i1 + 2;
			if(i2Max > endCol) {
				i2Max = endCol;
			}
		} else {
			i2Max = endCol;
		}
		rowSum.setVal(0.0);
		k11 = i3*colDim + i2Min;
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			tmp.setVal(mat[k11]);
			tmp.mult(xVec[i2]);
			rowSum.add(tmp);
			k11++;
		}
		tmp.setVal(bVec[i3]);
		tmp.sub(rowSum);
		k11 = i3*colDim + i1;
		tmp.dvd(mat[k11]);
		xVec[i1].setVal(tmp);
	}
	
	return;
}

void getDetInv(DiffDoub1& det, DiffDoub1 inv[], DiffDoub1 mat[], int colDim, int triDiag, DiffDoub1 xVec[], DiffDoub1 bVec[]) {
	qRFactor(mat, colDim, 0, (colDim-1), 0, (colDim-1), triDiag);
	int i1;
    int i2;
	int i3;
	det.setVal(1.0);
	for (i1 = 0; i1 < colDim; i1++) {
		for (i2 = 0; i2 < colDim; i2++) {
			if(i1 == i2) {
				bVec[i2].setVal(1.0);
			} else {
			    bVec[i2].setVal(0.0);
			}
		}
		solveqRxEqb(xVec,mat,bVec,colDim,0,(colDim-1),0,(colDim-1),triDiag);
		for (i2 = 0; i2 < colDim; i2++) {
			i3 = i2*colDim + i1;
			inv[i3].setVal(xVec[i2]);
		}
		i3 = i1*colDim + i1;
		det.mult(mat[i3]);
	}
	return;
}

//end dup
 
//end skip 
 
 
 
 
 
 
//dup2
void matMul(DiffDoub0 prod[], DiffDoub0 mat1[], DiffDoub0 mat2[], int m1Rows, int m1Cols, int m2Cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub0 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1Rows; i1++) {
		for (i2 = 0; i2 < m2Cols; i2++) {
			i6 = i2;
			prod[i4].setVal(0.0);
			for (i3 = 0; i3 < m1Cols; i3++) {
				tmp.setVal(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2Cols;
			}
			i5 -= m1Cols;
			i4++;
		}
		i5 += m1Cols;
	}
	return;
}

void transpose(DiffDoub0 matT[], DiffDoub0 mat[], int rowDim, int colDim) {
	DiffDoub0 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < rowDim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < colDim; i2++) {
			matT[i4].setVal(mat[i3]);
			i3++;
			i4+= rowDim;
		}
	}
	return;
}

void crossProd(DiffDoub0 prod[], DiffDoub0 v1[], DiffDoub0 v2[]) {
	DiffDoub0 tmp;
	
	prod[0].setVal(v1[1]);
	prod[0].mult(v2[2]);
	tmp.setVal(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].setVal(v1[2]);
	prod[1].mult(v2[0]);
	tmp.setVal(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].setVal(v1[0]);
	prod[2].mult(v2[1]);
	tmp.setVal(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotateOrient(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[]) {
	DiffDoub0 mag;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 oneHalf;
	DiffDoub0 a1[9];
	DiffDoub0 a2[9];
	DiffDoub0 a3[9];
	int i1;
	int i2;
	int i3;
	
	oneHalf.setVal(0.5);
	
	mag.setVal(rot[0]);
	mag.sqr();
	tmp.setVal(rot[1]);
	tmp.sqr();
	tmp2.setVal(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub0 locRot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			locRot[i1].setVal(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.setVal(locOri[i3]);
				tmp.mult(rot[i2]);
				locRot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].setVal(1.0);
		tmp.setVal(locRot[1]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[0].sub(tmp);
		
		a1[1].setVal(0.5);
		a1[1].mult(locRot[0]);
		a1[1].mult(locRot[1]);
		a1[1].add(locRot[2]);
		
		a1[2].setVal(0.5);
		a1[2].mult(locRot[0]);
		a1[2].mult(locRot[2]);
		a1[2].sub(locRot[1]);
		
		a1[3].setVal(0.5);
		a1[3].mult(locRot[0]);
		a1[3].mult(locRot[1]);
		a1[3].sub(locRot[2]);
		
		a1[4].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[4].sub(tmp);
		
		a1[5].setVal(0.5);
		a1[5].mult(locRot[1]);
		a1[5].mult(locRot[2]);
		a1[5].add(locRot[0]);
		
		a1[6].setVal(0.5);
		a1[6].mult(locRot[0]);
		a1[6].mult(locRot[2]);
		a1[6].add(locRot[1]);
		
		a1[7].setVal(0.5);
		a1[7].mult(locRot[1]);
		a1[7].mult(locRot[2]);
		a1[7].sub(locRot[0]);
		
		a1[8].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[8].sub(tmp);
		
		matMul(instOri, a1, locOri, 3, 3, 3);
	} else {
		DiffDoub0 sth;
		DiffDoub0 cth;
		
		sth.setVal(mag);
		sth.sn();
		cth.setVal(mag);
		cth.cs();
		
		DiffDoub0 unitRot[3];
		tmp.setVal(1.0);
		tmp.dvd(mag);
		unitRot[0].setVal(rot[0]);
		unitRot[0].mult(tmp);
		unitRot[1].setVal(rot[1]);
		unitRot[1].mult(tmp);
		unitRot[2].setVal(rot[2]);
		unitRot[2].mult(tmp);
		
		a1[0].setVal(unitRot[0]);
		a1[1].setVal(unitRot[1]);
		a1[2].setVal(unitRot[2]);
		
		i1 = 0;
		if(abs(unitRot[1].val) < abs(unitRot[0].val)) {
			i1 = 1;
		}
		if(abs(unitRot[2].val) < abs(unitRot[i1].val)) {
			i1 = 2;
		}
		tmp.setVal(1.0);
		tmp2.setVal(unitRot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unitRot[i1]^2)
		a1[3+i1].setVal(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].setVal(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		crossProd(&a1[6],&a1[0],&a1[3]);
		
		a2[0].setVal(a1[0]);
		a2[1].setVal(a1[1]);
		a2[2].setVal(a1[2]);
		
		a2[3].setVal(cth);
		a2[3].mult(a1[3]);
		tmp.setVal(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].setVal(cth);
		a2[4].mult(a1[4]);
		tmp.setVal(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].setVal(cth);
		a2[5].mult(a1[5]);
		tmp.setVal(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].setVal(sth);
		a2[6].mult(a1[3]);
		tmp.setVal(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].setVal(sth);
		a2[7].mult(a1[4]);
		tmp.setVal(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].setVal(sth);
		a2[8].mult(a1[5]);
		tmp.setVal(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose(a3,a2,3,3); // a3 = a2^T
		matMul(a2,a3,a1,3,3,3); //a2 = a2^T*a1
		matMul(instOri,locOri,a2,3,3,3);
	}
	return;
}
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup2
void matMul(DiffDoub1 prod[], DiffDoub1 mat1[], DiffDoub1 mat2[], int m1Rows, int m1Cols, int m2Cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub1 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1Rows; i1++) {
		for (i2 = 0; i2 < m2Cols; i2++) {
			i6 = i2;
			prod[i4].setVal(0.0);
			for (i3 = 0; i3 < m1Cols; i3++) {
				tmp.setVal(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2Cols;
			}
			i5 -= m1Cols;
			i4++;
		}
		i5 += m1Cols;
	}
	return;
}

void transpose(DiffDoub1 matT[], DiffDoub1 mat[], int rowDim, int colDim) {
	DiffDoub1 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < rowDim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < colDim; i2++) {
			matT[i4].setVal(mat[i3]);
			i3++;
			i4+= rowDim;
		}
	}
	return;
}

void crossProd(DiffDoub1 prod[], DiffDoub1 v1[], DiffDoub1 v2[]) {
	DiffDoub1 tmp;
	
	prod[0].setVal(v1[1]);
	prod[0].mult(v2[2]);
	tmp.setVal(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].setVal(v1[2]);
	prod[1].mult(v2[0]);
	tmp.setVal(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].setVal(v1[0]);
	prod[2].mult(v2[1]);
	tmp.setVal(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotateOrient(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[]) {
	DiffDoub1 mag;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 oneHalf;
	DiffDoub1 a1[9];
	DiffDoub1 a2[9];
	DiffDoub1 a3[9];
	int i1;
	int i2;
	int i3;
	
	oneHalf.setVal(0.5);
	
	mag.setVal(rot[0]);
	mag.sqr();
	tmp.setVal(rot[1]);
	tmp.sqr();
	tmp2.setVal(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub1 locRot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			locRot[i1].setVal(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.setVal(locOri[i3]);
				tmp.mult(rot[i2]);
				locRot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].setVal(1.0);
		tmp.setVal(locRot[1]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[0].sub(tmp);
		
		a1[1].setVal(0.5);
		a1[1].mult(locRot[0]);
		a1[1].mult(locRot[1]);
		a1[1].add(locRot[2]);
		
		a1[2].setVal(0.5);
		a1[2].mult(locRot[0]);
		a1[2].mult(locRot[2]);
		a1[2].sub(locRot[1]);
		
		a1[3].setVal(0.5);
		a1[3].mult(locRot[0]);
		a1[3].mult(locRot[1]);
		a1[3].sub(locRot[2]);
		
		a1[4].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[4].sub(tmp);
		
		a1[5].setVal(0.5);
		a1[5].mult(locRot[1]);
		a1[5].mult(locRot[2]);
		a1[5].add(locRot[0]);
		
		a1[6].setVal(0.5);
		a1[6].mult(locRot[0]);
		a1[6].mult(locRot[2]);
		a1[6].add(locRot[1]);
		
		a1[7].setVal(0.5);
		a1[7].mult(locRot[1]);
		a1[7].mult(locRot[2]);
		a1[7].sub(locRot[0]);
		
		a1[8].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[8].sub(tmp);
		
		matMul(instOri, a1, locOri, 3, 3, 3);
	} else {
		DiffDoub1 sth;
		DiffDoub1 cth;
		
		sth.setVal(mag);
		sth.sn();
		cth.setVal(mag);
		cth.cs();
		
		DiffDoub1 unitRot[3];
		tmp.setVal(1.0);
		tmp.dvd(mag);
		unitRot[0].setVal(rot[0]);
		unitRot[0].mult(tmp);
		unitRot[1].setVal(rot[1]);
		unitRot[1].mult(tmp);
		unitRot[2].setVal(rot[2]);
		unitRot[2].mult(tmp);
		
		a1[0].setVal(unitRot[0]);
		a1[1].setVal(unitRot[1]);
		a1[2].setVal(unitRot[2]);
		
		i1 = 0;
		if(abs(unitRot[1].val) < abs(unitRot[0].val)) {
			i1 = 1;
		}
		if(abs(unitRot[2].val) < abs(unitRot[i1].val)) {
			i1 = 2;
		}
		tmp.setVal(1.0);
		tmp2.setVal(unitRot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unitRot[i1]^2)
		a1[3+i1].setVal(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].setVal(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		crossProd(&a1[6],&a1[0],&a1[3]);
		
		a2[0].setVal(a1[0]);
		a2[1].setVal(a1[1]);
		a2[2].setVal(a1[2]);
		
		a2[3].setVal(cth);
		a2[3].mult(a1[3]);
		tmp.setVal(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].setVal(cth);
		a2[4].mult(a1[4]);
		tmp.setVal(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].setVal(cth);
		a2[5].mult(a1[5]);
		tmp.setVal(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].setVal(sth);
		a2[6].mult(a1[3]);
		tmp.setVal(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].setVal(sth);
		a2[7].mult(a1[4]);
		tmp.setVal(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].setVal(sth);
		a2[8].mult(a1[5]);
		tmp.setVal(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose(a3,a2,3,3); // a3 = a2^T
		matMul(a2,a3,a1,3,3,3); //a2 = a2^T*a1
		matMul(instOri,locOri,a2,3,3,3);
	}
	return;
}
//end dup
 
//Diff2Doub versions: 
//dup2
void matMul(DiffDoub2 prod[], DiffDoub2 mat1[], DiffDoub2 mat2[], int m1Rows, int m1Cols, int m2Cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub2 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1Rows; i1++) {
		for (i2 = 0; i2 < m2Cols; i2++) {
			i6 = i2;
			prod[i4].setVal(0.0);
			for (i3 = 0; i3 < m1Cols; i3++) {
				tmp.setVal(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2Cols;
			}
			i5 -= m1Cols;
			i4++;
		}
		i5 += m1Cols;
	}
	return;
}

void transpose(DiffDoub2 matT[], DiffDoub2 mat[], int rowDim, int colDim) {
	DiffDoub2 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < rowDim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < colDim; i2++) {
			matT[i4].setVal(mat[i3]);
			i3++;
			i4+= rowDim;
		}
	}
	return;
}

void crossProd(DiffDoub2 prod[], DiffDoub2 v1[], DiffDoub2 v2[]) {
	DiffDoub2 tmp;
	
	prod[0].setVal(v1[1]);
	prod[0].mult(v2[2]);
	tmp.setVal(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].setVal(v1[2]);
	prod[1].mult(v2[0]);
	tmp.setVal(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].setVal(v1[0]);
	prod[2].mult(v2[1]);
	tmp.setVal(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotateOrient(DiffDoub2 instOri[], DiffDoub2 locOri[], DiffDoub2 rot[]) {
	DiffDoub2 mag;
	DiffDoub2 tmp;
	DiffDoub2 tmp2;
	DiffDoub2 oneHalf;
	DiffDoub2 a1[9];
	DiffDoub2 a2[9];
	DiffDoub2 a3[9];
	int i1;
	int i2;
	int i3;
	
	oneHalf.setVal(0.5);
	
	mag.setVal(rot[0]);
	mag.sqr();
	tmp.setVal(rot[1]);
	tmp.sqr();
	tmp2.setVal(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub2 locRot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			locRot[i1].setVal(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.setVal(locOri[i3]);
				tmp.mult(rot[i2]);
				locRot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].setVal(1.0);
		tmp.setVal(locRot[1]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[0].sub(tmp);
		
		a1[1].setVal(0.5);
		a1[1].mult(locRot[0]);
		a1[1].mult(locRot[1]);
		a1[1].add(locRot[2]);
		
		a1[2].setVal(0.5);
		a1[2].mult(locRot[0]);
		a1[2].mult(locRot[2]);
		a1[2].sub(locRot[1]);
		
		a1[3].setVal(0.5);
		a1[3].mult(locRot[0]);
		a1[3].mult(locRot[1]);
		a1[3].sub(locRot[2]);
		
		a1[4].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[4].sub(tmp);
		
		a1[5].setVal(0.5);
		a1[5].mult(locRot[1]);
		a1[5].mult(locRot[2]);
		a1[5].add(locRot[0]);
		
		a1[6].setVal(0.5);
		a1[6].mult(locRot[0]);
		a1[6].mult(locRot[2]);
		a1[6].add(locRot[1]);
		
		a1[7].setVal(0.5);
		a1[7].mult(locRot[1]);
		a1[7].mult(locRot[2]);
		a1[7].sub(locRot[0]);
		
		a1[8].setVal(1.0);
		tmp.setVal(locRot[0]);
		tmp.sqr();
		tmp2.setVal(locRot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(oneHalf);
		a1[8].sub(tmp);
		
		matMul(instOri, a1, locOri, 3, 3, 3);
	} else {
		DiffDoub2 sth;
		DiffDoub2 cth;
		
		sth.setVal(mag);
		sth.sn();
		cth.setVal(mag);
		cth.cs();
		
		DiffDoub2 unitRot[3];
		tmp.setVal(1.0);
		tmp.dvd(mag);
		unitRot[0].setVal(rot[0]);
		unitRot[0].mult(tmp);
		unitRot[1].setVal(rot[1]);
		unitRot[1].mult(tmp);
		unitRot[2].setVal(rot[2]);
		unitRot[2].mult(tmp);
		
		a1[0].setVal(unitRot[0]);
		a1[1].setVal(unitRot[1]);
		a1[2].setVal(unitRot[2]);
		
		i1 = 0;
		if(abs(unitRot[1].val) < abs(unitRot[0].val)) {
			i1 = 1;
		}
		if(abs(unitRot[2].val) < abs(unitRot[i1].val)) {
			i1 = 2;
		}
		tmp.setVal(1.0);
		tmp2.setVal(unitRot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unitRot[i1]^2)
		a1[3+i1].setVal(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].setVal(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		crossProd(&a1[6],&a1[0],&a1[3]);
		
		a2[0].setVal(a1[0]);
		a2[1].setVal(a1[1]);
		a2[2].setVal(a1[2]);
		
		a2[3].setVal(cth);
		a2[3].mult(a1[3]);
		tmp.setVal(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].setVal(cth);
		a2[4].mult(a1[4]);
		tmp.setVal(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].setVal(cth);
		a2[5].mult(a1[5]);
		tmp.setVal(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].setVal(sth);
		a2[6].mult(a1[3]);
		tmp.setVal(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].setVal(sth);
		a2[7].mult(a1[4]);
		tmp.setVal(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].setVal(sth);
		a2[8].mult(a1[5]);
		tmp.setVal(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose(a3,a2,3,3); // a3 = a2^T
		matMul(a2,a3,a1,3,3,3); //a2 = a2^T*a1
		matMul(instOri,locOri,a2,3,3,3);
	}
	return;
}
//end dup
 
//end skip 
 
 
 
 
 
 
//dup1
void dOridThet(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[], int v1, int v2) {
	if(v1 + v2 == 0) {
		rotateOrient(instOri, locOri, rot);
		return;
	}
//preserve
    DiffDoub1 dRot[3];
    DiffDoub2 d2Rot[3];
	DiffDoub1 dLocOri[9];
	DiffDoub2 d2LocOri[9];
	DiffDoub1 dInstOri[9];
	DiffDoub2 d2InstOri[9];
//end preserve
    int i1;
	bool isDiff = locOri[0].diffType();
	double param;
	double param2;
	if(!isDiff) {
		if(v1*v2 == 0) {
			for (i1 = 0; i1 < 9; i1++) {
				param = locOri[i1].val;
			    dLocOri[i1].setVal(param);
		    }
			v1 = v1 + v2;
			for (i1 = 1; i1 < 4; i1++) {
				param = rot[i1-1].val;
				if(i1 == v1) {
					dRot[i1-1].setVal(param,1.0);
				} else {
					dRot[i1-1].setVal(param,0.0);
				}
			}
			rotateOrient(dInstOri, dLocOri, dRot);
			for (i1 = 0; i1 < 9; i1++) {
				param = dInstOri[i1].dval;
				instOri[i1].setVal(param);
			}
		} else {
			for (i1 = 0; i1 < 9; i1++) {
				param = locOri[i1].val;
				param2 = locOri[i1].dval;
			    d2LocOri[i1].setVal(param,param2);
		    }
			if(v1 == v2) {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2Rot[i1-1].setVal(param,1.0,1.0);
					} else {
						d2Rot[i1-1].setVal(param);
					}
				}	
			} else {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2Rot[i1-1].setVal(param,1.0);
					} else if(i1 == v2) {
						d2Rot[i1-1].setVal(param,0.0,1.0);
					} else {
						d2Rot[i1-1].setVal(param,0.0);
					}
				}
			}
			rotateOrient(d2InstOri, d2LocOri, d2Rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d2InstOri[i1].dv12;
				instOri[i1].setVal(param);
			}
		}
	} else {
		for (i1 = 0; i1 < 9; i1++) {
			param = locOri[i1].val;
			param2 = locOri[i1].dval;
			d2LocOri[i1].setVal(param,param2);
		}
		v1 = v1 + v2;
		for (i1 = 1; i1 < 4; i1++) {
			param = rot[i1-1].val;
			if(i1 == v1) {
				d2Rot[i1-1].setVal(param,0.0,1.0);
			} else {
				d2Rot[i1-1].setVal(param,0.0);
			}
		}
		rotateOrient(d2InstOri, d2LocOri, d2Rot);
		for (i1 = 0; i1 < 9; i1++) {
			param = d2InstOri[i1].dv2;
			param2 = d2InstOri[i1].dv12;
			instOri[i1].setVal(param,param2);
		}
	}
	return;
}

//end dup  
 
//skip 
 
//DiffDoub versions: 
//dup1
void dOridThet(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[], int v1, int v2) {
	if(v1 + v2 == 0) {
		rotateOrient(instOri, locOri, rot);
		return;
	}
    DiffDoub1 dRot[3];
    DiffDoub2 d2Rot[3];
	DiffDoub1 dLocOri[9];
	DiffDoub2 d2LocOri[9];
	DiffDoub1 dInstOri[9];
	DiffDoub2 d2InstOri[9];
    int i1;
	bool isDiff = locOri[0].diffType();
	double param;
	double param2;
	if(!isDiff) {
		if(v1*v2 == 0) {
			for (i1 = 0; i1 < 9; i1++) {
				param = locOri[i1].val;
			    dLocOri[i1].setVal(param);
		    }
			v1 = v1 + v2;
			for (i1 = 1; i1 < 4; i1++) {
				param = rot[i1-1].val;
				if(i1 == v1) {
					dRot[i1-1].setVal(param,1.0);
				} else {
					dRot[i1-1].setVal(param,0.0);
				}
			}
			rotateOrient(dInstOri, dLocOri, dRot);
			for (i1 = 0; i1 < 9; i1++) {
				param = dInstOri[i1].dval;
				instOri[i1].setVal(param);
			}
		} else {
			for (i1 = 0; i1 < 9; i1++) {
				param = locOri[i1].val;
				param2 = locOri[i1].dval;
			    d2LocOri[i1].setVal(param,param2);
		    }
			if(v1 == v2) {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2Rot[i1-1].setVal(param,1.0,1.0);
					} else {
						d2Rot[i1-1].setVal(param);
					}
				}	
			} else {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2Rot[i1-1].setVal(param,1.0);
					} else if(i1 == v2) {
						d2Rot[i1-1].setVal(param,0.0,1.0);
					} else {
						d2Rot[i1-1].setVal(param,0.0);
					}
				}
			}
			rotateOrient(d2InstOri, d2LocOri, d2Rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d2InstOri[i1].dv12;
				instOri[i1].setVal(param);
			}
		}
	} else {
		for (i1 = 0; i1 < 9; i1++) {
			param = locOri[i1].val;
			param2 = locOri[i1].dval;
			d2LocOri[i1].setVal(param,param2);
		}
		v1 = v1 + v2;
		for (i1 = 1; i1 < 4; i1++) {
			param = rot[i1-1].val;
			if(i1 == v1) {
				d2Rot[i1-1].setVal(param,0.0,1.0);
			} else {
				d2Rot[i1-1].setVal(param,0.0);
			}
		}
		rotateOrient(d2InstOri, d2LocOri, d2Rot);
		for (i1 = 0; i1 < 9; i1++) {
			param = d2InstOri[i1].dv2;
			param2 = d2InstOri[i1].dv12;
			instOri[i1].setVal(param,param2);
		}
	}
	return;
}

//end dup  
 
//end skip 
 
 
 
 
 
 
