#include "matrixFunctions.h"
#include <cmath>
#include "DiffDoubClass.h"
#include "ListEntClass.h"

using namespace std;

const double pi_2 = 1.57079632679489662;
const double tol = 1.0e-12;
const double magtol = 1.0e-8;
const double maxMag = 1.0e+12;
const double minMag = 1.0e-12;

void qRFactor(double& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Max;
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
	double tmp;

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
			k11 = (i2-1)*colDim + i1;
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
				mat[k11] = p1;
				mat[k12] =p2;
			}
			mat[k12] = theta;
		}
	}
	return;
}

void solveqRxEqb(double& xVec, double& mat, double& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i2Min;
	int i2Max;
	int i3Max;
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
	double tmp;
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

void symFactor(double& mat, double& qMat, matDim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
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
		for (i2 = i1+2: i2 < matDim; i2++) {
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

void getCharFun(DiffDoub& cFun, DiffDoub& mat, int matDim, DoubList& eVals; double lam, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Max;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	DiffDoub theta;
	DiffDoub sth;
	DiffDoub cth;
	DiffDoub p1;
	DiffDoub p2;
	DiffDoub tmp;

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
			k11 = (i2-1)*matDim + i1;
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
				mat[k11].setVal(p1);
				mat[k12].setVal(p2);
			}
			mat[k12].setVal(theta);
		}
	}
	
	cFun.setVal(mat[0]);
	i2 = matDim + 1;
	for (i1 = 0; i1 < matDim; i1++) {
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
	
	DoubListEnt *thisEval = evals.getFirst();
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

void getEvals(DoubList& eVals, double& mat, int matDim, double lamInit, double convTol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int matSize = matDim*matDim;
	DiffDoub *matCopy = new DiffDoub[matSize];
	DiffDoub cFun;
	double lam = lamInit;
	double dLam;
	double dTrm;
	double backStep = 1.0e+6*covTol;
	int numFound = 0;
	int maxIt = 100*matDim;
	int it = 0;
	while(numFound < matDim && it < maxIt) {
		i3 = 0;
		for (i1 = 0; i1 < matDim; i1++) {
			for (i2 = 0; i2 < matDim; i2++ {
				if(i1 == i2) {
					dTrm = mat[i3] - lam;
					matCopy[i3].setVal(dTrm,-1.0);
				} else {
					matCopy[i3].setVal(mat[i3]);
				}
				i3++;
			}
		}
		getCharFun(&cFun,matCopy,matDim,eVals,lam,triDiag);
		dLam = -cFun.val/cFun.dval;
		lam+= dLam;
		if(abs(dLam) < convTol) {
			eVals.addEntry(lam);
			numFound++;
			lam-= backStep;
		}
		it++;
	}
	delete[] matCopy;
	return;
}

double getLowEval(double& mat, int matDim) {
	double& vec = new double[matDim];
	double& prevVec new double[matDim];
	double high;
	double tmp;
	double mag;
	double dp;
	int i1;
	int i2;
	int i3;
	
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
	while(!conv) {
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
	
	conv = false;
	while(!conv) {
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

void eigenSolve(DoubList& eVals, double& eVecs, double& mat, int matDim, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i4;
	double tmp;
	int matSize = matDim*matDim;
	newTd = triDiag;
	if(triDiag == 0) {
		double& qMat = new double[matSize];
		symFactor(mat,qMat,matDim);
		newTd = 1;
	}
	
    double low = getLowEval();
	
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
	
	DoubListEnt *thisVal = eVals.getFirst();
	double& matCopy = new double[matSize];
	double& vec = new double[matDim];
	double& prevVec = new double[matDim];
	double& bVec = new double[matDim];
	double shift;
	bool conv;
	double mag;
	double dp;
	int valNum = 0;
	while(thisVal) {
		shift = thisVal->value + 1.0e-8*matMag;
		i3 = 0;
		for (i1 = 0; i1 < matDim; i1++) {
			for (i2 = 0; i2 < matDim; i2++) {
				if(i1 == i2) {
					matCopy[i3] = mat[i3] - shift;
				} else {
					matCopy[i3] = mat[i3];
				}
				i3++;
			}
		}
		qRFactor(matCopy,matDim,0,(matDim-1),0,(matDim-1),newTd);
		mag = 0.0;
		for (i1 = 0; i1 < matDim; i1++) {
			tmp = sin(1.0*i1*(valNum+1));
			prevVec[i1] = tmp;
			mag+= tmp*tmp;
		}
		mag = 1.0/sqrt(mag);
		for (i1 = 0; i1 < matDim; i1++) {
			prevVec[i1] = mag*prevVec[i1];
		}
		conv = false;
		while(!conv) {
			for (i1 = 0; i1 < matDim; i1++) {
				bVec[i1] = prevVec[i1];
			}
			solveqRxEqb(vec,matCopy,bVec,matDim,0,(matDim-1),0,(matDim-1),newTd);
			mag = 0.0;
			dp = 0.0;
			for(i1 = 0; i1 < matDim; i1++) {
				mag+= vec[i1]*vec[i1];
				dp+= vec[i1]*prevVec[i1];
			}
			mag = sqrt(mag);
			if(abs(dp) > 0.999999*mag) {
				conv = true;
				if(triDiag == 0) {
					i3 = 0;
					i4 = valNum;
					for (i1 = 0; i1 < matDim; i1++) {
						eVecs[i4] = 0.0;
						for (i2 = 0; i2 < matDim; i2++) {
						    eVecs[i4]+= qMat[i3]*prevVec[i2];
							i3++;
						}
						i4+= matDim;
					}
				} else {
					i4 = valNum;
					for (i1 = 0; i1 < matDim; i1++) {
						eVecs[i4] = prevVec[i1];
						i4+= matDim;
					}
				}
			} else {
				mag = 1.0/mag;
				for (i1 = 0; i1 < matDim; i1++) {
					prevVec[i1] = mag*vec[i1];
				}
			}
		}
		valNum++;
		thisVal = thisVal->next;
	}
	
	if(triDiag == 0) {
		delete[] qMat;
	}
	
	delete[] matCopy;
	delete[] vec;
	delete[] prevVec;
	delete[] bVec;
	return;
}

//dup1
void qRFactor(Doub& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i2Min;
	int i2Max;
	int i3Max;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	Doub theta;
	Doub sth;
	Doub cth;
	Doub p1;
	Doub p2;
	Doub tmp;

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
			k11 = (i2-1)*colDim + i1;
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
				mat[k11].setVal(p1);
				mat[k12].setVal(p2);
			}
			mat[k12].setVal(theta);
		}
	}
	return;
}

void solveqRxEqb(Doub& xVec, Doub& mat, Doub& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i2Min;
	int i2Max;
	int i3Max;
	int i3Max;
	int k11;
	int k12;
	int k22;
	int k23;
	Doub theta;
	Doub sth;
	Doub cth;
	Doub p1;
	Doub p2;
	Doub tmp;
	Doub rowSum;

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

void getDetInv(Doub& det, Doub& inv, Doub& mat, int colDim, int triDiag, Doub& xVec, Doub& bVec) {
	qRFactor(mat, colDim, 0, (colDim-1), 0, (colDim-1), triDiag);
	int i1;
    int i2;
	int i3;
	det.setVal(1.0);
	for (i1 = 0; i1 < colDim; i1++) {
		for (i2 = 0; i2 < colDim; i2++) {
			if(i1 == i2) {
				bVec[i1].setVal(1.0);
			} else {
			    bVec[i1].setVal(0.0);
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

//dup2
void matMul(Doub& prod, Doub& mat1, Doub& mat2, int m1Rows, int m1Cols, int m2Cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	Doub tmp;
	
	i4 = 0;
	for (i1 = 0; i1 < m1Rows; i1++) {
		for (i2 = 0; i2 < m2Cols; i2++) {
			instOri[i4].setVal(0.0);
			i5 = i1*m1Cols;
			i6 = i2;
			prod[i4].setVal(0.0);
			for (i3 = 0; i3 < m1Cols; i3++) {
				tmp.setVal(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2Cols;
			}
			i4++;
		}
	}
	return;
}

void transpose(Doub& matT, Doub& mat, int rowDim, int colDim) {
	Doub tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < rowDim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < colDim; i2++) {
			matT[i4].setVal(mat[i3]);
			i4+= rowDim;
		}
	}
	return;
}

void crossProd(Doub& prod, Doub& v1, Doub& v2) {
	Doub tmp;
	
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

void rotateOrient(Doub& instOri, Doub& locOri, Doub& rot) {
	Doub mag;
	Doub tmp;
	Doub tmp2;
	Doub a1[9];
	Doub a2[9];
	Doub a3[9];
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
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
		Doub locRot[3];
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
		tmp.mult(0.5);
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
		tmp.mult(0.5);
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
		tmp.mult(0.5);
		a1[8].sub(tmp);
		
		matMul(instOri, &a1, locOri, 3, 3, 3);
	} else {
		Doub sth;
		Doub cth;
		
		sth.setVal(mag);
		sth.sn();
		cth.setVal(mag);
		cth.cs();
		
		Doub unitRot[3];
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
		tmp.setVal(1.0)
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
		
		transpose(&a3,&a2,3,3); // a3 = a2^T
		matMul(&a2,&a3,&a1,3,3,3); //a2 = a2^T*a1
		matMul(instOri,locOri,&a2,3,3,3);
	}
	return;
}
//end dup

//dup1
void dOridThet(Doub& instOri, Doub& locOri, Doub& rot, int v1, int v2) {
	if(v1 + v2 == 0) {
		rotateOrient(instOri, locOri, rot);
		return;
	}
//preserve
    DiffDoub dRot[3];
    Diff2Doub d2Rot[3];
	DiffDoub dLocOri[9];
	Diff2Doub d2LocOri[9];
	DiffDoub dInstOri[9];
	Diff2Doub d2InstOri[9];
//end preserve
    int i1;
	bool isDiff = locOri[0].diffType();
	if(!isDiff) {
		if(v1*v2 == 0) {
			for (i1 = 0; i1 < 9; i1++) {
			    dLocOri[i1].setVal(locOri[i1].val);
		    }
			v1 = v1 + v2;
			for (i1 = 1; i1 < 4; i1++) {
				if(i1 == v1) {
					dRot[i1-1].setVal(rot[i1-1].val,1.0);
				} else {
					dRot[i1-1].setVal(rot[i1-1].val,0.0);
				}
			}
			rotateOrient(&dInstOri, &dLocOri, &dRot);
			for (i1 = 0; i1 < 9; i1++) {
				instOri[i1].setVal(dInstOri[i1].dval);
			}
		} else {
			for (i1 = 0; i1 < 9; i1++) {
			    d2LocOri[i1].setVal(locOri[i1].val,locOri[i1].dval);
		    }
			if(v1 == v2) {
				for (i1 = 1; i1 < 4; i1++) {
					if(i1 == v1) {
						d2Rot[i1-1].setVal(rot[i1-1].val,1.0,1.0);
					} else {
						d2Rot[i1-1].setVal(rot[i1-1].val);
					}
				}	
			} else {
				for (i1 = 1; i1 < 4; i1++) {
					if(i1 == v1) {
						d2Rot[i1-1].setVal(rot[i1-1].val,1.0);
					} else if(i1 == v2) {
						d2Rot[i1-1].setVal(rot[i1-1].val,0.0,1.0);
					} else {
						d2Rot[i1-1].setVal(rot[i1-1].val,0.0);
					}
				}
			}
			rotateOrient(&d2InstOri, &d2LocOri, &d2Rot);
			for (i1 = 0; i1 < 9; i1++) {
				instOri[i1].setVal(d2InstOri[i1].dv12);
			}
		}
	} else {
		for (i1 = 0; i1 < 9; i1++) {
			d2LocOri[i1].setVal(locOri[i1].val,locOri[i1].dval);
		}
		v1 = v1 + v2;
		for (i1 = 1; i1 < 4; i1++) {
			if(i1 == v1) {
				d2Rot[i1-1].setVal(rot[i1-1].val,0.0,1.0);
			} else {
				d2Rot[i1-1].setVal(rot[i1-1].val,0.0);
			}
		}
		rotateOrient(&d2InstOri, &d2LocOri, &d2Rot);
		for (i1 = 0; i1 < 9; i1++) {
			instOri[i1].setVal(d2InstOri[i1].dv12);
		}
	}
	return;
}

//end dup