#include "utilities.h"
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

const double tol = 1.0e-12;

void crossProd(double prod[], double v1[], double v2[]) {
	prod[0] = v1[1] * v2[2] - v1[2] * v2[1];
	prod[1] = v1[2] * v2[0] - v1[0] * v2[2];
	prod[2] = v1[0] * v2[1] - v1[1] * v2[0];
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
		if (triDiag == 0) {
			i2Max = endRow;
		}
		else {
			i2Max = stRow + (i1 - stCol) + 1;
			if (i2Max > endRow) {
				i2Max = endRow;
			}
		}
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k11 = (i2Min - 1) * colDim + i1;
			k12 = i2 * colDim + i1;
			if (abs(mat[k11]) < tol) {
				mat[k11] = tol;
			}
			theta = atan(mat[k12] / mat[k11]);
			sth = sin(theta);
			cth = cos(theta);
			i3Min = i1;
			if (triDiag == 2) {
				i3Max = i1 + 2;
				if (i3Max > endCol) {
					i3Max = endCol;
				}
			}
			else {
				i3Max = endCol;
			}
			for (i3 = i3Min; i3 <= i3Max; i3++) {
				k22 = (i2Min - 1) * colDim + i3;
				k23 = i2 * colDim + i3;
				p1 = cth * mat[k22] + sth * mat[k23];
				p2 = -sth * mat[k22] + cth * mat[k23];
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
		if (triDiag == 0) {
			i2Max = endRow;
		}
		else {
			i2Max = stRow + (i1 - stCol) + 1;
			if (i2Max > endRow) {
				i2Max = endRow;
			}
		}
		i3 = i2Min - 1;
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			k12 = i2 * colDim + i1;
			theta = mat[k12];
			sth = sin(theta);
			cth = cos(theta);
			p1 = cth * bVec[i3] + sth * bVec[i2];
			p2 = -sth * bVec[i3] + cth * bVec[i2];
			bVec[i3] = p1;
			bVec[i2] = p2;
		}
		xVec[i1] = 0.0;
	}

	for (i1 = endCol; i1 >= stCol; i1--) {
		i3 = stRow + (i1 - stCol);
		i2Min = i1 + 1;
		if (triDiag == 2) {
			i2Max = i1 + 2;
			if (i2Max > endCol) {
				i2Max = endCol;
			}
		}
		else {
			i2Max = endCol;
		}
		rowSum = 0.0;
		k11 = i3 * colDim + i2Min;
		for (i2 = i2Min; i2 <= i2Max; i2++) {
			rowSum += mat[k11] * xVec[i2];
			k11++;
		}
		k11 = i3 * colDim + i1;
		xVec[i1] = (bVec[i3] - rowSum) / mat[k11];
	}

	return;
}

void readInputLine(ifstream& inFile, string& fileLine, string headings[], int hdLdSpace[], string data[], int& dataLen) {
	int i1;
	int i2;
	int lnLen;
	int wrdLen;
	getline(inFile, fileLine);
	i1 = fileLine.find("#");
	if (i1 > -1) {
		fileLine = fileLine.substr(0, i1);
	}
	fileLine = fileLine + " ";
	lnLen = fileLine.length();
	i1 = fileLine.find(":");
	dataLen = 0;
	if (i1 > -1) {
		i2 = fileLine.find_first_not_of(" -\n\t");
		wrdLen = i1 - i2;
		if (headings[0] == "" || hdLdSpace[0] == i2) {
			headings[0] = fileLine.substr(i2, wrdLen);
			hdLdSpace[0] = i2;
			headings[1] = "";
			hdLdSpace[1] = 0;
			headings[2] = "";
			hdLdSpace[2] = 0;
			headings[3] = "";
			hdLdSpace[3] = 0;
		}
		else if (headings[1] == "" || hdLdSpace[1] == i2) {
			headings[1] = fileLine.substr(i2, wrdLen);
			hdLdSpace[1] = i2;
			headings[2] = "";
			hdLdSpace[2] = 0;
			headings[3] = "";
			hdLdSpace[3] = 0;
		}
		else if (headings[2] == "" || hdLdSpace[2] == i2) {
			headings[2] = fileLine.substr(i2, wrdLen);
			hdLdSpace[2] = i2;
			headings[3] = "";
			hdLdSpace[3] = 0;
		}
		else {
			headings[3] = fileLine.substr(i2, wrdLen);
			hdLdSpace[3] = i2;
		}
		i1++;
		while (i1 < lnLen) {
			fileLine = fileLine.substr(i1);
			i1 = fileLine.find_first_not_of(" ,[]\t\n");
			if (i1 > -1) {
				fileLine = fileLine.substr(i1);
				lnLen = fileLine.length();
				i1 = fileLine.find_first_of(" ,[]\t\n");
				if (i1 > -1) {
					data[dataLen] = fileLine.substr(0, i1);
					dataLen++;
				}
				else {
					i1 = lnLen;
				}
			}
			else {
				i1 = lnLen;
			}
		}
	}
	else {
		i1 = fileLine.find("- ");
		if (i1 > -1) {
			i1++;
			while (i1 < lnLen) {
				fileLine = fileLine.substr(i1);
				i1 = fileLine.find_first_not_of(" ,[]\t\n");
				if (i1 > -1) {
					fileLine = fileLine.substr(i1);
					lnLen = fileLine.length();
					i1 = fileLine.find_first_of(" ,[]\t\n");
					if (i1 > -1) {
						data[dataLen] = fileLine.substr(0, i1);
						dataLen++;
					}
					else {
						i1 = lnLen;
					}
				}
				else {
					i1 = lnLen;
				}
			}
		}
	}

	return;
}