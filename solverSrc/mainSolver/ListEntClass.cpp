#include <cmath>
#include <cstddef>
#include <string>
#include <fstream>
#include <vector>
#include "ListEntClass.h"

using namespace std;

IDCapsule::IDCapsule() {
	intDat = 0;
	doubDat = 0.0;
}

MatrixEnt::MatrixEnt() {
	row = 0;
	col = 0;
	value = 0.0;
}

MatrixRow::MatrixRow() {
	rowVec.clear();
}

void MatrixRow::addEntry(int row, int col, double val) {
	auto iter = rowVec.begin();
	bool inserted = false;
	while (iter != rowVec.end() && !inserted) {
		if (iter->col == col) {
			iter->value += val;
			inserted = true;
		}
		++iter;
	}
	if (!inserted) {
		MatrixEnt newEnt;
		newEnt.row = row;
		newEnt.col = col;
		newEnt.value = val;
		rowVec.push_back(newEnt);
	}
	return;
}

SparseMat::SparseMat() {
	dim = 0;
	matrix.clear();
	return;
}

void SparseMat::setDim(int newDim) {
	int i1;
	dim = newDim;
	matrix = vector<MatrixRow>(newDim);
	return;
}

void SparseMat::zeroAll() {
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.rowVec) {
			i2.value = 0.0;
		}
	}
	
	return;
}

void SparseMat::addEntry(int row, int col, double val) {
	matrix[row].addEntry(row, col, val);
	return;
}

void SparseMat::addMatrix(SparseMat& inpMat) {
	for (auto& i1 : inpMat.matrix) {
		for (auto& i2 : i1.rowVec) {
			addEntry(i2.row, i2.col, i2.value);
		}
	}
	return;
}

void SparseMat::vectorMultiply(vector<double>& prod, vector<double>& inpVec, bool transpose) {
	if(transpose) {
		for (auto& i1 : matrix) {
			for (auto& i2 : i1.rowVec) {
				prod[i2.col] += i2.value * inpVec[i2.row];
			}
		}
	} else {
		for (auto& i1 : matrix) {
			for (auto& i2 : i1.rowVec) {
				prod[i2.row] += i2.value * inpVec[i2.col];
			}
		}
	}
	return;
}

double SparseMat::getMaxAbsVal() {
	double thisVal;
	double maxVal = 0.0;
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.rowVec) {
			thisVal = abs(i2.value);
			if (thisVal > maxVal) {
				maxVal = thisVal;
			}
		}
	}

	return maxVal;
}

void SparseMat::writeToFile(ofstream& outFile) {
	for (auto& i1 : matrix) {
		for (auto& i2 : i1.rowVec) {
			outFile << "    - [" << i2.row << ", " << i2.col << ", " << i2.value << "]\n";
		}
	}
	return;
}
