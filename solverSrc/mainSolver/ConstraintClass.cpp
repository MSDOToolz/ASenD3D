#include "ConstraintClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

const int max_int = 2000000000;

ConstraintTerm::ConstraintTerm() {
	nodeSet = "";
	nsPtr = max_int;
	int dof = max_int;
	coef = 0.0;
	return;
}


Constraint::Constraint() {
	terms.clear();
	rhs = 0.0;
	scaleFact = 0.0;
	return;
}

void Constraint::buildMat(vector<Node>& ndAr, vector<Set>& setAr) {
	int setLen = 1;
	int setiLen;
	int ndIndex;
	int row;
	int col;
	int dof;
	double coef;
	string setNm;
	bool setFound;
	int thisSet;

	for (auto& tm : terms) {
		thisSet = tm.nsPtr;
		setiLen = setAr[thisSet].labels.size();
		if (setiLen > setLen) {
			setLen = setiLen;
		}
	}

	mat.setDim(setLen);
	for (auto& tm : terms) {
		coef = tm.coef;
		dof = tm.dof;
		thisSet = tm.nsPtr;
		list<int>& setLabs = setAr[thisSet].labels;
		if (setLabs.size() == 1) {
			ndIndex = *setLabs.begin();
			if (type == "displacement") {
				col = ndAr[ndIndex].dofIndex[dof - 1];
			}
			else {
				col = ndAr[ndIndex].sortedRank;
			}
			for (row = 0; row < setLen; row++) {
				mat.addEntry(row, col, coef);
			}
		}
		else {
			row = 0;
			for (auto& nd : setLabs) {
				if (type == "displacement") {
					col = ndAr[nd].dofIndex[dof - 1];
				}
				else {
					col = ndAr[nd].sortedRank;
				}
				mat.addEntry(row, col, coef);
				row++;
			}
		}
	}
	return;
}

void Constraint::fullVecMultiply(vector<double>& prod, vector<double>& vec, vector<double>& tmpV) {
	// Compute scaleFact*[mat]^T*([mat]*vec + qVec)
	int i1;
	for (auto& tv : tmpV) {
		tv = 0.0;
	}
	mat.vectorMultiply(tmpV, vec, false);
	for (auto& tv : tmpV) {
		tv *= scaleFact;
	}
	mat.vectorMultiply(prod, tmpV, true);
	return;
}

void Constraint::getLoad(vector<double>& cLd, vector<double>& uVec, vector<double>& qVec, int resDim) {
	int i1;
	int dim = mat.dim;
	for (i1 = 0; i1 < dim; i1++) {
		qVec[i1] = -rhs;
	}
	mat.vectorMultiply(qVec, uVec, false);
	for (i1 = 0; i1 < dim; i1++) {
		qVec[i1] *= -scaleFact;
	}
	mat.vectorMultiply(cLd, qVec, true);

	return;
}

void Constraint::writeToFile(ofstream& outFile) {
	outFile << "type: " << type << "\n";
	outFile << "scale factor: " << scaleFact << "\n";
	outFile << "rhs: " << rhs << "\n";
	outFile << "matrix: \n";
	mat.writeToFile(outFile);
	return;
}

ConstraintList::ConstraintList() {
	constVec.clear();
	return;
}

void ConstraintList::setScaleFact(double newSF) {
	for (auto& cnst : constVec) {
		cnst.scaleFact = newSF;
	}
	return;
}

void ConstraintList::getTotalVecMult(vector<double>& prod, vector<double>& vec, vector<double>& tmpV) {
	for (auto& cnst : constVec) {
		cnst.fullVecMultiply(prod, vec, tmpV);
	}
	return;
}

void ConstraintList::getTotalLoad(vector<double>& cLd, vector<double>& uVec, vector<double>& qVec, int resDim) {
	for (auto& cnst : constVec) {
		cnst.getLoad(cLd,uVec,qVec,resDim);
	}
	return;
}

void ConstraintList::writeAllToFile(string fileName) {
	ofstream outFile;
	outFile.open(fileName);
	for (auto& cnst : constVec) {
		cnst.writeToFile(outFile);
		outFile << "##------------------------------------------------------\n";
	}
	outFile.close();
	return;
}