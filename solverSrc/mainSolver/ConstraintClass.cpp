#include "ConstraintClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <iostream>

using namespace std;


ConstraintTerm::ConstraintTerm(string newNSet) {
	nodeSet = newNSet;
	nsPtr = nullptr;
	nextTerm = nullptr;
	return;
}

void ConstraintTerm::setNsPtr(Set* newPtr) {
	nsPtr = newPtr;
	return;
}

void ConstraintTerm::setDof(int newDof) {
	dof = newDof;
	return;
}

void ConstraintTerm::setCoef(double newCoef) {
	coef = newCoef;
	return;
}

string ConstraintTerm::getSetName() {
	return nodeSet;
}

Set* ConstraintTerm::getSetPtr() {
	return nsPtr;
}

int ConstraintTerm::getDof() {
	return dof;
}

double ConstraintTerm::getCoef() {
	return coef;
}

void ConstraintTerm::setNext(ConstraintTerm *newNext) {
	nextTerm = newNext;
	return;
}

ConstraintTerm* ConstraintTerm::getNext() {
	return nextTerm;
}


Constraint::Constraint() {
	firstTerm = nullptr;
	lastTerm = nullptr;
	nextConst = nullptr;
	rhs = 0.0;
	scaleFact = 0.0;
	return;
}

void Constraint::setType(string newType) {
	type = newType;
	return;
}

void Constraint::addTerm(ConstraintTerm *newTerm) {
	if(!firstTerm) {
		firstTerm = newTerm;
		lastTerm = newTerm;
	} else {
		lastTerm->setNext(newTerm);
		lastTerm = newTerm;
	}
	return;
}

void Constraint::setRhs(double newRhs) {
	rhs = newRhs;
	return;
}

void Constraint::setScaleFact(double newSFact) {
	scaleFact = newSFact;
	return;
}

void Constraint::setNext(Constraint *newNext) {
	nextConst = newNext;
	return;
}

Constraint* Constraint::getNext() {
	return nextConst;
}

void Constraint::buildMat(Node* ndAr[]) {
	int setLen = 1;
	int setiLen;
	int ndIndex;
	int row;
	int col;
	int dof;
	double coef;
	string setNm;
	bool setFound;
	IntListEnt *thisNd;
	Set *thisSet;

	ConstraintTerm *thisTerm = firstTerm;
	while(thisTerm) {
		thisSet = thisTerm->getSetPtr();
		setiLen = thisSet->getLength();
		if (setiLen > setLen) {
			setLen = setiLen;
		}
		thisTerm = thisTerm->getNext();
	}

	mat.setDim(setLen);
	thisTerm = firstTerm;
	while(thisTerm) {
		coef = thisTerm->getCoef();
		dof = thisTerm->getDof();
		thisSet = thisTerm->getSetPtr();
		if (thisSet->getLength() == 1) {
			thisNd = thisSet->getFirstEntry();
			ndIndex = thisNd->value;
			if (type == "displacement") {
				col = ndAr[ndIndex]->getDofIndex(dof - 1);
			}
			else {
				col = ndAr[ndIndex]->getSortedRank();
			}
			for (row = 0; row < setLen; row++) {
				mat.addEntry(row, col, coef);
			}
		}
		else {
			thisNd = thisSet->getFirstEntry();
			row = 0;
			while (thisNd) {
				ndIndex = thisNd->value;
				if (type == "displacement") {
					col = ndAr[ndIndex]->getDofIndex(dof - 1);
				}
				else {
					col = ndAr[ndIndex]->getSortedRank();
				}
				mat.addEntry(row, col, coef);
				thisNd = thisNd->next;
				row++;
			}
		}

		thisTerm = thisTerm->getNext();
	}
	return;
}

ConstraintTerm* Constraint::getFirst() {
	return firstTerm;
}

int Constraint::getMatDim() {
	return mat.getDim();
}

double Constraint::getScaleFact() {
	return scaleFact;
}

MatrixEnt* Constraint::getMatFirst(int row) {
	return mat.getFirstEnt(row);
}

void Constraint::getLoad(double cLd[], double uVec[], double qVec[], int resDim) {
	int i1;
	int dim = mat.getDim();
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

void Constraint::destroy() {
	mat.destroy();
	ConstraintTerm *thisTerm = firstTerm;
	ConstraintTerm *nextTerm;
	while(thisTerm) {
		nextTerm = thisTerm->getNext();
		delete thisTerm;
		thisTerm = nextTerm;
	}
	return;
}

ConstraintList::ConstraintList() {
	firstConst = nullptr;
	lastConst = nullptr;
	return;
}

void ConstraintList::addConstraint(Constraint *newConst) {
	if(!firstConst) {
		firstConst = newConst;
		lastConst = newConst;
	} else {
		lastConst->setNext(newConst);
		lastConst = newConst;
	}
	return;
}

Constraint* ConstraintList::getFirst() {
	return firstConst;
}

void ConstraintList::setScaleFact(double newSF) {
	Constraint* thisConst = firstConst;
	while (thisConst) {
		thisConst->setScaleFact(newSF);
		thisConst = thisConst->getNext();
	}
	return;
}

void ConstraintList::getTotalLoad(double cLd[], double uVec[], double qVec[], int resDim) {
	int i1;
	Constraint* thisConst = firstConst;
	while (thisConst) {
		thisConst->getLoad(cLd,uVec,qVec,resDim);
		thisConst = thisConst->getNext();
	}
	return;
}

void ConstraintList::writeAllToFile(string fileName) {
	ofstream outFile;
	outFile.open(fileName);
	Constraint* thisConst = firstConst;
	while (thisConst) {
		thisConst->writeToFile(outFile);
		outFile << "##------------------------------------------------------\n";
		thisConst = thisConst->getNext();
	}
	outFile.close();
	return;
}

void ConstraintList::destroy() {
	Constraint *thisConst = firstConst;
	Constraint *nextConst;
	while(thisConst) {
		nextConst = thisConst->getNext();
		thisConst->destroy();
		delete thisConst;
		thisConst = nextConst;
	}
	return;
}