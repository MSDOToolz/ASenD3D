#include "ConstraintClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>

using namespace std;


ConstraintTerm::ConstraintTerm(string newNSet) {
	nodeSet = newNSet;
	nextTerm = NULL;
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
	firstTerm = NULL;
	lastTerm = NULL;
	nextConst = NULL;
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

void Constraint::buildMat(Set *firstSet, NdPt ndAr[]) {
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
		setNm = thisTerm->getSetName();
		try {
			ndIndex = stoi(setNm);
		} catch(...) {
			thisSet = firstSet;
			setFound = false;
			while(thisSet && !setFound) {
				if(thisSet->getName() == setNm) {
					setFound = true;
					setiLen = thisSet->getLength();
					if(setiLen > setLen) {
						setLen = setiLen;
					}
				}
				thisSet = thisSet->getNext();
			}
		}
		thisTerm = thisTerm->getNext();
	}
	mat.setDim(setLen);
	thisTerm = firstTerm;
	while(thisTerm) {
		setNm = thisTerm->getSetName();
		try {
			ndIndex = stoi(setNm);
			coef = thisTerm->getCoef();
			dof = thisTerm->getDof();
			if (type == "displacement") {
				col = ndAr[ndIndex].ptr->getDofIndex(dof-1);
			}
			else if(type == "temperature") {
				col = ndAr[ndIndex].ptr->getSortedRank();
			}
			for(row = 0; row < setLen; row++) {
				mat.addEntry(row,col,coef);
			}
		} catch(...) {
			thisSet = firstSet;
			setFound = false;
			while(thisSet && !setFound) {
				if(thisSet->getName() == setNm) {
					setFound = true;
					thisNd = thisSet->getFirstEntry();
					row = 0;
					while(thisNd) {
						ndIndex = thisNd->value;
						coef = thisTerm->getCoef();
						dof = thisTerm->getDof();
						if (type == "displacement") {
							col = ndAr[ndIndex].ptr->getDofIndex(dof - 1);
						}
						else if (type == "temperature") {
							col = ndAr[ndIndex].ptr->getSortedRank();
						}
						mat.addEntry(row,col,coef);
						thisNd = thisNd->next;
						row++;
					}
				}
				thisSet = thisSet->getNext();
			}
		}
		thisTerm = thisTerm->getNext();
	}
	return;
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
	mat.vectorMultiply(cLd, qVec, true);
	for (i1 = 0; i1 < resDim; i1++) {
		cLd[i1] *= -scaleFact;
	}

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
	firstConst = NULL;
	lastConst = NULL;
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

void ConstraintList::scaleElastic(double newSF) {
	Constraint* thisConst = firstConst;
	while (thisConst) {
		thisConst->setScaleFact(newSF);
		thisConst = thisConst->getNext();
	}
	return;
}

void ConstraintList::getTotalLoad(double cLd[], double uVec[], double qVec[], int resDim) {
	Constraint* thisConst = firstConst;
	while (thisConst) {
		thisConst->getLoad(cLd,uVec,qVec,resDim);
		thisConst = thisConst->getNext();
	}
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