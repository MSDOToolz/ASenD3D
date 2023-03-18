#include "ConstraintClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>

using namespace std;


ConstraintTerm::ConstraintTerm(string newNSet, int newDof, double newCoef) {
	nodeSet = newNSet;
	dof = newDof;
	coef = newCoef;
	nextTerm = NULL;
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


Constraint::Constraint() {
	firstTerm = NULL;
	lastTerm = NULL;
	nextConst = NULL;
	return;
}

void Constraint::addTerm(string nodeSet, int dof, double coef) {
	ConstraintTerm *newTerm = new ConstraintTerm(nodeSet,dof,coef);
	if(!firstTerm) {
		firstTerm = newTerm;
		lastTerm = newTerm;
	} else {
		lastTerm->setNext(newTerm);
		lastTerm = newTerm;
	}
	return;
}

Constraint* Constraint::getNext() {
	return nextConst;
}

void Constraint::buildMat(Set *firstSet, NdPt *ndAr) {
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
			ndIndex = stoi(setName);
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
			ndIndex = stoi(setName);
			coef = thisTerm->getCoef();
			dof = thisTerm->getDof();
			col = ndAr[ndIndex].ptr->getDofIndex(dof);
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
						col = ndAr[ndIndex].ptr->getDofIndex(dof);
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

MatrixEnt* Constraint::getMatFirst(int row) {
	return mat.getFirstEnt(row);
}

ConstraintList::ConstraintList() {
	firstConst = NULL;
	lastConst = NULL;
	return;
}

Constraint* ConstraintList::getFirst() {
	return firstConst;
}