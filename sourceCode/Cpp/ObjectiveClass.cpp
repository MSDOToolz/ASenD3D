#include <cmath>
#include "ObjectiveClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

using namespace std;

ObjectiveTerm::ObjectiveTerm(string newCat) {
	category = newCat;
	elSetPtr = NULL;
	ndSetPtr = NULL;
	qVec = NULL;
	elVolVec = NULL;
	tgtVec = NULL;
	errNormVec = NULL;
	next = NULL;
	qLen = 0;
}

void ObjectiveTerm::setOperator(string newOp) {
	optr = newOp;
	return;
}

void ObjectiveTerm::setActiveTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
	return;
}

void ObjectiveTerm::setComponent(int newComp) {
	component = newComp;
	return;
}

void ObjectiveTerm::setLayer(int newLay) {
	layer = newLay;
	return;
}

void ObjectiveTerm::setCoef(double newCoef) {
	coef = newCoef;
	return;
}

void ObjectiveTerm::setExponent(double newExp) {
	expnt = newExp;
	return;
}

void ObjectiveTerm::setElset(string newElset) {
	elSetName = newElset;
	return;
}

void ObjectiveTerm::setElsetPtr(Set* newPtr) {
	elSetPtr = newPtr;
	return;
}

void ObjectiveTerm::setNdset(string newNdset) {
	ndSetName = newNdset;
	return;
}

void ObjectiveTerm::setNdsetPtr(Set* newPtr) {
	ndSetPtr = newPtr;
	return;
}

void ObjectiveTerm::setTgtTag(string newTag) {
	tgtTag = newTag;
	return;
}

void ObjectiveTerm::setNext(ObjectiveTerm* newNext) {
	next = newNext;
	return;
}

void ObjectiveTerm::addTargetValue(double newTgt) {
	tgtVals.addEntry(newTgt);
	return;
}

string ObjectiveTerm::getElsetName() {
	return elSetName;
}

string ObjectiveTerm::getNdsetName() {
	return ndSetName;
}

ObjectiveTerm* ObjectiveTerm::getNext() {
	return next;
}

void ObjectiveTerm::allocateObj() {
	if (qLen == 0) {
		if (elSetPtr) {
			qLen = elSetPtr->getLength();
		}
		else if (ndSetPtr) {
			qLen = ndSetPtr->getLength();
		}
		if (qLen > 0) {
			qVec = new double[qLen];
			if (optr == "powerNorm") {
				tgtVec = new double[qLen];
			} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
				elVolVec = new double[qLen];
			}
		}
	}
	return;
}

void ObjectiveTerm::allocateObjGrad() {
	allocateObj();
	if (qLen > 0 && dQdU.getDim() == 0) {
		dQdU.setDim(qLen);
		dQdV.setDim(qLen);
		dQdA.setDim(qLen);
		dQdT.setDim(qLen);
		dQdTdot.setDim(qLen);
		dQdD.setDim(qLen);
		dVdD.setDim(qLen);
		if (optr == "powerNorm") {
			errNormVec = new double[qLen];
		}
	}
	return;
}

double ObjectiveTerm::getPowerNorm() {
	int i1;
	double qErr;
	double pSum = 0.0;
	for (i1 = 0; i1 < qLen; i1++) {
		qErr = qVec[i1] - tgtVec[i1];
		pSum += pow(qErr, expnt);
	}
	return coef * pSum;
}

void ObjectiveTerm::dPowerNormdU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]) {
	dQdU.vectorMultiply(dLdU, errNormVec, true);
	dQdV.vectorMultiply(dLdV, errNormVec, true);
	dQdA.vectorMultiply(dLdA, errNormVec, true);
	dQdT.vectorMultiply(dLdT, errNormVec, true);
	dQdTdot.vectorMultiply(dLdTdot, errNormVec, true);

	return;
}

void ObjectiveTerm::dPowerNormdD(double dLdD[]) {
	dQdD.vectorMultiply(dLdD, errNormVec, true);

	return;
}

double ObjectiveTerm::getVolIntegral() {
	int i1;
	double vInt = 0.0;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
	}
	vInt -= tgtVec[0];
	return coef * pow(vInt, expnt);
}

void ObjectiveTerm::dVolIntegraldU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]) {
	int i1;
	double vInt = 0.0;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
	}
	vInt -= tgtVec[0];
	vInt = coef*expnt*pow(vInt, expnt - 1.0);
	
	for (i1 = 0; i1 < qLen; i1++) {
		 elVolVec[i1]*= vInt;
	}

	dQdU.vectorMultiply(dLdU, elVolVec, true);
	dQdV.vectorMultiply(dLdV, elVolVec, true);
	dQdA.vectorMultiply(dLdA, elVolVec, true);
	dQdT.vectorMultiply(dLdT, elVolVec, true);
	dQdTdot.vectorMultiply(dLdTdot, elVolVec, true);

	vInt = 1.0 / vInt;
	for (i1 = 0; i1 < qLen; i1++) {
		elVolVec[i1] *= vInt;
	}

	return;
}

void ObjectiveTerm::dVolIntegraldD(double dLdD[]) {
	int i1;
	double vInt = 0.0;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
	}
	vInt -= tgtVec[0];
	vInt = coef * expnt * pow(vInt, expnt - 1.0);

	for (i1 = 0; i1 < qLen; i1++) {
		qVec[i1] *= vInt;
		elVolVec[i1] *= vInt;
	}

	dQdD.vectorMultiply(dLdD, elVolVec, true);
	dVdD.vectorMultiply(dLdD, qVec, true);

	vInt = 1.0 / vInt;
	for (i1 = 0; i1 < qLen; i1++) {
		qVec[i1] *= vInt;
		elVolVec[i1] *= vInt;
	}

	return;
}

double ObjectiveTerm::getVolAverage() {
	int i1;
	double vInt = 0.0;
	double totVol = 0.0;
	double volAvg;
	double vaErr;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
		totVol += elVolVec[i1];
	}
	volAvg = vInt/totVol;
	vaErr = volAvg - tgtVec[0];
	return coef * pow(vaErr, expnt);
}

void ObjectiveTerm::dVolAveragedU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]) {
	int i1;
	double vInt = 0.0;
	double totVol = 0.0;
	double volAvg;
	double vaErr;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
		totVol += elVolVec[i1];
	}
	volAvg = vInt / totVol;
	vaErr = volAvg - tgtVec[0];
	vaErr = coef * expnt * pow(vaErr, expnt - 1.0)/totVol;

	for (i1 = 0; i1 < qLen; i1++) {
		elVolVec[i1] *= vaErr;
	}

	dQdU.vectorMultiply(dLdU, elVolVec, true);
	dQdV.vectorMultiply(dLdV, elVolVec, true);
	dQdA.vectorMultiply(dLdA, elVolVec, true);
	dQdT.vectorMultiply(dLdT, elVolVec, true);
	dQdTdot.vectorMultiply(dLdTdot, elVolVec, true);

	vaErr = 1.0 / vaErr;
	for (i1 = 0; i1 < qLen; i1++) {
		elVolVec[i1] *= vaErr;
	}

	return;
}

void ObjectiveTerm::dVolAveragedD(double dLdD[]) {
	int i1;
	int col;
	double vInt = 0.0;
	double totVol = 0.0;
	double volAvg;
	double vaErr;
	for (i1 = 0; i1 < qLen; i1++) {
		vInt += qVec[i1] * elVolVec[i1];
		totVol += elVolVec[i1];
	}
	volAvg = vInt / totVol;
	vaErr = volAvg - tgtVec[0];
	vaErr = coef * expnt * pow(vaErr, expnt - 1.0) / totVol;

	for (i1 = 0; i1 < qLen; i1++) {
		elVolVec[i1] *= vaErr;
		qVec[i1] *= vaErr;
	}

	dQdD.vectorMultiply(dLdD, elVolVec, true);
	dVdD.vectorMultiply(dLdD, qVec, true);

	vaErr = 1.0 / vaErr;
	for (i1 = 0; i1 < qLen; i1++) {
		elVolVec[i1] *= vaErr;
		qVec[i1] *= vaErr;
	}

	vaErr = vaErr * volAvg;

	MatrixEnt* thisEnt;
	for (i1 = 0; i1 < qLen; i1++) {
		thisEnt = dVdD.getFirstEnt(i1);
		while (thisEnt) {
			col = thisEnt->col;
			dLdD[col] -= vaErr * thisEnt->value;
			thisEnt = thisEnt->nextEnt;
		}
	}

	return;
}

double ObjectiveTerm::getObjVal(NdPt ndAr[], ElPt elAr[], DVPt dvAr[]) {


}

Objective::Objective() {
	firstTerm = NULL;
	lastTerm = NULL;
	length = 0;
}

void Objective::addTerm(ObjectiveTerm* newTerm) {
	if(!firstTerm) {
		firstTerm = newTerm;
		lastTerm = newTerm;
	} else {
		lastTerm->setNext(newTerm);
		lastTerm = newTerm;
	}
	length++;
}

int Objective::getLength() {
	return length;
}

ObjectiveTerm* Objective::getFirst() {
	return firstTerm;
}