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

void ObjectiveTerm::setValue(double newVal) {
	value = newVal;
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

double ObjectiveTerm::getValue() {
	return value;
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
				tgtVec = new double;
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

void ObjectiveTerm::getObjVal(double time, bool nLGeom, NdPt ndAr[], ElPt elAr[], DVPt dvAr[]) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int tgtLen;
	double tgtVal;
	allocateObj();
	IntListEnt* thisEnt;
	DoubListEnt* thisDb;
	int qInd;
	int label;
	double ndData[6];
	Element* thisEl;
	DoubStressPrereq stPre;
	double spt[3] = {0.0,0.0,0.0};
	Doub strain[6];
	Doub stress[6];
	double seDen;
	Doub eVol;
	Doub eDen;

	string catList = "displacement velocity acceleration temperature tdot";
	if (catList.find(category) > -1) {
		if (!ndSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid node set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the node set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = ndSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			if (category == "displacement") {
				ndAr[thisEnt->value].ptr->getDisp(ndData);
			} else if (category == "velocity") {
				ndAr[thisEnt->value].ptr->getVel(ndData);
			} else if (category == "acceleration") {
				ndAr[thisEnt->value].ptr->getAcc(ndData);
			} else if (category == "temperature") {
				ndData[0] = ndAr[thisEnt->value].ptr->getTemperature();
				component = 1;
			} else if (category == "tdot") {
				ndData[0] = ndAr[thisEnt->value].ptr->getTdot();
				component = 1;
			}
			qVec[qInd] = ndData[component - 1];
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.getLength();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == 0.0;
				}
			} else if (tgtLen == 1) {
				tgtVal = tgtVals.getFirst()->value;
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == tgtVal;
				}
			} else {
				thisDb = tgtVals.getFirst();
				qInd = 0;
				while (thisDb) {
					tgtVec[qInd] = thisDb->value;
					thisDb = thisDb->next;
					qInd++;
				}
			}
			value+= getPowerNorm();
			return;
		} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
			for (i1 = 0; i1 < qLen; i1++) {
				elVolVec[i1] = 1.0;
			}
			tgtLen = tgtVals.getLength();
			if (tgtLen == 0) {
				*tgtVec = 0.0;
			} else {
				*tgtVec = tgtVals.getFirst()->value;
			}
			if (optr == "volumeIntegral") {
				stPre.destroy();
				value+= getVolIntegral();
				return;
			} else {
				stPre.destroy();
				value+= getVolAverage();
				return;
			}
		}
	}
	catList = "stress strain strainEnergyDen";
	if (catList.find(category) > -1) {
		if (!elSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = elSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			thisEl = elAr[thisEnt->value].ptr;
			thisEl->getStressPrereq(stPre, ndAr, dvAr);
			thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
			if (category == "stress") {
				qVec[qInd] = stress[component - 1].val;
			} else if (category == "strain") {
				qVec[qInd] = strain[component - 1].val;
			} else {
				seDen = 0.0;
				for (i1 = 0; i1 < 6; i1++) {
					seDen += stress[i1].val * strain[i1].val;
				}
				seDen *= 0.5;
				qVec[qInd] = seDen;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				thisEl->getVolume(eVol, stPre, layer);
				elVolVec[qInd] = eVol.val;
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.getLength();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == 0.0;
				}
			}
			else if (tgtLen == 1) {
				tgtVal = tgtVals.getFirst()->value;
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == tgtVal;
				}
			}
			else {
				thisDb = tgtVals.getFirst();
				qInd = 0;
				while (thisDb) {
					tgtVec[qInd] = thisDb->value;
					thisDb = thisDb->next;
					qInd++;
				}
			}
			stPre.destroy();
			value+= getPowerNorm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgtVals.getLength() == 0) {
				*tgtVec = 0.0;
			} else {
				*tgtVec = tgtVals.getFirst()->value;
			}
			if (optr == "volumeIntegral") {
				stPre.destroy();
				value+= getVolIntegral();
				return;
			} else {
				stPre.destroy();
				value+= getVolAverage();
				return;
			}
		}
	}
	// rem: insert flux and grad(T) section

	catList = "mass volume";
	if (catList.find(category) > -1) {
		if (!elSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = elSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			thisEl = elAr[thisEnt->value].ptr;
			thisEl->getStressPrereq(stPre, ndAr, dvAr);
			thisEl->getVolume(eVol,stPre,layer);
			elVolVec[qInd] = eVol.val;
			if (category == "volume") {
				qVec[qInd] = 1.0;
			} else {
				thisEl->getDensity(eDen, layer, dvAr);
				qVec[qInd] = eDen.val;
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (tgtVals.getLength() == 0) {
			*tgtVec = 0.0;
		} else {
			*tgtVec = tgtVals.getFirst()->value;
		}
		stPre.destroy();
		value+= getVolIntegral();
		return;
	}

	stPre.destroy();

	return;
}

void ObjectiveTerm::getdLdU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[], double time, bool nLGeom, NdPt ndAr[], ElPt elAr[], DVPt dvAr[]) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int i2;
	int i3;
	int tgtLen;
	double tgtVal;
	allocateObj();
	allocateObjGrad();
	IntListEnt* thisEnt;
	DoubListEnt* thisDb;
	int qInd;
	int label;
	int dofInd;
	int currRank;
	double ndData[6];
	Element* thisEl;
	DoubStressPrereq stPre;
	double spt[3] = { 0.0,0.0,0.0 };
	Doub strain[6];
	Doub stress[6];
	Doub dsdU[198];
	Doub dedU[198];
	Doub dseDendU[33];
	int elNumNds;
	int elDofPerNd;
	int elNumIntDof;
	int elTotDof;
	double seDen;
	Doub eVol;
	Doub eDen;

	string catList = "displacement velocity acceleration temperature tdot";
	if (catList.find(category) > -1) {
		if (!ndSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid node set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the node set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = ndSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			label = thisEnt->value;
			dofInd = ndAr[label].ptr->getDofIndex(component - 1);
			currRank = ndAr[label].ptr->getSortedRank();
			if (category == "displacement") {
				dQdU.addEntry(qInd, dofInd, 1.0);
			}
			else if (category == "velocity") {
				dQdV.addEntry(qInd, dofInd, 1.0);
			}
			else if (category == "acceleration") {
				dQdA.addEntry(qInd, dofInd, 1.0);
			}
			else if (category == "temperature") {
				dQdT.addEntry(qInd, currRank, 1.0);
			}
			else if (category == "tdot") {
				dQdTdot.addEntry(qInd, currRank, 1.0);
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
			return;
		}
		else if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				stPre.destroy();
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				stPre.destroy();
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}
	catList = "stress strain strainEnergyDen";
	if (catList.find(category) > -1) {
		if (!elSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = elSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			thisEl = elAr[thisEnt->value].ptr;
			thisEl->getStressPrereq(stPre, ndAr, dvAr);
			thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
			thisEl->dStressStraindU(dsdU, dedU, spt, layer, nLGeom, stPre);
			elNumNds = thisEl->getNumNds();
			elDofPerNd = thisEl->getDofPerNd();
			elNumIntDof = thisEl->getNumIntDof();
			elTotDof = elNumNds * elDofPerNd + elNumIntDof;
			i1 = elTotDof * (component - 1);
			if (category == "stress") {
				thisEl->putVecToGlobMat(dQdU, &dsdU[i1], qInd, ndAr);
			}
			else if (category == "strain") {
				thisEl->putVecToGlobMat(dQdU, &dedU[i1], qInd, ndAr);
			}
			else {
				for (i2 = 0; i2 < elTotDof; i2++) {
					dseDendU[i2].setVal(0.0);
					i3 = i2;
					for (i1 = 0; i1 < 6; i1++) {
						dseDendU[i2].val += stress[i1].val * dedU[i3].val + dsdU[i3].val * strain[i1].val;
						i3 += elTotDof;
					}
					dseDendU[i2].val *= 0.5;
				}
				thisEl->putVecToGlobMat(dQdU, dseDendU, qInd, ndAr);
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (optr == "powerNorm") {
			stPre.destroy();
			dPowerNormdU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				stPre.destroy();
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				stPre.destroy();
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}
	// rem: insert flux and grad(T) section

	stPre.destroy();

	return;
}

void ObjectiveTerm::getdLdD(double dLdD[], double time, bool nLGeom, NdPt ndAr[], ElPt elAr[], DVPt dvAr[]) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int tgtLen;
	double tgtVal;
	allocateObj();
	IntListEnt* thisEnt;
	DoubListEnt* thisDb;
	int qInd;
	int label;
	double ndData[6];
	Element* thisEl;
	DoubStressPrereq stPre;
	double spt[3] = { 0.0,0.0,0.0 };
	Doub strain[6];
	Doub stress[6];
	double seDen;
	Doub eVol;
	Doub eDen;


	string catList = "stress strain strainEnergyDen";
	if (catList.find(category) > -1) {
		if (!elSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = elSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			thisEl = elAr[thisEnt->value].ptr;
			thisEl->getStressPrereq(stPre, ndAr, dvAr);
			thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
			if (category == "stress") {
				qVec[qInd] = stress[component - 1].val;
			}
			else if (category == "strain") {
				qVec[qInd] = strain[component - 1].val;
			}
			else {
				seDen = 0.0;
				for (i1 = 0; i1 < 6; i1++) {
					seDen += stress[i1].val * strain[i1].val;
				}
				seDen *= 0.5;
				qVec[qInd] = seDen;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				thisEl->getVolume(eVol, stPre, layer);
				elVolVec[qInd] = eVol.val;
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.getLength();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == 0.0;
				}
			}
			else if (tgtLen == 1) {
				tgtVal = tgtVals.getFirst()->value;
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] == tgtVal;
				}
			}
			else {
				thisDb = tgtVals.getFirst();
				qInd = 0;
				while (thisDb) {
					tgtVec[qInd] = thisDb->value;
					thisDb = thisDb->next;
					qInd++;
				}
			}
			stPre.destroy();
			value += getPowerNorm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgtVals.getLength() == 0) {
				*tgtVec = 0.0;
			}
			else {
				*tgtVec = tgtVals.getFirst()->value;
			}
			if (optr == "volumeIntegral") {
				stPre.destroy();
				value += getVolIntegral();
				return;
			}
			else {
				stPre.destroy();
				value += getVolAverage();
				return;
			}
		}
	}
	// rem: insert flux and grad(T) section

	catList = "mass volume";
	if (catList.find(category) > -1) {
		if (!elSetPtr) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		thisEnt = elSetPtr->getFirstEntry();
		qInd = 0;
		while (thisEnt) {
			thisEl = elAr[thisEnt->value].ptr;
			thisEl->getStressPrereq(stPre, ndAr, dvAr);
			thisEl->getVolume(eVol, stPre, layer);
			elVolVec[qInd] = eVol.val;
			if (category == "volume") {
				qVec[qInd] = 1.0;
			}
			else {
				thisEl->getDensity(eDen, layer, dvAr);
				qVec[qInd] = eDen.val;
			}
			thisEnt = thisEnt->next;
			qInd++;
		}
		if (tgtVals.getLength() == 0) {
			*tgtVec = 0.0;
		}
		else {
			*tgtVec = tgtVals.getFirst()->value;
		}
		stPre.destroy();
		value += getVolIntegral();
		return;
	}

	stPre.destroy();

	return;
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