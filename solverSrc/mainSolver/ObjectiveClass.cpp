#include <cmath>
#include "ObjectiveClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

using namespace std;

ObjectiveTerm::ObjectiveTerm() {
	category = "";
	elSetPtr = -1;
	ndSetPtr = -1;
	qVec.clear();
	elVolVec.clear();
	tgtVec.clear();
	errNormVec.clear();
	activeTime[0] = -1.0;
	activeTime[1] = 1.0e+100;
	qLen = 0;
}

void ObjectiveTerm::setActiveTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
	return;
}

void ObjectiveTerm::allocateObj(vector<Set>& ndSets, vector<Set>& elSets) {
	int i1;
	if (qLen == 0) {
		if (elSetPtr > -1) {
			qLen = elSets[elSetPtr].labels.size();
		}
		else if (ndSetPtr > -1) {
			qLen = ndSets[ndSetPtr].labels.size();
		}
		if (qLen > 0) {
			qVec = vector<double>(qLen);
			if (optr == "powerNorm") {
				tgtVec = vector<double>(qLen);
			} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
				elVolVec = vector<double>(qLen);
				tgtVec = vector<double>(1);
			}
		}
	}

	for (i1 = 0; i1 < qLen; i1++) {
		qVec[i1] = 0.0;
		if (optr == "powerNorm") {
			tgtVec[i1] = 0.0;
		}
		else {
			elVolVec[i1] = 0.0;
		}
	}
	if (optr == "volumeIntegral" || optr == "volumeAverage") {
		tgtVec[0] = 0.0;
	}

	return;
}

void ObjectiveTerm::allocateObjGrad() {
	int i1;
	if (qLen > 0 && dQdU.dim == 0) {
		dQdU.setDim(qLen);
		dQdV.setDim(qLen);
		dQdA.setDim(qLen);
		dQdT.setDim(qLen);
		dQdTdot.setDim(qLen);
		dQdD.setDim(qLen);
		dVdD.setDim(qLen);
		if (optr == "powerNorm") {
			errNormVec = vector<double>(qLen);
		}
	}
	dQdU.zeroAll();
	dQdV.zeroAll();
	dQdA.zeroAll();
	dQdT.zeroAll();
	dQdTdot.zeroAll();
	dQdD.zeroAll();
	dVdD.zeroAll();
	for (i1 = 0; i1 < qLen; i1++) {
		errNormVec[i1] = 0.0;
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

void ObjectiveTerm::dPowerNormdU(vector<double>& dLdU, vector<double>& dLdV, vector<double>& dLdA, vector<double>& dLdT, vector<double>& dLdTdot) {
	dQdU.vectorMultiply(dLdU, errNormVec, true);
	dQdV.vectorMultiply(dLdV, errNormVec, true);
	dQdA.vectorMultiply(dLdA, errNormVec, true);
	dQdT.vectorMultiply(dLdT, errNormVec, true);
	dQdTdot.vectorMultiply(dLdTdot, errNormVec, true);

	return;
}

void ObjectiveTerm::dPowerNormdD(vector<double>& dLdD) {
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

void ObjectiveTerm::dVolIntegraldU(vector<double>& dLdU, vector<double>& dLdV, vector<double>& dLdA, vector<double>& dLdT, vector<double>& dLdTdot) {
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

void ObjectiveTerm::dVolIntegraldD(vector<double>& dLdD) {
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

void ObjectiveTerm::dVolAveragedU(vector<double>& dLdU, vector<double>& dLdV, vector<double>& dLdA, vector<double>& dLdT, vector<double>& dLdTdot) {
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

void ObjectiveTerm::dVolAveragedD(vector<double>& dLdD) {
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
	vaErr = 1.0 / vaErr;

	vaErr = vaErr * volAvg;

	for (i1 = 0; i1 < qLen; i1++) {
		MatrixRow& thisRow = dVdD.matrix[i1];
		for (auto& me : thisRow.rowVec) {
			col = me.col;
			dLdD[col] -= vaErr * me.value;
		}
	}

	return;
}

void ObjectiveTerm::getObjVal(double time, bool nLGeom, vector<Node>& ndAr, vector<Element>& elAr, vector<Set>& ndSets, vector<Set>& elSets, vector<Section>& secAr, vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int fi;
	int numLay;
	int tgtLen;
	double tgtVal;
	allocateObj(ndSets,elSets);
	int qInd;
	double ndData[6];
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	double seDen;
	DiffDoub0 def[9];
	DiffDoub0 frcMom[9];
	DiffDoub0 flux[3];
	DiffDoub0 tGrad[3];
	DiffDoub0 eVol;
	DiffDoub0 eDen;
	DiffDoub0 tmp;

	string catList = "displacement velocity acceleration temperature tdot";
	fi = catList.find(category);
	if (fi > -1) {
		if (ndSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid node set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the node set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& ndi : ndSets[ndSetPtr].labels) {
			if (category == "displacement") {
				qVec[qInd] = ndAr[ndi].displacement[component - 1];
			} else if (category == "velocity") {
				qVec[qInd] = ndAr[ndi].velocity[component - 1];
			} else if (category == "acceleration") {
				qVec[qInd] = ndAr[ndi].acceleration[component - 1];
			} else if (category == "temperature") {
				qVec[qInd] = ndAr[ndi].temperature;
			} else if (category == "tdot") {
				qVec[qInd] = ndAr[ndi].tempChangeRate;
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.size();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = 0.0;
				}
			} else if (tgtLen == 1) {
				tgtVal = *tgtVals.begin();
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = tgtVal;
				}
			} else {
				qInd = 0;
				for (auto& tv : tgtVals) {
					tgtVec[qInd] = tv;
					qInd++;
				}
			}
			value+= getPowerNorm();
			return;
		} else if (optr == "volumeIntegral" || optr == "volumeAverage") {
			for (i1 = 0; i1 < qLen; i1++) {
				elVolVec[i1] = 1.0;
			}
			tgtLen = tgtVals.size();
			if (tgtLen == 0) {
				tgtVec[0] = 0.0;
			} else {
				tgtVec[0] = *tgtVals.begin();
			}
			if (optr == "volumeIntegral") {
				value+= getVolIntegral();
				return;
			} else {
				value+= getVolAverage();
				return;
			}
		}
	}
	catList = "stress strain strainEnergyDen";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			thisEl.getStressStrain(stress, strain, thisEl.sCent, layer, nLGeom, stPre);
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
				thisEl.getVolume(eVol, stPre, layer, secAr, dvAr);
				elVolVec[qInd] = eVol.val;
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.size();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = 0.0;
				}
			}
			else if (tgtLen == 1) {
				tgtVal = *tgtVals.begin();
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = tgtVal;
				}
			}
			else {
				qInd = 0;
				for (auto& tv : tgtVals) {
					tgtVec[qInd] = tv;
					qInd++;
				}
			}
			value+= getPowerNorm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgtVals.size() == 0) {
				tgtVec[0] = 0.0;
			} else {
				tgtVec[0] = *tgtVals.begin();
			}
			if (optr == "volumeIntegral") {
				value+= getVolIntegral();
				return;
			} else {
				value+= getVolAverage();
				return;
			}
		}
	}

	catList = "sectionDef sectionFrcMom";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			//thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
			thisEl.getDefFrcMom(def, frcMom, thisEl.sCent, nLGeom, stPre);
			if (category == "sectionDef") {
				qVec[qInd] = def[component - 1].val;
			}
			else if (category == "sectionFrcMom") {
				qVec[qInd] = frcMom[component - 1].val;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				numLay = secAr[thisEl.sectPtr].layers.size();
				if (numLay > 1) {
					eVol.setVal(0.0);
					for (i1 = 0; i1 < numLay; i1++) {
						thisEl.getVolume(tmp, stPre, i1, secAr, dvAr);
						eVol.add(tmp);
					}
				}
				else {
					thisEl.getVolume(eVol, stPre, 0, secAr, dvAr);
				}
				elVolVec[qInd] = eVol.val;
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.size();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = 0.0;
				}
			}
			else if (tgtLen == 1) {
				tgtVal = *tgtVals.begin();
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = tgtVal;
				}
			}
			else {
				qInd = 0;
				for (auto& tv : tgtVals) {
					tgtVec[qInd] = tv;
					qInd++;
				}
			}
			value += getPowerNorm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgtVals.size() == 0) {
				tgtVec[0] = 0.0;
			}
			else {
				tgtVec[0] = *tgtVals.begin();
			}
			if (optr == "volumeIntegral") {
				value += getVolIntegral();
				return;
			}
			else {
				value += getVolAverage();
				return;
			}
		}
	}
	
	catList = "flux tempGradient";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			thisEl.getFluxTGrad(flux, tGrad, thisEl.sCent, layer, stPre);
			if (category == "flux") {
				qVec[qInd] = flux[component - 1].val;
			}
			else if (category == "tempGradient") {
				qVec[qInd] = tGrad[component - 1].val;
			}
			if (optr == "volumeIntegral" || optr == "volumeAverage") {
				thisEl.getVolume(eVol, stPre, layer, secAr, dvAr);
				elVolVec[qInd] = eVol.val;
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			tgtLen = tgtVals.size();
			if (tgtLen == 0) {
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = 0.0;
				}
			}
			else if (tgtLen == 1) {
				tgtVal = *tgtVals.begin();
				for (i1 = 0; i1 < qLen; i1++) {
					tgtVec[i1] = tgtVal;
				}
			}
			else {
				qInd = 0;
				for (auto& tv : tgtVals) {
					tgtVec[qInd] = tv;
					qInd++;
				}
			}
			value += getPowerNorm();
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (tgtVals.size() == 0) {
				tgtVec[0] = 0.0;
			}
			else {
				tgtVec[0] = *tgtVals.begin();
			}
			if (optr == "volumeIntegral") {
				value += getVolIntegral();
				return;
			}
			else {
				value += getVolAverage();
				return;
			}
		}
	}

	catList = "mass volume";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			thisEl.getVolume(eVol,stPre,layer,secAr,dvAr);
			elVolVec[qInd] = eVol.val;
			if (category == "volume") {
				qVec[qInd] = 1.0;
			} else {
				thisEl.getDensity(eDen, layer, secAr, matAr, dvAr);
				qVec[qInd] = eDen.val;
			}
			qInd++;
		}
		if (tgtVals.size() == 0) {
			tgtVec[0] = 0.0;
		} else {
			tgtVec[0] = *tgtVals.begin();
		}
		value+= getVolIntegral();
		return;
	}

	return;
}

void ObjectiveTerm::getdLdU(vector<double>& dLdU, vector<double>& dLdV, vector<double>& dLdA, vector<double>& dLdT, vector<double>& dLdTdot, double time, bool nLGeom, vector<Node>& ndAr, vector<Element>& elAr, vector<Set>& ndSets, vector<Set>& elSets, vector<Section>& secAr, vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int i2;
	int i3;
	int fi;
	allocateObjGrad();
	int qInd;
	int dofInd;
	int currRank;
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	vector<DiffDoub0> dsdU(288);
	vector<DiffDoub0> dedU(288);
	vector<DiffDoub0> dsdT(90);
	DiffDoub0 dseDendU[33];
	DiffDoub0 dseDendT[10];
	DiffDoub0 def[9];
	DiffDoub0 frcMom[9];
	DiffDoub0 flux[3];
	DiffDoub0 tGrad[3];
	vector<DiffDoub0> dFdT(30);
	vector<DiffDoub0> dTGdT(30);
	int elNumNds;
	int elDofPerNd;
	int elNumIntDof;
	int elTotDof;
	DiffDoub0 eVol;
	DiffDoub0 eDen;

	string catList = "displacement velocity acceleration temperature tdot";
	fi = catList.find(category);
	if (fi > -1) {
		if (ndSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid node set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the node set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& ndi : ndSets[ndSetPtr].labels) {
			dofInd = ndAr[ndi].dofIndex[component - 1];
			currRank = ndAr[ndi].sortedRank;
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
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}

	catList = "stress strain strainEnergyDen";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			thisEl.getStressStrain(stress, strain, thisEl.sCent, layer, nLGeom, stPre);
			thisEl.dStressStraindU(dsdU, dedU, dsdT, thisEl.sCent, layer, nLGeom, stPre);
			elNumNds = thisEl.numNds;
			elDofPerNd = thisEl.dofPerNd;
			elNumIntDof = thisEl.numIntDof;
			elTotDof = elNumNds * elDofPerNd + elNumIntDof;
			i1 = elTotDof * (component - 1);
			i2 = elNumNds * (component - 1);
			if (category == "stress") {
				subVec(stPre.scrVec1, dsdU, i1, i1 + elTotDof);
				thisEl.putVecToGlobMat(dQdU, stPre.scrVec1, false, qInd, ndAr);
				subVec(stPre.scrVec1, dsdT, i2, i2 + elNumNds);
				thisEl.putVecToGlobMat(dQdT, stPre.scrVec1, true, qInd, ndAr);
			}
			else if (category == "strain") {
				subVec(stPre.scrVec1, dedU, i1, i1 + elTotDof);
				thisEl.putVecToGlobMat(dQdU, stPre.scrVec1, false, qInd, ndAr);
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
				for (i2 = 0; i2 < elNumNds; i2++) {
					dseDendT[i2].setVal(0.0);
					i3 = i2;
					for (i1 = 0; i1 < 6; i1++) {
						dseDendT[i2].val += dsdT[i3].val * strain[i1].val;
						i3 += elNumNds;
					}
					dseDendT[i2].val *= 0.5;
				}
				arToVec(dseDendU, stPre.scrVec1, 0, 33);
				thisEl.putVecToGlobMat(dQdU, stPre.scrVec1, false, qInd, ndAr);
				arToVec(dseDendT, stPre.scrVec1, 0, 10);
				thisEl.putVecToGlobMat(dQdT, stPre.scrVec1, true, qInd, ndAr);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}
	
	catList = "sectionDef sectionFrcMom";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			//thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
			//thisEl->dStressStraindU(dsdU, dedU, dsdT, spt, layer, nLGeom, stPre);
			thisEl.getDefFrcMom(def, frcMom, thisEl.sCent, nLGeom, stPre);
			thisEl.dDefFrcMomdU(dedU, dsdU, dsdT, thisEl.sCent, nLGeom, stPre);
			elNumNds = thisEl.numNds;
			elDofPerNd = thisEl.dofPerNd;
			elNumIntDof = thisEl.numIntDof;
			elTotDof = elNumNds * elDofPerNd + elNumIntDof;
			i1 = elTotDof * (component - 1);
			i2 = elNumNds * (component - 1);
			if (category == "sectionFrcMom") {
				subVec(stPre.scrVec1, dsdU, i1, i1 + elTotDof);
				thisEl.putVecToGlobMat(dQdU, stPre.scrVec1, false, qInd, ndAr);
				subVec(stPre.scrVec1, dsdT, i2, i2 + elNumNds);
				thisEl.putVecToGlobMat(dQdT, stPre.scrVec1, true, qInd, ndAr);
			}
			else if (category == "sectionDef") {
				subVec(stPre.scrVec1, dedU, i1, i1 + elTotDof);
				thisEl.putVecToGlobMat(dQdU, stPre.scrVec1, false, qInd, ndAr);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}

	catList = "flux tempGradient";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
			thisEl.getFluxTGrad(flux, tGrad, thisEl.sCent, layer, stPre);
			thisEl.dFluxTGraddT(dFdT, dTGdT, thisEl.sCent, layer, stPre);
			elNumNds = thisEl.numNds;
			i1 = elNumNds * (component - 1);
			if (category == "flux") {
				subVec(stPre.scrVec1, dFdT, i1, i1 + elNumNds);
				thisEl.putVecToGlobMat(dQdT, stPre.scrVec1, true, qInd, ndAr);
			}
			else if (category == "tempGradient") {
				subVec(stPre.scrVec1, dTGdT, i1, i1 + elNumNds);
				thisEl.putVecToGlobMat(dQdT, stPre.scrVec1, true, qInd, ndAr);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
			else {
				dVolAveragedU(dLdU, dLdV, dLdA, dLdT, dLdTdot);
				return;
			}
		}
	}

	return;
}

void ObjectiveTerm::getdLdD(vector<double>& dLdD, double time, bool nLGeom, vector<Node>& ndAr,  vector<Element>& elAr, vector<Set>& ndSets, vector<Set>& elSets, vector<Section>& secAr, vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub1StressPrereq& stPre) {
	if (time < activeTime[0] || time > activeTime[1]) {
		return;
	}

	int i1;
	int numLay;
	int fi;
	int qInd;
	DiffDoub0 dvVal;
	DiffDoub1 strain[6];
	DiffDoub1 stress[6];
	double seDen;
	DiffDoub1 def[9];
	DiffDoub1 frcMom[9];
	DiffDoub1 flux[3];
	DiffDoub1 tGrad[3];
	DiffDoub1 eVol;
	DiffDoub1 eDen;
	DiffDoub1 tmp;

	string catList = "stress strain strainEnergyDen";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			for (auto& dvi : thisEl.compDVars) {
				DesignVariable& thisDV = dvAr[dvi];
				thisDV.getValue(dvVal);
				thisDV.diffVal.setVal(dvVal.val, 1.0);
				thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
				thisEl.getStressStrain(stress, strain, thisEl.sCent, layer, nLGeom, stPre);
				if (category == "stress") {
					dQdD.addEntry(qInd, dvi, stress[component - 1].dval);
				}
				else if (category == "strain") {
					dQdD.addEntry(qInd, dvi, strain[component - 1].dval);
				}
				else {
					seDen = 0.0;
					for (i1 = 0; i1 < 6; i1++) {
						seDen += stress[i1].val * strain[i1].dval + stress[i1].dval * strain[i1].val;
					}
					seDen *= 0.5;
					dQdD.addEntry(qInd, dvi, seDen);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					thisEl.getVolume(eVol, stPre, layer, secAr, dvAr);
					dVdD.addEntry(qInd, dvi, eVol.dval);
				}
				thisDV.diffVal.setVal(dvVal.val, 0.0);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdD(dLdD);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldD(dLdD);
				return;
			}
			else {
				dVolAveragedD(dLdD);
				return;
			}
		}
	}

	catList = "sectionDef sectionFrcMom";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			for (auto& dvi : thisEl.compDVars) {
				DesignVariable& thisDV = dvAr[dvi];
				thisDV.getValue(dvVal);
				thisDV.diffVal.setVal(dvVal.val, 1.0);
				thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
				//thisEl->getStressStrain(stress, strain, spt, layer, nLGeom, stPre);
				thisEl.getDefFrcMom(def, frcMom, thisEl.sCent, nLGeom, stPre);
				if (category == "sectionFrcMom") {
					dQdD.addEntry(qInd, dvi, frcMom[component - 1].dval);
				}
				else if (category == "sectionDef") {
					dQdD.addEntry(qInd, dvi, def[component - 1].dval);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					numLay = secAr[thisEl.sectPtr].layers.size();
					if (numLay > 1) {
						eVol.setVal(0.0);
						for (i1 = 0; i1 < numLay; i1++) {
							thisEl.getVolume(tmp, stPre, i1, secAr, dvAr);
							eVol.add(tmp);
						}
					}
					else {
						thisEl.getVolume(eVol, stPre, 0, secAr, dvAr);
					}
					dVdD.addEntry(qInd, dvi, eVol.dval);
				}
				thisDV.diffVal.setVal(dvVal.val, 0.0);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdD(dLdD);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldD(dLdD);
				return;
			}
			else {
				dVolAveragedD(dLdD);
				return;
			}
		}
	}
	
	catList = "flux tempGradient";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			for (auto& dvi : thisEl.compDVars) {
				DesignVariable& thisDV = dvAr[dvi];
				thisDV.getValue(dvVal);
				thisDV.diffVal.setVal(dvVal.val, 1.0);
				thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
				thisEl.getFluxTGrad(flux, tGrad, thisEl.sCent, layer, stPre);
				if (category == "flux") {
					dQdD.addEntry(qInd, dvi, flux[component - 1].dval);
				}
				else if (category == "tempGradient") {
					dQdD.addEntry(qInd, dvi, tGrad[component - 1].dval);
				}
				if (optr == "volumeIntegral" || optr == "volumeAverage") {
					thisEl.getVolume(eVol, stPre, layer, secAr, dvAr);
					dVdD.addEntry(qInd, dvi, eVol.dval);
				}
				thisDV.diffVal.setVal(dvVal.val, 0.0);
			}
			qInd++;
		}
		if (optr == "powerNorm") {
			for (i1 = 0; i1 < qLen; i1++) {
				errNormVec[i1] = coef * expnt * pow((qVec[i1] - tgtVec[i1]), (expnt - 1.0));
			}
			dPowerNormdD(dLdD);
			return;
		}
		if (optr == "volumeIntegral" || optr == "volumeAverage") {
			if (optr == "volumeIntegral") {
				dVolIntegraldD(dLdD);
				return;
			}
			else {
				dVolAveragedD(dLdD);
				return;
			}
		}
	}

	catList = "mass volume";
	fi = catList.find(category);
	if (fi > -1) {
		if (elSetPtr == -1) {
			string errStr = "Error: objective terms of category '" + category + "' must have a valid element set specified.\n";
			errStr = errStr + "Check the objective input file to make sure the element set name is correct and defined in the model input file.";
			throw runtime_error(errStr);
		}
		qInd = 0;
		for (auto& eli : elSets[elSetPtr].labels) {
			Element& thisEl = elAr[eli];
			for (auto& dvi : thisEl.compDVars) {
				DesignVariable& thisDV = dvAr[dvi];
				thisDV.getValue(dvVal);
				thisDV.diffVal.setVal(dvVal.val, 1.0);
				thisEl.getStressPrereq(stPre, secAr, matAr, ndAr, dvAr);
				thisEl.getVolume(eVol, stPre, layer, secAr, dvAr);
				dVdD.addEntry(qInd, dvi, eVol.dval);
				if(category == "mass") {
					thisEl.getDensity(eDen, layer, secAr, matAr, dvAr);
					dQdD.addEntry(qInd, dvi, eDen.dval);
				}
				thisDV.diffVal.setVal(dvVal.val, 0.0);
			}
			qInd++;
		}
		dVolIntegraldD(dLdD);
		return;
	}

	return;
}

Objective::Objective() {
	terms.clear();
}

void Objective::clearValues() {
	for (auto& tm : terms) {
		tm.value = 0.0;
	}
	return;
}

void Objective::calculateTerms(double time, bool nLGeom, vector<Node>& ndAr,  vector<Element>& elAr, vector<Set>& ndSets, vector<Set>& elSets, vector<Section>& secAr, vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre) {
	for (auto& tm : terms) {
		tm.getObjVal(time, nLGeom, ndAr, elAr, ndSets, elSets, secAr, matAr, dvAr, stPre);
	}
	return;
}

void Objective::calculatedLdU(vector<double>& dLdU, vector<double>& dLdV, vector<double>& dLdA, vector<double>& dLdT, vector<double>& dLdTdot, double time, bool nLGeom, vector<Node>& ndAr, vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre) {
	for (auto& tm : terms) {
		tm.getdLdU(dLdU, dLdV, dLdA, dLdT, dLdTdot, time, nLGeom, ndAr, elAr, ndSets, elSets, secAr, matAr, dvAr, stPre);
	}
	return;
}

void Objective::calculatedLdD(vector<double>& dLdD, double time, bool nLGeom, vector<Node>& ndAr,  vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, vector<DesignVariable>& dvAr, DiffDoub1StressPrereq& stPre) {
	for (auto& tm : terms) {
		tm.getdLdD(dLdD, time, nLGeom, ndAr, elAr, ndSets, elSets, secAr, matAr, dvAr, stPre);
	}
	return;
}