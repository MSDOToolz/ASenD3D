#include "NodeClass.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

Node::Node() {
	label = 0;
	fluid = false;
	numDof = 3;
	sortedRank = 0;
	coord[0] = 0.0;
	coord[1] = 0.0;
	coord[2] = 0.0;
	int i1;
	for (i1=0; i1 < 6; i1++) {
		displacement[i1] = 0.0;
		velocity[i1] = 0.0;
		acceleration[i1] = 0.0;
		prevDisp[i1] = 0.0;
		prevVel[i1] = 0.0;
		prevAcc[i1] = 0.0;
		initialDisp[i1] = 0.0;
		initialVel[i1] = 0.0;
		initialAcc[i1] = 0.0;
		dofIndex[i1] = 0;
	}
	for (i1 = 0; i1 < 3; i1++) {
		flVel[i1] = 0.0;
		flVelDot[i1] = 0.0;
		prevFlVel[i1] = 0.0;
		pFlVelLF[i1] = 0.0;
		prevFlVelDot[i1] = 0.0;
		initialFlVel[i1] = 0.0;
		initialFlVelDot[i1] = 0.0;
	}
	temperature = 0.0;
	tempChangeRate = 0.0;
	flDen = 0.0;
	flDenDot = 0.0;
	turbE = 0.0;
	turbEDot = 0.0;
	prevTemp = 0.0;
	pTempLF = 0.0;
	prevTdot = 0.0;
	prevFlDen = 0.0;
	pFlDenLF = 0.0;
	prevFlDenDot = 0.0;
	prevTurbE = 0.0;
	prevTurbEDot = 0.0;
	pTurbELF = 0.0;
	initialTemp = 0.0;
	initialTdot = 0.0;
	initialFlDen = 0.0;
	initialFlDenDot = 0.0;
	initialTurbE = 0.0;
	initialTurbEDot = 0.0;
	dVarLst.clear();
}

void Node::setCrd(double newCrd[]) {
	coord[0] = newCrd[0];
	coord[1] = newCrd[1];
	coord[2] = newCrd[2];
	return;
}

void Node::setDisplacement(double newDisp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] = newDisp[i1];
	}
	return;
}

void Node::setVelocity(double newVel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		velocity[i1] = newVel[i1];
	}
	return;
}

void Node::setAcceleration(double newAcc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		acceleration[i1] = newAcc[i1];
	}
	return;
}

void Node::setFlVel(double newVel[]) {
	flVel[0] = newVel[0];
	flVel[1] = newVel[1];
	flVel[2] = newVel[2];
	return;
}

void Node::setFlVelDot(double newFlVDot[]) {
	flVelDot[0] = newFlVDot[0];
	flVelDot[1] = newFlVDot[1];
	flVelDot[2] = newFlVDot[2];
	return;
}

void Node::addToDisplacement(double delDisp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] += delDisp[i1];
	}	
	return;
}

void Node::addToFlVel(double delFlVel[]) {
	flVel[0] += delFlVel[0];
	flVel[1] += delFlVel[1];
	flVel[2] += delFlVel[2];
	return;
}

void Node::setInitialDisp(double newDisp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initialDisp[i1] = newDisp[i1];
	}
	return;
}

void Node::setInitialVel(double newVel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initialVel[i1] = newVel[i1];
	}
	return;
}

void Node::setInitialAcc(double newAcc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initialAcc[i1] = newAcc[i1];
	}
	return;
}

void Node::setInitialFlVel(double newFlVel[]) {
	initialFlVel[0] = newFlVel[0];
	initialFlVel[1] = newFlVel[1];
	initialFlVel[2] = newFlVel[2];
	return;
}

void Node::setInitialFlVDot(double newFlVDot[]) {
	initialFlVelDot[0] = newFlVDot[0];
	initialFlVelDot[1] = newFlVDot[1];
	initialFlVelDot[2] = newFlVDot[2];
	return;
}

void Node::setPrevDisp(double newDisp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prevDisp[i1] = newDisp[i1];
	}
	return;
}

void Node::setPrevVel(double newVel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prevVel[i1] = newVel[i1];
	}
	return;
}

void Node::setPrevAcc(double newAcc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prevAcc[i1] = newAcc[i1];
	}
	return;
}

void Node::setPrevFlVel(double newFlVel[]) {
	prevFlVel[0] = newFlVel[0];
	prevFlVel[1] = newFlVel[1];
	prevFlVel[2] = newFlVel[2];
	return;
}

void Node::setPFlVelLF(double newVel[]) {
	pFlVelLF[0] = newVel[0];
	pFlVelLF[1] = newVel[1];
	pFlVelLF[2] = newVel[2];
	return;
}

void Node::setPrevFlVDot(double newFlVDot[]) {
	prevFlVelDot[0] = newFlVDot[0];
	prevFlVelDot[1] = newFlVDot[1];
	prevFlVelDot[2] = newFlVDot[2];
	return;
}

void Node::initializeDisp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prevDisp[i1] = initialDisp[i1];
		prevVel[i1] = initialVel[i1];
		prevAcc[i1] = initialAcc[i1];
		displacement[i1] = initialDisp[i1];
	}
	return;
}

void Node::initializeFlVel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		prevFlVel[i1] = initialFlVel[i1];
		prevFlVelDot[i1] = initialFlVelDot[i1];
		flVel[i1] = initialFlVel[i1];
	}
	return;
}

void Node::initializeTemp() {
	prevTemp = initialTemp;
	prevTdot = initialTdot;
	temperature = initialTemp;
	return;
}

void Node::initializeFlDen() {
	prevFlDen = initialFlDen;
	prevFlDenDot = initialFlDenDot;
	return;
}

void Node::initializeTurbE() {
	prevTurbE = initialTurbE;
	prevTurbEDot = initialTurbEDot;
	return;
}

void Node::updateVelAcc(double nmBeta, double nmGamma, double delT) {
	double c1;
	double c2;
	int i1;
	c1 = 1.0 / (delT * delT * (nmBeta - nmGamma));
	c2 = delT * delT * (0.5 + nmBeta - nmGamma);
	for (i1 = 0; i1 < 6; i1++) {
		acceleration[i1] = c1 * (prevDisp[i1] - displacement[i1] + delT * prevVel[i1] + c2 * prevAcc[i1]);
		velocity[i1] = prevVel[i1] + delT * ((1.0 - nmGamma) * prevAcc[i1] + nmGamma * acceleration[i1]);
	}
	
	return;
}

void Node::updateFlVelDot(double nmGamma, double delT) {
	int i1;
	double c1;
	double c2;
	c1 = 1.0 / nmGamma;
	c2 = 1.0 / delT;
	for (i1 = 0; i1 < 3; i1++) {
		flVelDot[i1] = c1 * (c2 * (flVel[i1] - pFlVelLF[i1]) - (1.0 - nmGamma) * prevFlVelDot[i1]);
	}
	return;
}

void Node::updateTdot(double nmGamma, double delT) {
	double c1;
	double c2;
	c1 = 1.0 / nmGamma;
	c2 = 1.0 / delT;
	if (fluid) {
		tempChangeRate = c1 * (c2 * (temperature - pTempLF) - (1.0 - nmGamma) * prevTdot);
	}
	else {
		tempChangeRate = c1 * (c2 * (temperature - prevTemp) - (1.0 - nmGamma) * prevTdot);
	}
	return;
}

void Node::updateFlDenDot(double nmGamma, double delT) {
	double c1;
	double c2;
	c1 = 1.0 / nmGamma;
	c2 = 1.0 / delT;
	if (fluid) {
		flDenDot = c1 * (c2 * (flDen - pFlDenLF) - (1.0 - nmGamma) * prevFlDenDot);
	}
	else {
	    flDenDot = c1 * (c2 * (flDen - prevFlDen) - (1.0 - nmGamma) * prevFlDenDot);
	}
	return;
}

void Node::updateTurbEDot(double nmGamma, double delT) {
	double c1;
	double c2;
	c1 = 1.0 / nmGamma;
	c2 = 1.0 / delT;
	if (fluid) {
		turbEDot = c1 * (c2 * (turbE - pTurbELF) - (1.0 - nmGamma) * prevTurbEDot);
	}
	return;
}

void Node::advanceDisp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prevDisp[i1] = displacement[i1];
		prevVel[i1] = velocity[i1];
		prevAcc[i1] = acceleration[i1];
	}
	return;	
}

void Node::advanceFlVel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		prevFlVel[i1] = flVel[i1];
		prevFlVelDot[i1] = flVelDot[i1];
	}
	return;
}

void Node::advanceTemp() {
	prevTemp = temperature;
	prevTdot = tempChangeRate;
	return;
}

void Node::advanceFlDen() {
	prevFlDen = flDen;
	prevFlDenDot = flDenDot;
	return;
}

void Node::advanceTurbE() {
	prevTurbE = turbE;
	prevTurbEDot = turbEDot;
	return;
}

void Node::backstepDisp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] = prevDisp[i1];
		velocity[i1] = prevVel[i1];
		acceleration[i1] = prevAcc[i1];
	}
	return;
}

void Node::backstepFlVel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		flVel[i1] = prevFlVel[i1];
		flVelDot[i1] = prevFlVelDot[i1];
	}
	return;
}

void Node::backstepTemp() {
	temperature = prevTemp;
	tempChangeRate = prevTdot;
	return;
}

void Node::backstepFlDen() {
	flDen = prevFlDen;
	flDenDot = prevFlDenDot;
	return;
}

void Node::backstepTurbE() {
	turbE = prevTurbE;
	turbEDot = prevTurbEDot;
	return;
}

void Node::addDesignVariable(int dIndex, double coef) {
	IDCapsule dv;
	dv.intDat = dIndex;
	dv.doubDat = coef;
	dVarLst.push_back(dv);
	return;
}

//dup1
void Node::getCrd(DiffDoub0 crdOut[], vector<DesignVariable>& dvAr) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	int dIndex;
	DiffDoub0 dVal;
	int comp;
	string cat;
	DiffDoub0 coef;
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		comp = thisDV.component - 1;
		if(cat == "nodeCoord") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			crdOut[comp].add(coef);
		}
	}
	return;
}

void Node::getDisp(DiffDoub0 disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

void Node::getElasticDVLoad(DiffDoub0 ld[], vector<DesignVariable>& dvAr) {
	int i1;
	int dIndex;
	DiffDoub0 dVal;
	DiffDoub0 coef;
	string cat;
	int comp;
	
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		comp = thisDV.component - 1;
		if(cat == "elasticLoad") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			ld[comp].add(coef);
		}
	}
	return;
}

void Node::getThermalDVLoad(DiffDoub0& ld, vector<DesignVariable>& dvAr) {
	int dIndex;
	DiffDoub0 dVal;
	DiffDoub0 coef;
	string cat;
	
	ld.setVal(0.0);
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		if (cat == "thermalLoad") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			ld.add(coef);
		}
	}
	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void Node::getCrd(DiffDoub1 crdOut[], vector<DesignVariable>& dvAr) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	int dIndex;
	DiffDoub1 dVal;
	int comp;
	string cat;
	DiffDoub1 coef;
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		comp = thisDV.component - 1;
		if(cat == "nodeCoord") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			crdOut[comp].add(coef);
		}
	}
	return;
}

void Node::getDisp(DiffDoub1 disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

void Node::getElasticDVLoad(DiffDoub1 ld[], vector<DesignVariable>& dvAr) {
	int i1;
	int dIndex;
	DiffDoub1 dVal;
	DiffDoub1 coef;
	string cat;
	int comp;
	
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		comp = thisDV.component - 1;
		if(cat == "elasticLoad") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			ld[comp].add(coef);
		}
	}
	return;
}

void Node::getThermalDVLoad(DiffDoub1& ld, vector<DesignVariable>& dvAr) {
	int dIndex;
	DiffDoub1 dVal;
	DiffDoub1 coef;
	string cat;
	
	ld.setVal(0.0);
	for (auto& dv : dVarLst) {
		dIndex = dv.intDat;
		DesignVariable& thisDV = dvAr[dIndex];
		thisDV.getValue(dVal);
		cat = thisDV.category;
		if (cat == "thermalLoad") {
			coef.setVal(dv.doubDat);
			coef.mult(dVal);
			ld.add(coef);
		}
	}
	return;
}

//end dup
 
//end skip 
 
 
 
 