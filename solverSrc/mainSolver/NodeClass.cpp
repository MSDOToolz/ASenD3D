#include "NodeClass.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"
#include <string>
#include <iostream>

using namespace std;


Node::Node(int newLab) {
	label = newLab;
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
	prevTemp = 0.0;
	pTempLF = 0.0;
	prevTdot = 0.0;
	prevFlDen = 0.0;
	pFlDenLF = 0.0;
	prevFlDenDot = 0.0;
	initialTemp = 0.0;
	initialTdot = 0.0;
	initialFlDen = 0.0;
	initialFlDenDot = 0.0;
	nextNd = nullptr;
}

void Node::setFluid(bool newFluid) {
	fluid = newFluid;
	return;
}

void Node::setNumDof(int nDof) {
	numDof = nDof;
	return;
}

void Node::setCrd(double newCrd[]) {
	coord[0] = newCrd[0];
	coord[1] = newCrd[1];
	coord[2] = newCrd[2];
	return;
}

void Node::setDofIndex(int dof, int index) {
	dofIndex[dof] = index;
	return;
}

void Node::setSortedRank(int newRank) {
	sortedRank = newRank;
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

void Node::setTemperature(double newTemp) {
	temperature = newTemp;
	return;
}

void Node::setTdot(double newTdot) {
	tempChangeRate = newTdot;
	return;
}

void Node::setFlDen(double newFlDen) {
	flDen = newFlDen;
	return;
}

void Node::setFlDenDot(double newFlDDot) {
	flDenDot = newFlDDot;
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

void Node::addToTemperature(double delTemp) {
	temperature += delTemp;
	return;
}

void Node::addToFlDen(double delFlDen) {
	flDen += delFlDen;
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

void Node::setInitialTemp(double newTemp) {
	initialTemp = newTemp;
	return;
}

void Node::setInitialTdot(double newTdot) {
	initialTdot = newTdot;
	return;
}

void Node::setInitialFlDen(double newFlDen) {
	initialFlDen = newFlDen;
	return;
}

void Node::setInitialFlDDot(double newFlDen) {
	initialFlDenDot = newFlDen;
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

void Node::setPrevTemp(double newTemp) {
	prevTemp = newTemp;
	return;
}

void Node::setPTempLF(double newTemp) {
	pTempLF = newTemp;
	return;
}

void Node::setPrevTdot(double newTdot) {
	prevTdot = newTdot;
	return;
}

void Node::setPrevFlDen(double newFlDen) {
	prevFlDen = newFlDen;
	return;
}

void Node::setPFlDenLF(double newDen) {
	pFlDenLF = newDen;
	return;
}

void Node::setPrevFlDDot(double newFlDDot) {
	prevFlDenDot = newFlDDot;
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

void Node::addDesignVariable(int dIndex, double coef) {
	dVars.addEntry(dIndex);
	coefs.addEntry(coef);
	return;
}

IntList* Node::getDesignVars() {
	return &dVars;
}

void Node::getCrd(double crdOut[]) {
	crdOut[0] = coord[0];
	crdOut[1] = coord[1];
	crdOut[2] = coord[2];
	return;
}

int Node::getLabel() {
	return label;
}

bool Node::isFluid() {
	return fluid;
}

int Node::getNumDof() {
	return numDof;
}

void Node::getDisp(double dispOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		dispOut[i1] = displacement[i1];
	}
	return;
}

void Node::getVel(double velOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		velOut[i1] = velocity[i1];
	}
	return;
}

void Node::getAcc(double accOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		accOut[i1] = acceleration[i1];
	}
	return;
}

void Node::getFlVel(double velOut[]) {
	velOut[0] = flVel[0];
	velOut[1] = flVel[1];
	velOut[2] = flVel[2];
	return;
}

void Node::getFlVelDot(double velDotOut[]) {
	velDotOut[0] = flVelDot[0];
	velDotOut[1] = flVelDot[1];
	velDotOut[2] = flVelDot[2];
	return;
}

double Node::getTemperature() {
	return temperature;
}

double Node::getTdot() {
	return tempChangeRate;
}

double Node::getFlDen() {
	return flDen;
}

double Node::getFlDenDot() {
	return flDenDot;
}

void Node::getPrevDisp(double dispOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		dispOut[i1] = prevDisp[i1];
	}
	return;
}

void Node::getPrevVel(double velOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		velOut[i1] = prevVel[i1];
	}
	return;
}

void Node::getPrevAcc(double accOut[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		accOut[i1] = prevAcc[i1];
	}
	return;
}

void Node::getPrevFlVel(double flVelOut[]) {
	flVelOut[0] = prevFlVel[0];
	flVelOut[1] = prevFlVel[1];
	flVelOut[2] = prevFlVel[2];
	return;
}

void Node::getPrevFlVelDot(double flVDOut[]) {
	flVDOut[0] = prevFlVelDot[0];
	flVDOut[1] = prevFlVelDot[1];
	flVDOut[2] = prevFlVelDot[2];
	return;
}

double Node::getPrevTemp() {
	return prevTemp;
}

double Node::getPrevTdot() {
	return prevTdot;
}

double Node::getPrevFlDen() {
	return prevFlDen;
}

double Node::getPrevFlDenDot() {
	return prevFlDenDot;
}

//dup1
void Node::getCrd(DiffDoub0 crdOut[], DesignVariable* dvAr[]) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars.getFirst();
	DoubListEnt *currCoef = coefs.getFirst();
	int dIndex;
	DesignVariable *dPtr;
	DiffDoub0 dVal;
	int comp;
	string cat;
	DiffDoub0 coef;
	while(currD) {
		dIndex = currD->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		comp = dPtr->getComponent() - 1;
		if(cat == "nodeCoord") {
			coef.setVal(currCoef->value);
			coef.mult(dVal);
			crdOut[comp].add(coef);
		}
		currD = currD->next;
		currCoef = currCoef->next;
	}
	return;
}

void Node::getDefCrd(DiffDoub0 crdOut[], DesignVariable* dvAr[]) {
	DiffDoub0 dsp[6];
	getCrd(crdOut, dvAr);
	getDisp(dsp);
	crdOut[0].add(dsp[0]);
	crdOut[1].add(dsp[1]);
	crdOut[2].add(dsp[2]);
	return;
}

void Node::getDisp(DiffDoub0 disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

void Node::getElasticDVLoad(DiffDoub0 ld[], DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *dPtr;
	DiffDoub0 dVal;
	DiffDoub0 coef;
	string cat;
	int comp;
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	thisDV = dVars.getFirst();
	thisCoef = coefs.getFirst();
	while(thisDV) {
		dIndex = thisDV->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		comp = dPtr->getComponent() - 1;
		if(cat == "elasticLoad") {
			coef.setVal(thisCoef->value);
			coef.mult(dVal);
			ld[comp].add(coef);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	return;
}

void Node::getThermalDVLoad(DiffDoub0& ld, DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* dPtr;
	DiffDoub0 dVal;
	DiffDoub0 coef;
	string cat;
	int comp;
	
	ld.setVal(0.0);
	thisDV = dVars.getFirst();
	thisCoef = coefs.getFirst();
	while (thisDV) {
		dIndex = thisDV->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		if (cat == "thermalLoad") {
			coef.setVal(thisCoef->value);
			coef.mult(dVal);
			ld.add(coef);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void Node::getCrd(DiffDoub1 crdOut[], DesignVariable* dvAr[]) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars.getFirst();
	DoubListEnt *currCoef = coefs.getFirst();
	int dIndex;
	DesignVariable *dPtr;
	DiffDoub1 dVal;
	int comp;
	string cat;
	DiffDoub1 coef;
	while(currD) {
		dIndex = currD->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		comp = dPtr->getComponent() - 1;
		if(cat == "nodeCoord") {
			coef.setVal(currCoef->value);
			coef.mult(dVal);
			crdOut[comp].add(coef);
		}
		currD = currD->next;
		currCoef = currCoef->next;
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

void Node::getElasticDVLoad(DiffDoub1 ld[], DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *dPtr;
	DiffDoub1 dVal;
	DiffDoub1 coef;
	string cat;
	int comp;
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	thisDV = dVars.getFirst();
	thisCoef = coefs.getFirst();
	while(thisDV) {
		dIndex = thisDV->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		comp = dPtr->getComponent() - 1;
		if(cat == "elasticLoad") {
			coef.setVal(thisCoef->value);
			coef.mult(dVal);
			ld[comp].add(coef);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	return;
}

void Node::getThermalDVLoad(DiffDoub1& ld, DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* dPtr;
	DiffDoub1 dVal;
	DiffDoub1 coef;
	string cat;
	int comp;
	
	ld.setVal(0.0);
	thisDV = dVars.getFirst();
	thisCoef = coefs.getFirst();
	while (thisDV) {
		dIndex = thisDV->value;
		dPtr = dvAr[dIndex];
		dPtr->getValue(dVal);
		cat = dPtr->getCategory();
		if (cat == "thermalLoad") {
			coef.setVal(thisCoef->value);
			coef.mult(dVal);
			ld.add(coef);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	return;
}

//end dup
 
//end skip 
 
 
int Node::getDofIndex(int dof) {
	return dofIndex[dof];
}

int Node::getSortedRank() {
	return sortedRank;
}

Node* Node::getNext() {
	return nextNd;
}

void Node::setNext(Node *newNext) {
	nextNd = newNext;
}


NodeList::NodeList() {
	firstNode = nullptr;
	lastNode = nullptr;
	length = 0;
}

void NodeList::addNode(Node *newNd) {
	if(!firstNode) {
		firstNode = newNd;
		lastNode = newNd;
	} else {
		lastNode->setNext(newNd);
		lastNode = newNd;
	}
	length++;
}

int NodeList::getLength() {
	return length;
}

Node* NodeList::getFirst() {
	return firstNode;
}

NodeList::~NodeList() {
	Node* thisNd = firstNode;
	Node* nextNd;
	while (thisNd) {
		nextNd = thisNd->getNext();
		delete thisNd;
		thisNd = nextNd;
	}
	firstNode = nullptr;
	lastNode = nullptr;
	length = 0;
	return;
}