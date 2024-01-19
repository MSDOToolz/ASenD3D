#include "NodeClass.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"
#include <string>
#include <iostream>

using namespace std;


Node::Node(int newLab) {
	label = newLab;
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
	temperature = 0.0;
	tempChangeRate = 0.0;
	prevTemp = 0.0;
	prevTdot = 0.0;
	initialTemp = 0.0;
	initialTdot = 0.0;
	dVars = new IntList;
	coefs = new DoubList;
	nextNd = nullptr;
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

void Node::setTemperature(double newTemp) {
	temperature = newTemp;
	return;
}

void Node::setTdot(double newTdot) {
	tempChangeRate = newTdot;
	return;
}

void Node::addToDisplacement(double delDisp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] += delDisp[i1];
	}	
	return;
}

void Node::addToTemperature(double delTemp) {
	temperature += delTemp;
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

void Node::setInitialTemp(double newTemp) {
	initialTemp = newTemp;
	return;
}

void Node::setInitialTdot(double newTdot) {
	initialTdot = newTdot;
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

void Node::setPrevTemp(double newTemp) {
	prevTemp = newTemp;
	return;
}

void Node::setPrevTdot(double newTdot) {
	prevTdot = newTdot;
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

void Node::initializeTemp() {
	prevTemp = initialTemp;
	prevTdot = initialTdot;
	temperature = initialTemp;
	return;
}

void Node::updateVelAcc(double nmBeta, double nmGamma, double delT) {
	double c1 = 1.0/(delT*delT*(nmBeta-nmGamma));
	double c2 = delT*delT*(0.5 + nmBeta - nmGamma);
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		acceleration[i1] = c1*(prevDisp[i1] - displacement[i1] + delT*prevVel[i1] + c2*prevAcc[i1]);
		velocity[i1] = prevVel[i1] + delT*((1.0-nmGamma)*prevAcc[i1] + nmGamma*acceleration[i1]);
	}
	
	return;
}

void Node::updateTdot(double nmGamma, double delT) {
	double c1 = 1.0 / nmGamma;
	double c2 = -1.0 / delT;
	tempChangeRate = c1 * (c2 * (prevTemp - temperature) - (1.0 - nmGamma) * prevTdot);
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

void Node::advanceTemp() {
	prevTemp = temperature;
	prevTdot = tempChangeRate;
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

void Node::backstepTemp() {
	temperature = prevTemp;
	tempChangeRate = prevTdot;
	return;
}

void Node::addDesignVariable(int dIndex, double coef) {
	dVars->addEntry(dIndex);
	coefs->addEntry(coef);
	return;
}

IntList* Node::getDesignVars() {
	return dVars;
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

double Node::getTemperature() {
	return temperature;
}

double Node::getTdot() {
	return tempChangeRate;
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

double Node::getPrevTemp() {
	return prevTemp;
}

double Node::getPrevTdot() {
	return prevTdot;
}

//dup1
void Node::getCrd(Doub crdOut[], DesignVariable* dvAr[]) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars->getFirst();
	DoubListEnt *currCoef = coefs->getFirst();
	int dIndex;
	DesignVariable *dPtr;
	Doub dVal;
	int comp;
	string cat;
	Doub coef;
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

void Node::getDisp(Doub disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

void Node::getElasticDVLoad(Doub ld[], DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *dPtr;
	Doub dVal;
	Doub coef;
	string cat;
	int comp;
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	thisDV = dVars->getFirst();
	thisCoef = coefs->getFirst();
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

void Node::getThermalDVLoad(Doub& ld, DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* dPtr;
	Doub dVal;
	Doub coef;
	string cat;
	int comp;
	
	ld.setVal(0.0);
	thisDV = dVars->getFirst();
	thisCoef = coefs->getFirst();
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
 
//DiffDoub versions: 
//dup1
void Node::getCrd(DiffDoub crdOut[], DesignVariable* dvAr[]) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars->getFirst();
	DoubListEnt *currCoef = coefs->getFirst();
	int dIndex;
	DesignVariable *dPtr;
	DiffDoub dVal;
	int comp;
	string cat;
	DiffDoub coef;
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

void Node::getDisp(DiffDoub disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

void Node::getElasticDVLoad(DiffDoub ld[], DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *dPtr;
	DiffDoub dVal;
	DiffDoub coef;
	string cat;
	int comp;
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].setVal(0.0);
	}
	thisDV = dVars->getFirst();
	thisCoef = coefs->getFirst();
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

void Node::getThermalDVLoad(DiffDoub& ld, DesignVariable* dvAr[]) {
	int i1;
	int dIndex;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* dPtr;
	DiffDoub dVal;
	DiffDoub coef;
	string cat;
	int comp;
	
	ld.setVal(0.0);
	thisDV = dVars->getFirst();
	thisCoef = coefs->getFirst();
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

void Node::destroy() {
	dVars->destroy();
	delete dVars;
	coefs->destroy();
	delete coefs;
	return;
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

void NodeList::destroy() {
	Node* thisNd = firstNode;
	Node* nextNd;
	while (thisNd) {
		nextNd = thisNd->getNext();
		thisNd->destroy();
		delete thisNd;
		thisNd = nextNd;
	}
	firstNode = nullptr;
	lastNode = nullptr;
	length = 0;
	return;
}