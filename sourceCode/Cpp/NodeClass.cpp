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
	int i1;
	for (i1=0; i1 < 6; i1++) {
		displacement[i1] = 0.0;
		velocity[i1] = 0.0;
		acceleration[i1] = 0.0;
		dofIndex[i1] = 0;
	}
	nextNd = NULL;
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

void Node::setInitialDisp(double newDisp[]) {
	initialDisp[0] = newDisp[0];
	initialDisp[1] = newDisp[1];
	initialDisp[2] = newDisp[2];
	initialDisp[3] = newDisp[3];
	initialDisp[4] = newDisp[4];
	initialDisp[5] = newDisp[5];
	return;
}

void Node::setInitialVel(double newVel[]) {
	initialVel[0] = newVel[0];
	initialVel[1] = newVel[1];
	initialVel[2] = newVel[2];
	initialVel[3] = newVel[3];
	initialVel[4] = newVel[4];
	initialVel[5] = newVel[5];
	return;
}

void Node::setInitialAcc(double newAcc[]) {
	initialAcc[0] = newAcc[0];
	initialAcc[1] = newAcc[1];
	initialAcc[2] = newAcc[2];
	initialAcc[3] = newAcc[3];
	initialAcc[4] = newAcc[4];
	initialAcc[5] = newAcc[5];
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

int Node::getNumDof() {
	return numDof;
}

//dup1
void Node::getCrd(Doub crdOut[], DVPt dvAr[]) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars.getFirst();
	DoubListEnt *currCoef = coefs.getFirst();
	int dIndex;
	DesignVariable *dPtr;
	Doub dVal;
	int comp;
	string cat;
	Doub coef;
	while(currD) {
		dIndex = currD->value;
		dPtr = dvAr[dIndex].ptr;
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

//end dup

int Node::getDofIndex(int dof) {
	return dofIndex[dof];
}

Node* Node::getNext() {
	return nextNd;
}

void Node::setNext(Node *newNext) {
	nextNd = newNext;
}

void Node::destroy() {
	dVars.destroy();
	coefs.destroy();
	return;
}

NdPt::NdPt() {
	ptr = NULL;
}


NodeList::NodeList() {
	firstNode = NULL;
	lastNode = NULL;
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