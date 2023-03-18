#include "NodeClass.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"
#include <string>
#include <iostream>

using namespace std;


Node::Node(int newLab, double& newCrd) {
	label = newLab;
	coord[0] = newCrd[0];
	coord[1] = newCrd[1];
	coord[2] = newCrd[2];
	int i1;
	for (i1=0; i1 < 6; i1++) {
		displacement[i1] = 0.0;
		velocity[i1] = 0.0;
		acceleration[i1] = 0.0;
		dofIndex[i1] = 0;
	}
	dVars = IntList();
	coefs = DoubList();
}

void Node::addDesignVariable(int dIndex, double coef) {
	dVars.addEntry(dIndex);
	coefs.addEntry(coef);
}

//dup1
void Node::getCrd(Doub& crdOut, DVPt& dvAr) {
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

void Node::getDisp(Doub& disp) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].setVal(displacement[i1]);
	}
	return;
}

//end dup

//skip
void Node::getCrd(DiffDoub& crdOut, DVPt& dvAr) {
	crdOut[0].setVal(coord[0]);
	crdOut[1].setVal(coord[1]);
	crdOut[2].setVal(coord[2]);
	IntListEnt *currD = dVars.getFirst();
	DoubListEnt *currCoef = coefs.getFirst();
	int dIndex;
	DesignVariable *dPtr;
	DiffDoub dVal;
	int comp;
	string cat;
	DiffDoub coef;
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
//end skip

int Node::getDofIndex(int dof) {
	return dofIndex[dof-1];
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

void NodeList::addNode(int newLab, double newCrd[]) {
	Node *newNd = new Node(newLab,newCrd);
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