#include "DesignVariableClass.h"
#include <string>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
using namespace std;


DesignVariable::DesignVariable(string newCat) {
	category = newCat;
	component = 1;
	layer = 0;
	elSetName = "";
	ndSetName = "";
	activeTime[0] = 0.0;
	activeTime[1] = 1.0e+100;
	value.setVal(0.0);
	diffVal.setVal(0.0);
	nextDV = NULL;
}

void DesignVariable::setActiveTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
}

void DesignVariable::setLayer(int newLay) {
	layer = newLay;
	return;
}

void DesignVariable::setElset(string newElset) {
	elSetName = newElset;
	return;
}

void DesignVariable::setNdset(std::string newNdset) {
	ndSetName = newNdset;
	return;
}

void DesignVariable::setValue(double newVal) {
	value.setVal(newVal);
	return;
}

void DesignVariable::setDiffVal(double newVal, double dNewVal) {
	diffVal.setVal(newVal,dNewVal);
	return;
}

void DesignVariable::setComponent(int newComp) {
	component = newComp;
	return;
}

void DesignVariable::setNext(DesignVariable* newNext) {
	nextDV = newNext;
}

void DesignVariable::addCoefficient(double newCoef) {
	coefs.addEntry(newCoef);
	return;
}

string DesignVariable::getCategory() {
	return category;
}

void DesignVariable::getValue(Doub& inp) {
	inp.setVal(value);
	return;
}

void DesignVariable::getValue(DiffDoub& inp) {
	inp.setVal(diffVal);
	return;
}

int DesignVariable::getComponent() {
	return component;
}

DesignVariable* DesignVariable::getNext() {
	return nextDV;
}

void DesignVariable::destroy() {
	coefs.destroy();
	return;
}


DVPt::DVPt() {
	ptr = NULL;
	return;
}


DesVarList::DesVarList() {
	firstDVar = NULL;
	lastDVar = NULL;
	length = 0;
}

void DesVarList::addDVar(DesignVariable* newDVar) {
	if(!firstDVar) {
		firstDVar = newDVar;
		lastDVar = newDVar;
	} else {
		lastDVar->setNext(newDVar);
		lastDVar = newDVar;
	}
	length++;
}

int DesVarList::getLength() {
	return length;
}

DesignVariable* DesVarList::getFirst() {
	return firstDVar;
}