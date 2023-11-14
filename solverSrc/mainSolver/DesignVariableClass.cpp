#include "DesignVariableClass.h"
#include <string>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"
using namespace std;


DesignVariable::DesignVariable(string newCat) {
	category = newCat;
	component = 1;
	layer = 0;
	elSetName = "";
	elSetPtr = nullptr;
	ndSetName = "";
	ndSetPtr = nullptr;
	activeTime[0] = 0.0;
	activeTime[1] = 1.0e+100;
	value.setVal(0.0);
	diffVal.setVal(0.0);
	nextDV = nullptr;
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

void DesignVariable::setElsetPtr(Set* newPtr) {
	elSetPtr = newPtr;
	return;
}

void DesignVariable::setNdset(std::string newNdset) {
	ndSetName = newNdset;
	return;
}

void DesignVariable::setNdsetPtr(Set* newSet) {
	ndSetPtr = newSet;
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

void DesignVariable::addCompEl(int newEl) {
	compElList.addIfAbsent(newEl);
	return;
}

string DesignVariable::getCategory() {
	return category;
}

string DesignVariable::getElSet() {
	return elSetName;
}

Set* DesignVariable::getElsetPtr() {
	return elSetPtr;
}

string DesignVariable::getNdSet() {
	return ndSetName;
}

Set* DesignVariable::getNdsetPtr() {
	return ndSetPtr;
}

DoubList* DesignVariable::getCoefs() {
	return &coefs;
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

int DesignVariable::getLayer() {
	return layer;
}

IntListEnt* DesignVariable::getFirstEl() {
	return compElList.getFirst();
}

IntListEnt* DesignVariable::getFirstNd() {
	if (ndSetPtr) {
		return ndSetPtr->getFirstEntry();
	}
	return nullptr;
}

DesignVariable* DesignVariable::getNext() {
	return nextDV;
}

void DesignVariable::destroy() {
	coefs.destroy();
	compElList.destroy();
	return;
}


DVPt::DVPt() {
	ptr = nullptr;
	return;
}


DesVarList::DesVarList() {
	firstDVar = nullptr;
	lastDVar = nullptr;
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

void DesVarList::destroy() {
	DesignVariable* thisDV = firstDVar;
	DesignVariable* nextDV;
	while (thisDV) {
		nextDV = thisDV->getNext();
		thisDV->destroy();
		delete thisDV;
		thisDV = nextDV;
	}
	return;
}