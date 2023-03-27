#include "ObjectiveClass.h"
#include "ListEntClass.h"

using namespace std;

ObjectiveTerm::ObjectiveTerm(string newCat) {
	category = newCat;
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

void ObjectiveTerm::setNdset(string newNdset) {
	ndSetName = newNdset;
	return;
}

void ObjectiveTerm::setTgtTag(string newTag) {
	tgtTag = newTag;
	return;
}

void ObjectiveTerm::addTargetValue(double newTgt) {
	tgtVals.addEntry(newTgt);
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