#include "ObjectiveClass.h"


ObjectiveTerm::ObjectiveTerm(string newCat) {
	category = newCat;
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