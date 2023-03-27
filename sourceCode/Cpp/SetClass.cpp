#include "SetClass.h"
#include <string>
#include "ListEntClass.h"

using namespace std;

Set::Set(string newNm) {
    name = newNm;
    nextSet = NULL;
}

void Set::addEntry(int newLabel) {
    labels.addEntry(newLabel);
}

string Set::getName() {
    return name;
}

int Set::getLength() {
	return labels.getLength();
}

IntListEnt* Set::getFirstEntry() {
    return labels.getFirst();
}

Set* Set::getNext() {
    return nextSet;
}

void Set::setNext(Set* newNext) {
    nextSet = newNext;
    return;
}


SetList::SetList() {
	firstSet = NULL;
	lastSet = NULL;
	length = 0;
}

void SetList::addSet(Set *newSet) {
	if(!firstSet) {
		firstSet = newSet;
		lastSet = newSet;
	} else {
		lastSet->setNext(newSet);
		lastSet = newSet;
	}
	length++;
}

int SetList::getLength() {
	return length;
}

Set* SetList::getFirst() {
	return firstSet;
}