#include "SetClass.h"
#include <string>
#include "ListEntClass.h"

using namespace std;

Set::Set() {
    name = "";
    nextSet = nullptr;
	return;
}

void Set::setName(string newNm) {
	name = newNm;
	return;
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

void Set::destroy() {
	labels.destroy();
	return;
}

SetPt::SetPt() {
	ptr = nullptr;
	return;
}

SetList::SetList() {
	firstSet = nullptr;
	lastSet = nullptr;
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

Set* SetList::getLast() {
	return lastSet;
}

void SetList::destroy() {
	Set* thisSet = firstSet;
	Set* nextSet;
	while (thisSet) {
		nextSet = thisSet->getNext();
		thisSet->destroy();
		delete thisSet;
		thisSet = nextSet;
	}
	firstSet = nullptr;
	lastSet = nullptr;
	length = 0;
	return;
}