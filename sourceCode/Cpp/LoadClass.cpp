#include <string>
#include <cmath>
#include "LoadClass.h"
#include "SetClass.h"

using namespace std;


Load::Load(string newType) {
	type = newType;
	activeTime[0] = 0.0;
	activeTime[1] = 1.0e+100;
	ndSetPtr = NULL;
	elSetPtr = NULL;
	return;
}

void Load::setActTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
	return;
}

void Load::setNodeSet(string newSet) {
	nodeSet = newSet;
	return;
}

void Load::setElSet(string newSet) {
	elementSet = newSet;
	return;
}

void Load::setLoad(double newLd[]) {
	load[0] = newLd[0];
	load[1] = newLd[1];
	load[2] = newLd[2];
	load[3] = newLd[3];
	load[4] = newLd[4];
	load[5] = newLd[5];
	return;
}

void Load::setNormDir(double newNDir[]) {
	normalDir[0] = newLd[0];
	normalDir[1] = newLd[1];
	normalDir[2] = newLd[2];
	double mag = normalDir[0]*normalDir[0] + normalDir[1]*normalDir[1] + normalDir[2]*normalDir[2];
	mag = 1.0/sqrt(mag);
	normalDir[0] = mag*normalDir[0];
	normalDir[1] = mag*normalDir[1];
	normalDir[2] = mag*normalDir[2];
	return;
}

void Load::setNormTol(double newTol) {
	normTol = newTol;
	return;
}

void Load::setCenter(double newCent[]) {
	center[0] = newCent[0];
	center[1] = newCent[1];
	center[2] = newCent[2];
	return;
}

void Load::setAxis(double newAxis[]) {
	axis[0] = newAxis[0];
	axis[1] = newAxis[1];
	axis[2] = newAxis[2];
	double mag = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
	mag = 1.0/sqrt(mag);
	axis[0] = mag*axis[0];
	axis[1] = mag*axis[1];
	axis[2] = mag*axis[2];
	return;
}

void Load::setAngVel(double newAngVel) {
	angularVel = newAngVel;
	return;
}

void Load::setNext(Load *newNext) {
	nextLd = newNext;
	return;
}

// get

string Load::getType() {
	return type;
}

void Load::getActTime(double actTm[]) {
	actTm[0] = activeTime[0];
	actTm[1] = activeTime[1];
	return;
}

string Load::getNodeSet() {
	return nodeSet;
}

Set* Load::getNdSetPtr() {
	return ndSetPtr;
}

string Load::getElSet() {
	return elementSet;
}

Set* Load::getElSetPtr() {
	return elSetPtr;
}

void Load::getLoad(double ldOut[]) {
	ldOut[0] = load[0];
	ldOut[1] = load[1];
	ldOut[2] = load[2];
	ldOut[3] = load[3];
	ldOut[4] = load[4];
	ldOut[5] = load[5];
	return;
}

void Load::getNormDir(double ndOut[]) {
	ndOut[0] = normalDir[0];
	ndOut[1] = normalDir[1];
	ndOut[2] = normalDir[2];
	return;
}

double Load::getNormTol() {
	return normTol;
}

void Load::getCenter(double centOut[]) {
	centOut[0] = center[0];
	centOut[1] = center[1];
	centOut[2] = center[2];
	return;
}

void Load::getAxis(double axisOut[]) {
	axisOut[0] = axis[0];
	axisOut[1] = axis[1];
	axisOut[2] = axis[2];
	return;
}

double Load::getAngVel() {
	return angularVel;
}

Load* Load::getNext() {
	return nextLd;
}


//LoadList

LoadList::LoadList() {
	firstLoad = NULL;
	lastLoad = NULL;
	length = 0;
}

void LoadList::addLoad(Load *newLd) {
	if(!firstLoad) {
		firstLoad = newLd;
		lastLoad = newLd;
	} else {
		lastLoad->setNext(newLd);
		lastLoad = newLd;
	}
	length++;
}

int LoadList::getLength() {
	return length;
}

Load* LoadList::getFirst() {
	return firstLoad;
}