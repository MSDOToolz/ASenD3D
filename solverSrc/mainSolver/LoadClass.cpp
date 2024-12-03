#include <string>
#include <cmath>
#include "LoadClass.h"
#include "SetClass.h"

using namespace std;


Load::Load() {
	type = "";
	activeTime[0] = 0.0;
	activeTime[1] = 1.0e+100;
	ndSetPtr = -1;
	elSetPtr = -1;
	return;
}

void Load::setActTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
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
	normalDir[0] = newNDir[0];
	normalDir[1] = newNDir[1];
	normalDir[2] = newNDir[2];
	double mag = normalDir[0]*normalDir[0] + normalDir[1]*normalDir[1] + normalDir[2]*normalDir[2];
	mag = 1.0/sqrt(mag);
	normalDir[0] = mag*normalDir[0];
	normalDir[1] = mag*normalDir[1];
	normalDir[2] = mag*normalDir[2];
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