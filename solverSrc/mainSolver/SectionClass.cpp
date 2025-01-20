#include "SectionClass.h"
#include "matrixFunctions.h"
#include <string>
#include <cmath>
#include <vector>

using namespace std;

const int max_int = 2000000000;

Material::Material() {
	name = "";
	int i1;
	for(i1 = 0; i1 < 36; i1++) {
		stiffness[i1] = 0.0;
		damping[i1] = 0.0;
	}
	return;
}

void Material::setStiffness(int row, int col, double val) {
	stiffness[(6 * row + col)] = val;
	stiffness[(6 * col + row)] = val;
	return;
}

void Material::setDamping(int row, int col, double val) {
	damping[(6 * row + col)] = val;
	damping[(6 * col + row)] = val;
	return;
}

Fluid::Fluid() {
	name = "";
	return;
}

Layer::Layer() {
	matName = "";
	matPtr = max_int;
	return;
}

Section::Section() {
	type = "";
	int i1;
	for (i1 = 0; i1 < 36; i1++) {
		stiffness[i1] = 0.0;
		mass[i1] = 0.0;
		damping[i1] = 0.0;
	}
	for (i1 = 1; i1 < 8; i1++) {
		orientation[i1] = 0.0;
	}
	orientation[0] = 1.0;
	orientation[4] = 1.0;
	orientation[8] = 1.0;
	zOffset = 0.0;
	area = 0.0;
	for (i1 = 0; i1 < 5; i1++) {
		areaMoment[i1] = 0.0;
	}
	for (i1 = 0; i1 < 6; i1++) {
		expLoadCoef[i1] = 0.0;
	}
	conductivity = 0.0;
	specHeat = 0.0;
	potCoef = 0.0;
	potExp = 1.0;
	dampCoef = 0.0;
	dampExp = 1.0;
	condCoef = 0.0;
	radCoef = 0.0;
	denVisCoef = 0.0;
	tempVisCoef = 0.0;
	turbVisCoef = 0.0;
	gradVTurbCoef = 0.0;
	dissTurbCoef = 0.0;
	enthCoef = 0.0;
	enthExp = 1.0;
	presCoef = 0.0;
	presExp = 1.0;
	refTemp = 0.0;
	refDen = 0.0;
	refTurbE = 0.0;
	refEnth = 0.0;
	matPtr = max_int;
	flPtr = max_int;
	return;
}

void Section::setOrientation(double newOri[]) {
	orientation[0] = newOri[0];
	orientation[1] = newOri[1];
	orientation[2] = newOri[2];
	orientation[3] = newOri[3];
	orientation[4] = newOri[4];
	orientation[5] = newOri[5];
	
	double mag = orientation[0]*orientation[0] + orientation[1]*orientation[1] + orientation[2]*orientation[2];
	mag = 1.0/sqrt(mag);
	orientation[0] = mag*orientation[0];
	orientation[1] = mag*orientation[1];
	orientation[2] = mag*orientation[2];
	
	crossProd(&orientation[6],&orientation[0],&orientation[3]);
	
	mag = orientation[6]*orientation[6] + orientation[7]*orientation[7] + orientation[8]*orientation[8];
	mag = 1.0/sqrt(mag);
	orientation[6] = mag*orientation[6];
	orientation[7] = mag*orientation[7];
	orientation[8] = mag*orientation[8];
	
	crossProd(&orientation[3],&orientation[6],&orientation[0]);
	
	return;
}

void Section::setAreaMoment(double newI[]) {
    areaMoment[0] = newI[0];
	areaMoment[1] = newI[1];
	areaMoment[2] = newI[2];
	areaMoment[3] = newI[3];
	areaMoment[4] = newI[4];
	return;
}

void Section::setStiffness(int row, int col, double val) {
	stiffness[(6*row+col)] = val;
	stiffness[(6*col+row)] = val;
	return;
}

void Section::setMass(int row, int col, double val) {
	mass[(6*row+col)] = val;
	mass[(6*col+row)] = val;
	return;
}

void Section::setDamping(int row, int col, double val) {
	damping[(6 * row + col)] = val;
	damping[(6 * col + row)] = val;
	return;
}

void Section::setExpLd(double newExpLd[]) {
	expLoadCoef[0] = newExpLd[0];
	expLoadCoef[1] = newExpLd[1];
	expLoadCoef[2] = newExpLd[2];
	expLoadCoef[3] = newExpLd[3];
	expLoadCoef[4] = newExpLd[4];
	expLoadCoef[5] = newExpLd[5];
	return;
}

int Section::getLayerMatPtr(int layi) {
	int i1 = 0;
	for (auto& lay : layers) {
		if (i1 == layi) {
			return lay.matPtr;
		}
		i1++;
	}
}