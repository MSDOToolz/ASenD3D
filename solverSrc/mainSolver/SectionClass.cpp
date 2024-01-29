#include "SectionClass.h"
#include "matrixFunctions.h"
#include <string>
#include <cmath>

using namespace std;


Material::Material(string newName) {
	name = newName;
	int i1;
	for(i1 = 0; i1 < 36; i1++) {
		stiffness[i1] = 0.0;
		damping[i1] = 0.0;
	}
	return;
}

void Material::setDensity(double newDen) {
	density = newDen;
	return;
}

void Material::setModulus(double newMod[]) {
	modulus[0] = newMod[0];
	modulus[1] = newMod[1];
	modulus[2] = newMod[2];
	return;
}

void Material::setPoissonRatio(double newPR[]) {
	poissonRatio[0] = newPR[0];
	poissonRatio[1] = newPR[1];
	poissonRatio[2] = newPR[2];
	return;
}

void Material::setShearMod(double newMod[]) {
	shearMod[0] = newMod[0];
	shearMod[1] = newMod[1];
	shearMod[2] = newMod[2];
	return;
}

void Material::setStiffness(int row, int col, double val) {
	stiffness[(6*row+col)] = val;
	stiffness[(6*col+row)] = val;
	return;
}

void Material::setConductivity(double newCond[]) {
	conductivity[0] = newCond[0];
	conductivity[1] = newCond[1];
	conductivity[2] = newCond[2];
	conductivity[3] = newCond[3];
	conductivity[4] = newCond[4];
	conductivity[5] = newCond[5];
	return;
}

void Material::setExpansion(double newExp[]) {
	expansion[0] = newExp[0];
	expansion[1] = newExp[1];
	expansion[2] = newExp[2];
	expansion[3] = newExp[3];
	expansion[4] = newExp[4];
	expansion[5] = newExp[5];
	return;
}

void Material::setSpecHeat(double newSpecHeat) {
	specHeat = newSpecHeat;
	return;
}

void Material::setDamping(int row, int col, double val) {
	damping[6 * row + col] = val;
	damping[6 * row + col] = val;
	return;
}

void Material::setMaxTenStress(double newMaxStr[]) {
	maxTenStress[0] = newMaxStr[0];
	maxTenStress[1] = newMaxStr[1];
	maxTenStress[2] = newMaxStr[2];
	return;
}

void Material::setMaxCompStress(double newMaxStr[]) {
	maxCompStress[0] = newMaxStr[0];
	maxCompStress[1] = newMaxStr[1];
	maxCompStress[2] = newMaxStr[2];
	return;
}

void Material::setMaxShearStress(double newMaxStr[]) {
	maxShearStress[0] = newMaxStr[0];
	maxShearStress[1] = newMaxStr[1];
	maxShearStress[2] = newMaxStr[2];
	return;
}

void Material::setMaxTenStrain(double newMaxStr[]) {
	maxTenStrain[0] = newMaxStr[0];
	maxTenStrain[1] = newMaxStr[1];
	maxTenStrain[2] = newMaxStr[2];
	return;
}

void Material::setMaxCompStrain(double newMaxStr[]) {
	maxCompStrain[0] = newMaxStr[0];
	maxCompStrain[1] = newMaxStr[1];
	maxCompStrain[2] = newMaxStr[2];
	return;
}

void Material::setMaxShearStrain(double newMaxStr[]) {
	maxShearStrain[0] = newMaxStr[0];
	maxShearStrain[1] = newMaxStr[1];
	maxShearStrain[2] = newMaxStr[2];
	return;
}

void Material::setMaxStrEng(double newMax) {
	maxStrEng = newMax;
	return;
}

void Material::setMaxMises(double newMax) {
	maxMises = newMax;
	return;
}

void Material::setNext(Material *newMat) {
	nextMat = newMat;
	return;
}

string Material::getName() {
	return name;
}

double Material::getDensity() {
	return density;
}

double* Material::getModulus() {
	return &modulus[0];
}

double* Material::getPoissonRatio() {
	return &poissonRatio[0];
}

double* Material::getShearMod() {
	return &shearMod[0];
}

double* Material::getStiffMat() {
	return &stiffness[0];
}

double* Material::getThermExp() {
	return &expansion[0];
}

double* Material::getConductivity() {
	return &conductivity[0];
}

double Material::getSpecificHeat() {
	return specHeat;
}

double* Material::getDamping() {
	return &damping[0];
}

Material* Material::getNext() {
	return nextMat;
}

MaterialList::MaterialList() {
	firstMat = nullptr;
	lastMat = nullptr;
	length = 0;
}

void MaterialList::addMaterial(Material *newMat) {
	if(!firstMat) {
		firstMat = newMat;
		lastMat = newMat;
	} else {
		lastMat->setNext(newMat);
		lastMat = newMat;
	}
	length++;
}

int MaterialList::getLength() {
	return length;
}

Material* MaterialList::getFirst() {
	return firstMat;
}

MaterialList::~MaterialList() {
	Material* thisMat = firstMat;
	Material* nextMat;
	while (thisMat) {
		nextMat = thisMat->getNext();
		delete thisMat;
		thisMat = nextMat;
	}
	firstMat = nullptr;
	lastMat = nullptr;
	length = 0;
	return;
}

Layer::Layer(string newNm) {
	matName = newNm;
	return;
}

void Layer::setThickness(double newThk) {
	thickness = newThk;
	return;
}
		
void Layer::setAngle(double newAngle) {
	angle = newAngle;
	return;
}

string Layer::getMatName() {
	return matName;
}

Material* Layer::getMatPt() {
	return matPtr;
}

double Layer::getThickness() {
	return thickness;
}

double Layer::getAngle() {
	return angle;
}

Layer* Layer::getNext() {
	return nextLay;
}

void Layer::setNext(Layer *newNext) {
	nextLay = newNext;
	return;
}

void Layer::setMatPtr(Material *newPtr) {
	matPtr = newPtr;
	return;
}


LayerList::LayerList() {
	firstLay = nullptr;
	lastLay = nullptr;
	length = 0;
}

void LayerList::addLayer(Layer *newLay) {
	if(!firstLay) {
		firstLay = newLay;
		lastLay = newLay;
	} else {
		lastLay->setNext(newLay);
		lastLay = newLay;
	}
	length++;
	return;
}

int LayerList::getLength() {
	return length;
}

Layer* LayerList::getFirst() {
	return firstLay;
}

LayerList::~LayerList() {
	Layer* thisLay = firstLay;
	Layer* nextLay;
	while (thisLay) {
		nextLay = thisLay->getNext();
		delete thisLay;
		thisLay = nextLay;
	}
	firstLay = nullptr;
	lastLay = nullptr;
	length = 0;
	return;
}

Section::Section(string newType) {
	type = newType;
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
	layers = new LayerList;
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
	matPtr = nullptr;
	nextSection = nullptr;
	return;
}

void Section::setElset(string newSet) {
	elSetName = newSet;
	return;
}

void Section::setMaterial(string newMat) {
	matName = newMat;
	return;
}

void Section::setMatPtr(Material *newMat) {
	matPtr = newMat;
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

void Section::setZOffset(double newZOff) {
	zOffset = newZOff;
	return;
}

void Section::addLayer(Layer *newLay) {
	layers->addLayer(newLay);
	return;
}

void Section::setArea(double newArea) {
	area = newArea;
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

void Section::setPolarMoment(double newJ) {
	polarMoment = newJ;
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

void Section::setConductivity(double newCond) {
	conductivity = newCond;
	return;
}

void Section::setSpecHeat(double specHeat) {
	specHeat = specHeat;
	return;
}

void Section::setMassPerEl(double newMass) {
	massPerEl = newMass;
	return;
}

void Section::setPotCoef(double newCoef) {
	potCoef = newCoef;
	return;
}

void Section::setPotExp(double newExp) {
	potExp = newExp;
	return;
}

void Section::setDampCoef(double newCoef) {
	dampCoef = newCoef;
	return;
}

void Section::setDampExp(double newExp) {
	dampExp = newExp;
	return;
}

string Section::getElset() {
	return elSetName;
}

string Section::getMaterial() {
	return matName;
}

Material* Section::getMatPtr() {
	return matPtr;
}

double* Section::getOrientation() {
	return &orientation[0];
}

double Section::getZOffset() {
	return zOffset;
}

int Section::getNumLayers() {
	return layers->getLength();
}

Layer* Section::getFirstLayer() {
	return layers->getFirst();
}

double Section::getArea() {
	return area;
}

double* Section::getAreaMoment() {
	return &areaMoment[0];
}

double Section::getPolarMoment() {
	return polarMoment;
}

double* Section::getStiffMat() {
	return &stiffness[0];
}

double* Section::getMassMat() {
	return &mass[0];
}

double* Section::getDampMat() {
	return &damping[0];
}

double* Section::getExpLoad() {
	return &expLoadCoef[0];
}

double Section::getConductivity() {
	return conductivity;
}

double Section::getSpecificHeat() {
	return specHeat;
}

double Section::getMassPerEl() {
	return massPerEl;
}

double Section::getPotCoef() {
	return potCoef;
}

double Section::getPotExp() {
	return potExp;
}

double Section::getDampCoef() {
	return dampCoef;
}

double Section::getDampExp() {
	return dampExp;
}

Section* Section::getNext() {
	return nextSection;
}

void Section::setNext(Section *newNext) {
	nextSection = newNext;
	return;
}

Section::~Section() {
	delete layers;
	return;
}

SectionList::SectionList() {
	firstSec = nullptr;
	lastSec = nullptr;
	length = 0;
	return;
}

void SectionList::addSection(Section *newSec) {
	if(!firstSec) {
		firstSec = newSec;
		lastSec = newSec;
	} else {
		lastSec->setNext(newSec);
		lastSec = newSec;
	}
	length++;
	return;
}

int SectionList::getLength() {
	return length;
}

Section* SectionList::getFirst() {
	return firstSec;
}

SectionList::~SectionList() {
	Section* thisSec = firstSec;
	Section* nextSec;
	while (thisSec) {
		nextSec = thisSec->getNext();
		delete thisSec;
		thisSec = nextSec;
	}
	firstSec = nullptr;
	lastSec = nullptr;
	length = 0;
	return;
}