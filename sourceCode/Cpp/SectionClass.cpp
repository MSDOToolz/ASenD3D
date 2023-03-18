#include "SectionClass.h"
#include <string>

using namespace std;


Material::Material(string newName) {
	name = newName;
}

double Material::getDensity() {
	return density;
}

double* Material::getElastic() {
	return &elastic[0];
}

double* Material::getStiffMat() {
	return &stiffness[0];
}


MaterialList::MaterialList() {
	firstMat = NULL;
	lastMat = NULL;
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


Layer::Layer(string newNm, double newThk, double newAng) {
	matName = newNm;
	thickness = newThk;
	angle = newAng;
}

string Layer::getMatName() {
	return matName;
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
}

void Layer::setMatPtr(Material *newPtr) {
	matPtr = newPtr;
	return;
}


LayerList::LayerList() {
	firstLay = NULL;
	lastLay = NULL;
	length = 0;
}

void LayerList::addLayer(string newNm, double newThk, double newAng) {
	Layer *newLay = new Layer(newNm,newThk,newAng);
	if(!firstLay) {
		firstLay = newLay;
		lastLay = newLay;
	} else {
		lastLay->setNext(newLay);
		lastLay = newLay;
	}
	length++;
}

int LayerList::getLength() {
	return length;
}

Layer* LayerList::getFirst() {
	return firstLay;
}


Section::Section(string newType) {
	type = newType;
	return;
}

Material* Section::getMaterial() {
	return matPtr;
}

void Section::setNext(Section *newNext) {
	nextSection = newNext;
	return;
}


SectionList::SectionList() {
	firstSec = NULL;
	lastSec = NULL;
	length = 0;
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
}

int SectionList::getLength() {
	return length;
}

Section* SectionList::getFirst() {
	return firstSec;
}