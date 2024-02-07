#include "FaceClass.h"
#include "DiffDoubClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"
#include "matrixFunctions.h"

using namespace std;

Face::Face(int newNumNds) {
	numNds = newNumNds;
	onSurf = true;
	next = nullptr;
	return;
}

void Face::setNode(int place, int locNd, int globNd) {
	locNodes[place] = locNd;
	globNodes[place] = globNd;
	return;
}

void Face::setOnSurface(bool newOnSurf) {
	onSurf = newOnSurf;
	return;
}

void Face::setNext(Face *newNext) {
	next = newNext;
}

void Face::sortedNodes(int srtNds[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int swap;
	for (i1 = 0; i1 < numNds; i1++) {
		srtNds[i1] = globNodes[i1];
	}
	i3 = numNds - 1;
	for (i1 = 0; i1 < i3; i1++) {
		for (i2 = 0; i2 < i3; i2++) {
			i4 = i2 + 1;
			if(srtNds[i4] < srtNds[i2]) {
				swap = srtNds[i2];
				srtNds[i2] = srtNds[i4];
				srtNds[i4] = swap;
			}
		}
	}
	return;
}

int Face::getNumNds() {
	return numNds;
}

int* Face::getLocNds() {
	return locNodes;
}

bool Face::onSurface() {
	return onSurf;
}

int Face::getLowNd() {
	int i1;
	int lowNd = globNodes[0];
	for (i1 = 1; i1 < numNds; i1++) {
		if(globNodes[i1] < lowNd) {
			lowNd = globNodes[i1];
		}
	}
	return lowNd;
}

Face* Face::getNext() {
	return next;
}

//dup1

void Face::getAreaNormal(Doub& area, Doub norm[], Node* ndAr[], DesignVariable* dvAr[]) {
	Doub v1[3];
	Doub v2[3];
	Doub tmpV[3];
	Doub tmp;

	if (numNds == 4) {
		ndAr[globNodes[2]]->getCrd(v1, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[3]]->getCrd(v2, dvAr);
		ndAr[globNodes[1]]->getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}
	else {
		ndAr[globNodes[1]]->getCrd(v1, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[2]]->getCrd(v2, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}

	crossProd(norm, v1, v2);
	area.setVal(norm[0]);
	area.sqr();
	tmp.setVal(norm[1]);
	tmp.sqr();
	area.add(tmp);
	tmp.setVal(norm[2]);
	tmp.sqr();
	area.add(tmp);
	area.sqt();

	tmp.setVal(1.0);
	tmp.dvd(area);
	norm[0].mult(tmp);
	norm[1].mult(tmp);
	norm[2].mult(tmp);

	tmp.setVal(0.5);
	area.mult(tmp);

	return;
}

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

void Face::getAreaNormal(DiffDoub& area, DiffDoub norm[], Node* ndAr[], DesignVariable* dvAr[]) {
	DiffDoub v1[3];
	DiffDoub v2[3];
	DiffDoub tmpV[3];
	DiffDoub tmp;

	if (numNds == 4) {
		ndAr[globNodes[2]]->getCrd(v1, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[3]]->getCrd(v2, dvAr);
		ndAr[globNodes[1]]->getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}
	else {
		ndAr[globNodes[1]]->getCrd(v1, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[2]]->getCrd(v2, dvAr);
		ndAr[globNodes[0]]->getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}

	crossProd(norm, v1, v2);
	area.setVal(norm[0]);
	area.sqr();
	tmp.setVal(norm[1]);
	tmp.sqr();
	area.add(tmp);
	tmp.setVal(norm[2]);
	tmp.sqr();
	area.add(tmp);
	area.sqt();

	tmp.setVal(1.0);
	tmp.dvd(area);
	norm[0].mult(tmp);
	norm[1].mult(tmp);
	norm[2].mult(tmp);

	tmp.setVal(0.5);
	area.mult(tmp);

	return;
}

//end dup
 
//end skip 

// FacePt

FacePt::FacePt(Face* newFc) {
	fc = newFc;
	next = nullptr;
	return;
}

// Begin FaceList
FaceList::FaceList() {
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}

void FaceList::addFace(Face *newFc) {
	if(!first) {
		first = newFc;
		last = newFc;
	} else {
		last->setNext(newFc);
		last = newFc;
	}
	len++;
	return;
}

void FaceList::addIfAbsent(Face *newFc) {
	int i1;
	int newNumNds;
	int newSrtd[8];
	int thisNumNds;
	int thisSrtd[8];
	bool allMatch;
	Face *thisFc;
	
	newNumNds = newFc->getNumNds();
	newFc->sortedNodes(newSrtd);
	thisFc = first;
	while(thisFc) {
		thisNumNds = thisFc->getNumNds();
		if(thisNumNds == newNumNds) {
			thisFc->sortedNodes(thisSrtd);
			allMatch = true;
			for (i1 = 0; i1 < thisNumNds; i1++) {
				if(thisSrtd[i1] != newSrtd[i1]) {
					allMatch = false;
				}
			}
			if(allMatch) {
				newFc->setOnSurface(false);
				thisFc->setOnSurface(false);
				return;
			}
		}
		thisFc = thisFc->getNext();
	}
	
	addFace(newFc);
	
	return;
}

Face* FaceList::getFirst() {
	return first;
}

FaceList::~FaceList() {
	Face *thisFc;
	Face *nextFc;
	thisFc = first;
	while(thisFc) {
		nextFc = thisFc->getNext();
		delete thisFc;
		thisFc = nextFc;
	}
	return;
}

//FacePtList

FacePtList::FacePtList() {
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}

void FacePtList::addFace(FacePt* newFpt) {
	if (!first) {
		first = newFpt;
		last = newFpt;
	}
	else {
		last->next = newFpt;
		last = newFpt;
	}
	len++;
	return;
}

bool FacePtList::addIfAbsent(Face* newFc) {
	int i1;
	int newNumNds;
	int newSrtd[8];
	int thisNumNds;
	int thisSrtd[8];
	bool allMatch;
	FacePt* thisFc;

	newNumNds = newFc->getNumNds();
	newFc->sortedNodes(newSrtd);
	thisFc = first;
	while (thisFc) {
		thisNumNds = thisFc->fc->getNumNds();
		if (thisNumNds == newNumNds) {
			thisFc->fc->sortedNodes(thisSrtd);
			allMatch = true;
			for (i1 = 0; i1 < thisNumNds; i1++) {
				if (thisSrtd[i1] != newSrtd[i1]) {
					allMatch = false;
				}
			}
			if (allMatch) {
				newFc->setOnSurface(false);
				thisFc->fc->setOnSurface(false);
				return false;
			}
		}
		thisFc = thisFc->next;
	}

	FacePt* newFpt = new FacePt(newFc);
	addFace(newFpt);
	return true;
}

FacePtList::~FacePtList() {
	FacePt* thisFpt = first;
	FacePt* nextFpt;
	while (thisFpt) {
		nextFpt = thisFpt->next;
		delete thisFpt;
		thisFpt = nextFpt;
	}
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}