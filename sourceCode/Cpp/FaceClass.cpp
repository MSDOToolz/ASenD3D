#include "FaceClass.h"

using namespace std;
Face::Face(int newNumNds) {
	numNds = newNumNds;
	locNodes = new int[numNds];
	globNodes = new int[numNds];
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

void Face::destroy() {
	delete[] locNodes;
	delete[] globNodes;
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

void FaceList::destroy() {
	Face *thisFc;
	Face *nextFc;
	thisFc = first;
	while(thisFc) {
		nextFc = thisFc->getNext();
		thisFc->destroy();
		delete thisFc;
		thisFc = nextFc;
	}
	return;
}