#include "SpatialGrid.h"

MeshEnt::MeshEnt() {
	nd = nullptr;
	el = nullptr;
	fc = nullptr;
	next = nullptr;
	return;
}

void MeshEnt::setPt(MeshNode* newNd) {
	nd = newNd;
	return;
}

void MeshEnt::setPt(MeshElement* newEl) {
	el = newEl;
	return;
}

void MeshEnt::setPt(MeshFace* newFc) {
	fc = newFc;
	return;
}

void MeshEnt::setNext(MeshEnt* newNext) {
	next = newNext;
}

MeshNode* MeshEnt::getPt(MeshNode* pt) {
	return nd;
}

MeshElement* MeshEnt::getPt(MeshElement* pt) {
	return el;
}

MeshFace* MeshEnt::getPt(MeshFace* pt) {
	return fc;
}

MeshEnt* MeshEnt::getNext() {
	return next;
}

EntityList::EntityList() {
	first = nullptr;
	last = nullptr;
	length = 0;
	return;
}

void EntityList::addEntry(MeshNode* newNd) {
	MeshEnt* newEnt = new MeshEnt;
	newEnt->setPt(newNd);
	if (!first) {
		first = newEnt;
		last = newEnt;
	}
	else {
		last->setNext(newEnt);
		last = newEnt;
	}
	length++;
}

void EntityList::addEntry(MeshElement* newEl) {
	MeshEnt* newEnt = new MeshEnt;
	newEnt->setPt(newEl);
	if (!first) {
		first = newEnt;
		last = newEnt;
	}
	else {
		last->setNext(newEnt);
		last = newEnt;
	}
	length++;
}

void EntityList::addEntry(MeshFace* newFc) {
	MeshEnt* newEnt = new MeshEnt;
	newEnt->setPt(newFc);
	if (!first) {
		first = newEnt;
		last = newEnt;
	}
	else {
		last->setNext(newEnt);
		last = newEnt;
	}
	length++;
}

int EntityList::copyToArray(MeshEnt* outLst[], int maxLen) {
	int lstLen = 0;
	MeshEnt* thisEnt = first;
	while (thisEnt && lstLen < maxLen) {
		outLst[lstLen] = thisEnt;
		thisEnt = thisEnt->getNext();
		lstLen++;
	}
	return lstLen;
}

EntityList::~EntityList() {
	MeshEnt* thisEnt = first;
	MeshEnt* nextEnt = nullptr;
	while (thisEnt) {
		nextEnt = thisEnt->getNext();
		delete thisEnt;
		thisEnt = nextEnt;
	}
	return;
}

SpatialGrid::SpatialGrid() {
	listAr = nullptr;
	return;
}

void SpatialGrid::initialize(double xRange[], double xSpacing, double yRange[], double ySpacing, double zRange[], double zSpacing) {
	xMin = xRange[0];
	xSp = xSpacing;
	yMin = yRange[0];
	ySp = ySpacing;
	zMin = zRange[0];
	zSp = zSpacing;
	xBins = (xRange[1] - xMin) / xSp + 1;
	yBins = (yRange[1] - yMin) / ySp + 1;
	zBins = (zRange[1] - zMin) / zSp + 1;
	int totBins = xBins * yBins * zBins;
	listAr = new EntityList[totBins];
	return;
}

void SpatialGrid::addEnt(MeshNode* newNd, double crd[]) {
	int xB = (crd[0] - xMin) / xSp;
	if (xB < 0) {
		xB = 0;
	}
	if (xB >= xBins) {
		xB = xBins - 1;
	}
	int yB = (crd[1] - yMin) / ySp;
	if (yB < 0) {
		yB = 0;
	}
	if (yB >= yBins) {
		yB = yBins - 1;
	}
	int zB = (crd[2] - zMin) / zSp;
	if (zB < 0) {
		zB = 0;
	}
	if (zB >= zBins) {
		zB = zBins - 1;
	}
	int ind = (zB*yBins + yB)*xBins + xB;
	listAr[ind].addEntry(newNd);
}

void SpatialGrid::addEnt(MeshElement* newEl, double crd[]) {
	int xB = (crd[0] - xMin) / xSp;
	if (xB < 0) {
		xB = 0;
	}
	if (xB >= xBins) {
		xB = xBins - 1;
	}
	int yB = (crd[1] - yMin) / ySp;
	if (yB < 0) {
		yB = 0;
	}
	if (yB >= yBins) {
		yB = yBins - 1;
	}
	int zB = (crd[2] - zMin) / zSp;
	if (zB < 0) {
		zB = 0;
	}
	if (zB >= zBins) {
		zB = zBins - 1;
	}
	int ind = (zB * yBins + yB) * xBins + xB;
	listAr[ind].addEntry(newEl);
}

void SpatialGrid::addEnt(MeshFace* newFc, double crd[]) {
	int xB = (crd[0] - xMin) / xSp;
	if (xB < 0) {
		xB = 0;
	}
	if (xB >= xBins) {
		xB = xBins - 1;
	}
	int yB = (crd[1] - yMin) / ySp;
	if (yB < 0) {
		yB = 0;
	}
	if (yB >= yBins) {
		yB = yBins - 1;
	}
	int zB = (crd[2] - zMin) / zSp;
	if (zB < 0) {
		zB = 0;
	}
	if (zB >= zBins) {
		zB = zBins - 1;
	}
	int ind = (zB * yBins + yB) * xBins + xB;
	listAr[ind].addEntry(newFc);
}

int SpatialGrid::getInXYZRange(MeshEnt* outLst[], int maxLen, double xRange[], double yRange[], double zRange[]) {
	int iMin;
	int iMax;
	int jMin;
	int jMax;
	int kMin;
	int kMax;

	if (xRange[0] < xRange[1]) {
		iMin = (xRange[0] - xMin) / xSp;
		if (iMin < 0) {
			iMin = 0;
		}
		iMax = (xRange[1] - xMin) / xSp;
		if (iMax >= xBins) {
			iMax = xBins - 1;
		}
	}
	else {
		iMin = 0;
		iMax = xBins - 1;
	}

	if (yRange[0] < yRange[1]) {
		jMin = (yRange[0] - yMin) / ySp;
		if (jMin < 0) {
			jMin = 0;
		}
		jMax = (yRange[1] - yMin) / ySp;
		if (jMax >= yBins) {
			jMax = yBins - 1;
		}
	}
	else {
		jMin = 0;
		jMax = yBins - 1;
	}

	if (zRange[0] < zRange[1]) {
		kMin = (zRange[0] - zMin) / zSp;
		if (kMin < 0) {
			kMin = 0;
		}
		kMax = (zRange[1] - zMin) / zSp;
		if (kMax >= zBins) {
			kMax = zBins - 1;
		}
	}
	else {
		kMin = 0;
		kMax = zBins - 1;
	}

	int i;
	int j;
	int k;
	int ind;
	int lstLen = 0;
	int numAdded;
	for (k = kMin; k <= kMax; k++) {
		for (j = jMin; j <= jMax; j++) {
			for (i = iMin; i <= iMax; i++) {
				ind = (k * yBins + j)*xBins + i;
				numAdded = listAr[ind].copyToArray(&outLst[lstLen], (maxLen - lstLen));
				lstLen += numAdded;
			}
		}
	}

	return lstLen;
}

int SpatialGrid::getInRadius(MeshEnt* outList[], int maxLen, double pt[], double rad) {
	double range[6];
	range[0] = pt[0] - rad;
	range[1] = pt[0] + rad;
	range[2] = pt[1] - rad;
	range[3] = pt[1] + rad;
	range[4] = pt[2] - rad;
	range[5] = pt[2] + rad;

	return getInXYZRange(outList, maxLen, &range[0], &range[2], &range[4]);
}

SpatialGrid::~SpatialGrid() {
	if (listAr) {
		delete[] listAr;
	}
	return;
}