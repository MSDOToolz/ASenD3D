#include "MeshFace.h"
#include "MeshNode.h"
#include "MeshElement.h"
#include "utilities.h"
#include <cmath>

MeshFace::MeshFace() {
	nodes[0] = nullptr;
	nodes[1] = nullptr;
	nodes[2] = nullptr;
	elements[0] = nullptr;
	elements[1] = nullptr;
	next = nullptr;
	return;
}

void MeshFace::setNodeLabs(int newLabs[]) {
	nodeLabs[0] = newLabs[0];
	nodeLabs[1] = newLabs[1];
	nodeLabs[2] = newLabs[2];
	return;
}

void MeshFace::setNodePt(MeshNode* newNds[]) {
	nodes[0] = newNds[0];
	nodes[1] = newNds[1];
	nodes[2] = newNds[2];
	return;
}

void MeshFace::setPtFromLabs(MeshNode* ndAr[]) {
	nodes[0] = ndAr[nodeLabs[0]];
	nodes[1] = ndAr[nodeLabs[1]];
	nodes[2] = ndAr[nodeLabs[2]];
	return;
}

void MeshFace::setElPt(MeshElement* newEls[]) {
	elements[0] = newEls[0];
	elements[1] = newEls[1];
	return;
}

void MeshFace::initNormDir() {
	double v1[3];
	double v2[3];
	double cp[3];
	double mag;
	double area;
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		v1[i1] = n2[i1] - n1[i1];
		v2[i1] = n3[i1] - n1[i1];
	}
	crossProd(cp, v1, v2);
	mag = sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]);
	area = 0.5 * mag;
	mag = 1.0 / mag;
	normDir[0] = mag * cp[0];
	normDir[1] = mag * cp[1];
	normDir[2] = mag * cp[2];
	projDist = 1.2408064788027997 * sqrt(area); // (2^1.5/3^0.75)*sqrt(area)
	return;
}

void MeshFace::normDirFromElCent(double cent[]) {
	double fcCent[3];
	double dVec[3];
	getCentroid(fcCent);
	dVec[0] = fcCent[0] - cent[0];
	dVec[1] = fcCent[1] - cent[1];
	dVec[2] = fcCent[2] - cent[2];
	double dp = dVec[0] * normDir[0] + dVec[1] * normDir[1] + dVec[2] * normDir[2];
	if (dp < 0.0) {
		normDir[0] *= -1.0;
		normDir[1] *= -1.0;
		normDir[2] *= -1.0;
	}
	return;
}

void MeshFace::setNext(MeshFace* newNext) {
	next = newNext;
	return;
}

MeshNode** MeshFace::getNdPt() {
	return &nodes[0];
}

MeshElement** MeshFace::getElPt() {
	return &elements[0];
}

double* MeshFace::getNormDir() {
	return &normDir[0];
}

void MeshFace::getCentroid(double cent[]) {
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	cent[0] = 0.333333333333333 * (n1[0] + n2[0] + n3[0]);
	cent[1] = 0.333333333333333 * (n1[1] + n2[1] + n3[1]);
	cent[2] = 0.333333333333333 * (n1[2] + n2[2] + n3[2]);
	return;
}

double MeshFace::getProjDist() {
	return projDist; 
}

bool MeshFace::getIntersection(double outParam[], double pt[], double vec[]) {
	double mat[9];
	double rhs[3];
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	mat[0] = vec[0];
	mat[3] = vec[1];
	mat[6] = vec[2];
	mat[1] = n1[0] - n2[0];
	mat[4] = n1[1] - n2[1];
	mat[7] = n1[2] - n2[2];
	mat[2] = n1[0] - n3[0];
	mat[5] = n1[1] - n3[1];
	mat[8] = n1[2] - n3[2];
	qRFactor(mat, 3, 0, 2, 0, 2, 0);
	double det = mat[0] * mat[4] * mat[8];
	if (det == 0.0) {
		return false;
	}
	else {
		rhs[0] = n1[0] - pt[0];
		rhs[1] = n1[1] - pt[1];
		rhs[2] = n1[2] - pt[2];
		solveqRxEqb(outParam, mat, rhs, 3, 0, 2, 0, 2, 0);
		if (outParam[1] > 0.0000000001 && outParam[2] > 0.0000000001) {
			if ((outParam[1] + outParam[2]) < 0.9999999999) {
				return true;
			}
		}
		return false;
	}
}

int MeshFace::getSharedNodes(MeshNode* ndPts[], bool shared[], MeshFace* fc) {
	int i1;
	int i2;
	MeshNode** fcNds = fc->getNdPt();
	ndPts[0] = nodes[0];
	ndPts[1] = nodes[1];
	ndPts[2] = nodes[2];
	shared[0] = false;
	shared[1] = false;
	shared[2] = false;
	int numShared = 0;
	for (i1 = 0; i1 < 3; i1++) {
		for (i2 = 0; i2 < 3; i2++) {
			if (ndPts[i1] == fcNds[i2]) {
				shared[i1] = true;
				numShared++;
			}
		}
	}
	return numShared;
}

MeshFace* MeshFace::getNext() {
	return next;
}

MFList::MFList() {
	first = nullptr;
	last = nullptr;
	length = 0;
	return;
}

void MFList::addEnt(MeshFace* newFc) {
	if (!first) {
		first = newFc;
		last = newFc;
	}
	else {
		last->setNext(newFc);
		last = newFc;
	}
	length++;
	return;
}

MeshFace* MFList::getFirst() {
	return first;
}

int MFList::getLength() {
	return length;
}

MFList::~MFList() {
	MeshFace* thisEnt = first;
	MeshFace* nextEnt = nullptr;
	while (thisEnt) {
		nextEnt = thisEnt->getNext();
		delete thisEnt;
		thisEnt = nextEnt;
	}
	return;
}