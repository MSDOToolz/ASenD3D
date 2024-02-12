#include "MeshFace.h"
#include "MeshNode.h"
#include "MeshElement.h"
#include "utilities.h"
#include <iostream>
#include <cmath>

using namespace std;

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

double MeshFace::getLongestEdgeLen() {
	double dVec[3];
	double dist;
	double maxDist = 0.0;
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	dVec[0] = n2[0] - n1[0];
	dVec[1] = n2[1] - n1[1];
	dVec[2] = n2[2] - n1[2];
	dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
	if (dist > maxDist) {
		maxDist = dist;
	}
	dVec[0] = n3[0] - n2[0];
	dVec[1] = n3[1] - n2[1];
	dVec[2] = n3[2] - n2[2];
	dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
	if (dist > maxDist) {
		maxDist = dist;
	}
	dVec[0] = n1[0] - n3[0];
	dVec[1] = n1[1] - n3[1];
	dVec[2] = n1[2] - n3[2];
	dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
	if (dist > maxDist) {
		maxDist = dist;
	}
	return maxDist;
}

bool MeshFace::getIntersection(double outParam[], double pt[], double vec[]) {
	int i1;
	double mat[9];
	double matMag;
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
	matMag = 0.0;
	for (i1 = 0; i1 < 9; i1++) {
		matMag += (mat[i1] * mat[i1]);
	}
	matMag = sqrt(matMag);
	qRFactor(mat, 3, 0, 2, 0, 2, 0);
	double det = mat[0] * mat[4] * mat[8];
	if (abs(det) < 1.0e-6*matMag*matMag*matMag) {
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

bool MeshFace::edgesIntersect(MeshFace* fc, double distTol) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double mat[6];
	double rhs[3];
	double soln[2];
	double det;
	double dVec[3];
	double dist;
	double v1[3];
	double v2[3];
	double matMag;
	double* n11;
	double* n21;
	double* n12;
	double* n22;
	MeshNode** fcNodes = fc->getNdPt();
	
	for (i1 = 0; i1 < 3; i1++) {
		n11 = nodes[i1]->getCrd();
		i3 = i1 + 1;
		if (i3 > 2) {
			i3 = 0;
		}
		n21 = nodes[i3]->getCrd();
		for (i2 = 0; i2 < 3; i2++) {
			n12 = fcNodes[i2]->getCrd();
			i4 = i2 + 1;
			if (i4 > 2) {
				i4 = 0;
			}
			n22 = fcNodes[i4]->getCrd();
			v1[0] = n21[0] - n11[0];
			v1[1] = n21[1] - n11[1];
			v1[2] = n21[2] - n11[2];
			v2[0] = n22[0] - n12[0];
			v2[1] = n22[1] - n12[1];
			v2[2] = n22[2] - n12[2];
			mat[0] = v1[0];
			mat[1] = -v2[0];
			mat[2] = v1[1];
			mat[3] = -v2[1];
			mat[4] = v1[2];
			mat[5] = -v2[2];
			rhs[0] = n12[0] - n11[0];
			rhs[1] = n12[1] - n11[1];
			rhs[2] = n12[2] - n11[2];
			matMag = 0.0;
			for (i5 = 0; i5 < 6; i5++) {
				matMag += (mat[i5] * mat[i5]);
			}
			matMag = sqrt(matMag);
			qRFactor(mat, 2, 0, 2, 0, 1, 0);
			det = mat[0] * mat[3];
			if (abs(det) > 1.0e-6*matMag*matMag) {
				solveqRxEqb(soln, mat, rhs, 2, 0, 2, 0, 1, 0);
				if (soln[0] > 0.0000000001 && soln[0] < 0.9999999999 && soln[1] > 0.0000000001 && soln[1] < 0.9999999999) {
					dVec[0] = n11[0] + soln[0] * v1[0] - n12[0] - soln[1] * v2[0];
					dVec[1] = n11[1] + soln[0] * v1[1] - n12[1] - soln[1] * v2[1];
					dVec[2] = n11[2] + soln[0] * v1[2] - n12[2] - soln[1] * v2[2];
					dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
					if (dist < distTol) {
						return true;
					}
				}
			}
		}
	}
	return false;
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

void MeshFace::printInfo() {
	int i1;
	double* crd;
	cout << "nodes:" << endl;
	for (i1 = 0; i1 < 3; i1++) {
		crd = nodes[i1]->getCrd();
		cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
	}
	return;
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