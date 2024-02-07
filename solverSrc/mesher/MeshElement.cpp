#include "MeshElement.h"
#include "MeshNode.h"
#include "utilities.h"

MeshElement::MeshElement() {
	nodes[0] = nullptr;
	nodes[1] = nullptr;
	nodes[2] = nullptr;
	nodes[3] = nullptr;
	next = nullptr;
	return;
}

void MeshElement::setNdPt(MeshNode* newNds[]) {
	int i1;
	for (i1 = 0; i1 < 4; i1++) {
		nodes[i1] = newNds[i1];
	}
	double mat[9];
	double det;
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	double* n4 = nodes[3]->getCrd();
	MeshNode* swap;
	mat[0] = n2[0] - n1[0];
	mat[3] = n2[1] - n1[1];
	mat[6] = n2[2] - n1[2];
	mat[1] = n3[0] - n1[0];
	mat[4] = n3[1] - n1[1];
	mat[7] = n3[2] - n1[2];
	mat[2] = n4[0] - n1[0];
	mat[5] = n4[1] - n1[1];
	mat[8] = n4[2] - n1[2];
	qRFactor(mat, 3, 0, 2, 0, 2, 0);
	det = mat[0] * mat[4] * mat[8];
	if (det < 0) {
		swap = nodes[1];
		nodes[1] = nodes[2];
		nodes[2] = swap;
	}
	return;
}

void MeshElement::setNext(MeshElement* newNext) {
	next = newNext;
	return;
}

MeshNode** MeshElement::getNodes() {
	return &nodes[0];
}

void MeshElement::getCentroid(double cent[]) {
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	double* n4 = nodes[3]->getCrd();
	cent[0] = 0.25 * (n1[0] + n2[0] + n3[0] + n4[0]);
	cent[1] = 0.25 * (n1[1] + n2[1] + n3[1] + n4[1]);
	cent[2] = 0.25 * (n1[2] + n2[2] + n3[2] + n4[2]);
	return;
}

double MeshElement::getVolume() {
	double mat[9];
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	double* n4 = nodes[3]->getCrd();
	mat[0] = n2[0] - n1[0];
	mat[3] = n2[1] - n1[1];
	mat[6] = n2[2] - n1[2];
	mat[1] = n3[0] - n1[0];
	mat[4] = n3[1] - n1[1];
	mat[7] = n3[2] - n1[2];
	mat[2] = n4[0] - n1[0];
	mat[5] = n4[1] - n1[1];
	mat[8] = n4[2] - n1[2];
	qRFactor(mat, 3, 0, 2, 0, 2, 0);
	double vol = mat[0] * mat[4] * mat[8];
	return vol;
}

bool MeshElement::pointIn(double pt[]) {
	double mat[9];
	double rhs[3];
	double soln[3];
	double solSum;
	double* n1 = nodes[0]->getCrd();
	double* n2 = nodes[1]->getCrd();
	double* n3 = nodes[2]->getCrd();
	double* n4 = nodes[3]->getCrd();
	mat[0] = n2[0] - n1[0];
	mat[3] = n2[1] - n1[1];
	mat[6] = n2[2] - n1[2];
	mat[1] = n3[0] - n1[0];
	mat[4] = n3[1] - n1[1];
	mat[7] = n3[2] - n1[2];
	mat[2] = n4[0] - n1[0];
	mat[5] = n4[1] - n1[1];
	mat[8] = n4[2] - n1[2];
	rhs[0] = pt[0] - n1[0];
	rhs[1] = pt[1] - n1[1];
	rhs[2] = pt[2] - n1[2];
	qRFactor(mat, 3, 0, 2, 0, 2, 0);
	solveqRxEqb(soln, mat, rhs, 3, 0, 2, 0, 2, 0);
	if (soln[0] > 0.0000000001 && soln[1] > 0.0000000001 && soln[2] > 0.0000000001) {
		solSum = soln[0] + soln[1] + soln[2];
		if (solSum < 0.9999999999) {
			return true;
		}
	}
	return false;
}

MeshElement* MeshElement::getNext() {
	return next;
}

MEList::MEList() {
	first = nullptr;
	last = nullptr;
	length = 0;
	return;
}

void MEList::addEnt(MeshElement* newEl) {
	if (!first) {
		first = newEl;
		last = newEl;
	}
	else {
		last->setNext(newEl);
		last = newEl;
	}
	length++;
	return;
}

MeshElement* MEList::getFirst() {
	return first;
}

int MEList::getLength() {
	return length;
}

MEList::~MEList() {
	MeshElement* thisEnt = first;
	MeshElement* nextEnt = nullptr;
	while (thisEnt) {
		nextEnt = thisEnt->getNext();
		delete thisEnt;
		thisEnt = nextEnt;
	}
	return;
}