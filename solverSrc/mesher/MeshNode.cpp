#include "MeshNode.h"

MeshNode::MeshNode() {
	next = nullptr;
	return;
}

void MeshNode::setLabel(int newLab) {
	label = newLab;
	return;
}

void MeshNode::setCrd(double newCrd[]) {
	coord[0] = newCrd[0];
	coord[1] = newCrd[1];
	coord[2] = newCrd[2];
	return;
}

void MeshNode::setNext(MeshNode* newNd) {
	next = newNd;
	return;
}

int MeshNode::getLabel() {
	return label;
}

double* MeshNode::getCrd() {
	return &coord[0];
}

MeshNode* MeshNode::getNext() {
	return next;
}

MNList::MNList() {
	first = nullptr;
	last = nullptr;
	length = 0;
	return;
}

void MNList::addEnt(MeshNode* newNd) {
	if (!first) {
		first = newNd;
		last = newNd;
	}
	else {
		last->setNext(newNd);
		last = newNd;
	}
	newNd->setLabel(length);
	length++;
	return;
}

MeshNode* MNList::getFirst() {
	return first;
}

int MNList::getLength() {
	return length;
}

MNList::~MNList() {
	MeshNode* thisEnt = first;
	MeshNode* nextEnt = nullptr;
	while (thisEnt) {
		nextEnt = thisEnt->getNext();
		delete thisEnt;
		thisEnt = nextEnt;
	}
	return;
}