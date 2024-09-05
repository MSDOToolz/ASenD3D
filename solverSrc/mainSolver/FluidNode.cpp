#include "FluidNode.h"

using namespace std;

FluidNode::FluidNode(int newLab) {
	label = newLab;
	nextNd = nullptr;
	return;
}

void FluidNode::setCrd(double newCrd[]) {
	coord[0] = newCrd[0];
	coord[1] = newCrd[1];
	coord[2] = newCrd[2];
	return;
}

void FluidNode::setDofIndex(int dof, int index) {
	dofIndex[dof] = index;
	return;
}

void FluidNode::setSortedRank(int newRank) {
	sortedRank = newRank;
	return;
}