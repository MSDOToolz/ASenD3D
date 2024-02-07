#include "Mesher.h"
#include "utilities.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

Mesher::Mesher() {
	gridOut1 = nullptr;
	gridOut2 = nullptr;
	globProjWt = 0.75;
	return;
}

void Mesher::readInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4];
	int hdLdSpace[4];
	string data[3];
	int dataLen;
	MeshNode* newNd;
	double ndCrd[3];
	MeshFace* newFc;
	int fcNd[3];

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "nodes" && dataLen == 3) {
				ndCrd[0] = stod(data[0]);
				ndCrd[1] = stod(data[1]);
				ndCrd[2] = stod(data[2]);
				newNd = new MeshNode(ndCrd);
				nodes.addEnt(newNd);
			}
			else if (headings[0] == "faces" && dataLen == 3) {
				fcNd[0] = stoi(data[0]);
				fcNd[1] = stoi(data[1]);
				fcNd[2] = stoi(data[2]);
				newFc = new MeshFace();
				newFc->setNodeLabs(fcNd);
				faces.addEnt(newFc);
			}
			else if (headings[0] == "globProjWt" && dataLen == 1) {
				globProjWt = stod(data[0]);
			}
		}
		inFile.close();
	}
	else {
		string erSt = "Error: Could not open input file " + fileName + " for unstructured 3D mesh generation.";
		throw invalid_argument(erSt);
	}
	return;
}

void Mesher::initBoundaryNormals() {
	int i1;
	int lstLen;
	double cent[3];
	double* normDir;
	int srchDir;
	double proj;
	double xRange[2];
	double yRange[2];
	double zRange[2];
	double vec[3];
	double intOut[3];
	double dp;
	bool intersects;
	int intCt;
	MeshFace* thisFc = faces.getFirst();
	MeshFace* thisFc2 = nullptr;
	while (thisFc) {
		thisFc->getCentroid(cent);
		normDir = thisFc->getNormDir();
		srchDir = 0;
		if (abs(normDir[1]) > abs(normDir[0])) {
			srchDir = 1;
		}
		if (abs(normDir[2]) > abs(normDir[1])) {
			srchDir = 2;
		}
		if (srchDir == 0) {
			vec[0] = 1.0;
			vec[1] = 0.0;
			vec[2] = 0.0;
			xRange[0] = 1;
			xRange[1] = 0;
			yRange[0] = cent[1] - avgProj;
			yRange[1] = cent[1] + avgProj;
			zRange[0] = cent[2] - avgProj;
			zRange[1] = cent[2] + avgProj;
		}
		else if (srchDir == 1) {
			vec[0] = 0.0;
			vec[1] = 1.0;
			vec[2] = 0.0;
			xRange[0] = cent[0] - avgProj;
			xRange[1] = cent[0] + avgProj;
			yRange[0] = 1;
			yRange[1] = 0;
			zRange[0] = cent[2] - avgProj;
			zRange[1] = cent[2] + avgProj;
		}
		else {
			vec[0] = 0.0;
			vec[1] = 0.0;
			vec[2] = 1.0;
			xRange[0] = cent[0] - avgProj;
			xRange[1] = cent[0] + avgProj;
			yRange[0] = cent[1] - avgProj;
			yRange[1] = cent[1] + avgProj;
			zRange[0] = 1;
			zRange[1] = 0;
		}
		lstLen = faceGrid.getInXYZRange(gridOut1, gOLen, xRange, yRange, zRange);
		intCt = 1;
		for (i1 = 0; i1 < lstLen; i1++) {
			thisFc2 = gridOut1[i1]->getPt(thisFc2);
			if (thisFc2 != thisFc) {
				intersects = thisFc2->getIntersection(intOut, cent, vec);
				if (intersects && intOut[0] > 0.0) {
					intCt *= -1;
				}
			}
		}
		// if inCt > 0, number of intersections is even, vec points out of the surface
		dp = vec[0] * normDir[0] + vec[1] * normDir[1] + vec[2] * normDir[2];
		if ((intCt > 0 && dp > 0.0) || (intCt < 0 && dp < 0.0)) {
			normDir[0] *= -1.0;
			normDir[1] *= -1.0;
			normDir[2] *= -1.0;
		}
		thisFc = thisFc->getNext();
	}
	return;
}

void Mesher::prep() {
	//Set Boundary face node pointers
	numBoundNds = nodes.getLength();
	MeshNode** ndAr = new MeshNode * [numBoundNds];
	MeshNode* thisNd = nodes.getFirst();
	int i1 = 0;
	while (thisNd) {
		ndAr[i1] = thisNd;
		thisNd = thisNd->getNext();
		i1++;
	}

	MeshFace* thisFc = faces.getFirst();
	while (thisFc) {
		thisFc->setPtFromLabs(ndAr);
		thisFc->initNormDir();
		thisFc = thisFc->getNext();
	}
	
	//Allocate the grid output list
	int numFaces = faces.getLength();
	if (numFaces > numBoundNds) {
		gOLen = 2.0 * numFaces;
	}
	else {
		gOLen = 2.0 * numBoundNds;
	}
	gridOut1 = new MeshEnt * [gOLen];
	gridOut2 = new MeshEnt * [gOLen];


	//Initialize grids
	double xRange[2] = { 1.0e+100,-1.0e+100 };
	double yRange[2] = { 1.0e+100,-1.0e+100 };
	double zRange[2] = { 1.0e+100,-1.0e+100 };
	double* crd;
	thisNd = nodes.getFirst();
	while (thisNd) {
		crd = thisNd->getCrd();
		if (crd[0] < xRange[0]) {
			xRange[0] = crd[0];
		}
		if (crd[0] > xRange[1]) {
			xRange[1] = crd[0];
		}
		if (crd[1] < yRange[0]) {
			yRange[0] = crd[1];
		}
		if (crd[1] > yRange[1]) {
			yRange[1] = crd[1];
		}
		if (crd[2] < zRange[0]) {
			zRange[0] = crd[2];
		}
		if (crd[2] > zRange[1]) {
			zRange[1] = crd[2];
		}
		thisNd = thisNd->getNext();
	}

	double spacing = 0;
	thisFc = faces.getFirst();
	while (thisFc) {
		spacing += thisFc->getProjDist();
		thisFc = thisFc->getNext();
	}
	spacing /= numFaces;
	avgProj = spacing;
	spacing *= 2.0;
	nodeGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);
	elementGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);
	faceGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);

	//Add boundary nodes and faces to their grids
	double cent[3];
	thisNd = nodes.getFirst();
	while (thisNd) {
		nodeGrid.addEnt(thisNd, thisNd->getCrd());
		thisNd = thisNd->getNext();
	}

	thisFc = faces.getFirst();
	while (thisFc) {
		thisFc->getCentroid(cent);
		faceGrid.addEnt(thisFc, cent);
		thisFc = thisFc->getNext();
	}

	// Check face normal directions

	initBoundaryNormals();

	delete[] ndAr;
	return;
}

bool Mesher::checkNewEl(MeshElement* newEl, MeshFace* newFaces[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int n1;
	int n2;
	int lstLen;
	double cent[3];
	int elEdges[12] = { 0,1,0,2,0,3,1,2,1,3,2,3 };
	int faceEdges[6] = { 0,1,0,2,1,2 };
	double* pt;
	double* pt2;
	double vec[3];
	double outP[3];
	bool intersects;
	MeshNode** elNds = newEl->getNodes();
	MeshNode** fcNds = nullptr;
	MeshFace* thisFc = nullptr;
	MeshElement** fcEls;
	MeshElement* thisEl = nullptr;
	MeshNode* thisNd = nullptr;
	MeshNode* fcNdOut[3];
	bool shared[3];

	newEl->getCentroid(cent);
	lstLen = faceGrid.getInRadius(gridOut2, gOLen, cent, 4.0 * avgProj);
	for (i1 = 0; i1 < lstLen; i1++) {
		thisFc = gridOut2[i1]->getPt(thisFc);
		// Any face of new element already exists and is closed

		for (i2 = 0; i2 < 4; i2++) {
			i3 = thisFc->getSharedNodes(fcNdOut, shared, newFaces[i2]);
			fcEls = thisFc->getElPt();
			if (i3 == 3 && fcEls[1]) {
				return false;
			}
		}

		// Any edge of new element intersects existing face
		i3 = 0;
		for (i2 = 0; i2 < 6; i2++) {
			n1 = elEdges[i3];
			n2 = elEdges[i3 + 1];
			pt = elNds[n1]->getCrd();
			pt2 = elNds[n2]->getCrd();
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			intersects = thisFc->getIntersection(outP, pt, vec);
			if (intersects && outP[0] < 0.9999999999 && outP[0] > 0.0000000001) {
				return false;
			}
			i3 += 2;
		}

		// Any edge of existing face intersects face of new element

		fcNds = thisFc->getNdPt();
		i3 = 0;
		for (i2 = 0; i2 < 3; i2++) {
			n1 = faceEdges[i3];
			n2 = faceEdges[i3 + 1];
			pt = fcNds[n1]->getCrd();
			pt2 = fcNds[n2]->getCrd();
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			for (i4 = 0; i4 < 4; i4++) {
				intersects = newFaces[i4]->getIntersection(outP, pt, vec);
				if (intersects && outP[0] < 0.9999999999 && outP[0] > 0.0000000001) {
					return false;
				}
			}
			i3 += 2;
		}
	}

	// Any node of new element in an existing element

	lstLen = elementGrid.getInRadius(gridOut2, gOLen, cent, 4.0 * avgProj);
	for (i1 = 0; i1 < lstLen; i1++) {
		thisEl = gridOut2[i1]->getPt(thisEl);
		for (i2 = 0; i2 < 4; i2++) {
			pt = elNds[i2]->getCrd();
			intersects = thisEl->pointIn(pt);
			if (intersects) {
				return false;
			}
		}
	}

	// Any existing node in the new element

	lstLen = nodeGrid.getInRadius(gridOut2, gOLen, cent, 4.0 * avgProj);
	for (i1 = 0; i1 < lstLen; i1++) {
		thisNd = gridOut2[i1]->getPt(thisNd);
		pt = thisNd->getCrd();
		intersects = newEl->pointIn(pt);
		if (intersects) {
			return false;
		}
	}

	return true;
}

bool Mesher::addFaceIfAbsent(MeshFace* newFace, MeshElement* newEl) {
	int i1;
	double cent[3];
	newFace->getCentroid(cent);
	MeshFace* thisFc = nullptr;
	MeshNode* thisNds[3];
	MeshElement** fcEls;
	MeshElement* newFcEls[2];
	bool shared[3];
	int numShared;

	int lstLen = faceGrid.getInRadius(gridOut2, gOLen, cent, avgProj);
	for (i1 = 0; i1 < lstLen; i1++) {
		thisFc = gridOut2[i1]->getPt(thisFc);
		numShared = thisFc->getSharedNodes(thisNds, shared, newFace);
		if (numShared == 3) {
			fcEls = thisFc->getElPt();
			newFcEls[0] = fcEls[0];
			newFcEls[1] = newEl;
			thisFc->setElPt(newFcEls);
			delete newFace;
			return false;
		}
	}

	faces.addEnt(newFace);
	faceGrid.addEnt(newFace, cent);

	return true;
}

bool Mesher::adoptConnectedNd(MeshFace* thisFc, double tgtPt[], double srchRad) {
	int i1;
	int i2;
	MeshFace* listFc = nullptr;
	MeshNode* faceNds[3];
	bool shared[3];
	int numShared;
	MeshNode* unShared;
	double* crd;
	double dVec[3];
	double dist;
	double dp;
	double cent[3];
	MeshElement* newEl;
	MeshNode* newElNds[4];
	MeshFace* newElFcs[4];
	MeshNode* newFcNds[3];
	MeshElement* newFcEls[2];
	bool elCheck;

	newEl = new MeshElement;
	newElFcs[0] = thisFc;
	newElFcs[1] = new MeshFace;
	newElFcs[2] = new MeshFace;
	newElFcs[3] = new MeshFace;
	newFcEls[0] = newEl;
	newFcEls[1] = nullptr;
	
	double* normDir = thisFc->getNormDir();
	MeshNode** thisFcNds = thisFc->getNdPt();
	MeshElement** thisFcEls = thisFc->getElPt();
	thisFc->getCentroid(cent);
	int lstLen = faceGrid.getInRadius(gridOut1, gOLen, tgtPt, 4.0*srchRad);
	for (i1 = 0; i1 < lstLen; i1++) {
		listFc = gridOut1[i1]->getPt(listFc);
		numShared = listFc->getSharedNodes(faceNds, shared, thisFc);
		if (numShared == 2) {
			unShared = nullptr;
			for (i2 = 0; i2 < 3; i2++) {
				if (!shared[i2]) {
					unShared = faceNds[i2];
				}
			}
			crd = unShared->getCrd();
			dVec[0] = crd[0] - tgtPt[0];
			dVec[1] = crd[1] - tgtPt[1];
			dVec[2] = crd[2] - tgtPt[2];
			dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
			if (dist < srchRad) {
				dVec[0] = crd[0] - cent[0];
				dVec[1] = crd[1] - cent[1];
				dVec[2] = crd[2] - cent[2];
				dp = dVec[0] * normDir[0] + dVec[1] * normDir[1] + dVec[2] * normDir[2];
				if (dp > 0) {
					newElNds[0] = thisFcNds[0];
					newElNds[1] = thisFcNds[1];
					newElNds[2] = thisFcNds[2];
					newElNds[3] = unShared;
					newEl->setNdPt(newElNds);

					newFcNds[0] = thisFcNds[0];
					newFcNds[1] = thisFcNds[1];
					newFcNds[2] = unShared;
					newElFcs[1]->setNodePt(newFcNds);
					newElFcs[1]->setElPt(newFcEls);
					newElFcs[1]->initNormDir();

					newFcNds[0] = thisFcNds[1];
					newFcNds[1] = thisFcNds[2];
					newFcNds[2] = unShared;
					newElFcs[2]->setNodePt(newFcNds);
					newElFcs[2]->setElPt(newFcEls);
					newElFcs[2]->initNormDir();

					newFcNds[0] = thisFcNds[2];
					newFcNds[1] = thisFcNds[0];
					newFcNds[2] = unShared;
					newElFcs[3]->setNodePt(newFcNds);
					newElFcs[3]->setElPt(newFcEls);
					newElFcs[3]->initNormDir();

					elCheck = checkNewEl(newEl, newElFcs);
					if (elCheck) {
						elements.addEnt(newEl);
						newEl->getCentroid(cent);
						elementGrid.addEnt(newEl, cent);

						newFcEls[0] = thisFcEls[0];
						newFcEls[1] = newEl;
						thisFc->setElPt(newFcEls);

						newElFcs[1]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[1], newEl);

						newElFcs[2]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[2], newEl);

						newElFcs[3]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[3], newEl);

						return true;
					}
				}
			}
		}
	}
	delete newElFcs[3];
	delete newElFcs[2];
	delete newElFcs[1];
	delete newEl;

	return false;
}

bool Mesher::adoptAnyNd(MeshFace* thisFc, double tgtPt[], double srchRad) {
	int i1;
	int i2;
	MeshNode* listNd = nullptr;
	MeshNode* faceNds[3];
	double* crd;
	double dVec[3];
	double dist;
	double dp;
	double cent[3];
	MeshElement* newEl;
	MeshNode* newElNds[4];
	MeshFace* newElFcs[4];
	MeshNode* newFcNds[3];
	MeshElement* newFcEls[2];
	bool elCheck;

	newEl = new MeshElement;
	newElFcs[0] = thisFc;
	newElFcs[1] = new MeshFace;
	newElFcs[2] = new MeshFace;
	newElFcs[3] = new MeshFace;
	newFcEls[0] = newEl;
	newFcEls[1] = nullptr;

	double* normDir = thisFc->getNormDir();
	MeshNode** thisFcNds = thisFc->getNdPt();
	MeshElement** thisFcEls = thisFc->getElPt();
	thisFc->getCentroid(cent);
	int lstLen = nodeGrid.getInRadius(gridOut1,gOLen,cent,srchRad);
	for (i1 = 0; i1 < lstLen; i1++) {
		listNd = gridOut1[i1]->getPt(listNd);
		if (listNd != thisFcNds[0] && listNd != thisFcNds[1] && listNd != thisFcNds[2]) {
			crd = listNd->getCrd();
			dVec[0] = crd[0] - tgtPt[0];
			dVec[1] = crd[1] - tgtPt[1];
			dVec[2] = crd[2] - tgtPt[2];
			dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
			if (dist < srchRad) {
				dVec[0] = crd[0] - cent[0];
				dVec[1] = crd[1] - cent[1];
				dVec[2] = crd[2] - cent[2];
				dp = dVec[0] * normDir[0] + dVec[1] * normDir[1] + dVec[2] * normDir[2];
				if (dp > 0) {
					newElNds[0] = thisFcNds[0];
					newElNds[1] = thisFcNds[1];
					newElNds[2] = thisFcNds[2];
					newElNds[3] = listNd;
					newEl->setNdPt(newElNds);

					newFcNds[0] = thisFcNds[0];
					newFcNds[1] = thisFcNds[1];
					newFcNds[2] = listNd;
					newElFcs[1]->setNodePt(newFcNds);
					newElFcs[1]->setElPt(newFcEls);
					newElFcs[1]->initNormDir();

					newFcNds[0] = thisFcNds[1];
					newFcNds[1] = thisFcNds[2];
					newFcNds[2] = listNd;
					newElFcs[2]->setNodePt(newFcNds);
					newElFcs[2]->setElPt(newFcEls);
					newElFcs[2]->initNormDir();

					newFcNds[0] = thisFcNds[2];
					newFcNds[1] = thisFcNds[0];
					newFcNds[2] = listNd;
					newElFcs[3]->setNodePt(newFcNds);
					newElFcs[3]->setElPt(newFcEls);
					newElFcs[3]->initNormDir();

					elCheck = checkNewEl(newEl, newElFcs);
					if (elCheck) {
						elements.addEnt(newEl);
						newEl->getCentroid(cent);
						elementGrid.addEnt(newEl, cent);

						newFcEls[0] = thisFcEls[0];
						newFcEls[1] = newEl;
						thisFc->setElPt(newFcEls);

						newElFcs[1]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[1], newEl);

						newElFcs[2]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[2], newEl);

						newElFcs[3]->normDirFromElCent(cent);
						addFaceIfAbsent(newElFcs[3], newEl);

						return true;
					}
				}
			}
		}
	}
	delete newElFcs[3];
	delete newElFcs[2];
	delete newElFcs[1];
	delete newEl;

	return false;
}

bool Mesher::createNewNd(MeshFace* thisFc, double tgtPt[]) {
	MeshNode* newNd = new MeshNode(tgtPt);
	MeshElement* newEl;
	MeshNode* newElNds[4];
	MeshFace* newElFcs[4];
	MeshNode* newFcNds[3];
	MeshElement* newFcEls[2];
	bool elCheck;
	double cent[3];

	newEl = new MeshElement;
	MeshNode** thisFcNds = thisFc->getNdPt();
	MeshElement** thisFcEls = thisFc->getElPt();
	newElNds[0] = thisFcNds[0];
	newElNds[1] = thisFcNds[1];
	newElNds[2] = thisFcNds[2];
	newElNds[3] = newNd;
	newEl->setNdPt(newElNds);
	newEl->getCentroid(cent);

	newElFcs[0] = thisFc;

	newFcEls[0] = newEl;
	newFcEls[1] = nullptr;

	newElFcs[1] = new MeshFace;
	newFcNds[0] = thisFcNds[0];
	newFcNds[1] = thisFcNds[1];
	newFcNds[2] = newNd;
	newElFcs[1]->setNodePt(newFcNds);
	newElFcs[1]->setElPt(newFcEls);
	newElFcs[1]->initNormDir();
	newElFcs[1]->normDirFromElCent(cent);

	newElFcs[2] = new MeshFace;
	newFcNds[0] = thisFcNds[1];
	newFcNds[1] = thisFcNds[2];
	newFcNds[2] = newNd;
	newElFcs[2]->setNodePt(newFcNds);
	newElFcs[2]->setElPt(newFcEls);
	newElFcs[2]->initNormDir();
	newElFcs[2]->normDirFromElCent(cent);

	newElFcs[3] = new MeshFace;
	newFcNds[0] = thisFcNds[2];
	newFcNds[1] = thisFcNds[0];
	newFcNds[2] = newNd;
	newElFcs[3]->setNodePt(newFcNds);
	newElFcs[3]->setElPt(newFcEls);
	newElFcs[3]->initNormDir();
	newElFcs[3]->normDirFromElCent(cent);

	elCheck = checkNewEl(newEl, newElFcs);
	if (elCheck) {
		nodes.addEnt(newNd);
		nodeGrid.addEnt(newNd, newNd->getCrd());

		elements.addEnt(newEl);
		newEl->getCentroid(cent);
		elementGrid.addEnt(newEl, cent);

		newFcEls[0] = thisFcEls[0];
		newFcEls[1] = newEl;
		thisFc->setElPt(newFcEls);

		addFaceIfAbsent(newElFcs[1], newEl);

		addFaceIfAbsent(newElFcs[2], newEl);

		addFaceIfAbsent(newElFcs[3], newEl);

		return true;
	}
	delete newElFcs[3];
	delete newElFcs[2];
	delete newElFcs[1];
	delete newEl;
	delete newNd;

	return false;
}

void Mesher::generateMesh() {
	MeshFace* thisFc;
	MeshElement** fcEls;
	double fcCent[3];
	double fcProj;
	double* fcNorm;
	double proj;
	double tgtPt[3];
	double srchRad;
	bool edgeClosed;
	
	bool elAdded = true;
	while (elAdded) {
		elAdded = false;
		thisFc = faces.getFirst();
		while (thisFc) {
			fcEls = thisFc->getElPt();
			if (!fcEls[1]) {
				thisFc->getCentroid(fcCent);
				fcProj = thisFc->getProjDist();
				proj = globProjWt * avgProj + (1.0 - globProjWt) * fcProj;
				fcNorm = thisFc->getNormDir();
				tgtPt[0] = fcCent[0] + 0.5 * proj * fcNorm[0];
				tgtPt[1] = fcCent[1] + 0.5 * proj * fcNorm[1];
				tgtPt[2] = fcCent[2] + 0.5 * proj * fcNorm[2];
				srchRad = 0.75 * proj;
				edgeClosed = adoptConnectedNd(thisFc, tgtPt, srchRad);
				if (!edgeClosed) {
					edgeClosed = adoptAnyNd(thisFc, tgtPt, srchRad);
				}
				if (!edgeClosed) {
					tgtPt[0] = fcCent[0] + proj * fcNorm[0];
					tgtPt[1] = fcCent[1] + proj * fcNorm[1];
					tgtPt[2] = fcCent[2] + proj * fcNorm[2];
					edgeClosed = createNewNd(thisFc, tgtPt);
				}
				if (edgeClosed) {
					elAdded = true;
				}
			}
			thisFc = thisFc->getNext();
		}
	}
	return;
}

void Mesher::distributeNodes() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int dim = nodes.getLength() * 3;
	MeshNode* thisNd;
	MeshElement* thisEl;
	MeshNode** elNds;
	double eVol;
	double avgWt;
	double* crd;
	double res;
	double alpha;
	double beta;
	double rNext;
	double dp;

	double* elMat = new double[144];
	double* xVec = new double[dim];
	double* hVec = new double[dim];
	double* zVec = new double[dim];
	double* gVec = new double[dim];
	double* wVec = new double[dim];
	double* dMat = new double[dim];
	double* pMat = new double[dim];
	double* pInv = new double[dim];
	
	i1 = elements.getLength();
	double* elWt = new double[i1];

	// Make static element matrix
	for (i1 = 0; i1 < 144; i1++) {
		elMat[i1] = 0.0;
	}

	i3 = 0;
	for (i1 = 0; i1 < 12; i1++) {
		i4 = i3 - 9;
		i5 = i1 * 12;
		i6 = (i1 + 1) * 12;
		while (i4 <= (i3 + 9)) {
			if (i4 >= i5 && i4 < i6) {
				elMat[i4] = -1.0;
			}
			i4 += 3;
		}
		elMat[i3] = 3.0;
		i3 += 13;
	}

	// Determine element weights
	i1 = 0;
	avgWt = 0.0;
	thisEl = elements.getFirst();
	while (thisEl) {
		eVol = thisEl->getVolume();
		elWt[i1] = eVol;
		avgWt += eVol;
		thisEl = thisEl->getNext();
		i1++;
	}
	avgWt /= i1;
	for (i2 = 0; i2 < i1; i2++) {
		elWt[i2] /= avgWt;
	}

	// Initialize matrices
	for (i1 = 0; i1 < dim; i1++) {
		dMat[i1] = 0.0;
		pMat[i1] = 30.0;
		gVec[i1] = 0.0;
	}

	thisNd = nodes.getFirst();
	i1 = 0;
	while (thisNd && i1 < numBoundNds) {
		i2 = thisNd->getLabel();
		crd = thisNd->getCrd();
		for (i3 = 0; i3 < 3; i3++) {
			i4 = i2 * 3 + i3;
			dMat[i4] = 100000.0;
			pMat[i4] += 100000.0;
			gVec[i4] = -100000.0 * crd[i3];
		}
		thisNd = thisNd->getNext();
		i1++;
	}

	for (i1 = 0; i1 < dim; i1++) {
		pInv[i1] = 1.0 / pMat[i1];
	}

	// Intialize Vectors

	res = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		xVec[i1] = 0.0;
		wVec[i1] = pInv[i1] * gVec[i1];
		hVec[i1] = -wVec[i1];
		res += wVec[i1] * gVec[i1];
	}

	i1 = 0;
	while (i1 < dim && res > 1.0e-12) {
		for (i2 = 0; i2 < dim; i2++) {
			zVec[i2] = 0.0;
		}
		thisEl = elements.getFirst();
		i2 = 0;
		while (thisEl) {
			elNds = thisEl->getNodes();
			i7 = 0; //Index in elMat
			for (i3 = 0; i3 < 4; i3++) {
				for (i4 = 0; i4 < 3; i4++) {
					i8 = 3*elNds[i3]->getLabel() + i4; // global row
					for (i5 = 0; i5 < 4; i5++) {
						for (i6 = 0; i6 < 3; i6++) {
							i9 = 3 * elNds[i5]->getLabel() + i6; //global col
							zVec[i8] += (elWt[i2] * elMat[i7] * hVec[i9]);
						}
					}
				}
			}
			thisEl = thisEl->getNext();
			i2++;
		}
		for (i2 = 0; i2 < dim; i2++) {
			zVec[i2] += (dMat[i2] * hVec[i2]);
		}
		dp = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			dp += (zVec[i2] * hVec[i2]);
		}
		alpha = res / dp;
		rNext = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			xVec[i2] += (alpha * hVec[i2]);
			gVec[i2] += (alpha * zVec[i2]);
			wVec[i2] = pInv[i2] * gVec[i2];
			rNext += (gVec[i2] * wVec[i2]);
		}
		beta = rNext / res;
		for (i2 = 0; i2 < dim; i2++) {
			hVec[i2] = -wVec[i2] + beta * hVec[i2];
		}
		res = rNext;
		i1++;
	}

	i1 = 0;
	thisNd = nodes.getFirst();
	while (thisNd) {
		if (i1 >= numBoundNds) {
			crd = thisNd->getCrd();
			i2 = thisNd->getLabel();
			for (i3 = 0; i3 < 3; i3++) {
				i4 = 3 * i2 + i3;
				crd[i3] = xVec[i4];
			}
		}
		thisNd = thisNd->getNext();
		i1++;
	}

	delete[] elMat;
	delete[] xVec;
	delete[] hVec;
	delete[] zVec;
	delete[] gVec;
	delete[] wVec;
	delete[] dMat;
	delete[] pMat;
	delete[] pInv;
	delete[] elWt;

	return;
}

void Mesher::writeOutput(string fileName) {
	int i1;
	ofstream outFile;
	MeshNode* thisNd;
	MeshElement* thisEl;
	MeshNode** elNds;
	double* crd;

	outFile.open(fileName);
	
	outFile << "nodes:\n";
	thisNd = nodes.getFirst();
	while (thisNd) {
		crd = thisNd->getCrd();
		outFile << "- [" << crd[0] << ", " << crd[1] << ", " << crd[2] << "]\n";
		thisNd = thisNd->getNext();
	}

	outFile << "elements:\n";
	thisEl = elements.getFirst();
	while (thisEl) {
		elNds = thisEl->getNodes();
		outFile << "- [" << elNds[0]->getLabel();
		for (i1 = 1; i1 < 4; i1++) {
			outFile << ", " << elNds[i1]->getLabel();
		}
		outFile << "]\n";
		thisEl = thisEl->getNext();
	}

	outFile.close();
	return;
}

void Mesher::printCurrentMesh() {
	MeshNode* thisNd;
	MeshElement* thisEl;
	MeshFace* thisFc;
	double* crd;
	MeshNode** elNds;
	MeshElement** fcEls;
	cout << "Current Mesh:" << endl;
	cout << "Nodes:" << endl;
	thisNd = nodes.getFirst();
	while (thisNd) {
		crd = thisNd->getCrd();
		cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
		thisNd = thisNd->getNext();
	}
	cout << "Elements:" << endl;
	thisEl = elements.getFirst();
	while (thisEl) {
		elNds = thisEl->getNodes();
		cout << elNds[0]->getLabel() << ", ";
		cout << elNds[1]->getLabel() << ", ";
		cout << elNds[2]->getLabel() << ", ";
		cout << elNds[3]->getLabel() << endl;
		thisEl = thisEl->getNext();
	}

	cout << "Faces:" << endl;
	thisFc = faces.getFirst();
	while (thisFc) {
		cout << "  nodes: ";
		elNds = thisFc->getNdPt();
		cout << elNds[0]->getLabel() << ", ";
		cout << elNds[1]->getLabel() << ", ";
		cout << elNds[2]->getLabel() << ", ";
		fcEls = thisFc->getElPt();
		cout << "element pt: ";
		cout << fcEls[0] << ", " << fcEls[1] << endl;
		thisFc = thisFc->getNext();
	}

	return;
}

Mesher::~Mesher() {
	if (gridOut1) {
		delete[] gridOut1;
		delete[] gridOut2;
	}
	return;
}