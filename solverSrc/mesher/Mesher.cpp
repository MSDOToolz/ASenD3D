#include "Mesher.h"
#include "constants.h"
#include "utilities.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

Mesher::Mesher() {
	nd_ct = 0;
	fc_ct = 0;
	el_ct = 0;
	gridOut1.clear();
	gridOut2.clear();
	globProjWt = 0.75;
	maxNumEls = max_int;
	newEl = max_int;
	newElFcs = vector<MeshFace>(4);
	newNd = max_int;
	return;
}

void Mesher::readInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4];
	int hdLdSpace[4];
	string data[3];
	int dataLen;
	int newNd;
	double ndCrd[3];
	int newFc;
	int fcNd[3];

	nd_ct = 0;
	fc_ct = 0;

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			getline(inFile, fileLine);
			readInputLine(fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "nodes" && dataLen == 3) {
				nd_ct++;
			}
			else if (headings[0] == "faces" && dataLen == 3) {
				fc_ct++;
			}
			else if (headings[0] == "globProjWt" && dataLen == 1) {
				globProjWt = stod(data[0]);
			}
			else if (headings[0] == "maxNumEls" && dataLen == 1) {
				maxNumEls = stoi(data[0]);
			}
		}
		inFile.close();
	}
	else {
		string erSt = "Error: Could not open input file " + fileName + " for unstructured 3D mesh generation.";
		throw invalid_argument(erSt);
	}

	int rad = 1;
	while (12*rad*rad < nd_ct) {
		rad++;
	}

	nd_cap = 15 * rad * rad * rad;
	el_cap = 6 * nd_cap;
	if (maxNumEls == max_int) {
		maxNumEls = el_cap;
	}
	else if (el_cap > maxNumEls) {
		el_cap = maxNumEls;
	}
	fc_cap = 2 * el_cap;

	nodes = vector<MeshNode>(nd_cap);
	elements = vector<MeshElement>(el_cap);
	faces = vector<MeshFace>(fc_cap);

	nd_ct = 0;
	fc_ct = 0;

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			getline(inFile, fileLine);
			readInputLine(fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "nodes" && dataLen == 3) {
				nodes[nd_ct].coord[0] = stod(data[0]);
				nodes[nd_ct].coord[1] = stod(data[1]);
				nodes[nd_ct].coord[2] = stod(data[2]);
				nd_ct++;
			}
			else if (headings[0] == "faces" && dataLen == 3) {
				faces[fc_ct].nodes[0] = stoi(data[0]);
				faces[fc_ct].nodes[1] = stoi(data[1]);
				faces[fc_ct].nodes[2] = stoi(data[2]);
				fc_ct++;
			}
		}
		inFile.close();
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
	int fci;
	int fci2;
	for (fci = 0; fci < fc_ct; fci++) {
		MeshFace& thisFc = faces[fci];
		thisFc.getCentroid(cent,nodes);
		normDir = &thisFc.normDir[0];
		srchDir = 0;
		if (abs(normDir[1]) > abs(normDir[0])) {
			srchDir = 1;
		}
		if (abs(normDir[2]) > abs(normDir[1])) {
			srchDir = 2;
		}
		if (srchDir == 0) {
			vec[0] = 1.0;
			vec[1] = 7.52893558402e-4;
			vec[2] = 5.39890329009e-4;
			xRange[0] = 1;
			xRange[1] = 0;
			yRange[0] = cent[1] - avgProj;
			yRange[1] = cent[1] + avgProj;
			zRange[0] = cent[2] - avgProj;
			zRange[1] = cent[2] + avgProj;
		}
		else if (srchDir == 1) {
			vec[0] = 7.52893558402e-4;
			vec[1] = 1.0;
			vec[2] = 5.39890329009e-4;
			xRange[0] = cent[0] - avgProj;
			xRange[1] = cent[0] + avgProj;
			yRange[0] = 1;
			yRange[1] = 0;
			zRange[0] = cent[2] - avgProj;
			zRange[1] = cent[2] + avgProj;
		}
		else {
			vec[0] = 7.52893558402e-4;
			vec[1] = 5.39890329009e-4;
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
			fci2 = gridOut1[i1];
			MeshFace& thisFc2 = faces[fci2];
			if (fci != fci2) {
				intersects = thisFc2.getIntersection(intOut, cent, vec, nodes);
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
		//
		if (cent[0] > 9.99) {
			dp = dp;
		}
		//
	}
	return;
}

void Mesher::prep() {
	int i1;
	int i2;
	numBoundNds = nd_ct;

	//Allocate the grid output list
	int numFaces = fc_ct;
	if (numFaces > numBoundNds) {
		gOLen = 2 * numFaces;
	}
	else {
		gOLen = 2 * numBoundNds;
	}
	gridOut1 = vector<int>(gOLen);
	gridOut2 = vector<int>(gOLen);

	//Initialize grids
	double xRange[2] = { 1.0e+100,-1.0e+100 };
	double yRange[2] = { 1.0e+100,-1.0e+100 };
	double zRange[2] = { 1.0e+100,-1.0e+100 };
	double* crd;
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& thisNd = nodes[i1];
		crd = &thisNd.coord[0];
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
	}

	double spacing = 0.0;
	avgProj = 0.0;
	maxProj = 0.0;
	maxEdgeLen = 0.0;
	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& thisFc = faces[i1];
		spacing = thisFc.projDist;
		avgProj += spacing;
		if (spacing > maxProj) {
			maxProj = spacing;
		}
		spacing = thisFc.getLongestEdgeLen(nodes);
		if (spacing > maxEdgeLen) {
			maxEdgeLen = spacing;
		}
	}
	avgProj /= numFaces;
	spacing = avgProj;
	nodeGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);
	elementGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);
	faceGrid.initialize(xRange, spacing, yRange, spacing, zRange, spacing);

	//Add boundary nodes and faces to their grids
	double cent[3];
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& thisNd = nodes[i1];
		nodeGrid.addEnt(i1, thisNd.coord);
	}

	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& thisFc = faces[i1];
		thisFc.getCentroid(cent,nodes);
		faceGrid.addEnt(i1, cent);
	}

	// Check boundary completeness

	int lstLen;
	int lstNd[3];
	bool shared[3];
	int numShared;
	int neighbCt;
	for (i2 = 0; i2 < fc_ct; i2++) {
		MeshFace& thisFc = faces[i2];
		thisFc.getCentroid(cent,nodes);
		lstLen = faceGrid.getInRadius(gridOut1, gOLen, cent, 1.01*maxEdgeLen);
		neighbCt = 0;
		for (i1 = 0; i1 < lstLen; i1++) {
			MeshFace& lstFc = faces[gridOut1[i1]];
			//fcNds = lstFc->getNdPt();
			//cout << fcNds[0]->getLabel() << ", " << fcNds[1]->getLabel() << ", " << fcNds[2]->getLabel() << endl;
			numShared = lstFc.getSharedNodes(lstNd, shared, i2, nodes, faces);
			if (numShared == 2) {
				neighbCt++;
			}
		}
		if (neighbCt != 3) {
			//cout << "neighbCt " << neighbCt;
			//fcNds = thisFc->getNdPt();
			//cout << fcNds[0]->getLabel() << ", " << fcNds[1]->getLabel() << ", " << fcNds[2]->getLabel() << endl;
			string erSt = "Error: Invalid boundary surface input to 3D unstructured mesh generator.\n";
			erSt = erSt + "Boundary must be a completely closed surface of triangular faces.";
			throw invalid_argument(erSt);
		}
	}

	// Check face normal directions

	initBoundaryNormals();

	return;
}

bool Mesher::checkNewEl(MeshElement& newEl, std::vector<MeshFace>& newFaces) {
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
	int* elNds = &newEl.nodes[0];
	int* fcNds = nullptr;
	int* fcEls = nullptr;
	int fcNdOut[3];
	bool shared[3];

	newEl.getCentroid(cent,nodes);
	lstLen = faceGrid.getInRadius(gridOut2, gOLen, cent, 1.01 * maxEdgeLen);
	for (i1 = 0; i1 < lstLen; i1++) {
		MeshFace& thisFc = faces[gridOut2[i1]];

		for (i2 = 0; i2 < 4; i2++) {
			// Any face of new element already exists and is closed
			i3 = thisFc.getSharedNodes(fcNdOut, shared, i2, nodes, newFaces);
			fcEls = &thisFc.elements[0];
			if (i3 == 3 && fcEls[1]) {
				return false;
			}
			// Any edge of new element intersects existing edge
			intersects = thisFc.edgesIntersect(i2, 1.0e-6 * avgProj, nodes, newFaces);
			if (intersects) {
				return false;
			}
		}

		// Any edge of new element intersects existing face
		i3 = 0;
		for (i2 = 0; i2 < 6; i2++) {
			n1 = elEdges[i3];
			n2 = elEdges[i3 + 1];
			pt = &nodes[elNds[n1]].coord[0];
			pt2 = &nodes[elNds[n2]].coord[0];
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			intersects = thisFc.getIntersection(outP, pt, vec, nodes);
			if (intersects && outP[0] < 0.9999999999 && outP[0] > 0.0000000001) {
				return false;
			}
			i3 += 2;
		}

		// Any edge of existing face intersects face of new element

		fcNds = &thisFc.nodes[0];
		i3 = 0;
		for (i2 = 0; i2 < 3; i2++) {
			n1 = faceEdges[i3];
			n2 = faceEdges[i3 + 1];
			pt = &nodes[fcNds[n1]].coord[0];
			pt2 = &nodes[fcNds[n2]].coord[0];
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			for (i4 = 0; i4 < 4; i4++) {
				intersects = newFaces[i4].getIntersection(outP, pt, vec, nodes);
				if (intersects && outP[0] < 0.9999999999 && outP[0] > 0.0000000001) {
					return false;
				}
			}
			i3 += 2;
		}
	}

	// Any node of new element in an existing element

	lstLen = elementGrid.getInRadius(gridOut2, gOLen, cent, 1.01 * maxEdgeLen);
	for (i1 = 0; i1 < lstLen; i1++) {
		MeshElement& thisEl = elements[gridOut2[i1]];
		for (i2 = 0; i2 < 4; i2++) {
			pt = &nodes[elNds[i2]].coord[0];
			intersects = thisEl.pointIn(pt,nodes);
			if (intersects) {
				return false;
			}
		}
	}

	// Any existing node in the new element

	lstLen = nodeGrid.getInRadius(gridOut2, gOLen, cent, 1.01 * maxEdgeLen);
	for (i1 = 0; i1 < lstLen; i1++) {
		MeshNode& thisNd = nodes[gridOut2[i1]];
		pt = &thisNd.coord[0];
		intersects = newEl.pointIn(pt, nodes);
		if (intersects) {
			return false;
		}
	}

	return true;
}

bool Mesher::addFaceIfAbsent(int newEl) {
	int i1;
	double cent[3];
	MeshFace& newFace = faces[fc_ct];
	newFace.getCentroid(cent, nodes);
	MeshFace* thisFc = nullptr;
	int thisNds[3];
	int* fcEls;
	int newFcEls[2];
	bool shared[3];
	int numShared;

	int lstLen = faceGrid.getInRadius(gridOut2, gOLen, cent, 0.1*avgProj);
	for (i1 = 0; i1 < lstLen; i1++) {
		MeshFace& thisFc = faces[gridOut2[i1]];
		numShared = thisFc.getSharedNodes(thisNds, shared, fc_ct, nodes, faces);
		if (numShared == 3) {
			thisFc.elements[1] = newEl;
			return false;
		}
	}

	faceGrid.addEnt(fc_ct, cent);
	fc_ct++;

	return true;
}

bool Mesher::adoptConnectedNd(int fc_i, double tgtPt[], double srchRad) {
	int i1;
	int i2;
	MeshFace& thisFc = faces[fc_i];
	int faceNds[3];
	bool shared[3];
	int numShared;
	int unShared;
	double* crd;
	double dVec[3];
	double dist;
	double dp;
	double cent[3];
	int* newElNds = &elements[el_ct].nodes[0];
	bool elCheck;
	bool fcAdded;
	MeshFace* newFc = nullptr;

	newElFcs[0].copy_data(thisFc);
	
	double* normDir = &thisFc.normDir[0];
	int* thisFcNds = &thisFc.nodes[0];
	int* thisFcEls = &thisFc.elements[0];
	thisFc.getCentroid(cent, nodes);
	int lstLen = faceGrid.getInRadius(gridOut1, gOLen, tgtPt, 1.01 * maxEdgeLen);
	for (i1 = 0; i1 < lstLen; i1++) {
		MeshFace& listFc = faces[gridOut1[i1]];
		numShared = listFc.getSharedNodes(faceNds, shared, fc_i, nodes, faces);
		if (numShared == 2) {
			unShared = max_int;
			for (i2 = 0; i2 < 3; i2++) {
				if (!shared[i2]) {
					unShared = faceNds[i2];
				}
			}
			crd = &nodes[unShared].coord[0];
			dVec[0] = crd[0] - tgtPt[0];
			dVec[1] = crd[1] - tgtPt[1];
			dVec[2] = crd[2] - tgtPt[2];
			dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
			if (dist < srchRad) {
				dVec[0] = crd[0] - cent[0];
				dVec[1] = crd[1] - cent[1];
				dVec[2] = crd[2] - cent[2];
				dp = dVec[0] * normDir[0] + dVec[1] * normDir[1] + dVec[2] * normDir[2];
				if (dp > 1.0e-3*avgProj) {
					newElNds[0] = thisFcNds[0];
					newElNds[1] = thisFcNds[1];
					newElNds[2] = thisFcNds[2];
					newElNds[3] = unShared;

					newFc = &newElFcs[1];
					newFc->nodes[0] = thisFcNds[0];
					newFc->nodes[1] = thisFcNds[1];
					newFc->nodes[2] = unShared;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					newFc = &newElFcs[2];
					newFc->nodes[0] = thisFcNds[1];
					newFc->nodes[1] = thisFcNds[2];
					newFc->nodes[2] = unShared;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					newFc = &newElFcs[3];
					newFc->nodes[0] = thisFcNds[2];
					newFc->nodes[1] = thisFcNds[0];
					newFc->nodes[2] = unShared;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					elCheck = checkNewEl(elements[el_ct], newElFcs);
					if (elCheck) {
						elements[el_ct].getCentroid(cent, nodes);
						elementGrid.addEnt(el_ct, cent);

						thisFc.elements[1] = el_ct;

						newElFcs[1].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[1]);
						fcAdded = addFaceIfAbsent(el_ct);

						newElFcs[2].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[2]);
						fcAdded = addFaceIfAbsent(el_ct);

						newElFcs[3].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[3]);
						fcAdded = addFaceIfAbsent(el_ct);

						el_ct++;
						return true;
					}
				}
			}
		}
	}

	return false;
}

bool Mesher::adoptAnyNd(int fc_i, double tgtPt[], double srchRad) {
	int i1;
	int i2;
	MeshFace& thisFc = faces[fc_i];
	MeshFace* newFc = nullptr;
	MeshNode* listNd = nullptr;
	int ndi;
	MeshNode* faceNds[3];
	double* crd;
	double dVec[3];
	double dist;
	double dp;
	double cent[3];
	int* newElNds = &elements[el_ct].nodes[0];
	bool elCheck;
	bool fcAdded;

	newElFcs[0] = thisFc;

	double* normDir = &thisFc.normDir[0];
	int* thisFcNds = &thisFc.nodes[0];
	int* thisFcEls = &thisFc.elements[0];
	thisFc.getCentroid(cent, nodes);
	int lstLen = nodeGrid.getInRadius(gridOut1,gOLen,tgtPt,srchRad);
	for (i1 = 0; i1 < lstLen; i1++) {
		//listNd = gridOut1[i1]->getPt(listNd);
		ndi = gridOut1[i1];
		MeshNode& listNd = nodes[ndi];
		if (ndi != thisFcNds[0] && ndi != thisFcNds[1] && ndi != thisFcNds[2]) {
			crd = &listNd.coord[0];
			dVec[0] = crd[0] - tgtPt[0];
			dVec[1] = crd[1] - tgtPt[1];
			dVec[2] = crd[2] - tgtPt[2];
			dist = sqrt(dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2]);
			if (dist < srchRad) {
				dVec[0] = crd[0] - cent[0];
				dVec[1] = crd[1] - cent[1];
				dVec[2] = crd[2] - cent[2];
				dp = dVec[0] * normDir[0] + dVec[1] * normDir[1] + dVec[2] * normDir[2];
				if (dp > 1.0e-3*avgProj) {
					newElNds[0] = thisFcNds[0];
					newElNds[1] = thisFcNds[1];
					newElNds[2] = thisFcNds[2];
					newElNds[3] = ndi;

					newFc = &newElFcs[1];
					newFc->nodes[0] = thisFcNds[0];
					newFc->nodes[1] = thisFcNds[1];
					newFc->nodes[2] = ndi;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					newFc = &newElFcs[2];
					newFc->nodes[0] = thisFcNds[1];
					newFc->nodes[1] = thisFcNds[2];
					newFc->nodes[2] = ndi;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					newFc = &newElFcs[3];
					newFc->nodes[0] = thisFcNds[2];
					newFc->nodes[1] = thisFcNds[0];
					newFc->nodes[2] = ndi;
					newFc->elements[0] = el_ct;
					newFc->elements[1] = max_int;
					newFc->initNormDir(nodes);

					elCheck = checkNewEl(elements[el_ct], newElFcs);
					if (elCheck) {
						//elements.addEnt(newEl);
						elements[el_ct].getCentroid(cent, nodes);
						elementGrid.addEnt(el_ct, cent);

						thisFc.elements[1] = el_ct;

						newElFcs[1].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[1]);
						fcAdded = addFaceIfAbsent(el_ct);

						newElFcs[2].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[2]);
						fcAdded = addFaceIfAbsent(el_ct);

						newElFcs[3].normDirFromElCent(cent, nodes);
						faces[fc_ct].copy_data(newElFcs[3]);
						fcAdded = addFaceIfAbsent(el_ct);

						el_ct++;

						return true;
					}
				}
			}
		}
	}

	return false;
}

bool Mesher::createNewNd(int fc_i, double tgtPt[]) {
	nodes[nd_ct].coord[0] = tgtPt[0];
	nodes[nd_ct].coord[1] = tgtPt[1];
	nodes[nd_ct].coord[2] = tgtPt[2];
	MeshFace& thisFc = faces[fc_i];
	bool elCheck;
	bool fcAdded;
	double cent[3];

	int* thisFcNds = &thisFc.nodes[0];
	int* thisFcEls = &thisFc.elements[0];
	int* newElNds = &elements[el_ct].nodes[0];
	newElNds[0] = thisFcNds[0];
	newElNds[1] = thisFcNds[1];
	newElNds[2] = thisFcNds[2];
	newElNds[3] = nd_ct;
	elements[el_ct].getCentroid(cent, nodes);

	newElFcs[0] = thisFc;

	MeshFace* newFc = &newElFcs[1];
	newFc->nodes[0] = thisFcNds[0];
	newFc->nodes[1] = thisFcNds[1];
	newFc->nodes[2] = nd_ct;
	newFc->elements[0] = el_ct;
	newFc->elements[1] = max_int;
	newFc->initNormDir(nodes);
	newFc->normDirFromElCent(cent, nodes);

	newFc = &newElFcs[2];
	newFc->nodes[0] = thisFcNds[1];
	newFc->nodes[1] = thisFcNds[2];
	newFc->nodes[2] = nd_ct;
	newFc->elements[0] = el_ct;
	newFc->elements[1] = max_int;
	newFc->initNormDir(nodes);
	newFc->normDirFromElCent(cent, nodes);

	newFc = &newElFcs[3];
	newFc->nodes[0] = thisFcNds[2];
	newFc->nodes[1] = thisFcNds[0];
	newFc->nodes[2] = nd_ct;
	newFc->elements[0] = el_ct;
	newFc->elements[1] = max_int;
	newFc->initNormDir(nodes);
	newFc->normDirFromElCent(cent, nodes);

	elCheck = checkNewEl(elements[el_ct], newElFcs);
	if (elCheck) {
		nodeGrid.addEnt(nd_ct, tgtPt);

		elementGrid.addEnt(el_ct, cent);

		thisFc.elements[1] = el_ct;

		faces[fc_ct].copy_data(newElFcs[1]);
		fcAdded = addFaceIfAbsent(el_ct);

		faces[fc_ct].copy_data(newElFcs[2]);
		fcAdded = addFaceIfAbsent(el_ct);

		faces[fc_ct].copy_data(newElFcs[3]);
		fcAdded = addFaceIfAbsent(el_ct);

		nd_ct++;
		el_ct++;

		return true;
	}

	return false;
}

bool Mesher::generateMesh() {
	int fc_i;
	int* fcEls;
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
		fc_i = 0;
		while (fc_i < fc_ct && el_ct < maxNumEls) {
			MeshFace& thisFc = faces[fc_i];
			fcEls = &thisFc.elements[0];
			if (fcEls[1] == max_int) {
				thisFc.getCentroid(fcCent, nodes);
				fcProj = thisFc.projDist;
				proj = globProjWt * avgProj + (1.0 - globProjWt) * fcProj;
				fcNorm = &thisFc.normDir[0];
				tgtPt[0] = fcCent[0] + 0.5 * proj * fcNorm[0];
				tgtPt[1] = fcCent[1] + 0.5 * proj * fcNorm[1];
				tgtPt[2] = fcCent[2] + 0.5 * proj * fcNorm[2];
				srchRad = 0.75 * proj;
				edgeClosed = adoptConnectedNd(fc_i, tgtPt, srchRad);
				if (!edgeClosed) {
					edgeClosed = adoptAnyNd(fc_i, tgtPt, srchRad);
				}
				if (!edgeClosed) {
					tgtPt[0] = fcCent[0] + proj * fcNorm[0];
					tgtPt[1] = fcCent[1] + proj * fcNorm[1];
					tgtPt[2] = fcCent[2] + proj * fcNorm[2];
					edgeClosed = createNewNd(fc_i, tgtPt);
				}
				if (!edgeClosed) {
					tgtPt[0] = fcCent[0] + 0.5 * proj * fcNorm[0];
					tgtPt[1] = fcCent[1] + 0.5 * proj * fcNorm[1];
					tgtPt[2] = fcCent[2] + 0.5 * proj * fcNorm[2];
					edgeClosed = createNewNd(fc_i, tgtPt);
				}
				if (edgeClosed) {
					elAdded = true;
				}
			}
			fc_i++;
		}
	}
	
	if (el_ct >= maxNumEls) {
		return false;
	}

	return true;
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
	int dim = nd_ct * 3;
	int* elNds;
	double eVol;
	double avgWt;
	double* crd;
	double res;
	double alpha;
	double beta;
	double rNext;
	double dp;

	vector<double> elMat(144);
	vector<double> xVec(dim);
	vector<double> hVec(dim);
	vector<double> zVec(dim);
	vector<double> gVec(dim);
	vector<double> wVec(dim);
	vector<double> dMat(dim);
	vector<double> pMat(dim);
	vector<double> pInv(dim);
	
	i1 = el_ct;
	vector<double> elWt(i1);

	// Make static element matrix
	for (i1 = 0; i1 < 144; i1++) {
		elMat[i1] = 0.0;
	}

	i3 = 0;
	for (i1 = 0; i1 < 12; i1++) {
		if (i3 < 9) {
			i4 = 0;
		}
		else {
			i4 = i3 - 9;
		}
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
	for (i2 = 0; i2 < el_ct; i2++) {
		MeshElement& thisEl = elements[i2];
		eVol = thisEl.getVolume(nodes);
		elWt[i1] = eVol;
		avgWt += eVol;
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

	//thisNd = nodes.getFirst();
	i1 = 0;
	for (i1 = 0; i1 < numBoundNds; i1++) {
		MeshNode& thisNd = nodes[i1];
		crd = &thisNd.coord[0];
		for (i3 = 0; i3 < 3; i3++) {
			i4 = i1 * 3 + i3;
			dMat[i4] = 100000.0;
			pMat[i4] += 100000.0;
			gVec[i4] = -100000.0 * crd[i3];
		}
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
		i2 = 0;
		for (i2 = 0; i2 < el_ct; i2++) {
			MeshElement& thisEl = elements[i2];
			elNds = &thisEl.nodes[0];
			i7 = 0; //Index in elMat
			for (i3 = 0; i3 < 4; i3++) {
				for (i4 = 0; i4 < 3; i4++) {
					i8 = 3*elNds[i3] + i4; // global row
					for (i5 = 0; i5 < 4; i5++) {
						for (i6 = 0; i6 < 3; i6++) {
							i9 = 3 * elNds[i5] + i6; //global col
							zVec[i8] += (elWt[i2] * elMat[i7] * hVec[i9]);
							i7++;
						}
					}
				}
			}
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
	for (i1 = numBoundNds; i1 < nd_ct; i1++) {
		MeshNode& thisNd = nodes[i1];
		crd = &thisNd.coord[0];
		for (i3 = 0; i3 < 3; i3++) {
			i4 = 3 * i1 + i3;
			crd[i3] = xVec[i4];
		}
	}

	return;
}

void Mesher::writeOutput(string fileName) {
	int i1;
	ofstream outFile;
	int* elNds;
	double* crd;

	outFile.open(fileName);
	
	outFile << "nodes:\n";
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& thisNd = nodes[i1];
		crd = &thisNd.coord[0];
		outFile << "  - [" << crd[0] << ", " << crd[1] << ", " << crd[2] << "]\n";
	}

	outFile << "elements:\n";
	for (i1 = 0; i1 < el_ct; i1++) {
		MeshElement& thisEl = elements[i1];
		elNds = &thisEl.nodes[0];
		outFile << "  - [" << elNds[0];
		for (i1 = 1; i1 < 4; i1++) {
			outFile << ", " << elNds[i1];
		}
		outFile << "]\n";
	}

	outFile.close();
	return;
}

void Mesher::printCurrentMesh() {
	int i1;
	double* crd;
	int* elNds;
	int* fcEls;
	
	cout << "Current Mesh:" << endl;
	cout << "Nodes:" << endl;
	
	//thisNd = nodes.getFirst();
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& thisNd = nodes[i1];
		crd = &thisNd.coord[0];
		cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
	}
	cout << "Elements:" << endl;
	//thisEl = elements.getFirst();
	for (i1 = 0; i1 < el_ct; i1++) {
		MeshElement& thisEl = elements[i1];
		elNds = &thisEl.nodes[0];
		cout << elNds[0] << ", ";
		cout << elNds[1] << ", ";
		cout << elNds[2] << ", ";
		cout << elNds[3] << endl;
	}

	cout << "Faces:" << endl;
	//thisFc = faces.getFirst();
	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& thisFc = faces[i1];
		cout << "  nodes: ";
		elNds = &thisFc.nodes[0];
		cout << elNds[0] << ", ";
		cout << elNds[1] << ", ";
		cout << elNds[2] << ", ";
		fcEls = &thisFc.elements[0];
		cout << "element pt: ";
		cout << fcEls[0] << ", " << fcEls[1] << endl;
	}

	return;
}