#include <vector>
#include "FaceClass.h"
#include "DiffDoubClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"
#include "matrixFunctions.h"

using namespace std;

Face::Face() {
	numNds = 0;
	onSurf = true;
	return;
}

void Face::setNode(int place, int locNd, int globNd) {
	locNodes[place] = locNd;
	globNodes[place] = globNd;
	return;
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

//dup1

void Face::getAreaNormal(DiffDoub0& area, DiffDoub0 norm[], vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	DiffDoub0 v1[3];
	DiffDoub0 v2[3];
	DiffDoub0 tmpV[3];
	DiffDoub0 tmp;

	if (numNds == 4) {
		ndAr[globNodes[2]].getCrd(v1, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[3]].getCrd(v2, dvAr);
		ndAr[globNodes[1]].getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}
	else {
		ndAr[globNodes[1]].getCrd(v1, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[2]].getCrd(v2, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
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
 
//DiffDoub1 versions: 
//dup1

void Face::getAreaNormal(DiffDoub1& area, DiffDoub1 norm[], vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	DiffDoub1 v1[3];
	DiffDoub1 v2[3];
	DiffDoub1 tmpV[3];
	DiffDoub1 tmp;

	if (numNds == 4) {
		ndAr[globNodes[2]].getCrd(v1, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[3]].getCrd(v2, dvAr);
		ndAr[globNodes[1]].getCrd(tmpV, dvAr);
		v2[0].sub(tmpV[0]);
		v2[1].sub(tmpV[1]);
		v2[2].sub(tmpV[2]);
	}
	else {
		ndAr[globNodes[1]].getCrd(v1, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
		v1[0].sub(tmpV[0]);
		v1[1].sub(tmpV[1]);
		v1[2].sub(tmpV[2]);
		ndAr[globNodes[2]].getCrd(v2, dvAr);
		ndAr[globNodes[0]].getCrd(tmpV, dvAr);
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
 
 
//FacePtList

FacePtList::FacePtList() {
	return;
}

void FacePtList::addFace(int newI) {
	fcList.push_back(newI);
	return;
}

bool FacePtList::addIfAbsent(int newI, vector<Face>& globFaces) {
	int i1;
	int newNumNds;
	int newSrtd[8];
	int thisNumNds;
	int thisSrtd[8];
	bool allMatch;

	Face& newFc = globFaces[newI];
	newNumNds = newFc.numNds;
	newFc.sortedNodes(newSrtd);
	for (auto& fi : fcList) {
		Face& thisFc = globFaces[fi];
		thisNumNds = thisFc.numNds;
		if (thisNumNds == newNumNds) {
			thisFc.sortedNodes(thisSrtd);
			allMatch = true;
			for (i1 = 0; i1 < thisNumNds; i1++) {
				if (thisSrtd[i1] != newSrtd[i1]) {
					allMatch = false;
				}
			}
			if (allMatch) {
				newFc.onSurf = false;
				thisFc.onSurf = false;
				return false;
			}
		}
	}

	fcList.push_back(newI);
	return true;
}