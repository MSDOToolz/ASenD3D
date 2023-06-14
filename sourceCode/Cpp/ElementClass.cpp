#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"

using namespace std;

const double r_2o3 = 0.666666666666666667;
const double r_1o3 = 0.333333333333333333;
const double r_1o6 = 0.166666666666666667;
const double r_pi = 3.14159265358979324;
const double r_pio180 = 0.0174532925199432958;
const double r_1ort3 = 0.577350269189625765;


Element::Element(int newType) {
	type = newType;
	int i1;
	int i2;
	int i3;
	int i4;
	if(type == 4) {
		numNds = 4;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 4;
		defDim = 6;
		numIP = 1;
		numFaces = 4;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		intPts[0][0] = 0.25;
		intPts[0][1] = 0.25;
		intPts[0][2] = 0.25;
		ipWt[0] = r_1o6;
		faceNds[0][0] = 0;
		faceNds[0][1] = 2;
		faceNds[0][2] = 1;
		faceNds[1][0] = 0;
		faceNds[1][1] = 1;
		faceNds[1][2] = 3;
		faceNds[2][0] = 1;
		faceNds[2][1] = 2;
		faceNds[2][2] = 3;
		faceNds[3][0] = 0;
		faceNds[3][1] = 3;
		faceNds[3][2] = 2;
	} else if(type == 6) {
		numNds = 6;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 6;
		defDim = 6;
		numIP = 2;
		numFaces = 5;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		intPts[0][0] = r_1o3;
		intPts[0][1] = r_1o3;
		intPts[0][2] = -r_1ort3;
		intPts[1][0] = r_1o3;
		intPts[1][1] = r_1o3;
		intPts[1][2] = r_1ort3;
		ipWt[0] = 0.5;
		ipWt[1] = 0.5;
		faceNds[0][0] = 0;
		faceNds[0][1] = 2;
		faceNds[0][2] = 1;
		faceNds[1][0] = 3;
		faceNds[1][1] = 4;
		faceNds[1][2] = 5;
		faceNds[2][0] = 0;
		faceNds[2][1] = 1;
		faceNds[2][2] = 4;
		faceNds[2][3] = 3;
		faceNds[3][0] = 1;
		faceNds[3][1] = 2;
		faceNds[3][2] = 5;
		faceNds[3][3] = 4;
		faceNds[4][0] = 0;
		faceNds[4][1] = 3;
		faceNds[4][2] = 5;
		faceNds[4][3] = 2;
	} else if(type == 8 || type == 81) {
		numNds = 8;
		dofPerNd = 3;
		if(type == 8) {
		    numIntDof = 0;
			nDim = 8;
		} else {
			numIntDof = 9;
			nDim = 11;
		}
		defDim = 6;
		numIP = 8;
		numFaces = 6;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		int sVal[2] = {-r_1ort3,r_1ort3};
		i4 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 2; i2++) {
				for (i3 = 0; i3 < 2; i3++) {
					intPts[i4][0] = sVal[i3];
					intPts[i4][1] = sVal[i2];
					intPts[i4][2] = sVal[i1];
					ipWt[i4] = 1.0;
				}
			}
		}
		faceNds[0][0] = 3;
		faceNds[0][1] = 2;
		faceNds[0][2] = 1;
		faceNds[0][3] = 0;
		faceNds[1][0] = 4;
		faceNds[1][1] = 5;
		faceNds[1][2] = 6;
		faceNds[1][3] = 7;
		faceNds[2][0] = 0;
		faceNds[2][1] = 1;
		faceNds[2][2] = 5;
		faceNds[2][3] = 4;
		faceNds[3][0] = 1;
		faceNds[3][1] = 2;
		faceNds[3][2] = 6;
		faceNds[3][3] = 5;
		faceNds[4][0] = 2;
		faceNds[4][1] = 3;
		faceNds[4][2] = 7;
		faceNds[4][3] = 6;
		faceNds[5][0] = 3;
		faceNds[5][1] = 0;
		faceNds[5][2] = 4;
		faceNds[5][3] = 7;
		if(type == 81) {
			i3 = 24;
			for (i1 = 8; i1 < 11; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					dofTable[i3][0] = i1;
					dofTable[i3][1] = i2;
					i3++;
				}
			}
		}
	} else if(type == 3) {
		numNds = 3;
		dofPerNd = 6;
		numIntDof = 3;
		nDim = 6;
		defDim = 9;
		numIP = 3;
		numFaces = 2;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		intPts[0][0] = r_1o6;
		intPts[0][1] = r_1o6;
		intPts[0][2] = 0.0;
		intPts[1][0] = r_2o3;
		intPts[1][1] = r_1o6;
		intPts[1][2] = 0.0;
		intPts[2][0] = r_1o6;
		intPts[2][1] = r_2o3;
		intPts[2][2] = 0.0;
		ipWt[0] = r_1o6;
		ipWt[1] = r_1o6;
		ipWt[2] = r_1o6;
		faceNds[0][0] = 0;
		faceNds[0][1] = 1;
		faceNds[0][2] = 2;
		faceNds[1][0] = 0;
		faceNds[1][1] = 2;
		faceNds[1][2] = 1;
		dofTable[18][0] = 3;
		dofTable[18][1] = 2;
		dofTable[19][0] = 4;
		dofTable[19][1] = 2;
		dofTable[20][0] = 5;
		dofTable[20][1] = 2;
	} else if(type == 41) {
		numNds = 4;
		dofPerNd = 6;
		numIntDof = 8;
		nDim = 10;
		defDef = 9;
		numIP = 4;
		numFaces = 2;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		intPts[0][0] = -r_1ort3;
		intPts[0][1] = -r_1ort3;
		intPts[0][2] = 0.0;
		intPts[1][0] = r_1ort3;
		intPts[1][1] = -r_1ort3;
		intPts[1][2] = 0.0;
		intPts[2][0] = -r_1ort3;
		intPts[2][1] = r_1ort3;
		intPts[2][2] = 0.0;
		intPts[3][0] = r_1ort3;
		intPts[3][1] = r_1ort3;
		intPts[3][2] = 0.0;
		ipWt[0] = 1.0;
		ipWt[1] = 1.0;
		ipWt[2] = 1.0;
		ipWt[3] = 1.0;
	    faceNds[0][0] = 0;
		faceNds[0][1] = 1;
		faceNds[0][2] = 2;
		faceNds[0][3] = 3;
		faceNds[1][0] = 0;
		faceNds[1][1] = 3;
		faceNds[1][2] = 2;
		faceNds[1][3] = 1;
		dofTable[24][0] = 4;
		dofTable[24][1] = 0;
		dofTable[25][0] = 5;
		dofTable[25][1] = 0;
		dofTable[26][0] = 4;
		dofTable[26][1] = 1;
		dofTable[27][0] = 5;
		dofTable[27][1] = 1;
		dofTable[28][0] = 6;
		dofTable[28][1] = 2;
		dofTable[29][0] = 7;
		dofTable[29][1] = 2;
		dofTable[30][0] = 8;
		dofTable[30][1] = 2;
		dofTable[31][0] = 9;
		dofTable[31][1] = 2;
	} else if(type == 2) {
		numNds = 2;
		dofPerNd = 6;
		numIntDof = 2;
		nDim = 3;
		defDim = 6;
		numIP = 2;
		numFaces = 1;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1][2];
		intPts = new double[numIP][3];
		ipWt = new double[numIP];
		faceNds = new int[numFaces][6];
		for (i1 = 0; i1 < numFaces; i1++) {
			for (i2 = 0; i2 < 6; i2++) {
				faceNds[i1][i2] = -1;
			}
		}
		intPts[0][0] = -r_1ort3;
		intPts[0][1] = 0.0;
		intPts[0][2] = 0.0;
		intPts[1][0] = r_1ort3;
		intPts[1][1] = 0.0;
		intPts[1][2] = 0.0;
		ipWt[0] = 1.0;
		ipWt[1] = 1.0;
		dofTable[12][0] = 2;
		dofTable[12][1] = 1;
		dofTable[13][0] = 2;
		dofTable[13][1] = 2;
	}
	
	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		for (i2 = 0; i2 < dofPerNd; i2++) {
			dofTable[i3][0] = i1;
			dofTable[i3][1] = i2;
			i3++;
		}
	}
	
	if(numIntDof == 0) {
		internalDisp = NULL;
		internaldLdu = NULL;
		internalAdj = NULL;
		internalRu = NULL;
		internalMat = NULL;
	} else {
		internalDisp = new double[numIntDof];
		internaldLdu = new double[numIntDof];
		internalAdj = new double[numIntDof];
		internalRu = new Doub[numIntDof];
		i1 = (numNds*dofPerNd + numIntDof)*numIntDof;
		internalMat = new double[i1];
	}
	
	sectPtr = NULL;
	nextEl = NULL;
	
	return;
}

void Element::setLabel(int newLab) {
	label = newLab;
	return;
}

void Element::setNodes(int newNds[]) {
	int i1;
	for (i1 = 0; i1 < numNds; i1++) {
		nodes[i1] = newNds[i1];
	}
	return;
}

void Element::setNext(Element *newEl) {
	nextEl = newEl;
	return;
}

Element* Element::getNext() {
	return nextEl;
}

ElPt::ElPt() {
	ptr = NULL;
}

//ElementList begin
ElementList::ElementList() {
	firstEl = NULL;
	lastEl = NULL;
	length = 0;
	return;
}

void ElementList::addElement(Element *newEl) {
	if(!firstEl) {
		firstEl = newEl;
		lastEl = newEl;
	} else {
		lastEl->setNext(newEl);
		lastEl = newEl;
	}
	length++;
	return;
}

int ElementList::getLength() {
	return length;
}

Element* ElementList::getFirst() {
	return firstEl;
}