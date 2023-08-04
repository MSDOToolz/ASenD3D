#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "FaceClass.h"
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
	int i5;
	if(type == 4) {
		numNds = 4;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 4;
		defDim = 6;
		numIP = 1;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1 * 2]; //int[i1][2];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		intPts[0] = 0.25;
		intPts[1] = 0.25;
		intPts[2] = 0.25;
		ipWt[0] = r_1o6;
	} else if(type == 6) {
		numNds = 6;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 6;
		defDim = 6;
		numIP = 2;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1 * 2];// int[i1][2];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		intPts[0] = r_1o3;
		intPts[1] = r_1o3;
		intPts[2] = -r_1ort3;
		intPts[3] = r_1o3;
		intPts[4] = r_1o3;
		intPts[5] = r_1ort3;
		ipWt[0] = 0.5;
		ipWt[1] = 0.5;
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
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1 * 2];//int[i1][2];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		int sVal[2] = {-r_1ort3,r_1ort3};
		i4 = 0;
		i5 = 0;
		for (i1 = 0; i1 < 2; i1++) {
			for (i2 = 0; i2 < 2; i2++) {
				for (i3 = 0; i3 < 2; i3++) {
					intPts[i5] = sVal[i3];
					i5++;
					intPts[i5] = sVal[i2];
					i5++;
					intPts[i5] = sVal[i1];
					i5++;
					ipWt[i4] = 1.0;
					i4++;
				}
			}
		}
		if(type == 81) {
			i3 = 24;
			for (i1 = 8; i1 < 11; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					dofTable[2*i3] = i1;
					dofTable[2*i3+1] = i2;
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
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[2*i1];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		intPts[0] = r_1o6;
		intPts[1] = r_1o6;
		intPts[2] = 0.0;
		intPts[3] = r_2o3;
		intPts[4] = r_1o6;
		intPts[5] = 0.0;
		intPts[6] = r_1o6;
		intPts[7] = r_2o3;
		intPts[8] = 0.0;
		ipWt[0] = r_1o6;
		ipWt[1] = r_1o6;
		ipWt[2] = r_1o6;
		dofTable[36] = 3;
		dofTable[37] = 2;
		dofTable[38] = 4;
		dofTable[39] = 2;
		dofTable[40] = 5;
		dofTable[41] = 2;
	} else if(type == 41) {
		numNds = 4;
		dofPerNd = 6;
		numIntDof = 8;
		nDim = 10;
		defDim = 9;
		numIP = 4;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[2*i1];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		intPts[0] = -r_1ort3;
		intPts[1] = -r_1ort3;
		intPts[2] = 0.0;
		intPts[3] = r_1ort3;
		intPts[4] = -r_1ort3;
		intPts[5] = 0.0;
		intPts[6] = -r_1ort3;
		intPts[7] = r_1ort3;
		intPts[8] = 0.0;
		intPts[9] = r_1ort3;
		intPts[10] = r_1ort3;
		intPts[11] = 0.0;
		ipWt[0] = 1.0;
		ipWt[1] = 1.0;
		ipWt[2] = 1.0;
		ipWt[3] = 1.0;
		dofTable[48] = 4;
		dofTable[49] = 0;
		dofTable[50] = 5;
		dofTable[51] = 0;
		dofTable[52] = 4;
		dofTable[53] = 1;
		dofTable[54] = 5;
		dofTable[55] = 1;
		dofTable[56] = 6;
		dofTable[57] = 2;
		dofTable[58] = 7;
		dofTable[59] = 2;
		dofTable[60] = 8;
		dofTable[61] = 2;
		dofTable[62] = 9;
		dofTable[63] = 2;
	} else if(type == 2) {
		numNds = 2;
		dofPerNd = 6;
		numIntDof = 2;
		nDim = 3;
		defDim = 6;
		numIP = 2;
		nodes = new int[numNds];
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = new int[i1*2];
		intPts = new double[numIP*3];
		ipWt = new double[numIP];
		intPts[0] = -r_1ort3;
		intPts[1] = 0.0;
		intPts[2] = 0.0;
		intPts[3] = r_1ort3;
		intPts[4] = 0.0;
		intPts[5] = 0.0;
		ipWt[0] = 1.0;
		ipWt[1] = 1.0;
		dofTable[24] = 2;
		dofTable[25] = 1;
		dofTable[26] = 2;
		dofTable[27] = 2;
	}
	
	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		for (i2 = 0; i2 < dofPerNd; i2++) {
			dofTable[2*i3] = i1;
			dofTable[2*i3+1] = i2;
			i3++;
		}
	}
	
	if(numIntDof == 0) {
		internalDisp = nullptr;
		internaldLdu = nullptr;
		internalAdj = nullptr;
		internalRu = nullptr;
		internalMat = nullptr;
	} else {
		internalDisp = new double[numIntDof];
		internaldLdu = new double[numIntDof];
		internalAdj = new double[numIntDof];
		internalRu = new DiffDoub[numIntDof];
		i1 = (numNds*dofPerNd + numIntDof)*numIntDof;
		internalMat = new double[i1];
	}

	intDofIndex = 0;
	
	sectPtr = nullptr;
	nextEl = nullptr;
	
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

void Element::setSectPtr(Section *newSec) {
	sectPtr = newSec;
	return;
}

void Element::initializeFaces() {
	Face *newFc;
	
	if(type == 4) {
		newFc = new Face(3);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		faces.addFace(newFc);
		newFc = new Face(3);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,3,nodes[3]);
		faces.addFace(newFc);
		newFc = new Face(3);
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,3,nodes[3]);
		faces.addFace(newFc);
		newFc = new Face(3);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,2,nodes[2]);
		faces.addFace(newFc);
	} else if(type == 6) {
        newFc = new Face(3);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		faces.addFace(newFc);
		newFc = new Face(3);
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,4,nodes[4]);
		newFc->setNode(2,5,nodes[5]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,4,nodes[4]);
		newFc->setNode(3,3,nodes[3]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,4,nodes[4]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,2,nodes[2]);
		faces.addFace(newFc);
	} else if(type == 8 || type == 81) {
		newFc = new Face(4);
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		newFc->setNode(3,0,nodes[0]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,4,nodes[4]);
		newFc->setNode(1,5,nodes[5]);
		newFc->setNode(2,6,nodes[6]);
		newFc->setNode(3,7,nodes[7]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,4,nodes[4]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,6,nodes[6]);
		newFc->setNode(3,5,nodes[5]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,2,nodes[2]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,7,nodes[7]);
		newFc->setNode(3,6,nodes[6]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,0,nodes[0]);
		newFc->setNode(2,4,nodes[4]);
		newFc->setNode(3,7,nodes[7]);
		faces.addFace(newFc);
	} else if(type == 3) {
		newFc = new Face(3);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,2,nodes[2]);
		faces.addFace(newFc);
		newFc = new Face(3);
		newFc->setNode(0,2,nodes[2]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,0,nodes[0]);
		faces.addFace(newFc);
	} else if(type == 41) {
		newFc = new Face(4);
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,2,nodes[2]);
		newFc->setNode(3,3,nodes[3]);
		faces.addFace(newFc);
		newFc = new Face(4);
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		newFc->setNode(3,0,nodes[0]);
		faces.addFace(newFc);
	} 
	
	return;
}

void Element::setIntDofIndex(int newInd) {
	intDofIndex = newInd;
	return;
}

void Element::setIntdLdU(double globdLdU[]) {
	int i1;
	int i2 = intDofIndex;
	for (i1 = 0; i1 < numIntDof; i1++) {
		internaldLdu[i1] = globdLdU[i2];
		i2++;
	}
	return;
}

void Element::setNext(Element *newEl) {
	nextEl = newEl;
	return;
}

int Element::getType() {
	return type;
}

int Element::getLabel() {
	return label;
}

int Element::getNumNds() {
	return numNds;
}

int Element::getDofPerNd() {
	return dofPerNd;
}

int Element::getNumIntDof() {
	return numIntDof;
}

int Element::getNumIP() {
	return numIP;
}

double* Element::getIP() {
	return intPts;
}

int Element::getNumLayers() {
	return sectPtr->getNumLayers();
}

int* Element::getNodes() {
	return nodes;
}

IntList* Element::getDesignVars() {
	return &designVars;
}

IntListEnt* Element::getFirstCompDV() {
	return compDVars.getFirst();
}

Face* Element::getFirstFace() {
	return faces.getFirst();
}

Element* Element::getNext() {
	return nextEl;
}

void Element::addDesignVariable(int dIndex, double coef) {
	designVars.addEntry(dIndex);
	dvCoef.addEntry(coef);
	return;
}

void Element::addCompDVar(int dIndex) {
	compDVars.addIfAbsent(dIndex);
	return;
}

void Element::destroy() {
	delete[] nodes;
	delete[] dofTable;
	delete[] intPts;
	delete[] ipWt;
	if (numIntDof > 0) {
		delete[] internalDisp;
		delete[] internaldLdu;
		delete[] internalAdj;
		delete[] internalRu;
		delete[] internalMat;
	}
	faces.destroy();
	designVars.destroy();
	dvCoef.destroy();
	compDVars.destroy();
	return;
}

ElPt::ElPt() {
	ptr = nullptr;
}

//ElementList begin
ElementList::ElementList() {
	firstEl = nullptr;
	lastEl = nullptr;
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

void ElementList::destroy() {
	Element* thisEl = firstEl;
	Element* nextEl;
	while (thisEl) {
		nextEl = thisEl->getNext();
		thisEl->destroy();
		delete thisEl;
		thisEl = nextEl;
	}
	return;
}

// Stress prerequisite classes

//dup1
 DoubStressPrereq::DoubStressPrereq() {
	layerZ = nullptr;
	layerThk = nullptr;
	layerAng = nullptr;
	layerQ = nullptr;
	currentLayLen = 0;
	return;
}

void DoubStressPrereq::destroy() {
	if(layerZ) {
	    delete[] layerZ;
	    delete[] layerThk;
	    delete[] layerAng;
	    delete[] layerQ;
	}
	layerZ = nullptr;
	layerThk = nullptr;
	layerAng = nullptr;
	layerQ = nullptr;
	currentLayLen = 0;

	return;
}
//end dup
//skip 
 
//DiffDoub versions: 
 DiffDoubStressPrereq::DiffDoubStressPrereq() {
	layerZ = nullptr;
	layerThk = nullptr;
	layerAng = nullptr;
	layerQ = nullptr;
	currentLayLen = 0;
	return;
}

void DiffDoubStressPrereq::destroy() {
	if(layerZ) {
	    delete[] layerZ;
	    delete[] layerThk;
	    delete[] layerAng;
	    delete[] layerQ;
	}
	layerZ = nullptr;
	layerThk = nullptr;
	layerAng = nullptr;
	layerQ = nullptr;
	currentLayLen = 0;

	return;
}
 
//end skip 
