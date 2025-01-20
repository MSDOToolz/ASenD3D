#include <vector>
#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "FaceClass.h"
#include "matrixFunctions.h"

using namespace std;

const int max_int = 2000000000;
const double r_1ort3 = 0.577350269189625765; //0.577350269189625765;
const double r_2o3 = 0.666666666666666667;
const double r_1o3 = 0.333333333333333333;
const double r_1o6 = 0.166666666666666667;
const double r_1o24 = 0.041666666666666667;
const double r_pi = 3.14159265358979324;
const double r_pio180 = 0.0174532925199432958;
const double r_tet1 = 0.1381966011250071;
const double r_tet2 = 0.585410196624929;

// Stress prerequisite classes

//dup1
DiffDoub0StressPrereq::DiffDoub0StressPrereq() {
	globNds = vector<DiffDoub0>(30);
	locNds = vector<DiffDoub0>(30);
	locOri = vector<DiffDoub0>(9);
	instOri = vector<DiffDoub0>(720);
	globDisp = vector<DiffDoub0>(60);
	globVel = vector<DiffDoub0>(30);
	globAcc = vector<DiffDoub0>(30);
	globTemp = vector<DiffDoub0>(10);
	globTdot = vector<DiffDoub0>(10);
	Cmat = vector<DiffDoub0>(81);
	Mmat = vector<DiffDoub0>(36);
	Dmat = vector<DiffDoub0>(81);
	thermExp = vector<DiffDoub0>(6);
	Einit = vector<DiffDoub0>(6);
	TCmat = vector<DiffDoub0>(9);
	BMat = vector<DiffDoub0>(288);
	CBMat = vector<DiffDoub0>(288);
	frcFldCoef = vector<DiffDoub0>(2);
	frcFldExp = vector<DiffDoub0>(2);
	thrmFldCoef = vector<DiffDoub0>(2);
	scrMat1 = vector<double>(3600);
	scrMat2 = vector<double>(3600);
	scrMat3 = vector<double>(3600);
	scrMat4 = vector<double>(3600);
	scrMat5 = vector<double>(3600);
	scrVec1 = vector<DiffDoub0>(60);
	scrVec2 = vector<DiffDoub0>(60);
	scrVec3 = vector<DiffDoub0>(60);
	scrVec4 = vector<DiffDoub0>(60);
	scrVec5 = vector<DiffDoub0>(60);
	currentLayLen = 0;
	return;
}

void DiffDoub0StressPrereq::allocateLayers(int numLayers) {
	if (numLayers != 0) {
		layerZ = vector<DiffDoub0>(numLayers);
		layerThk = vector<DiffDoub0>(numLayers);
		layerAng = vector<DiffDoub0>(numLayers);
		layerQ = vector<DiffDoub0>(9 * numLayers);
		layerD = vector<DiffDoub0>(9 * numLayers);
		layerTE = vector<DiffDoub0>(3 * numLayers);
		layerE0 = vector<DiffDoub0>(3 * numLayers);
		layerDen = vector<DiffDoub0>(numLayers);
		layerTC = vector<DiffDoub0>(9 * numLayers);
		layerSH = vector<DiffDoub0>(numLayers);
	}
	currentLayLen = numLayers;
	return;
}

DiffDoub0FlPrereq::DiffDoub0FlPrereq() {
	globNds = vector<DiffDoub0>(30);
	globDisp = vector<DiffDoub0>(30);
	globVel = vector<DiffDoub0>(30);
	flDen = vector<DiffDoub0>(10);
	flVel = vector<DiffDoub0>(30);
	flTemp = vector<DiffDoub0>(10);
	flTurbE = vector<DiffDoub0>(10);
	flDenDot = vector<DiffDoub0>(10);
	flVelDot = vector<DiffDoub0>(30);
	flTDot = vector<DiffDoub0>(10);
	flTurbEDot = vector<DiffDoub0>(10);
	scratch = vector<double>(3600);

	return;
}

//end dup 
 
//skip 
 
//DiffDoub1 versions: 
//dup1
DiffDoub1StressPrereq::DiffDoub1StressPrereq() {
	globNds = vector<DiffDoub1>(30);
	locNds = vector<DiffDoub1>(30);
	locOri = vector<DiffDoub1>(9);
	instOri = vector<DiffDoub1>(720);
	globDisp = vector<DiffDoub1>(60);
	globVel = vector<DiffDoub1>(30);
	globAcc = vector<DiffDoub1>(30);
	globTemp = vector<DiffDoub1>(10);
	globTdot = vector<DiffDoub1>(10);
	Cmat = vector<DiffDoub1>(81);
	Mmat = vector<DiffDoub1>(36);
	Dmat = vector<DiffDoub1>(81);
	thermExp = vector<DiffDoub1>(6);
	Einit = vector<DiffDoub1>(6);
	TCmat = vector<DiffDoub1>(9);
	BMat = vector<DiffDoub1>(288);
	CBMat = vector<DiffDoub1>(288);
	frcFldCoef = vector<DiffDoub1>(2);
	frcFldExp = vector<DiffDoub1>(2);
	thrmFldCoef = vector<DiffDoub1>(2);
	scrMat1 = vector<double>(3600);
	scrMat2 = vector<double>(3600);
	scrMat3 = vector<double>(3600);
	scrMat4 = vector<double>(3600);
	scrMat5 = vector<double>(3600);
	scrVec1 = vector<DiffDoub1>(60);
	scrVec2 = vector<DiffDoub1>(60);
	scrVec3 = vector<DiffDoub1>(60);
	scrVec4 = vector<DiffDoub1>(60);
	scrVec5 = vector<DiffDoub1>(60);
	currentLayLen = 0;
	return;
}

void DiffDoub1StressPrereq::allocateLayers(int numLayers) {
	if (numLayers != 0) {
		layerZ = vector<DiffDoub1>(numLayers);
		layerThk = vector<DiffDoub1>(numLayers);
		layerAng = vector<DiffDoub1>(numLayers);
		layerQ = vector<DiffDoub1>(9 * numLayers);
		layerD = vector<DiffDoub1>(9 * numLayers);
		layerTE = vector<DiffDoub1>(3 * numLayers);
		layerE0 = vector<DiffDoub1>(3 * numLayers);
		layerDen = vector<DiffDoub1>(numLayers);
		layerTC = vector<DiffDoub1>(9 * numLayers);
		layerSH = vector<DiffDoub1>(numLayers);
	}
	currentLayLen = numLayers;
	return;
}

DiffDoub1FlPrereq::DiffDoub1FlPrereq() {
	globNds = vector<DiffDoub1>(30);
	globDisp = vector<DiffDoub1>(30);
	globVel = vector<DiffDoub1>(30);
	flDen = vector<DiffDoub1>(10);
	flVel = vector<DiffDoub1>(30);
	flTemp = vector<DiffDoub1>(10);
	flTurbE = vector<DiffDoub1>(10);
	flDenDot = vector<DiffDoub1>(10);
	flVelDot = vector<DiffDoub1>(30);
	flTDot = vector<DiffDoub1>(10);
	flTurbEDot = vector<DiffDoub1>(10);
	scratch = vector<double>(3600);

	return;
}

//end dup 
 
//end skip 
 
 
Element::Element() {
	return;
}

void Element::initializeType(int newType) {
	type = newType;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;

	sCent[0] = 0.0;
	sCent[1] = 0.0;
	sCent[2] = 0.0;
	if(type == 4 || type == 400) {
		numNds = 4;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 4;
		defDim = 6;
		numIP = 1;
		numFaces = 4;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(i1 * 2); //int[i1][2];
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
		intPts[0] = 0.25;
		intPts[1] = 0.25;
		intPts[2] = 0.25;
		ipWt[0] = r_1o6;
		ndSpts = vector<double>(12);
		ndSpts[0] = 0.0;
		ndSpts[1] = 0.0;
		ndSpts[2] = 0.0;
		ndSpts[3] = 1.0;
		ndSpts[4] = 0.0;
		ndSpts[5] = 0.0;
		ndSpts[6] = 0.0;
		ndSpts[7] = 1.0;
		ndSpts[8] = 0.0;
		ndSpts[9] = 0.0;
		ndSpts[10] = 0.0;
		ndSpts[11] = 1.0;
	} else if(type == 6 || type == 600) {
		numNds = 6;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 6;
		defDim = 6;
		numIP = 2;
		numFaces = 5;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(i1 * 2);// int[i1][2];
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
		intPts[0] = r_1o3;
		intPts[1] = r_1o3;
		intPts[2] = -r_1ort3;
		intPts[3] = r_1o3;
		intPts[4] = r_1o3;
		intPts[5] = r_1ort3;
		ipWt[0] = 0.5;
		ipWt[1] = 0.5;
		sCent[0] = r_1o3;
		sCent[1] = r_1o3;
		sCent[2] = 0.0;
		ndSpts = vector<double>(18);
		ndSpts[0] = 0.0;
		ndSpts[1] = 0.0;
		ndSpts[2] = -1.0;
		ndSpts[3] = 1.0;
		ndSpts[4] = 0.0;
		ndSpts[5] = -1.0;
		ndSpts[6] = 0.0;
		ndSpts[7] = 1.0;
		ndSpts[8] = -1.0;
		ndSpts[9] = 0.0;
		ndSpts[10] = 0.0;
		ndSpts[11] = 1.0;
		ndSpts[12] = 1.0;
		ndSpts[13] = 0.0;
		ndSpts[14] = 1.0;
		ndSpts[15] = 0.0;
		ndSpts[16] = 1.0;
		ndSpts[17] = 1.0;
	} else if(type == 8 || type == 81 || type == 800) {
		numNds = 8;
		dofPerNd = 3;
		if(type == 8 || type == 800) {
		    numIntDof = 0;
			nDim = 8;
		} else {
			numIntDof = 9;
			nDim = 11;
		}
		defDim = 6;
		numIP = 8;
		numFaces = 6;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(i1 * 2);//int[i1][2];
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
		double sVal[2] = {-r_1ort3,r_1ort3};
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
		ndSpts = vector<double>(24);
		ndSpts[0] = -1.0;
		ndSpts[1] = -1.0;
		ndSpts[2] = -1.0;
		ndSpts[3] = 1.0;
		ndSpts[4] = -1.0;
		ndSpts[5] = -1.0;
		ndSpts[6] = 1.0;
		ndSpts[7] = 1.0;
		ndSpts[8] = -1.0;
		ndSpts[9] = -1.0;
		ndSpts[10] = 1.0;
		ndSpts[11] = -1.0;
		ndSpts[12] = -1.0;
		ndSpts[13] = -1.0;
		ndSpts[14] = 1.0;
		ndSpts[15] = 1.0;
		ndSpts[16] = -1.0;
		ndSpts[17] = 1.0;
		ndSpts[18] = 1.0;
		ndSpts[19] = 1.0;
		ndSpts[20] = 1.0;
		ndSpts[21] = -1.0;
		ndSpts[22] = 1.0;
		ndSpts[23] = 1.0;
	}
	else if (type == 10 || type == 1000) {
		numNds = 10;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 10;
		defDim = 6;
		numIP = 4;
		numFaces = 4;
		nodes = vector<int>(numNds);
		i1 = numNds * dofPerNd + numIntDof;
		dofTable = vector<int>(2 * i1);
		intPts = vector<double>(numIP * 3);
		ipWt = vector<double>(numIP);
		intPts[0] = r_tet1;
		intPts[1] = r_tet1;
		intPts[2] = r_tet1;
		intPts[3] = r_tet2;
		intPts[4] = r_tet1;
		intPts[5] = r_tet1;
		intPts[6] = r_tet1;
		intPts[7] = r_tet2;
		intPts[8] = r_tet1;
		intPts[9] = r_tet1;
		intPts[10] = r_tet1;
		intPts[11] = r_tet2;
		ipWt[0] = r_1o24;
		ipWt[1] = r_1o24;
		ipWt[2] = r_1o24;
		ipWt[3] = r_1o24;
		sCent[0] = 0.25;
		sCent[1] = 0.25;
		sCent[2] = 0.25;
		ndSpts = vector<double>(30);
		ndSpts[0] = 0.0;
		ndSpts[1] = 0.0;
		ndSpts[2] = 0.0;
		ndSpts[3] = 1.0;
		ndSpts[4] = 0.0;
		ndSpts[5] = 0.0;
		ndSpts[6] = 0.0;
		ndSpts[7] = 1.0;
		ndSpts[8] = 0.0;
		ndSpts[9] = 0.0;
		ndSpts[10] = 0.0;
		ndSpts[11] = 1.0;
		ndSpts[12] = 0.5;
		ndSpts[13] = 0.0;
		ndSpts[14] = 0.0;
		ndSpts[15] = 0.5;
		ndSpts[16] = 0.5;
		ndSpts[17] = 0.0;
		ndSpts[18] = 0.0;
		ndSpts[19] = 0.5;
		ndSpts[20] = 0.0;
		ndSpts[21] = 0.0;
		ndSpts[22] = 0.0;
		ndSpts[23] = 0.5;
		ndSpts[24] = 0.5;
		ndSpts[25] = 0.0;
		ndSpts[26] = 0.5;
		ndSpts[27] = 0.0;
		ndSpts[28] = 0.5;
		ndSpts[29] = 0.5;
	}
	else if (type == 3) {
		numNds = 3;
		dofPerNd = 6;
		numIntDof = 3;
		nDim = 6;
		defDim = 9;
		numIP = 3;
		numFaces = 2;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(2*i1);
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
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
		sCent[0] = r_1o3;
		sCent[1] = r_1o3;
		sCent[2] = r_1o3;
	} else if(type == 41) {
		numNds = 4;
		dofPerNd = 6;
		numIntDof = 8;
		nDim = 10;
		defDim = 9;
		numIP = 4;
		numFaces = 2;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(2*i1);
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
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
		numFaces = 0;
		nodes = vector<int>(numNds);
		i1 = numNds*dofPerNd + numIntDof;
		dofTable = vector<int>(i1*2);
		intPts = vector<double>(numIP*3);
		ipWt = vector<double>(numIP);
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
	else if (type == 21) { //Force field
		numNds = 2;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 2;
		defDim = 0;
		numIP = 0;
		numFaces = 0;
		nodes = vector<int>(numNds);
		dofTable = vector<int>(12);
		intPts = vector<double>(1);
		ipWt = vector<double>(1);
	}
	else if (type == 1) { // Mass
		numNds = 1;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 1;
		defDim = 0;
		numIP = 0;
		numFaces = 0;
		nodes = vector<int>(numNds);
		dofTable = vector<int>(6);
		intPts = vector<double>(1);
		intPts = vector<double>(1);
	}
	
	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		for (i2 = 0; i2 < dofPerNd; i2++) {
			dofTable[i3] = i1;
			dofTable[i3+1] = i2;
			i3+= 2;
		}
	}
	
	if(numIntDof != 0) {
		internalDisp = vector<double>(numIntDof);
		intPrevDisp = vector<double>(numIntDof);
		internaldLdu = vector<double>(numIntDof);
		internalAdj = vector<double>(numIntDof);
		internalRu = vector<DiffDoub1>(numIntDof);
		i1 = (numNds*dofPerNd + numIntDof)*numIntDof;
		internalMat = vector<double>(i1);
	}

	intDofIndex = 0;
	
	sectPtr = max_int;
	
	return;
}

void Element::setNodes(int newNds[]) {
	int i1;
	for (i1 = 0; i1 < numNds; i1++) {
		nodes[i1] = newNds[i1];
	}
	return;
}

void Element::initializeFaces(vector<Face>& globFcLst, int& fi) {
	//fi = the current number of faces that have been written into globFcLst
	Face* newFc = nullptr;
	if(type == 4 || type == 400) {
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,2,nodes[2]);
		faces.push_back(fi);
		fi++;
	} else if(type == 6 || type == 600) {
        globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,4,nodes[4]);
		newFc->setNode(2,5,nodes[5]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,4,nodes[4]);
		newFc->setNode(3,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,4,nodes[4]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,2,nodes[2]);
		faces.push_back(fi);
		fi++;
	} else if(type == 8 || type == 81 || type == 800) {
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		newFc->setNode(3,0,nodes[0]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,4,nodes[4]);
		newFc->setNode(1,5,nodes[5]);
		newFc->setNode(2,6,nodes[6]);
		newFc->setNode(3,7,nodes[7]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,5,nodes[5]);
		newFc->setNode(3,4,nodes[4]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,1,nodes[1]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,6,nodes[6]);
		newFc->setNode(3,5,nodes[5]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,2,nodes[2]);
		newFc->setNode(1,3,nodes[3]);
		newFc->setNode(2,7,nodes[7]);
		newFc->setNode(3,6,nodes[6]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,0,nodes[0]);
		newFc->setNode(2,4,nodes[4]);
		newFc->setNode(3,7,nodes[7]);
		faces.push_back(fi);
		fi++;
	}
	else if (type == 10 || type == 1000) {
		globFcLst[fi].numNds = 6;
		newFc = &globFcLst[fi];
		newFc->setNode(0, 0, nodes[0]);
		newFc->setNode(1, 2, nodes[2]);
		newFc->setNode(2, 1, nodes[1]);
		newFc->setNode(3, 6, nodes[6]);
		newFc->setNode(4, 5, nodes[5]);
		newFc->setNode(5, 4, nodes[4]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 6;
		newFc = &globFcLst[fi];
		newFc->setNode(0, 0, nodes[0]);
		newFc->setNode(1, 1, nodes[1]);
		newFc->setNode(2, 3, nodes[3]);
		newFc->setNode(3, 4, nodes[4]);
		newFc->setNode(4, 8, nodes[8]);
		newFc->setNode(5, 7, nodes[7]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 6;
		newFc = &globFcLst[fi];
		newFc->setNode(0, 1, nodes[1]);
		newFc->setNode(1, 2, nodes[2]);
		newFc->setNode(2, 3, nodes[3]);
		newFc->setNode(3, 5, nodes[5]);
		newFc->setNode(4, 9, nodes[9]);
		newFc->setNode(5, 8, nodes[8]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 6;
		newFc = &globFcLst[fi];
		newFc->setNode(0, 0, nodes[0]);
		newFc->setNode(1, 3, nodes[3]);
		newFc->setNode(2, 2, nodes[2]);
		newFc->setNode(3, 7, nodes[7]);
		newFc->setNode(4, 9, nodes[9]);
		newFc->setNode(5, 6, nodes[6]);
		faces.push_back(fi);
		fi++;
	}
	else if (type == 3) {
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,2,nodes[2]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 3;
		newFc = &globFcLst[fi];
		newFc->setNode(0,2,nodes[2]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,0,nodes[0]);
		faces.push_back(fi);
		fi++;
	} else if(type == 41) {
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,0,nodes[0]);
		newFc->setNode(1,1,nodes[1]);
		newFc->setNode(2,2,nodes[2]);
		newFc->setNode(3,3,nodes[3]);
		faces.push_back(fi);
		fi++;
		globFcLst[fi].numNds = 4;
		newFc = &globFcLst[fi];
		newFc->setNode(0,3,nodes[3]);
		newFc->setNode(1,2,nodes[2]);
		newFc->setNode(2,1,nodes[1]);
		newFc->setNode(3,0,nodes[0]);
		faces.push_back(fi);
		fi++;
	} 
	
	return;
}

void Element::setIntDisp(double newDisp[]) {
	int i1;
	for (i1 = 0; i1 < numIntDof; i1++) {
		internalDisp[i1] = newDisp[i1];
	}
	return;
}

void Element::setIntPrevDisp(double newDisp[]) {
	int i1;
	for (i1 = 0; i1 < numIntDof; i1++) {
		intPrevDisp[i1] = newDisp[i1];
	}
	return;
}

void Element::advanceIntDisp() {
	int i1;
	for (i1 = 0; i1 < numIntDof; i1++) {
		intPrevDisp[i1] = internalDisp[i1];
	}
	return;
}

void Element::backstepIntDisp() {
	int i1;
	for (i1 = 0; i1 < numIntDof; i1++) {
		internalDisp[i1] = intPrevDisp[i1];
	}
	return;
}

void Element::setIntdLdU(vector<double>& globdLdU) {
	int i1;
	int i2 = intDofIndex;
	for (i1 = 0; i1 < numIntDof; i1++) {
		internaldLdu[i1] = globdLdU[i2];
		i2++;
	}
	return;
}

int Element::getNumLayers(vector<Section>& secLst) {
	return secLst[sectPtr].layers.size();
}


void Element::addDesignVariable(int dIndex, double coef) {
	int i1 = designVars.size();
	IDCapsule dv;
	dv.intDat = dIndex;
	dv.doubDat = coef;
	designVars.push_back(dv);
	return;
}

void Element::addCompDVar(int dIndex) {
	for (auto& dv : compDVars) {
		if (dv == dIndex) {
			return;
		}
	}
	compDVars.push_back(dIndex);
	return;
}


