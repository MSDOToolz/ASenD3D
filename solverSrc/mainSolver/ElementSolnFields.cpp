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

//dup1
void Element::getNdDisp(DiffDoub0 globDisp[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub0 ndDisp[6];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getDisp(ndDisp);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globDisp[i3].setVal(ndDisp[i2]);
			i3+= nDim;
		}
	}
	
	if(numIntDof > 0) {
		i2 = 2*numNds*dofPerNd;
		for (i3 = 0; i3 < numIntDof; i3++) {
			i4 = dofTable[i2];
			i5 = dofTable[i2+1];
			i6 = nDim*i5 + i4;
			globDisp[i6].setVal(internalDisp[i3]);
			i2+= 2;
		}
	}
	
	return;
}

void Element::getNdVel(DiffDoub0 globVel[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	Node* nPtr;
	double ndVel[6];

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getVel(ndVel);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globVel[i3].setVal(ndVel[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdAcc(DiffDoub0 globAcc[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	Node* nPtr;
	double ndAcc[6];

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getAcc(ndAcc);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globAcc[i3].setVal(ndAcc[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdTemp(DiffDoub0 globTemp[], Node* ndAr[]) {
	int i1;
	double temp;
	Node* nPtr;

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		temp = nPtr->getTemperature();
		globTemp[i1].setVal(temp);
	}
	return;
}

void Element::getNdTdot(DiffDoub0 globTdot[], Node* ndAr[]) {
	int i1;
	double tdot;
	Node* nPtr;

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		tdot = nPtr->getTdot();
		globTdot[i1].setVal(tdot);
	}
	return;
}

void Element::evalN(DiffDoub0 nVec[], DiffDoub0 dNds[], double spt[]) {
	if(type == 4) {
		nVec[0].setVal(1.0-spt[0]-spt[1]-spt[2]);
		dNds[0].setVal(-1.0);
		dNds[1].setVal(-1.0);
		dNds[2].setVal(-1.0);
        
		nVec[1].setVal(spt[0]);
		dNds[3].setVal(1.0);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
        
	    nVec[2].setVal(spt[1]);
		dNds[6].setVal(0.0);
		dNds[7].setVal(1.0);
		dNds[8].setVal(0.0);
		
        nVec[3].setVal(spt[2]);
		dNds[9].setVal(0.0);
		dNds[10].setVal(0.0);
		dNds[11].setVal(1.0);
	} else if(type == 6) {
		nVec[0].setVal(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
		dNds[0].setVal(-0.5*(1.0-spt[2]));
	    dNds[1].setVal(-0.5*(1.0-spt[2]));
		dNds[2].setVal(-0.5*(1.0-spt[0]-spt[1]));
		
	    nVec[1].setVal(0.5*spt[0]*(1.0-spt[2]));
		dNds[3].setVal(0.5*(1.0-spt[2]));
		dNds[4].setVal(0.0);
		dNds[5].setVal(-0.5*spt[0]);
		
	    nVec[2].setVal(0.5*spt[1]*(1.0-spt[2]));
		dNds[6].setVal(0.0);
		dNds[7].setVal(0.5*(1.0-spt[2]));
		dNds[8].setVal(-0.5*spt[1]);
		
		nVec[3].setVal(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
		dNds[9].setVal(-0.5*(1.0+spt[2]));
		dNds[10].setVal(-0.5*(1.0+spt[2]));
		dNds[11].setVal(0.5*(1.0-spt[0]-spt[1]));
		
	    nVec[4].setVal(0.5*spt[0]*(1.0+spt[2]));
		dNds[12].setVal(0.5*(1.0+spt[2]));
		dNds[13].setVal(0.0);
		dNds[14].setVal(0.5*spt[0]);
		
	    nVec[5].setVal(0.5*spt[1]*(1.0+spt[2]));
		dNds[15].setVal(0.0);
		dNds[16].setVal(0.5*(1.0+spt[2]));
		dNds[17].setVal(0.5*spt[1]);
	} else if(type == 8 || type == 81) {
		nVec[0].setVal(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
		dNds[0].setVal(-0.125*(1.0-spt[1])*(1.0-spt[2]));
		dNds[1].setVal(-0.125*(1.0-spt[0])*(1.0-spt[2]));
		dNds[2].setVal(-0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        nVec[1].setVal(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
		dNds[3].setVal(0.125*(1.0-spt[1])*(1.0-spt[2]));
		dNds[4].setVal(-0.125*(1.0+spt[0])*(1.0-spt[2]));
		dNds[5].setVal(-0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        nVec[2].setVal(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		dNds[6].setVal(0.125*(1.0+spt[1])*(1.0-spt[2]));
		dNds[7].setVal(0.125*(1.0+spt[0])*(1.0-spt[2]));
		dNds[8].setVal(-0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        nVec[3].setVal(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		dNds[9].setVal(-0.125*(1.0+spt[1])*(1.0-spt[2]));
		dNds[10].setVal(0.125*(1.0-spt[0])*(1.0-spt[2]));
		dNds[11].setVal(-0.125*(1.0-spt[0])*(1.0+spt[1]));
		
        nVec[4].setVal(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		dNds[12].setVal(-0.125*(1.0-spt[1])*(1.0+spt[2]));
		dNds[13].setVal(-0.125*(1.0-spt[0])*(1.0+spt[2]));
		dNds[14].setVal(0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        nVec[5].setVal(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		dNds[15].setVal(0.125*(1.0-spt[1])*(1.0+spt[2]));
		dNds[16].setVal(-0.125*(1.0+spt[0])*(1.0+spt[2]));
		dNds[17].setVal(0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        nVec[6].setVal(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		dNds[18].setVal(0.125*(1.0+spt[1])*(1.0+spt[2]));
		dNds[19].setVal(0.125*(1.0+spt[0])*(1.0+spt[2]));
		dNds[20].setVal(0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        nVec[7].setVal(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		dNds[21].setVal(-0.125*(1.0+spt[1])*(1.0+spt[2]));
		dNds[22].setVal(0.125*(1.0-spt[0])*(1.0+spt[2]));
		dNds[23].setVal(0.125*(1.0-spt[0])*(1.0+spt[1]));
		
		if(type == 81) {
			nVec[8].setVal(1.0-spt[0]*spt[0]);
			dNds[24].setVal(-2.0*spt[0]);
		    dNds[25].setVal(0.0);
		    dNds[26].setVal(0.0);
			
            nVec[9].setVal(1.0-spt[1]*spt[1]);
			dNds[27].setVal(0.0);
		    dNds[28].setVal(-2.0*spt[1]);
		    dNds[29].setVal(0.0);
			
            nVec[10].setVal(1.0-spt[2]*spt[2]);
			dNds[30].setVal(0.0);
		    dNds[31].setVal(0.0);
		    dNds[32].setVal(-2.0*spt[2]);
		}
	}
	else if (type == 10) {
		double p1 = 1.0 - spt[0] - spt[1] - spt[2];
		double p2 = p1 - 0.5;
		nVec[0].setVal(2.0 * p1 * p2);
		dNds[0].setVal(2.0 * (-p2 - p1));
		dNds[1].setVal(2.0 * (-p2 - p1));
		dNds[2].setVal(2.0 * (-p2 - p1));

		nVec[1].setVal(-2.0 * spt[0] * (0.5 - spt[0]));
		dNds[3].setVal(-2.0 * (0.5 - spt[0] - spt[0]));
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
		nVec[2].setVal(-2.0 * spt[1] * (0.5 - spt[1]));
		dNds[6].setVal(0.0);
		dNds[7].setVal(-2.0 * (0.5 - spt[1] - spt[1]));
		dNds[8].setVal(0.0);
		
		nVec[3].setVal(-2.0 * spt[2] * (0.5 - spt[2]));
		dNds[9].setVal(0.0);
		dNds[10].setVal(0.0);
		dNds[11].setVal(-2.0 * (0.5 - spt[2] - spt[2]));
		
		nVec[4].setVal(4.0 * spt[0] * p1);
		dNds[12].setVal(4.0 * (p1 - spt[0]));
		dNds[13].setVal(-4.0 * spt[0]);
		dNds[14].setVal(-4.0 * spt[0]);
		
		nVec[5].setVal(4.0 * spt[0] * spt[1]);
		dNds[15].setVal(4.0 * spt[1]);
		dNds[16].setVal(4.0 * spt[0]);
		dNds[17].setVal(0.0);
		
		nVec[6].setVal(4.0 * spt[1] * p1);
		dNds[18].setVal(-4.0 * spt[1]);
		dNds[19].setVal(4.0 * (p1 - spt[1]));
		dNds[20].setVal(-4.0 * spt[1]);
		
		nVec[7].setVal(4.0 * spt[2] * p1);
		dNds[21].setVal(-4.0 * spt[2]);
		dNds[22].setVal(-4.0 * spt[2]);
		dNds[23].setVal(4.0 * (p1 - spt[2]));
		
		nVec[8].setVal(4.0 * spt[0] * spt[2]);
		dNds[24].setVal(4.0 * spt[2]);
		dNds[25].setVal(0.0);
		dNds[26].setVal(4.0 * spt[0]);
		
		nVec[9].setVal(4.0 * spt[1] * spt[2]);
		dNds[27].setVal(0.0);
		dNds[28].setVal(4.0 * spt[2]);
		dNds[29].setVal(4.0 * spt[1]);
	}
	else if (type == 3) {
		nVec[0].setVal(1.0-spt[0]-spt[1]);
		dNds[0].setVal(-1.0);
		dNds[1].setVal(-1.0);
		dNds[2].setVal(0.0);
		
		nVec[1].setVal(spt[0]);
		dNds[3].setVal(1.0);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
		nVec[2].setVal(spt[1]);
		dNds[6].setVal(0.0);
		dNds[7].setVal(1.0);
		dNds[8].setVal(0.0);
		
		nVec[3].setVal(spt[0]*(1.0-spt[0]-spt[1]));
		dNds[9].setVal(1.0-spt[0]-spt[1] - spt[0]);
		dNds[10].setVal(-spt[0]);
		dNds[11].setVal(0.0);
		
		nVec[4].setVal(spt[0]*spt[1]);
		dNds[12].setVal(spt[1]);
		dNds[13].setVal(spt[0]);
		dNds[14].setVal(0.0);
		
		nVec[5].setVal(spt[1]*(1.0-spt[0]-spt[1]));
		dNds[15].setVal(-spt[1]);
		dNds[16].setVal(1.0-spt[0]-spt[1]-spt[1]);
		dNds[17].setVal(0.0);
		
	} else if(type == 41) {
		nVec[0].setVal(0.25*(1.0-spt[0])*(1.0-spt[1]));
		dNds[0].setVal(-0.25*(1.0-spt[1]));
		dNds[1].setVal(-0.25*(1.0-spt[0]));
		dNds[2].setVal(0.0);
		
        nVec[1].setVal(0.25*(1.0+spt[0])*(1.0-spt[1]));
		dNds[3].setVal(0.25*(1.0-spt[1]));
		dNds[4].setVal(-0.25*(1.0+spt[0]));
		dNds[5].setVal(0.0);
		
        nVec[2].setVal(0.25*(1.0+spt[0])*(1.0+spt[1]));
		dNds[6].setVal(0.25*(1.0+spt[1]));
		dNds[7].setVal(0.25*(1.0+spt[0]));
		dNds[8].setVal(0.0);
		
        nVec[3].setVal(0.25*(1.0-spt[0])*(1.0+spt[1]));
		dNds[9].setVal(-0.25*(1.0+spt[1]));
		dNds[10].setVal(0.25*(1.0-spt[0]));
		dNds[11].setVal(0.0);
		
		nVec[4].setVal(1.0-spt[0]*spt[0]);
		dNds[12].setVal(-2.0*spt[0]);
		dNds[13].setVal(0.0);
		dNds[14].setVal(0.0);
		
        nVec[5].setVal(1.0-spt[1]*spt[1]);
		dNds[15].setVal(0.0);
		dNds[16].setVal(-2.0*spt[1]);
		dNds[17].setVal(0.0);
		
	    nVec[6].setVal((1.0-spt[0]*spt[0])*(1.0-spt[1]));
		dNds[18].setVal(-2.0*spt[0]*(1.0-spt[1]));
		dNds[19].setVal(-1.0+spt[0]*spt[0]);
		dNds[20].setVal(0.0);
		
	    nVec[7].setVal((1.0-spt[1]*spt[1])*(1.0+spt[0]));
		dNds[21].setVal(1.0-spt[1]*spt[1]);
		dNds[22].setVal(-2.0*spt[1]*(1.0+spt[0]));
		dNds[23].setVal(0.0);
		
	    nVec[8].setVal((1.0-spt[0]*spt[0])*(1.0+spt[1]));
		dNds[24].setVal(-2.0*spt[0]*(1.0+spt[1]));
		dNds[25].setVal(1.0-spt[0]*spt[0]);
		dNds[26].setVal(0.0);
		
	    nVec[9].setVal((1.0-spt[1]*spt[1])*(1.0-spt[0]));
		dNds[27].setVal(-1.0+spt[1]*spt[1]);
		dNds[28].setVal(-2.0*spt[1]*(1.0-spt[0]));
		dNds[29].setVal(0.0);
		
	} else if(type == 2) {
		nVec[0].setVal(0.5*(1.0 - spt[0]));
		dNds[0].setVal(-0.5);
		dNds[1].setVal(0.0);
		dNds[2].setVal(0.0);
		
        nVec[1].setVal(0.5*(1.0 + spt[0]));
		dNds[3].setVal(0.5);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
        nVec[2].setVal(1.0 - spt[0]*spt[0]);
		dNds[0].setVal(-2.0*spt[0]);
		dNds[1].setVal(0.0);
		dNds[2].setVal(0.0);
	}
	return;
}

void Element::getIpData(DiffDoub0 nVec[], DiffDoub0 dNdx[], DiffDoub0& detJ, DiffDoub0 locNds[], double spt[]) {
	int i1;
	int i2;
	DiffDoub0 nCent[11];
	DiffDoub0 dNds[33];
	DiffDoub0 dNdsCent[33];
	DiffDoub0 jMat[9];
	DiffDoub0 jCent[9];
	DiffDoub0 detCent;
	DiffDoub0 jInv[9];
	DiffDoub0 jInvCent[9];
	double sCent[3] = { 0.0,0.0,0.0 };
	DiffDoub0 xVec[3];
	DiffDoub0 bVec[3];
	DiffDoub0 zDir;
	DiffDoub0 tmp;
	
	evalN(nVec, dNds, spt);
	matMul(jMat,locNds,dNds,3,numNds,3);

	if(type == 41 || type == 3) {
		zDir.setVal(jMat[0]);
		zDir.mult(jMat[4]);
		tmp.setVal(jMat[3]);
		tmp.mult(jMat[1]);
		zDir.sub(tmp);
		if(zDir.val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0);
		}
	} else if(type == 2) {
		jMat[4].setVal(1.0);
		if(jMat[0].val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0);
		}
	}
	
	getDetInv(detJ,jInv,jMat,3,0,xVec,bVec);
	
	// matMul(dNds,dNds,jInv,nDim,3,3);
	matMul(dNdx,dNds,jInv,numNds,3,3);

	if (nDim > numNds) {
		evalN(nCent, dNdsCent, sCent);
		matMul(jCent, locNds, dNdsCent, 3, numNds, 3);

		if (type == 41 || type == 3) {
			zDir.setVal(jCent[0]);
			zDir.mult(jCent[4]);
			tmp.setVal(jCent[3]);
			tmp.mult(jCent[1]);
			zDir.sub(tmp);
			if (zDir.val > 0.0) {
				jCent[8].setVal(1.0);
			}
			else {
				jCent[8].setVal(-1.0);
			}
		}
		else if (type == 2) {
			jCent[4].setVal(1.0);
			if (jCent[0].val > 0.0) {
				jCent[8].setVal(1.0);
			}
			else {
				jCent[8].setVal(-1.0);
			}
		}

		getDetInv(detCent, jInvCent, jCent, 3, 0, xVec, bVec);

		i1 = 3 * numNds;
		i2 = nDim - numNds;
		matMul(&dNdx[i1], &dNds[i1], jInvCent, i2, 3, 3);
	}
	
	return;
}

void Element::getInstOri(DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 globDisp[], int stat) {
	// stat = 1: nonlinear geometry, element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub0 version, 1st order for DiffDoub1 version
	DiffDoub0 rot[3];
	DiffDoub0 nnds;
	DiffDoub0 one;
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int stIndex;
	int iOriSize = (numNds+1)*144;
	for (i1 = 0; i1 < iOriSize; i1++) {
		instOriMat[i1].setVal(0.0);
	}
	bool isDiff = locOri[0].diffType();
	
	nnds.setVal(0.0);
	one.setVal(1.0);
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
	rot[2].setVal(0.0);
	i2 = 3 * nDim;
	for (i1 = 0; i1 < numNds; i1++) {
		rot[0].add(globDisp[i2]);
		rot[1].add(globDisp[i2 + nDim]);
		rot[2].add(globDisp[i2 + 2 * nDim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);
	
	if(stat == 1) {
		dOridThet(&instOriMat[0], locOri, rot, 0, 0);
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			stIndex = 144 * i1;
			dOridThet(&instOriMat[stIndex], locOri, rot, 0, 0);
			for (i2 = 1; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else if(isDiff) {
		for (i2 = 0; i2 < 4; i2++) {
			stIndex = 36*i2;
			dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
		}
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else {
		for (i2 = 0; i2 < 4; i2++) {
			stIndex = 36*i2;
			dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				stIndex = 36*i2 + 9*i6;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,i6);
				i3 = stIndex;
				i4 = 36*i6 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					stIndex = 144*i1 + 36*i2 + 9*i6;
					dOridThet(&instOriMat[stIndex],locOri,rot,i2,i6);
					i3 = stIndex;
					i4 = 144*i1 + 36*i6 + 9*i2;
					for (i5 = 0; i5 < 9; i5++) {
						instOriMat[i4].setVal(instOriMat[i3]);
						i3++;
						i4++;
					}
				}
			}
		}
	}

	return;
}

void Element::getInstDisp(DiffDoub0 instDisp[], DiffDoub0 globDisp[], DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 xGlob[], bool nLGeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int ndOInd;
	int dof;
	int dofOInd;
	int nd2;
	int dof2;
	int dof2OInd;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 nnInv;
	DiffDoub0 nnInv2;
	
	i2 = 6*nDim;
	for (i1 = 0; i1 < i2; i1++) {
		instDisp[i1].setVal(0.0);
	}
	
	if (!nLGeom) {
		if (dv1 + dv2 == -2) {
			i7 = 3 * nDim;
			for (i1 = 0; i1 < 3; i1++) {
				i4 = i1 * nDim;
				for (i2 = 0; i2 < numNds; i2++) {
					i5 = i1 * 3;
					i6 = i2;
					for (i3 = 0; i3 < 3; i3++) {
						//i4 = i1 * nDim + i2;
						//i5 = i1 * 3 + i3;
						//i6 = i3 * nDim + i2;
						tmp.setVal(locOri[i5]);
						tmp.mult(globDisp[i6]);
						instDisp[i4].add(tmp);
						tmp.setVal(locOri[i5]);
						tmp.mult(globDisp[i6 + i7]);
						instDisp[i4 + i7].add(tmp);
						i5++;
						i6 += nDim;
					}
					i4++;
				}
			}
			i2 = 2 * numNds * dofPerNd;
			for (i1 = 0; i1 < numIntDof; i1++) {
				nd = dofTable[i2];
				dof = dofTable[i2 + 1];
				i3 = dof * nDim + nd;
				instDisp[i3].setVal(globDisp[i3]);
				i2 += 2;
			}
		}
		else if ((dv1+1)*(dv2+1) == 0) {
			dv1 = dv1 + dv2 + 1;
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < numNds) {
					i2 = nd;
					i3 = dof;
					for (i1 = 0; i1 < 3; i1++) {
						//i2 = i1 * nDim + nd;
						//i3 = i1 * 3 + dof;
						instDisp[i2].setVal(locOri[i3]);
						i2 += nDim;
						i3 += 3;
					}
				}
				else {
					i1 = dof * nDim + nd;
					instDisp[i1].setVal(1.0);
				}
			}
			else {
				i2 = 3 * nDim + nd;
				i3 = dof - 3;
				for (i1 = 0; i1 < 3; i1++) {
					//i2 = (i1 + 3) * nDim + nd;
					//i3 = i1 * 3 + (dof - 3);
					instDisp[i2].setVal(locOri[i3]);
					i2 += nDim;
					i3 += 3;
				}
			}
		}
	}
	else {
		if (dv1 + dv2 == -2) {
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = i1 * nDim + i2;
					i5 = i2;
					i6 = i2;
					i7 = i1 * 3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp2.setVal(locOri[i7]);
						tmp2.mult(xGlob[i6]);
						tmp.sub(tmp2);
						instDisp[i4].add(tmp);
						i5 += nDim;
						i6 += numNds;
						i7++;
					}
				}
			}

			i2 = numNds * dofPerNd;
			i3 = i2 + numIntDof;
			for (i1 = i2; i1 < i3; i1++) {
				i4 = dofTable[2 * i1];
				i5 = dofTable[2 * i1 + 1];
				i6 = i5 * nDim + i4;
				instDisp[i6].setVal(globDisp[i6]);
			}

			for (i1 = 0; i1 < numNds; i1++) {
				i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
				i4 = 4 * nDim + i1;
				i5 = 5 * nDim + i1;
				ndOInd = 144 * (i1 + 1);
				for (i2 = 0; i2 < 3; i2++) {
					i6 = ndOInd + 3 + i2;
					i7 = 6 + i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i3].add(tmp);
					i6 = ndOInd + 6 + i2;
					i7 = i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i4].add(tmp);
					i6 = ndOInd + i2;
					i7 = 3 + i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i5].add(tmp);
				}
			}
		}
		else if ((dv1 + 1) * (dv2 + 1) == 0) {
			dv1 = dv1 + dv2 + 1;
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < numNds) {
					i1 = nd; // Index in instDisp
					i2 = dof; //Index in instOriMat
					instDisp[i1].setVal(instOriMat[i2]);
					i1 = nDim + nd;
					i2 = 3 + dof;
					instDisp[i1].setVal(instOriMat[i2]);
					i1 = 2 * nDim + nd;
					i2 = 6 + dof;
					instDisp[i1].setVal(instOriMat[i2]);
				}
				else {
					i1 = dof * nDim + nd;
					instDisp[i1].setVal(1.0);
				}
			}
			else { // dof is rotation
				nnInv.setVal(1.0 / numNds);
				dofOInd = 36 * (dof - 2);
				for (i1 = 0; i1 < 3; i1++) {
					for (i2 = 0; i2 < numNds; i2++) {
						i4 = i1 * nDim + i2;
						i5 = i2;
						i6 = i2;
						i7 = dofOInd + i1 * 3;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.setVal(globDisp[i5]);
							tmp.add(xGlob[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i5 += nDim;
							i6 += numNds;
							i7++;
						}
					}
				}

				for (i1 = 0; i1 < numNds; i1++) {
					i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
					i4 = 4 * nDim + i1;
					i5 = 5 * nDim + i1;
					ndOInd = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = ndOInd + 3 + i2;
						i7 = dofOInd + 6 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i3].add(tmp);
						i6 = ndOInd + 6 + i2;
						i7 = dofOInd + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i4].add(tmp);
						i6 = ndOInd + i2;
						i7 = dofOInd + 3 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i5].add(tmp);
						if (i1 == nd) {
							i6 = ndOInd + dofOInd + 3 + i2;
							i7 = 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + 6 + i2;
							i7 = i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + i2;
							i7 = 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i5].add(tmp);
						}
					}
				}
			}
		}
		else {
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			nd2 = dofTable[2 * dv2];
			dof2 = dofTable[2 * dv2 + 1];
			nnInv.setVal(1.0 / numNds);
			nnInv2.setVal(nnInv);
			nnInv2.sqr();
			if (dof > 2 && dof2 > 2) {
				for (i1 = 0; i1 < 3; i1++) {
					dofOInd = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
					for (i2 = 0; i2 < numNds; i2++) {
						i4 = nDim * i1 + i2;
						i5 = i2;
						i6 = i2;
						i7 = dofOInd;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.setVal(globDisp[i5]);
							tmp.add(xGlob[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv2);
							instDisp[i4].add(tmp);
							i5 += nDim;
							i6 += numNds;
							i7++;
						}
					}
				}

				dofOInd = 36 * (dof - 2);
				dof2OInd = 9 * (dof2 - 2);
				for (i1 = 0; i1 < numNds; i1++) {
					i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
					i4 = 4 * nDim + i1;
					i5 = 5 * nDim + i1;
					ndOInd = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = ndOInd + 3 + i2;
						i7 = dofOInd + dof2OInd + 6 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i3].add(tmp);
						i6 = ndOInd + 6 + i2;
						i7 = dofOInd + dof2OInd + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i4].add(tmp);
						i6 = ndOInd + i2;
						i7 = dofOInd + dof2OInd + 3 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i5].add(tmp);
						if (i1 == nd) {
							i6 = ndOInd + dofOInd + 3 + i2;
							i7 = dof2OInd + 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + 6 + i2;
							i7 = dof2OInd + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + i2;
							i7 = dof2OInd + 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i5].add(tmp);
						}
						if (i1 == nd2) {
							i6 = ndOInd + dof2OInd + 3 + i2;
							i7 = dofOInd + 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dof2OInd + 6 + i2;
							i7 = dofOInd + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dof2OInd + i2;
							i7 = dofOInd + 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i5].add(tmp);
						}
						if (i1 == nd && i1 == nd2) {
							i6 = ndOInd + dofOInd + dof2OInd + 3 + i2;
							i7 = 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + dof2OInd + 6 + i2;
							i7 = i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + dof2OInd + i2;
							i7 = 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i5].add(tmp);
						}
					}
				}
			}
			else if (dof < 3 && dof2 < 3) {
				return;
			}
			else {
				if (dof > 2) {
					i1 = dof;
					dof = dof2;
					dof2 = i1;
					i1 = nd;
					nd = nd2;
					nd2 = i1;
				}
				i1 = 36 * (dof2 - 2) + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = nd;
				instDisp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 3 + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = nDim + nd;
				instDisp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 6 + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = 2 * nDim + nd;
				instDisp[i2].add(tmp);
			}
		}
	}

	
	return;
}

void Element::getStressPrereq(DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int numLay;
	DiffDoub0 offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	getNdVel(pre.globVel, ndAr);
	getNdAcc(pre.globAcc, ndAr);
	getNdTemp(pre.globTemp, ndAr);
	getNdTdot(pre.globTdot, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		if (type != 2) {
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, dvAr);
			getLayerAngle(pre.layerAng, dvAr);
			getLayerQ(pre.layerQ, dvAr);
			getLayerD(pre.layerD, dvAr);
			getLayerThExp(pre.layerTE, dvAr);
			getLayerEinit(pre.layerE0, dvAr);
			getLayerDen(pre.layerDen, dvAr);
			getLayerCond(pre.layerTC, dvAr);
			getLayerSpecHeat(pre.layerSH, dvAr);
			getABD(pre.Cmat, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerAng);
			getShellDamp(pre.Dmat, pre.layerThk, pre.layerZ, pre.layerD, pre.layerAng);
			getShellExpLoad(pre.thermExp, pre.Einit, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerTE, pre.layerE0, pre.layerAng);
			getShellMass(pre.Mmat, pre.layerThk, pre.layerZ, pre.layerDen, dvAr);
			getShellCond(pre.TCmat, pre.layerThk, pre.layerAng, pre.layerTC, dvAr);
			getShellSpecHeat(pre.SpecHeat, pre.layerThk, pre.layerSH, pre.layerDen);
		}
		else {
			getBeamStiff(pre.Cmat, dvAr);
			getBeamDamp(pre.Dmat, dvAr);
			getBeamExpLoad(pre.thermExp, pre.Einit, dvAr);
			getBeamMass(pre.Mmat, dvAr);
			getBeamCond(pre.TCmat, dvAr);
			getBeamSpecHeat(pre.SpecHeat, dvAr);
		}
	}
	else if (type == 21) {
		getFrcFldConst(pre.frcFldCoef, pre.frcFldExp, dvAr);
	}
	else if (type == 1) {
		getMassPerEl(pre.massPerEl, dvAr);
	}
	else {
		getSolidStiff(pre.Cmat, dvAr);
		getSolidDamp(pre.Dmat, dvAr);
		getThermalExp(pre.thermExp, pre.Einit, dvAr);
		getDensity(pre.Mmat[0], 0, dvAr);
		getConductivity(pre.TCmat, dvAr);
		getSpecificHeat(pre.SpecHeat, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getVolume(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer) {
	int i1;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 thk;
	DiffDoub0 tmp;

	if (type == 2) {
		thk.setVal(sectPtr->getArea()); //rem: update to factor in dvars;
	} else if (type == 3 || type == 41) {
		thk.setVal(pre.layerThk[layer]);
	} else {
		thk.setVal(1.0);
	}

	vol.setVal(0.0);
	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		tmp.mult(detJ);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::getSectionDef(DiffDoub0 secDef[], DiffDoub0 globDisp[],  DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 xGlob[], DiffDoub0 dNdx[], DiffDoub0 nVec[], bool nLGeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub0 instDisp[60];
	DiffDoub0 ux[9];
	DiffDoub0 rx[9];
	DiffDoub0 rot[3];
	DiffDoub0 tmp;

	if (dofPerNd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			secDef[i1].setVal(0.0);
		}
		return;
	}
	
	getInstDisp(instDisp, globDisp, instOriMat, locOri, xGlob, nLGeom, dv1, dv2);
	
	if(dv1 + dv2 == -2) {
		matMul(ux,instDisp,dNdx,3,nDim,3);
		i1 = 3*nDim;
		matMul(rx,&instDisp[i1],dNdx,3,nDim,3);
		matMul(rot,&instDisp[i1],nVec,3,nDim,1);
	} else if((dv1+1)*(dv2+1) == 0) {
		dv1 = dv1 + dv2 + 1;
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		if(dof < 3) {
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].setVal(instDisp[i4]);
					ux[i3].mult(dNdx[i5]);
					rx[i3].setVal(0.0);
					i3++;
					i5++;
				}
				rot[i1].setVal(0.0);
				i4+= nDim;
			}
		} else {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = nDim*i1;
					i6 = i2;
					i7 = nDim*(i1+3);
					ux[i4].setVal(0.0);
					rx[i4].setVal(0.0);
					for (i3 = 0; i3 < numNds; i3++) {
						tmp.setVal(instDisp[i5]);
						tmp.mult(dNdx[i6]);
						ux[i4].add(tmp);
						tmp.setVal(instDisp[i7]);
						tmp.mult(dNdx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].setVal(0.0);
				i3 = nDim*(i1+3);
				for (i2 = 0; i2 < numNds; i2++) {
					tmp.setVal(instDisp[i3]);
					tmp.mult(nVec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		}
	} else {
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		nd2 = dofTable[2*dv2];
		dof2 = dofTable[2*dv2+1];
		if(dof < 3 && dof2 < 3) {
			for (i1 = 0; i1 < 9; i1++) {
				secDef[i1].setVal(0.0);
			}
			return;
		} else if(dof > 2 && dof2 > 2) {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = nDim*i1;
					i6 = i2;
					i7 = nDim*(i1+3);
					ux[i4].setVal(0.0);
					for (i3 = 0; i3 < numNds; i3++) {
						tmp.setVal(instDisp[i5]);
						tmp.mult(dNdx[i6]);
						ux[i4].add(tmp);
						tmp.setVal(instDisp[i7]);
						tmp.mult(dNdx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].setVal(0.0);
				i3 = nDim*(i1+3);
				for (i2 = 0; i2 < numNds; i2++) {
					tmp.setVal(instDisp[i3]);
					tmp.mult(nVec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		} else {
			if(dof > dof2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].setVal(instDisp[i4]);
					ux[i3].mult(dNdx[i5]);
					rx[i3].setVal(0.0);
					i3++;
					i5++;
				}
				rot[i1].setVal(0.0);
				i4+= nDim;
			}
		}
	}
	
	if(type == 2) {
		secDef[0].setVal(ux[0]);
		secDef[1].setVal(ux[3]);
		secDef[1].sub(rot[2]);
		secDef[2].setVal(ux[6]);
		secDef[2].add(rot[1]);
		secDef[3].setVal(rx[0]);
		secDef[4].setVal(rx[3]);
		secDef[5].setVal(rx[6]);
	} else {
		secDef[0].setVal(ux[0]);
		secDef[1].setVal(ux[4]);
		secDef[2].setVal(ux[1]);
		secDef[2].add(ux[3]);
		secDef[3].setVal(rx[3]);
		secDef[4].setVal(rx[1]);
		secDef[4].neg();
		secDef[5].setVal(rx[4]);
		secDef[5].sub(rx[0]);
		secDef[6].setVal(ux[6]);
		secDef[6].add(rot[1]);
		secDef[7].setVal(ux[7]);
		secDef[7].sub(rot[0]);
		secDef[8].setVal(2.0);
		secDef[8].mult(rot[2]);
		secDef[8].sub(ux[3]);
		secDef[8].add(ux[1]);
	}
	
	return;
}

void Element::getSolidStrain(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 dNdx[], DiffDoub0 locOri[], int dv1, int dv2, bool nLGeom) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub0 uxL[9];
	DiffDoub0 strnMat[9];
	DiffDoub0 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strnMat[i1].setVal(0.0);
	}
	if(dv1 + dv2 == -2) {
		matMul(uxL,locOri,ux,3,3,3);
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			i5 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				strnMat[i4].add(uxL[i4]);
				strnMat[i4].add(uxL[i5]);
				if(nLGeom) {
					i6 = i1;
					i7 = i2;
				    for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(ux[i6]);
						tmp.mult(ux[i7]);
						strnMat[i4].add(tmp);
						i6+= 3;
						i7+= 3;
				    }
				}
				i4++;
				i5+= 3;
			}
		}
	} else if((dv1+1)*(dv2+1) == 0) {
		if(dv1 > -1) {
			nd = dofTable[2*dv1];
			dof = dofTable[2*dv1+1];
		} else {
			nd = dofTable[2*dv2];
			dof = dofTable[2*dv2+1];
		}
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				i5 = 3*i1 + dof;
				i6 = 3*nd + i2;
				tmp.setVal(locOri[i5]);
				tmp.mult(dNdx[i6]);
				strnMat[i4].add(tmp);
				i5 = 3*i2 + dof;
				i6 = 3*nd + i1;
				tmp.setVal(locOri[i5]);
				tmp.mult(dNdx[i6]);
				strnMat[i4].add(tmp);
				if(nLGeom) {
					i5 = 3*nd + i1;
					i6 = 3*dof + i2;
					tmp.setVal(dNdx[i5]);
					tmp.mult(ux[i6]);
					strnMat[i4].add(tmp);
					i5 = 3*nd + i2;
					i6 = 3*dof + i1;
					tmp.setVal(dNdx[i5]);
					tmp.mult(ux[i6]);
					strnMat[i4].add(tmp);
				}
				i4++;
			}
		}
	} else {
		if(nLGeom) {
			nd = dofTable[2*dv1];
			dof = dofTable[2*dv1+1];
			nd2 = dofTable[2*dv2];
			dof2 = dofTable[2*dv2+1];
			if(dof == dof2) {
				for (i1 = 0; i1 < 3; i1++) {
					i4 = 4*i1;
					for (i2 = i1; i2 < 3; i2++) {
						i5 = 3*nd + i1;
						i6 = 3*nd2 + i2;
						tmp.setVal(dNdx[i5]);
						tmp.mult(dNdx[i6]);
						strnMat[i4].add(tmp);
						i5 = 3*nd2 + i1;
						i6 = 3*nd + i2;
						tmp.setVal(dNdx[i5]);
						tmp.mult(dNdx[i6]);
						strnMat[i4].add(tmp);
						i4++;
					}
				}
			}
		}
	}
	
	tmp.setVal(0.5);
	strain[0].setVal(strnMat[0]);
	strain[0].mult(tmp);
	strain[1].setVal(strnMat[4]);
	strain[1].mult(tmp);
	strain[2].setVal(strnMat[8]);
	strain[2].mult(tmp);
	strain[3].setVal(strnMat[1]);
	strain[4].setVal(strnMat[2]);
	strain[5].setVal(strnMat[5]);
	
	return;
}

void Element::getStressStrain(DiffDoub0 stress[], DiffDoub0 strain[], double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 ux[9];
	DiffDoub0 secDef[9];
	DiffDoub0 sectStrn[3];
	DiffDoub0 adjStn[6];
	DiffDoub0 tmp;
	DiffDoub0 ipTemp;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);

	ipTemp.setVal(0.0);
	for (i1 = 0; i1 < numNds; i1++) {
		tmp.setVal(pre.globTemp[i1]);
		tmp.mult(nVec[i1]);
		ipTemp.add(tmp);
	}

	if (type == 41 || type == 3) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}

		getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, -1, -1);
		
		sectStrn[0].setVal(secDef[0]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[3]);
		sectStrn[0].add(tmp);
		
		sectStrn[1].setVal(secDef[1]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[4]);
		sectStrn[1].sub(tmp);

		sectStrn[2].setVal(secDef[2]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[5]);
		sectStrn[2].add(tmp);

		tmp.setVal(pre.layerAng[layer]);
		tmp.neg();
		transformStrain(strain, sectStrn, tmp);
		i2 = 3 * layer;
		for (i1 = 0; i1 < 3; i1++) {
			tmp.setVal(pre.layerTE[i2]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.layerE0[i2]);
			i2++;
		}
		matMul(stress, &pre.layerQ[9 * layer], adjStn, 3, 3, 1);

		strain[3].setVal(strain[2]);
		strain[2].setVal(0.0);
		strain[4].setVal(0.0);
		strain[5].setVal(0.0);

		stress[3].setVal(stress[2]);
		stress[2].setVal(0.0);
		stress[4].setVal(0.0);
		stress[5].setVal(0.0);

	} else if(type != 2) {
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.setVal(pre.thermExp[i1]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.Einit[i1]);
		}
		matMul(stress, pre.Cmat, adjStn, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(DiffDoub0 dsdU[], DiffDoub0 dedU[], DiffDoub0 dsdT[], double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 ux[9];
	DiffDoub0 secDef[9];
	DiffDoub0 sectStrn[3];
	DiffDoub0 dStrain[6];
	DiffDoub0 dStress[6];
	DiffDoub0 CTE[6];
	DiffDoub0 CTEN[60];
	DiffDoub0 tmp;

	int totDof = numNds * dofPerNd + numIntDof;
	i2 = 6 * totDof;
	for (i1 = 0; i1 < i2; i1++) {
		dsdU[i1].setVal(0.0);
		dedU[i1].setVal(0.0);
	}
	i2 = 6 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		dsdT[i1].setVal(0.0);
	}

	if (type == 41 || type == 3) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}

		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		for (i1 = 0; i1 < totDof; i1++) {
			getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, i1, -1);

			sectStrn[0].setVal(secDef[0]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[3]);
			sectStrn[0].add(tmp);

			sectStrn[1].setVal(secDef[1]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[4]);
			sectStrn[1].sub(tmp);

			sectStrn[2].setVal(secDef[2]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[5]);
			sectStrn[2].add(tmp);

			tmp.setVal(pre.layerAng[layer]);
			tmp.neg();
			transformStrain(dStrain, sectStrn, tmp);
			matMul(dStress, &pre.layerQ[9 * layer], dStrain, 3, 3, 1);

			dedU[i1].setVal(dStrain[0]);
			dedU[i1 + totDof].setVal(dStrain[1]);
			dedU[i1 + 3 * totDof].setVal(dStrain[2]);

			dsdU[i1].setVal(dStress[0]);
			dsdU[totDof + i1].setVal(dStress[1]);
			dsdU[3 * totDof + i1].setVal(dStress[2]);
		}
		matMul(CTE, &pre.layerQ[9 * layer], &pre.layerTE[3 * layer], 3, 3, 1);
		CTE[0].neg();
		CTE[1].neg();
		CTE[2].neg();
		matMul(CTEN, CTE, nVec, 3, 1, numNds);
		for (i1 = 0; i1 < numNds; i1++) {
			dsdT[i1].setVal(CTEN[i1]);
			dsdT[i1 + numNds].setVal(CTEN[i1 + numNds]);
			dsdT[i1 + 3 * numNds].setVal(CTEN[i1 + 2 * numNds]);
		}
	}
	else if (type != 2) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		for (i1 = 0; i1 < totDof; i1++) {
			getSolidStrain(dStrain, ux, dNdx, pre.locOri, i1, -1, nLGeom);
			matMul(dStress, pre.Cmat, dStrain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				dedU[i3].setVal(dStrain[i2]);
				dsdU[i3].setVal(dStress[i2]);
				i3 += totDof;
			}
		}
		matMul(CTE, pre.Cmat, pre.thermExp, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			CTE[i1].neg();
		}
		matMul(dsdT, CTE, nVec, 6, 1, numNds);
	}

	return;
}

void Element::getDefFrcMom(DiffDoub0 def[], DiffDoub0 frcMom[], double spt[], bool nLGeom, DiffDoub0StressPrereq& pre) {
	int i1;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 ptTemp;
	DiffDoub0 tmp;

	if (dofPerNd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			def[i1].setVal(0.0);
		}
		return;
	}

	if (nLGeom) {
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
	}

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	getSectionDef(def, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, -1, -1);

	ptTemp.setVal(0.0);
	for (i1 = 0; i1 < numNds; i1++) {
		tmp.setVal(pre.globTemp[i1]);
		tmp.mult(nVec[i1]);
		ptTemp.add(tmp);
	}

	matMul(frcMom, pre.Cmat, def, defDim, defDim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frcMom[i1].sub(pre.Einit[i1]);
		tmp.setVal(ptTemp);
		tmp.mult(pre.thermExp[i1]);
		frcMom[i1].sub(tmp);
	}

	return;
}

void Element::dDefFrcMomdU(DiffDoub0 dDefdU[], DiffDoub0 dFrcMomdU[], DiffDoub0 dFrcMomdT[], double spt[], bool nLGeom, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	int totDof;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 ptTemp;
	DiffDoub0 tmp;
	DiffDoub0 def[9];

	if (nLGeom) {
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
	}


	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	totDof = numNds * dofPerNd + numIntDof;
	for (i1 = 0; i1 < totDof; i1++) {
		getSectionDef(def, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, i1, -1);
		i2 = i1;
		for (i3 = 0; i3 < defDim; i3++) {
			dDefdU[i2].setVal(def[i3]);
			i2 += totDof;
		}
	}

	matMul(dFrcMomdU, pre.Cmat, dDefdU, defDim, defDim, totDof);

	matMul(dFrcMomdT, pre.thermExp, nVec, 6, 1, numNds);

	return;
}

void Element::getFluxTGrad(DiffDoub0 flux[], DiffDoub0 tGrad[], double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	matMul(tGrad, pre.globTemp, dNdx, 1, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		matMul(flux, &pre.layerTC[i1], tGrad, 3, 3, 1);
	}
	else {
		matMul(flux, pre.TCmat, tGrad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::dFluxTGraddT(DiffDoub0 dFdT[], DiffDoub0 dTG[], double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	transpose(dTG, dNdx, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		matMul(dFdT, &pre.layerTC[i1], dTG, 3, 3, numNds);
	}
	else {
		matMul(dFdT, pre.TCmat, dTG, 3, 3, numNds);
	}
	i2 = 3 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		dFdT[i1].neg();
	}

	return;
}

void Element::putVecToGlobMat(SparseMat& qMat, DiffDoub0 elQVec[], bool forTherm, int matRow, Node* ndAr[]) {
	int i1;
	int ndDof = numNds*dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[i1]->getSortedRank();
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < totDof; i1++) {
			if (i1 < ndDof) {
				nd = nodes[dofTable[2 * i1]];
				dof = dofTable[2 * i1 + 1];
				globInd = ndAr[nd]->getDofIndex(dof);
				qMat.addEntry(matRow, globInd, elQVec[i1].val);
			}
			else {
				dof = i1 - ndDof;
				globInd = intDofIndex + dof;
				qMat.addEntry(matRow, globInd, elQVec[i1].val);
			}
		}
	}

	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void Element::getNdDisp(DiffDoub1 globDisp[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub1 ndDisp[6];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getDisp(ndDisp);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globDisp[i3].setVal(ndDisp[i2]);
			i3+= nDim;
		}
	}
	
	if(numIntDof > 0) {
		i2 = 2*numNds*dofPerNd;
		for (i3 = 0; i3 < numIntDof; i3++) {
			i4 = dofTable[i2];
			i5 = dofTable[i2+1];
			i6 = nDim*i5 + i4;
			globDisp[i6].setVal(internalDisp[i3]);
			i2+= 2;
		}
	}
	
	return;
}

void Element::getNdVel(DiffDoub1 globVel[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	Node* nPtr;
	double ndVel[6];

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getVel(ndVel);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globVel[i3].setVal(ndVel[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdAcc(DiffDoub1 globAcc[], Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	Node* nPtr;
	double ndAcc[6];

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getAcc(ndAcc);
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globAcc[i3].setVal(ndAcc[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdTemp(DiffDoub1 globTemp[], Node* ndAr[]) {
	int i1;
	double temp;
	Node* nPtr;

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		temp = nPtr->getTemperature();
		globTemp[i1].setVal(temp);
	}
	return;
}

void Element::getNdTdot(DiffDoub1 globTdot[], Node* ndAr[]) {
	int i1;
	double tdot;
	Node* nPtr;

	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		tdot = nPtr->getTdot();
		globTdot[i1].setVal(tdot);
	}
	return;
}

void Element::evalN(DiffDoub1 nVec[], DiffDoub1 dNds[], double spt[]) {
	if(type == 4) {
		nVec[0].setVal(1.0-spt[0]-spt[1]-spt[2]);
		dNds[0].setVal(-1.0);
		dNds[1].setVal(-1.0);
		dNds[2].setVal(-1.0);
        
		nVec[1].setVal(spt[0]);
		dNds[3].setVal(1.0);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
        
	    nVec[2].setVal(spt[1]);
		dNds[6].setVal(0.0);
		dNds[7].setVal(1.0);
		dNds[8].setVal(0.0);
		
        nVec[3].setVal(spt[2]);
		dNds[9].setVal(0.0);
		dNds[10].setVal(0.0);
		dNds[11].setVal(1.0);
	} else if(type == 6) {
		nVec[0].setVal(0.5*(1.0-spt[0]-spt[1])*(1.0-spt[2]));
		dNds[0].setVal(-0.5*(1.0-spt[2]));
	    dNds[1].setVal(-0.5*(1.0-spt[2]));
		dNds[2].setVal(-0.5*(1.0-spt[0]-spt[1]));
		
	    nVec[1].setVal(0.5*spt[0]*(1.0-spt[2]));
		dNds[3].setVal(0.5*(1.0-spt[2]));
		dNds[4].setVal(0.0);
		dNds[5].setVal(-0.5*spt[0]);
		
	    nVec[2].setVal(0.5*spt[1]*(1.0-spt[2]));
		dNds[6].setVal(0.0);
		dNds[7].setVal(0.5*(1.0-spt[2]));
		dNds[8].setVal(-0.5*spt[1]);
		
		nVec[3].setVal(0.5 * (1.0 - spt[0] - spt[1]) * (1.0 + spt[2]));
		dNds[9].setVal(-0.5*(1.0+spt[2]));
		dNds[10].setVal(-0.5*(1.0+spt[2]));
		dNds[11].setVal(0.5*(1.0-spt[0]-spt[1]));
		
	    nVec[4].setVal(0.5*spt[0]*(1.0+spt[2]));
		dNds[12].setVal(0.5*(1.0+spt[2]));
		dNds[13].setVal(0.0);
		dNds[14].setVal(0.5*spt[0]);
		
	    nVec[5].setVal(0.5*spt[1]*(1.0+spt[2]));
		dNds[15].setVal(0.0);
		dNds[16].setVal(0.5*(1.0+spt[2]));
		dNds[17].setVal(0.5*spt[1]);
	} else if(type == 8 || type == 81) {
		nVec[0].setVal(0.125 * (1.0 - spt[0]) * (1.0 - spt[1]) * (1.0 - spt[2]));
		dNds[0].setVal(-0.125*(1.0-spt[1])*(1.0-spt[2]));
		dNds[1].setVal(-0.125*(1.0-spt[0])*(1.0-spt[2]));
		dNds[2].setVal(-0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        nVec[1].setVal(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0-spt[2]));
		dNds[3].setVal(0.125*(1.0-spt[1])*(1.0-spt[2]));
		dNds[4].setVal(-0.125*(1.0+spt[0])*(1.0-spt[2]));
		dNds[5].setVal(-0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        nVec[2].setVal(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		dNds[6].setVal(0.125*(1.0+spt[1])*(1.0-spt[2]));
		dNds[7].setVal(0.125*(1.0+spt[0])*(1.0-spt[2]));
		dNds[8].setVal(-0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        nVec[3].setVal(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0-spt[2]));
		dNds[9].setVal(-0.125*(1.0+spt[1])*(1.0-spt[2]));
		dNds[10].setVal(0.125*(1.0-spt[0])*(1.0-spt[2]));
		dNds[11].setVal(-0.125*(1.0-spt[0])*(1.0+spt[1]));
		
        nVec[4].setVal(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		dNds[12].setVal(-0.125*(1.0-spt[1])*(1.0+spt[2]));
		dNds[13].setVal(-0.125*(1.0-spt[0])*(1.0+spt[2]));
		dNds[14].setVal(0.125*(1.0-spt[0])*(1.0-spt[1]));
		
        nVec[5].setVal(0.125*(1.0+spt[0])*(1.0-spt[1])*(1.0+spt[2]));
		dNds[15].setVal(0.125*(1.0-spt[1])*(1.0+spt[2]));
		dNds[16].setVal(-0.125*(1.0+spt[0])*(1.0+spt[2]));
		dNds[17].setVal(0.125*(1.0+spt[0])*(1.0-spt[1]));
		
        nVec[6].setVal(0.125*(1.0+spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		dNds[18].setVal(0.125*(1.0+spt[1])*(1.0+spt[2]));
		dNds[19].setVal(0.125*(1.0+spt[0])*(1.0+spt[2]));
		dNds[20].setVal(0.125*(1.0+spt[0])*(1.0+spt[1]));
		
        nVec[7].setVal(0.125*(1.0-spt[0])*(1.0+spt[1])*(1.0+spt[2]));
		dNds[21].setVal(-0.125*(1.0+spt[1])*(1.0+spt[2]));
		dNds[22].setVal(0.125*(1.0-spt[0])*(1.0+spt[2]));
		dNds[23].setVal(0.125*(1.0-spt[0])*(1.0+spt[1]));
		
		if(type == 81) {
			nVec[8].setVal(1.0-spt[0]*spt[0]);
			dNds[24].setVal(-2.0*spt[0]);
		    dNds[25].setVal(0.0);
		    dNds[26].setVal(0.0);
			
            nVec[9].setVal(1.0-spt[1]*spt[1]);
			dNds[27].setVal(0.0);
		    dNds[28].setVal(-2.0*spt[1]);
		    dNds[29].setVal(0.0);
			
            nVec[10].setVal(1.0-spt[2]*spt[2]);
			dNds[30].setVal(0.0);
		    dNds[31].setVal(0.0);
		    dNds[32].setVal(-2.0*spt[2]);
		}
	}
	else if (type == 10) {
		double p1 = 1.0 - spt[0] - spt[1] - spt[2];
		double p2 = p1 - 0.5;
		nVec[0].setVal(2.0 * p1 * p2);
		dNds[0].setVal(2.0 * (-p2 - p1));
		dNds[1].setVal(2.0 * (-p2 - p1));
		dNds[2].setVal(2.0 * (-p2 - p1));

		nVec[1].setVal(-2.0 * spt[0] * (0.5 - spt[0]));
		dNds[3].setVal(-2.0 * (0.5 - spt[0] - spt[0]));
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
		nVec[2].setVal(-2.0 * spt[1] * (0.5 - spt[1]));
		dNds[6].setVal(0.0);
		dNds[7].setVal(-2.0 * (0.5 - spt[1] - spt[1]));
		dNds[8].setVal(0.0);
		
		nVec[3].setVal(-2.0 * spt[2] * (0.5 - spt[2]));
		dNds[9].setVal(0.0);
		dNds[10].setVal(0.0);
		dNds[11].setVal(-2.0 * (0.5 - spt[2] - spt[2]));
		
		nVec[4].setVal(4.0 * spt[0] * p1);
		dNds[12].setVal(4.0 * (p1 - spt[0]));
		dNds[13].setVal(-4.0 * spt[0]);
		dNds[14].setVal(-4.0 * spt[0]);
		
		nVec[5].setVal(4.0 * spt[0] * spt[1]);
		dNds[15].setVal(4.0 * spt[1]);
		dNds[16].setVal(4.0 * spt[0]);
		dNds[17].setVal(0.0);
		
		nVec[6].setVal(4.0 * spt[1] * p1);
		dNds[18].setVal(-4.0 * spt[1]);
		dNds[19].setVal(4.0 * (p1 - spt[1]));
		dNds[20].setVal(-4.0 * spt[1]);
		
		nVec[7].setVal(4.0 * spt[2] * p1);
		dNds[21].setVal(-4.0 * spt[2]);
		dNds[22].setVal(-4.0 * spt[2]);
		dNds[23].setVal(4.0 * (p1 - spt[2]));
		
		nVec[8].setVal(4.0 * spt[0] * spt[2]);
		dNds[24].setVal(4.0 * spt[2]);
		dNds[25].setVal(0.0);
		dNds[26].setVal(4.0 * spt[0]);
		
		nVec[9].setVal(4.0 * spt[1] * spt[2]);
		dNds[27].setVal(0.0);
		dNds[28].setVal(4.0 * spt[2]);
		dNds[29].setVal(4.0 * spt[1]);
	}
	else if (type == 3) {
		nVec[0].setVal(1.0-spt[0]-spt[1]);
		dNds[0].setVal(-1.0);
		dNds[1].setVal(-1.0);
		dNds[2].setVal(0.0);
		
		nVec[1].setVal(spt[0]);
		dNds[3].setVal(1.0);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
		nVec[2].setVal(spt[1]);
		dNds[6].setVal(0.0);
		dNds[7].setVal(1.0);
		dNds[8].setVal(0.0);
		
		nVec[3].setVal(spt[0]*(1.0-spt[0]-spt[1]));
		dNds[9].setVal(1.0-spt[0]-spt[1] - spt[0]);
		dNds[10].setVal(-spt[0]);
		dNds[11].setVal(0.0);
		
		nVec[4].setVal(spt[0]*spt[1]);
		dNds[12].setVal(spt[1]);
		dNds[13].setVal(spt[0]);
		dNds[14].setVal(0.0);
		
		nVec[5].setVal(spt[1]*(1.0-spt[0]-spt[1]));
		dNds[15].setVal(-spt[1]);
		dNds[16].setVal(1.0-spt[0]-spt[1]-spt[1]);
		dNds[17].setVal(0.0);
		
	} else if(type == 41) {
		nVec[0].setVal(0.25*(1.0-spt[0])*(1.0-spt[1]));
		dNds[0].setVal(-0.25*(1.0-spt[1]));
		dNds[1].setVal(-0.25*(1.0-spt[0]));
		dNds[2].setVal(0.0);
		
        nVec[1].setVal(0.25*(1.0+spt[0])*(1.0-spt[1]));
		dNds[3].setVal(0.25*(1.0-spt[1]));
		dNds[4].setVal(-0.25*(1.0+spt[0]));
		dNds[5].setVal(0.0);
		
        nVec[2].setVal(0.25*(1.0+spt[0])*(1.0+spt[1]));
		dNds[6].setVal(0.25*(1.0+spt[1]));
		dNds[7].setVal(0.25*(1.0+spt[0]));
		dNds[8].setVal(0.0);
		
        nVec[3].setVal(0.25*(1.0-spt[0])*(1.0+spt[1]));
		dNds[9].setVal(-0.25*(1.0+spt[1]));
		dNds[10].setVal(0.25*(1.0-spt[0]));
		dNds[11].setVal(0.0);
		
		nVec[4].setVal(1.0-spt[0]*spt[0]);
		dNds[12].setVal(-2.0*spt[0]);
		dNds[13].setVal(0.0);
		dNds[14].setVal(0.0);
		
        nVec[5].setVal(1.0-spt[1]*spt[1]);
		dNds[15].setVal(0.0);
		dNds[16].setVal(-2.0*spt[1]);
		dNds[17].setVal(0.0);
		
	    nVec[6].setVal((1.0-spt[0]*spt[0])*(1.0-spt[1]));
		dNds[18].setVal(-2.0*spt[0]*(1.0-spt[1]));
		dNds[19].setVal(-1.0+spt[0]*spt[0]);
		dNds[20].setVal(0.0);
		
	    nVec[7].setVal((1.0-spt[1]*spt[1])*(1.0+spt[0]));
		dNds[21].setVal(1.0-spt[1]*spt[1]);
		dNds[22].setVal(-2.0*spt[1]*(1.0+spt[0]));
		dNds[23].setVal(0.0);
		
	    nVec[8].setVal((1.0-spt[0]*spt[0])*(1.0+spt[1]));
		dNds[24].setVal(-2.0*spt[0]*(1.0+spt[1]));
		dNds[25].setVal(1.0-spt[0]*spt[0]);
		dNds[26].setVal(0.0);
		
	    nVec[9].setVal((1.0-spt[1]*spt[1])*(1.0-spt[0]));
		dNds[27].setVal(-1.0+spt[1]*spt[1]);
		dNds[28].setVal(-2.0*spt[1]*(1.0-spt[0]));
		dNds[29].setVal(0.0);
		
	} else if(type == 2) {
		nVec[0].setVal(0.5*(1.0 - spt[0]));
		dNds[0].setVal(-0.5);
		dNds[1].setVal(0.0);
		dNds[2].setVal(0.0);
		
        nVec[1].setVal(0.5*(1.0 + spt[0]));
		dNds[3].setVal(0.5);
		dNds[4].setVal(0.0);
		dNds[5].setVal(0.0);
		
        nVec[2].setVal(1.0 - spt[0]*spt[0]);
		dNds[0].setVal(-2.0*spt[0]);
		dNds[1].setVal(0.0);
		dNds[2].setVal(0.0);
	}
	return;
}

void Element::getIpData(DiffDoub1 nVec[], DiffDoub1 dNdx[], DiffDoub1& detJ, DiffDoub1 locNds[], double spt[]) {
	int i1;
	int i2;
	DiffDoub1 nCent[11];
	DiffDoub1 dNds[33];
	DiffDoub1 dNdsCent[33];
	DiffDoub1 jMat[9];
	DiffDoub1 jCent[9];
	DiffDoub1 detCent;
	DiffDoub1 jInv[9];
	DiffDoub1 jInvCent[9];
	double sCent[3] = { 0.0,0.0,0.0 };
	DiffDoub1 xVec[3];
	DiffDoub1 bVec[3];
	DiffDoub1 zDir;
	DiffDoub1 tmp;
	
	evalN(nVec, dNds, spt);
	matMul(jMat,locNds,dNds,3,numNds,3);

	if(type == 41 || type == 3) {
		zDir.setVal(jMat[0]);
		zDir.mult(jMat[4]);
		tmp.setVal(jMat[3]);
		tmp.mult(jMat[1]);
		zDir.sub(tmp);
		if(zDir.val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0);
		}
	} else if(type == 2) {
		jMat[4].setVal(1.0);
		if(jMat[0].val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0);
		}
	}
	
	getDetInv(detJ,jInv,jMat,3,0,xVec,bVec);
	
	// matMul(dNds,dNds,jInv,nDim,3,3);
	matMul(dNdx,dNds,jInv,numNds,3,3);

	if (nDim > numNds) {
		evalN(nCent, dNdsCent, sCent);
		matMul(jCent, locNds, dNdsCent, 3, numNds, 3);

		if (type == 41 || type == 3) {
			zDir.setVal(jCent[0]);
			zDir.mult(jCent[4]);
			tmp.setVal(jCent[3]);
			tmp.mult(jCent[1]);
			zDir.sub(tmp);
			if (zDir.val > 0.0) {
				jCent[8].setVal(1.0);
			}
			else {
				jCent[8].setVal(-1.0);
			}
		}
		else if (type == 2) {
			jCent[4].setVal(1.0);
			if (jCent[0].val > 0.0) {
				jCent[8].setVal(1.0);
			}
			else {
				jCent[8].setVal(-1.0);
			}
		}

		getDetInv(detCent, jInvCent, jCent, 3, 0, xVec, bVec);

		i1 = 3 * numNds;
		i2 = nDim - numNds;
		matMul(&dNdx[i1], &dNds[i1], jInvCent, i2, 3, 3);
	}
	
	return;
}

void Element::getInstOri(DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 globDisp[], int stat) {
	// stat = 1: nonlinear geometry, element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub1 version, 1st order for DiffDoub1 version
	DiffDoub1 rot[3];
	DiffDoub1 nnds;
	DiffDoub1 one;
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int stIndex;
	int iOriSize = (numNds+1)*144;
	for (i1 = 0; i1 < iOriSize; i1++) {
		instOriMat[i1].setVal(0.0);
	}
	bool isDiff = locOri[0].diffType();
	
	nnds.setVal(0.0);
	one.setVal(1.0);
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
	rot[2].setVal(0.0);
	i2 = 3 * nDim;
	for (i1 = 0; i1 < numNds; i1++) {
		rot[0].add(globDisp[i2]);
		rot[1].add(globDisp[i2 + nDim]);
		rot[2].add(globDisp[i2 + 2 * nDim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);
	
	if(stat == 1) {
		dOridThet(&instOriMat[0], locOri, rot, 0, 0);
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			stIndex = 144 * i1;
			dOridThet(&instOriMat[stIndex], locOri, rot, 0, 0);
			for (i2 = 1; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else if(isDiff) {
		for (i2 = 0; i2 < 4; i2++) {
			stIndex = 36*i2;
			dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
		}
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
	} else {
		for (i2 = 0; i2 < 4; i2++) {
			stIndex = 36*i2;
			dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				stIndex = 36*i2 + 9*i6;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,i6);
				i3 = stIndex;
				i4 = 36*i6 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
			}
		}
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					stIndex = 144*i1 + 36*i2 + 9*i6;
					dOridThet(&instOriMat[stIndex],locOri,rot,i2,i6);
					i3 = stIndex;
					i4 = 144*i1 + 36*i6 + 9*i2;
					for (i5 = 0; i5 < 9; i5++) {
						instOriMat[i4].setVal(instOriMat[i3]);
						i3++;
						i4++;
					}
				}
			}
		}
	}

	return;
}

void Element::getInstDisp(DiffDoub1 instDisp[], DiffDoub1 globDisp[], DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 xGlob[], bool nLGeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int ndOInd;
	int dof;
	int dofOInd;
	int nd2;
	int dof2;
	int dof2OInd;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 nnInv;
	DiffDoub1 nnInv2;
	
	i2 = 6*nDim;
	for (i1 = 0; i1 < i2; i1++) {
		instDisp[i1].setVal(0.0);
	}
	
	if (!nLGeom) {
		if (dv1 + dv2 == -2) {
			i7 = 3 * nDim;
			for (i1 = 0; i1 < 3; i1++) {
				i4 = i1 * nDim;
				for (i2 = 0; i2 < numNds; i2++) {
					i5 = i1 * 3;
					i6 = i2;
					for (i3 = 0; i3 < 3; i3++) {
						//i4 = i1 * nDim + i2;
						//i5 = i1 * 3 + i3;
						//i6 = i3 * nDim + i2;
						tmp.setVal(locOri[i5]);
						tmp.mult(globDisp[i6]);
						instDisp[i4].add(tmp);
						tmp.setVal(locOri[i5]);
						tmp.mult(globDisp[i6 + i7]);
						instDisp[i4 + i7].add(tmp);
						i5++;
						i6 += nDim;
					}
					i4++;
				}
			}
			i2 = 2 * numNds * dofPerNd;
			for (i1 = 0; i1 < numIntDof; i1++) {
				nd = dofTable[i2];
				dof = dofTable[i2 + 1];
				i3 = dof * nDim + nd;
				instDisp[i3].setVal(globDisp[i3]);
				i2 += 2;
			}
		}
		else if ((dv1+1)*(dv2+1) == 0) {
			dv1 = dv1 + dv2 + 1;
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < numNds) {
					i2 = nd;
					i3 = dof;
					for (i1 = 0; i1 < 3; i1++) {
						//i2 = i1 * nDim + nd;
						//i3 = i1 * 3 + dof;
						instDisp[i2].setVal(locOri[i3]);
						i2 += nDim;
						i3 += 3;
					}
				}
				else {
					i1 = dof * nDim + nd;
					instDisp[i1].setVal(1.0);
				}
			}
			else {
				i2 = 3 * nDim + nd;
				i3 = dof - 3;
				for (i1 = 0; i1 < 3; i1++) {
					//i2 = (i1 + 3) * nDim + nd;
					//i3 = i1 * 3 + (dof - 3);
					instDisp[i2].setVal(locOri[i3]);
					i2 += nDim;
					i3 += 3;
				}
			}
		}
	}
	else {
		if (dv1 + dv2 == -2) {
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = i1 * nDim + i2;
					i5 = i2;
					i6 = i2;
					i7 = i1 * 3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp2.setVal(locOri[i7]);
						tmp2.mult(xGlob[i6]);
						tmp.sub(tmp2);
						instDisp[i4].add(tmp);
						i5 += nDim;
						i6 += numNds;
						i7++;
					}
				}
			}

			i2 = numNds * dofPerNd;
			i3 = i2 + numIntDof;
			for (i1 = i2; i1 < i3; i1++) {
				i4 = dofTable[2 * i1];
				i5 = dofTable[2 * i1 + 1];
				i6 = i5 * nDim + i4;
				instDisp[i6].setVal(globDisp[i6]);
			}

			for (i1 = 0; i1 < numNds; i1++) {
				i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
				i4 = 4 * nDim + i1;
				i5 = 5 * nDim + i1;
				ndOInd = 144 * (i1 + 1);
				for (i2 = 0; i2 < 3; i2++) {
					i6 = ndOInd + 3 + i2;
					i7 = 6 + i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i3].add(tmp);
					i6 = ndOInd + 6 + i2;
					i7 = i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i4].add(tmp);
					i6 = ndOInd + i2;
					i7 = 3 + i2;
					tmp.setVal(instOriMat[i6]);
					tmp.mult(instOriMat[i7]);
					instDisp[i5].add(tmp);
				}
			}
		}
		else if ((dv1 + 1) * (dv2 + 1) == 0) {
			dv1 = dv1 + dv2 + 1;
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			if (dof < 3) {
				if (nd < numNds) {
					i1 = nd; // Index in instDisp
					i2 = dof; //Index in instOriMat
					instDisp[i1].setVal(instOriMat[i2]);
					i1 = nDim + nd;
					i2 = 3 + dof;
					instDisp[i1].setVal(instOriMat[i2]);
					i1 = 2 * nDim + nd;
					i2 = 6 + dof;
					instDisp[i1].setVal(instOriMat[i2]);
				}
				else {
					i1 = dof * nDim + nd;
					instDisp[i1].setVal(1.0);
				}
			}
			else { // dof is rotation
				nnInv.setVal(1.0 / numNds);
				dofOInd = 36 * (dof - 2);
				for (i1 = 0; i1 < 3; i1++) {
					for (i2 = 0; i2 < numNds; i2++) {
						i4 = i1 * nDim + i2;
						i5 = i2;
						i6 = i2;
						i7 = dofOInd + i1 * 3;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.setVal(globDisp[i5]);
							tmp.add(xGlob[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i5 += nDim;
							i6 += numNds;
							i7++;
						}
					}
				}

				for (i1 = 0; i1 < numNds; i1++) {
					i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
					i4 = 4 * nDim + i1;
					i5 = 5 * nDim + i1;
					ndOInd = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = ndOInd + 3 + i2;
						i7 = dofOInd + 6 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i3].add(tmp);
						i6 = ndOInd + 6 + i2;
						i7 = dofOInd + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i4].add(tmp);
						i6 = ndOInd + i2;
						i7 = dofOInd + 3 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i5].add(tmp);
						if (i1 == nd) {
							i6 = ndOInd + dofOInd + 3 + i2;
							i7 = 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + 6 + i2;
							i7 = i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + i2;
							i7 = 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i5].add(tmp);
						}
					}
				}
			}
		}
		else {
			nd = dofTable[2 * dv1];
			dof = dofTable[2 * dv1 + 1];
			nd2 = dofTable[2 * dv2];
			dof2 = dofTable[2 * dv2 + 1];
			nnInv.setVal(1.0 / numNds);
			nnInv2.setVal(nnInv);
			nnInv2.sqr();
			if (dof > 2 && dof2 > 2) {
				for (i1 = 0; i1 < 3; i1++) {
					dofOInd = 36 * (dof - 2) + 9 * (dof2 - 2) + 3 * i1;
					for (i2 = 0; i2 < numNds; i2++) {
						i4 = nDim * i1 + i2;
						i5 = i2;
						i6 = i2;
						i7 = dofOInd;
						for (i3 = 0; i3 < 3; i3++) {
							tmp.setVal(globDisp[i5]);
							tmp.add(xGlob[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv2);
							instDisp[i4].add(tmp);
							i5 += nDim;
							i6 += numNds;
							i7++;
						}
					}
				}

				dofOInd = 36 * (dof - 2);
				dof2OInd = 9 * (dof2 - 2);
				for (i1 = 0; i1 < numNds; i1++) {
					i3 = 3 * nDim + i1; // indexes of thetax, y and z for node i1
					i4 = 4 * nDim + i1;
					i5 = 5 * nDim + i1;
					ndOInd = 144 * (i1 + 1);
					for (i2 = 0; i2 < 3; i2++) {
						i6 = ndOInd + 3 + i2;
						i7 = dofOInd + dof2OInd + 6 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i3].add(tmp);
						i6 = ndOInd + 6 + i2;
						i7 = dofOInd + dof2OInd + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i4].add(tmp);
						i6 = ndOInd + i2;
						i7 = dofOInd + dof2OInd + 3 + i2;
						tmp.setVal(instOriMat[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i5].add(tmp);
						if (i1 == nd) {
							i6 = ndOInd + dofOInd + 3 + i2;
							i7 = dof2OInd + 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + 6 + i2;
							i7 = dof2OInd + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + i2;
							i7 = dof2OInd + 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i5].add(tmp);
						}
						if (i1 == nd2) {
							i6 = ndOInd + dof2OInd + 3 + i2;
							i7 = dofOInd + 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dof2OInd + 6 + i2;
							i7 = dofOInd + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dof2OInd + i2;
							i7 = dofOInd + 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							tmp.mult(nnInv);
							instDisp[i5].add(tmp);
						}
						if (i1 == nd && i1 == nd2) {
							i6 = ndOInd + dofOInd + dof2OInd + 3 + i2;
							i7 = 6 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i3].add(tmp);
							i6 = ndOInd + dofOInd + dof2OInd + 6 + i2;
							i7 = i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i4].add(tmp);
							i6 = ndOInd + dofOInd + dof2OInd + i2;
							i7 = 3 + i2;
							tmp.setVal(instOriMat[i6]);
							tmp.mult(instOriMat[i7]);
							instDisp[i5].add(tmp);
						}
					}
				}
			}
			else if (dof < 3 && dof2 < 3) {
				return;
			}
			else {
				if (dof > 2) {
					i1 = dof;
					dof = dof2;
					dof2 = i1;
					i1 = nd;
					nd = nd2;
					nd2 = i1;
				}
				i1 = 36 * (dof2 - 2) + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = nd;
				instDisp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 3 + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = nDim + nd;
				instDisp[i2].add(tmp);
				i1 = 36 * (dof2 - 2) + 6 + dof;
				tmp.setVal(instOriMat[i1]);
				tmp.mult(nnInv);
				i2 = 2 * nDim + nd;
				instDisp[i2].add(tmp);
			}
		}
	}

	
	return;
}

void Element::getStressPrereq(DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int numLay;
	DiffDoub1 offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	getNdVel(pre.globVel, ndAr);
	getNdAcc(pre.globAcc, ndAr);
	getNdTemp(pre.globTemp, ndAr);
	getNdTdot(pre.globTdot, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		if (type != 2) {
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, dvAr);
			getLayerAngle(pre.layerAng, dvAr);
			getLayerQ(pre.layerQ, dvAr);
			getLayerD(pre.layerD, dvAr);
			getLayerThExp(pre.layerTE, dvAr);
			getLayerEinit(pre.layerE0, dvAr);
			getLayerDen(pre.layerDen, dvAr);
			getLayerCond(pre.layerTC, dvAr);
			getLayerSpecHeat(pre.layerSH, dvAr);
			getABD(pre.Cmat, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerAng);
			getShellDamp(pre.Dmat, pre.layerThk, pre.layerZ, pre.layerD, pre.layerAng);
			getShellExpLoad(pre.thermExp, pre.Einit, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerTE, pre.layerE0, pre.layerAng);
			getShellMass(pre.Mmat, pre.layerThk, pre.layerZ, pre.layerDen, dvAr);
			getShellCond(pre.TCmat, pre.layerThk, pre.layerAng, pre.layerTC, dvAr);
			getShellSpecHeat(pre.SpecHeat, pre.layerThk, pre.layerSH, pre.layerDen);
		}
		else {
			getBeamStiff(pre.Cmat, dvAr);
			getBeamDamp(pre.Dmat, dvAr);
			getBeamExpLoad(pre.thermExp, pre.Einit, dvAr);
			getBeamMass(pre.Mmat, dvAr);
			getBeamCond(pre.TCmat, dvAr);
			getBeamSpecHeat(pre.SpecHeat, dvAr);
		}
	}
	else if (type == 21) {
		getFrcFldConst(pre.frcFldCoef, pre.frcFldExp, dvAr);
	}
	else if (type == 1) {
		getMassPerEl(pre.massPerEl, dvAr);
	}
	else {
		getSolidStiff(pre.Cmat, dvAr);
		getSolidDamp(pre.Dmat, dvAr);
		getThermalExp(pre.thermExp, pre.Einit, dvAr);
		getDensity(pre.Mmat[0], 0, dvAr);
		getConductivity(pre.TCmat, dvAr);
		getSpecificHeat(pre.SpecHeat, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getVolume(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer) {
	int i1;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 thk;
	DiffDoub1 tmp;

	if (type == 2) {
		thk.setVal(sectPtr->getArea()); //rem: update to factor in dvars;
	} else if (type == 3 || type == 41) {
		thk.setVal(pre.layerThk[layer]);
	} else {
		thk.setVal(1.0);
	}

	vol.setVal(0.0);
	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		tmp.mult(detJ);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::getSectionDef(DiffDoub1 secDef[], DiffDoub1 globDisp[],  DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 xGlob[], DiffDoub1 dNdx[], DiffDoub1 nVec[], bool nLGeom, int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub1 instDisp[60];
	DiffDoub1 ux[9];
	DiffDoub1 rx[9];
	DiffDoub1 rot[3];
	DiffDoub1 tmp;

	if (dofPerNd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			secDef[i1].setVal(0.0);
		}
		return;
	}
	
	getInstDisp(instDisp, globDisp, instOriMat, locOri, xGlob, nLGeom, dv1, dv2);
	
	if(dv1 + dv2 == -2) {
		matMul(ux,instDisp,dNdx,3,nDim,3);
		i1 = 3*nDim;
		matMul(rx,&instDisp[i1],dNdx,3,nDim,3);
		matMul(rot,&instDisp[i1],nVec,3,nDim,1);
	} else if((dv1+1)*(dv2+1) == 0) {
		dv1 = dv1 + dv2 + 1;
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		if(dof < 3) {
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].setVal(instDisp[i4]);
					ux[i3].mult(dNdx[i5]);
					rx[i3].setVal(0.0);
					i3++;
					i5++;
				}
				rot[i1].setVal(0.0);
				i4+= nDim;
			}
		} else {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = nDim*i1;
					i6 = i2;
					i7 = nDim*(i1+3);
					ux[i4].setVal(0.0);
					rx[i4].setVal(0.0);
					for (i3 = 0; i3 < numNds; i3++) {
						tmp.setVal(instDisp[i5]);
						tmp.mult(dNdx[i6]);
						ux[i4].add(tmp);
						tmp.setVal(instDisp[i7]);
						tmp.mult(dNdx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].setVal(0.0);
				i3 = nDim*(i1+3);
				for (i2 = 0; i2 < numNds; i2++) {
					tmp.setVal(instDisp[i3]);
					tmp.mult(nVec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		}
	} else {
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		nd2 = dofTable[2*dv2];
		dof2 = dofTable[2*dv2+1];
		if(dof < 3 && dof2 < 3) {
			for (i1 = 0; i1 < 9; i1++) {
				secDef[i1].setVal(0.0);
			}
			return;
		} else if(dof > 2 && dof2 > 2) {
			i4 = 0;
			for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < 3; i2++) {
					i5 = nDim*i1;
					i6 = i2;
					i7 = nDim*(i1+3);
					ux[i4].setVal(0.0);
					for (i3 = 0; i3 < numNds; i3++) {
						tmp.setVal(instDisp[i5]);
						tmp.mult(dNdx[i6]);
						ux[i4].add(tmp);
						tmp.setVal(instDisp[i7]);
						tmp.mult(dNdx[i6]);
						rx[i4].add(tmp);
						i5++;
						i6+= 3;
						i7++;
					}
					i4++;
				}
			}
			for (i1 = 0; i1 < 3; i1++) {
				rot[i1].setVal(0.0);
				i3 = nDim*(i1+3);
				for (i2 = 0; i2 < numNds; i2++) {
					tmp.setVal(instDisp[i3]);
					tmp.mult(nVec[i2]);
					rot[i1].add(tmp);
					i3++;
				}
			}
		} else {
			if(dof > dof2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i3 = 0;
			i4 = nd;
			for (i1 = 0; i1 < 3; i1++) {
				i5 = 3*nd;
				for (i2 = 0; i2 < 3; i2++) {
					ux[i3].setVal(instDisp[i4]);
					ux[i3].mult(dNdx[i5]);
					rx[i3].setVal(0.0);
					i3++;
					i5++;
				}
				rot[i1].setVal(0.0);
				i4+= nDim;
			}
		}
	}
	
	if(type == 2) {
		secDef[0].setVal(ux[0]);
		secDef[1].setVal(ux[3]);
		secDef[1].sub(rot[2]);
		secDef[2].setVal(ux[6]);
		secDef[2].add(rot[1]);
		secDef[3].setVal(rx[0]);
		secDef[4].setVal(rx[3]);
		secDef[5].setVal(rx[6]);
	} else {
		secDef[0].setVal(ux[0]);
		secDef[1].setVal(ux[4]);
		secDef[2].setVal(ux[1]);
		secDef[2].add(ux[3]);
		secDef[3].setVal(rx[3]);
		secDef[4].setVal(rx[1]);
		secDef[4].neg();
		secDef[5].setVal(rx[4]);
		secDef[5].sub(rx[0]);
		secDef[6].setVal(ux[6]);
		secDef[6].add(rot[1]);
		secDef[7].setVal(ux[7]);
		secDef[7].sub(rot[0]);
		secDef[8].setVal(2.0);
		secDef[8].mult(rot[2]);
		secDef[8].sub(ux[3]);
		secDef[8].add(ux[1]);
	}
	
	return;
}

void Element::getSolidStrain(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 dNdx[], DiffDoub1 locOri[], int dv1, int dv2, bool nLGeom) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int nd2;
	int dof2;
	DiffDoub1 uxL[9];
	DiffDoub1 strnMat[9];
	DiffDoub1 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strnMat[i1].setVal(0.0);
	}
	if(dv1 + dv2 == -2) {
		matMul(uxL,locOri,ux,3,3,3);
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			i5 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				strnMat[i4].add(uxL[i4]);
				strnMat[i4].add(uxL[i5]);
				if(nLGeom) {
					i6 = i1;
					i7 = i2;
				    for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(ux[i6]);
						tmp.mult(ux[i7]);
						strnMat[i4].add(tmp);
						i6+= 3;
						i7+= 3;
				    }
				}
				i4++;
				i5+= 3;
			}
		}
	} else if((dv1+1)*(dv2+1) == 0) {
		if(dv1 > -1) {
			nd = dofTable[2*dv1];
			dof = dofTable[2*dv1+1];
		} else {
			nd = dofTable[2*dv2];
			dof = dofTable[2*dv2+1];
		}
		for (i1 = 0; i1 < 3; i1++) {
			i4 = 4*i1;
			for (i2 = i1; i2 < 3; i2++) {
				i5 = 3*i1 + dof;
				i6 = 3*nd + i2;
				tmp.setVal(locOri[i5]);
				tmp.mult(dNdx[i6]);
				strnMat[i4].add(tmp);
				i5 = 3*i2 + dof;
				i6 = 3*nd + i1;
				tmp.setVal(locOri[i5]);
				tmp.mult(dNdx[i6]);
				strnMat[i4].add(tmp);
				if(nLGeom) {
					i5 = 3*nd + i1;
					i6 = 3*dof + i2;
					tmp.setVal(dNdx[i5]);
					tmp.mult(ux[i6]);
					strnMat[i4].add(tmp);
					i5 = 3*nd + i2;
					i6 = 3*dof + i1;
					tmp.setVal(dNdx[i5]);
					tmp.mult(ux[i6]);
					strnMat[i4].add(tmp);
				}
				i4++;
			}
		}
	} else {
		if(nLGeom) {
			nd = dofTable[2*dv1];
			dof = dofTable[2*dv1+1];
			nd2 = dofTable[2*dv2];
			dof2 = dofTable[2*dv2+1];
			if(dof == dof2) {
				for (i1 = 0; i1 < 3; i1++) {
					i4 = 4*i1;
					for (i2 = i1; i2 < 3; i2++) {
						i5 = 3*nd + i1;
						i6 = 3*nd2 + i2;
						tmp.setVal(dNdx[i5]);
						tmp.mult(dNdx[i6]);
						strnMat[i4].add(tmp);
						i5 = 3*nd2 + i1;
						i6 = 3*nd + i2;
						tmp.setVal(dNdx[i5]);
						tmp.mult(dNdx[i6]);
						strnMat[i4].add(tmp);
						i4++;
					}
				}
			}
		}
	}
	
	tmp.setVal(0.5);
	strain[0].setVal(strnMat[0]);
	strain[0].mult(tmp);
	strain[1].setVal(strnMat[4]);
	strain[1].mult(tmp);
	strain[2].setVal(strnMat[8]);
	strain[2].mult(tmp);
	strain[3].setVal(strnMat[1]);
	strain[4].setVal(strnMat[2]);
	strain[5].setVal(strnMat[5]);
	
	return;
}

void Element::getStressStrain(DiffDoub1 stress[], DiffDoub1 strain[], double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 ux[9];
	DiffDoub1 secDef[9];
	DiffDoub1 sectStrn[3];
	DiffDoub1 adjStn[6];
	DiffDoub1 tmp;
	DiffDoub1 ipTemp;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);

	ipTemp.setVal(0.0);
	for (i1 = 0; i1 < numNds; i1++) {
		tmp.setVal(pre.globTemp[i1]);
		tmp.mult(nVec[i1]);
		ipTemp.add(tmp);
	}

	if (type == 41 || type == 3) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}

		getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, -1, -1);
		
		sectStrn[0].setVal(secDef[0]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[3]);
		sectStrn[0].add(tmp);
		
		sectStrn[1].setVal(secDef[1]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[4]);
		sectStrn[1].sub(tmp);

		sectStrn[2].setVal(secDef[2]);
		tmp.setVal(pre.layerZ[layer]);
		tmp.mult(secDef[5]);
		sectStrn[2].add(tmp);

		tmp.setVal(pre.layerAng[layer]);
		tmp.neg();
		transformStrain(strain, sectStrn, tmp);
		i2 = 3 * layer;
		for (i1 = 0; i1 < 3; i1++) {
			tmp.setVal(pre.layerTE[i2]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.layerE0[i2]);
			i2++;
		}
		matMul(stress, &pre.layerQ[9 * layer], adjStn, 3, 3, 1);

		strain[3].setVal(strain[2]);
		strain[2].setVal(0.0);
		strain[4].setVal(0.0);
		strain[5].setVal(0.0);

		stress[3].setVal(stress[2]);
		stress[2].setVal(0.0);
		stress[4].setVal(0.0);
		stress[5].setVal(0.0);

	} else if(type != 2) {
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.setVal(pre.thermExp[i1]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.Einit[i1]);
		}
		matMul(stress, pre.Cmat, adjStn, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(DiffDoub1 dsdU[], DiffDoub1 dedU[], DiffDoub1 dsdT[], double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 ux[9];
	DiffDoub1 secDef[9];
	DiffDoub1 sectStrn[3];
	DiffDoub1 dStrain[6];
	DiffDoub1 dStress[6];
	DiffDoub1 CTE[6];
	DiffDoub1 CTEN[60];
	DiffDoub1 tmp;

	int totDof = numNds * dofPerNd + numIntDof;
	i2 = 6 * totDof;
	for (i1 = 0; i1 < i2; i1++) {
		dsdU[i1].setVal(0.0);
		dedU[i1].setVal(0.0);
	}
	i2 = 6 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		dsdT[i1].setVal(0.0);
	}

	if (type == 41 || type == 3) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}

		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		for (i1 = 0; i1 < totDof; i1++) {
			getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, i1, -1);

			sectStrn[0].setVal(secDef[0]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[3]);
			sectStrn[0].add(tmp);

			sectStrn[1].setVal(secDef[1]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[4]);
			sectStrn[1].sub(tmp);

			sectStrn[2].setVal(secDef[2]);
			tmp.setVal(pre.layerZ[layer]);
			tmp.mult(secDef[5]);
			sectStrn[2].add(tmp);

			tmp.setVal(pre.layerAng[layer]);
			tmp.neg();
			transformStrain(dStrain, sectStrn, tmp);
			matMul(dStress, &pre.layerQ[9 * layer], dStrain, 3, 3, 1);

			dedU[i1].setVal(dStrain[0]);
			dedU[i1 + totDof].setVal(dStrain[1]);
			dedU[i1 + 3 * totDof].setVal(dStrain[2]);

			dsdU[i1].setVal(dStress[0]);
			dsdU[totDof + i1].setVal(dStress[1]);
			dsdU[3 * totDof + i1].setVal(dStress[2]);
		}
		matMul(CTE, &pre.layerQ[9 * layer], &pre.layerTE[3 * layer], 3, 3, 1);
		CTE[0].neg();
		CTE[1].neg();
		CTE[2].neg();
		matMul(CTEN, CTE, nVec, 3, 1, numNds);
		for (i1 = 0; i1 < numNds; i1++) {
			dsdT[i1].setVal(CTEN[i1]);
			dsdT[i1 + numNds].setVal(CTEN[i1 + numNds]);
			dsdT[i1 + 3 * numNds].setVal(CTEN[i1 + 2 * numNds]);
		}
	}
	else if (type != 2) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		for (i1 = 0; i1 < totDof; i1++) {
			getSolidStrain(dStrain, ux, dNdx, pre.locOri, i1, -1, nLGeom);
			matMul(dStress, pre.Cmat, dStrain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				dedU[i3].setVal(dStrain[i2]);
				dsdU[i3].setVal(dStress[i2]);
				i3 += totDof;
			}
		}
		matMul(CTE, pre.Cmat, pre.thermExp, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			CTE[i1].neg();
		}
		matMul(dsdT, CTE, nVec, 6, 1, numNds);
	}

	return;
}

void Element::getDefFrcMom(DiffDoub1 def[], DiffDoub1 frcMom[], double spt[], bool nLGeom, DiffDoub1StressPrereq& pre) {
	int i1;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 ptTemp;
	DiffDoub1 tmp;

	if (dofPerNd != 6) {
		for (i1 = 0; i1 < 9; i1++) {
			def[i1].setVal(0.0);
		}
		return;
	}

	if (nLGeom) {
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
	}

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	getSectionDef(def, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, -1, -1);

	ptTemp.setVal(0.0);
	for (i1 = 0; i1 < numNds; i1++) {
		tmp.setVal(pre.globTemp[i1]);
		tmp.mult(nVec[i1]);
		ptTemp.add(tmp);
	}

	matMul(frcMom, pre.Cmat, def, defDim, defDim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frcMom[i1].sub(pre.Einit[i1]);
		tmp.setVal(ptTemp);
		tmp.mult(pre.thermExp[i1]);
		frcMom[i1].sub(tmp);
	}

	return;
}

void Element::dDefFrcMomdU(DiffDoub1 dDefdU[], DiffDoub1 dFrcMomdU[], DiffDoub1 dFrcMomdT[], double spt[], bool nLGeom, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	int totDof;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 ptTemp;
	DiffDoub1 tmp;
	DiffDoub1 def[9];

	if (nLGeom) {
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
	}


	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	totDof = numNds * dofPerNd + numIntDof;
	for (i1 = 0; i1 < totDof; i1++) {
		getSectionDef(def, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, nLGeom, i1, -1);
		i2 = i1;
		for (i3 = 0; i3 < defDim; i3++) {
			dDefdU[i2].setVal(def[i3]);
			i2 += totDof;
		}
	}

	matMul(dFrcMomdU, pre.Cmat, dDefdU, defDim, defDim, totDof);

	matMul(dFrcMomdT, pre.thermExp, nVec, 6, 1, numNds);

	return;
}

void Element::getFluxTGrad(DiffDoub1 flux[], DiffDoub1 tGrad[], double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	matMul(tGrad, pre.globTemp, dNdx, 1, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		matMul(flux, &pre.layerTC[i1], tGrad, 3, 3, 1);
	}
	else {
		matMul(flux, pre.TCmat, tGrad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::dFluxTGraddT(DiffDoub1 dFdT[], DiffDoub1 dTG[], double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	transpose(dTG, dNdx, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		matMul(dFdT, &pre.layerTC[i1], dTG, 3, 3, numNds);
	}
	else {
		matMul(dFdT, pre.TCmat, dTG, 3, 3, numNds);
	}
	i2 = 3 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		dFdT[i1].neg();
	}

	return;
}

void Element::putVecToGlobMat(SparseMat& qMat, DiffDoub1 elQVec[], bool forTherm, int matRow, Node* ndAr[]) {
	int i1;
	int ndDof = numNds*dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[i1]->getSortedRank();
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < totDof; i1++) {
			if (i1 < ndDof) {
				nd = nodes[dofTable[2 * i1]];
				dof = dofTable[2 * i1 + 1];
				globInd = ndAr[nd]->getDofIndex(dof);
				qMat.addEntry(matRow, globInd, elQVec[i1].val);
			}
			else {
				dof = i1 - ndDof;
				globInd = intDofIndex + dof;
				qMat.addEntry(matRow, globInd, elQVec[i1].val);
			}
		}
	}

	return;
}

//end dup
 
//end skip 
 
 
void Element::getElVec(double elVec[], double globVec[], bool forTherm, bool intnl, Node* ndAr[]) {
	int i1;
	int i2;
	int nd;
	int dof;
	int globInd;
	int ndDof;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[nd]->getSortedRank();
			elVec[i1] = globVec[globInd];
		}
	}
	else {
		ndDof = numNds * dofPerNd;
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = nodes[dofTable[i2]];
			dof = dofTable[i2 + 1];
			globInd = ndAr[nd]->getDofIndex(dof);
			elVec[i1] = globVec[globInd];
			i2 += 2;
		}
		if (intnl) {
			for (i1 = 0; i1 < numIntDof; i1++) {
				elVec[i1 + ndDof] = globVec[i1 + intDofIndex];
			}
		}
	}

	return;
}

void Element::addToGlobVec(double elVec[], double globVec[], bool forTherm, bool intnl, Node* ndAr[]) {
	int i1;
	int i2;
	int nd;
	int dof;
	int globInd;
	int ndDof;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[nd]->getSortedRank();
			globVec[globInd] += elVec[i1];
		}
	}
	else {
		ndDof = numNds * dofPerNd;
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = nodes[dofTable[i2]];
			dof = dofTable[i2 + 1];
			globInd = ndAr[nd]->getDofIndex(dof);
			globVec[globInd] += elVec[i1];
			i2 += 2;
		}
		if (intnl) {
			for (i1 = 0; i1 < numIntDof; i1++) {
				globVec[i1 + intDofIndex] += elVec[i1 + ndDof];
			}
		}
	}

	return;
}
 
 
// 

