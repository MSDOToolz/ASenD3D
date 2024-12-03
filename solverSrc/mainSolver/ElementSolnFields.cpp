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
void Element::getNdDisp(vector<DiffDoub0>& globDisp, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globDisp[i3].setVal(thisNd.displacement[i2]);
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

void Element::getNdVel(vector<DiffDoub0>& globVel, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globVel[i3].setVal(thisNd.velocity[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdAcc(vector<DiffDoub0>& globAcc, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globAcc[i3].setVal(thisNd.acceleration[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdFlVel(vector<DiffDoub0>& flVel, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			flVel[i3].setVal(thisNd.flVel[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdFlVDot(vector<DiffDoub0>& flVDot, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			flVDot[i3].setVal(thisNd.flVelDot[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdTemp(vector<DiffDoub0>& globTemp, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		globTemp[i1].setVal(thisNd.temperature);
	}
	return;
}

void Element::getNdTdot(vector<DiffDoub0>& globTdot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		globTdot[i1].setVal(thisNd.tempChangeRate);
	}
	return;
}

void Element::getNdFlDen(vector<DiffDoub0>& flDen, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		flDen[i1].setVal(thisNd.flDen);
	}

	return;
}

void Element::getNdFlDenDot(vector<DiffDoub0>& flDenDot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		flDenDot[i1].setVal(thisNd.flDenDot);
	}

	return;
}

void Element::getNdTurbE(vector<DiffDoub0>& turbE, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		turbE[i1].setVal(thisNd.turbE);
	}

	return;
}

void Element::getNdTurbEDot(vector<DiffDoub0>& turbEDot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		turbEDot[i1].setVal(thisNd.turbEDot);
	}

	return;
}

void Element::evalN(DiffDoub0 nVec[], DiffDoub0 dNds[], double spt[]) {
	if(type == 4 || type == 400) {
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
	} else if(type == 6 || type == 600) {
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
	} else if(type == 8 || type == 81 || type == 800) {
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
	else if (type == 10 || type == 1000) {
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

void Element::getIpData(DiffDoub0 nVec[], DiffDoub0 dNdx[], DiffDoub0& detJ, vector<DiffDoub0>& locNds, double spt[]) {
	int i1;
	int i2;
	DiffDoub0 nCent[11];
	DiffDoub0 dNds[33];
	DiffDoub0 dNdsCent[33];
	DiffDoub0 jMat[9];
	DiffDoub0 tmpNds[30];
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
	vecToAr(tmpNds, locNds, 0, 30);
	matMul(jMat,tmpNds,dNds,3,numNds,3);

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
		matMul(jCent, tmpNds, dNdsCent, 3, numNds, 3);

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

void Element::getInstOri(vector<DiffDoub0>& instOriMat, vector<DiffDoub0>& locOri, vector<DiffDoub0>& globDisp, int stat) {
	// stat = 1: nonlinear geometry, element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub0 version, 1st order for DiffDoub1 version
	DiffDoub0 rot[3];
	DiffDoub0 nnds;
	DiffDoub0 one;
	DiffDoub0 tmpOri[9];
	DiffDoub0 tmpInst[9];
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

	vecToAr(tmpOri, locOri, 0, 9);
	
	if(stat == 1) {
		dOridThet(tmpInst, tmpOri, rot, 0, 0);
		arToVec(tmpInst, instOriMat, 0, 9);
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			stIndex = 144 * i1;
			dOridThet(tmpInst, tmpOri, rot, 0, 0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
			for (i2 = 1; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
			dOridThet(tmpInst,tmpOri,rot,i2,0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
			dOridThet(tmpInst,tmpOri,rot,i2,0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				stIndex = 36*i2 + 9*i6;
				dOridThet(tmpInst,tmpOri,rot,i2,i6);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					stIndex = 144*i1 + 36*i2 + 9*i6;
					dOridThet(tmpInst,tmpOri,rot,i2,i6);
					arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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

void Element::getInstDisp(DiffDoub0 instDisp[], vector<DiffDoub0>& globDisp, vector<DiffDoub0>& instOriMat, vector<DiffDoub0>& locOri, vector<DiffDoub0>& xGlob, bool nLGeom, int dv1, int dv2) {
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

void Element::getStressPrereq(DiffDoub0StressPrereq& pre, vector<Section>& secAr, vector<Material>& matAr, vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	int numLay;
	DiffDoub0 offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, secAr, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	getNdVel(pre.globVel, ndAr);
	getNdAcc(pre.globAcc, ndAr);
	getNdTemp(pre.globTemp, ndAr);
	getNdTdot(pre.globTdot, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		if (type != 2) {
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, secAr, dvAr);
			getLayerAngle(pre.layerAng, secAr, dvAr);
			getLayerQ(pre.layerQ, secAr, matAr, dvAr);
			getLayerD(pre.layerD, secAr, matAr, dvAr);
			getLayerThExp(pre.layerTE, secAr, matAr, dvAr);
			getLayerEinit(pre.layerE0, secAr, dvAr);
			getLayerDen(pre.layerDen, secAr, matAr, dvAr);
			getLayerCond(pre.layerTC, secAr, matAr, dvAr);
			getLayerSpecHeat(pre.layerSH, secAr, matAr, dvAr);
			getABD(pre.Cmat, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerAng, secAr);
			getShellDamp(pre.Dmat, pre.layerThk, pre.layerZ, pre.layerD, pre.layerAng, secAr);
			getShellExpLoad(pre.thermExp, pre.Einit, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerTE, pre.layerE0, pre.layerAng, secAr);
			getShellMass(pre.Mmat, pre.layerThk, pre.layerZ, pre.layerDen, secAr, dvAr);
			getShellCond(pre.TCmat, pre.layerThk, pre.layerAng, pre.layerTC, secAr, dvAr);
			getShellSpecHeat(pre.SpecHeat, pre.layerThk, pre.layerSH, pre.layerDen, secAr);
		}
		else {
			getBeamStiff(pre.Cmat, secAr, matAr, dvAr);
			getBeamDamp(pre.Dmat, secAr, matAr, dvAr);
			getBeamExpLoad(pre.thermExp, pre.Einit, secAr, matAr, dvAr);
			getBeamMass(pre.Mmat, secAr, matAr, dvAr);
			getBeamCond(pre.TCmat, secAr, matAr, dvAr);
			getBeamSpecHeat(pre.SpecHeat, secAr, matAr, dvAr);
		}
	}
	else if (type == 21) {
		getFrcFldConst(pre.frcFldCoef, pre.frcFldExp, secAr, dvAr);
		getThrmFldConst(pre.thrmFldCoef, pre.refTemp, secAr, dvAr);
	}
	else if (type == 1) {
		getMassPerEl(pre.massPerEl, secAr, dvAr);
		getSpecificHeat(pre.SpecHeat, secAr, matAr, dvAr);
	}
	else {
		getSolidStiff(pre.Cmat, secAr, matAr, dvAr);
		getSolidDamp(pre.Dmat, secAr, matAr, dvAr);
		getThermalExp(pre.thermExp, pre.Einit, secAr, matAr, dvAr);
		getDensity(pre.Mmat[0], 0, secAr, matAr, dvAr);
		getConductivity(pre.TCmat, secAr, matAr, dvAr);
		getSpecificHeat(pre.SpecHeat, secAr, matAr, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getFluidPrereq(DiffDoub0FlPrereq& pre, vector<Section>& secAr, vector<Fluid>& flAr, vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	int i1;
	int i2;
	double secProp;
	
	getNdCrds(pre.globNds, ndAr, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	i2 = 3 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		pre.globNds[i1].add(pre.globDisp[i1]);
	}
	getNdVel(pre.globVel, ndAr);
	getNdFlDen(pre.flDen, ndAr);
	getNdFlVel(pre.flVel, ndAr);
	getNdTemp(pre.flTemp, ndAr);
	getNdTurbE(pre.flTurbE, ndAr);
	getNdFlDenDot(pre.flDenDot, ndAr);
	getNdFlVDot(pre.flVelDot, ndAr);
	getNdTdot(pre.flTDot, ndAr);
	getNdTurbEDot(pre.flTurbEDot, ndAr);

	Section& thisSec = secAr[sectPtr];
	Fluid& thisFl = flAr[thisSec.flPtr];

	secProp = thisSec.refTurbE;
	pre.refTurbE.setVal(secProp);
	getGenProp(pre.refTurbE, "refTurbE", dvAr);

	secProp = thisSec.gradVTurbCoef;
	pre.gradVTurbCoef.setVal(secProp);
	getGenProp(pre.refTurbE, "gradVTurbCoef", dvAr);

	secProp = thisSec.dissTurbCoef;
	pre.dissTurbCoef.setVal(secProp);
	getGenProp(pre.refTurbE, "dissTurbCoef", dvAr);

	secProp = thisFl.viscosity;
	pre.refVisc.setVal(secProp);
	getGenProp(pre.refVisc, "viscosity", dvAr);

	secProp = thisSec.denVisCoef;
	pre.denVisCoef.setVal(secProp);
	getGenProp(pre.denVisCoef, "denVisCoef", dvAr);

	secProp = thisSec.tempVisCoef;
	pre.tempVisCoef.setVal(secProp);
	getGenProp(pre.tempVisCoef, "tempVisCoef", dvAr);

	secProp = thisSec.turbVisCoef;
	pre.turbVisCoef.setVal(secProp);
	getGenProp(pre.turbVisCoef, "turbVisCoef", dvAr);

	secProp = thisSec.refEnth;
	pre.refEnth.setVal(secProp);
	getGenProp(pre.refEnth, "refEnth", dvAr);

	secProp = thisSec.enthCoef;
	pre.denEnthCoef.setVal(secProp);
	getGenProp(pre.denEnthCoef, "denEnthCoef", dvAr);

	secProp = thisSec.enthExp;
	pre.denEnthExp.setVal(secProp);
	getGenProp(pre.denEnthExp, "denEnthExp", dvAr);

	secProp = thisSec.presCoef;
	pre.denPresCoef.setVal(secProp);
	getGenProp(pre.denPresCoef, "denPresCoef", dvAr);

	secProp = thisSec.presExp;
	pre.denPresExp.setVal(secProp);
	getGenProp(pre.denPresExp, "denPresExp", dvAr);

	secProp = thisSec.refDen;
	pre.refDen.setVal(secProp);
	getGenProp(pre.refDen, "refDen", dvAr);

	secProp = thisSec.refTemp;
	pre.refTemp.setVal(secProp);
	getGenProp(pre.refTemp, "refTemp", dvAr);

	secProp = thisFl.thermCond;
	pre.thermCond.setVal(secProp);
	getGenProp(pre.thermCond, "thermCond", dvAr);

	secProp = thisFl.specHeat;
	pre.specHeat.setVal(secProp);
	getGenProp(pre.specHeat, "specHeat", dvAr);

	secProp = thisFl.idealGas;
	pre.iGConst.setVal(secProp);
	getGenProp(pre.iGConst, "iGConst", dvAr);

	return;
}

void Element::getVolume(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer, vector<Section>& secAr, vector<DesignVariable>& dvAr) {
	int i1;
	int i2;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 thk;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;
    double sTmp[3];

	if (type == 2) {
		thk.setVal(secAr[sectPtr].area); //rem: update to factor in dvars;
		for (auto& dv : designVars) {
			DesignVariable& thisDV = dvAr[dv.intDat];
			if (thisDV.category == "area") {
				tmp.setVal(dv.doubDat);
				thisDV.getValue(dvVal);
				dvVal.mult(tmp);
				thk.add(dvVal);
			}
		}
	} else if (type == 3 || type == 41) {
		thk.setVal(pre.layerThk[layer]);
	} else {
		thk.setVal(1.0);
	}

	vol.setVal(0.0);
	for (i1 = 0; i1 < numIP; i1++) {
		i2 = 3 * i1;
		vecToAr(sTmp, intPts, i2, i2 + 3);
		getIpData(nVec, dNdx, detJ, pre.locNds, sTmp);
		tmp.setVal(ipWt[i1]);
		tmp.mult(detJ);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::getSectionDef(DiffDoub0 secDef[], vector<DiffDoub0>& globDisp,  vector<DiffDoub0>& instOriMat, vector<DiffDoub0>& locOri, vector<DiffDoub0>& xGlob, DiffDoub0 dNdx[], DiffDoub0 nVec[], bool nLGeom, int dv1, int dv2) {
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

void Element::getSolidStrain(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 dNdx[], vector<DiffDoub0>& locOri, int dv1, int dv2, bool nLGeom) {
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
	DiffDoub0 tmpOri[9];
	DiffDoub0 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strnMat[i1].setVal(0.0);
	}
	if(dv1 + dv2 == -2) {
		vecToAr(tmpOri, locOri, 0, 9);
		matMul(uxL,tmpOri,ux,3,3,3);
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
	DiffDoub0 tmpAr[60];

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
		vecToAr(tmpAr, pre.globDisp, 0, 60);
		matMul(ux, tmpAr, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.setVal(pre.thermExp[i1]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.Einit[i1]);
		}
		vecToAr(tmpAr, pre.Cmat, 0, 36);
		matMul(stress, tmpAr, adjStn, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(vector<DiffDoub0>& dsdU, vector<DiffDoub0>& dedU, vector<DiffDoub0>& dsdT, double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre) {
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
	DiffDoub0 tmpAr[60];
	DiffDoub0 tmpAr2[6];

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
		vecToAr(tmpAr, pre.globDisp, 0, 60);
		matMul(ux, tmpAr, dNdx, 3, nDim, 3);
		vecToAr(tmpAr, pre.Cmat, 0, 36);
		for (i1 = 0; i1 < totDof; i1++) {
			getSolidStrain(dStrain, ux, dNdx, pre.locOri, i1, -1, nLGeom);
			matMul(dStress, tmpAr, dStrain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				dedU[i3].setVal(dStrain[i2]);
				dsdU[i3].setVal(dStress[i2]);
				i3 += totDof;
			}
		}
		vecToAr(tmpAr2, pre.thermExp, 0, 6);
		matMul(CTE, tmpAr, tmpAr2, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			CTE[i1].neg();
		}
		matMul(tmpAr, CTE, nVec, 6, 1, numNds);
		arToVec(tmpAr, dsdT, 0, 60);
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
	DiffDoub0 tmpAr[36];

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

	vecToAr(tmpAr, pre.Cmat, 0, 36);
	matMul(frcMom, tmpAr, def, defDim, defDim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frcMom[i1].sub(pre.Einit[i1]);
		tmp.setVal(ptTemp);
		tmp.mult(pre.thermExp[i1]);
		frcMom[i1].sub(tmp);
	}

	return;
}

void Element::dDefFrcMomdU(vector<DiffDoub0>& dDefdU, vector<DiffDoub0>& dFrcMomdU, vector<DiffDoub0>& dFrcMomdT, double spt[], bool nLGeom, DiffDoub0StressPrereq& pre) {
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
	DiffDoub0 tmpAr[60];
	DiffDoub0 tmpAr2[6];

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

	vecToAr(tmpAr2, pre.thermExp, 0, 6);
	matMul(tmpAr, tmpAr2, nVec, 6, 1, numNds);
	arToVec(tmpAr, dFrcMomdT, 0, 60);

	return;
}

void Element::getFluxTGrad(DiffDoub0 flux[], DiffDoub0 tGrad[], double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 tmpAr[10];

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	vecToAr(tmpAr, pre.globTemp, 0, 10);
	matMul(tGrad, tmpAr, dNdx, 1, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vecToAr(tmpAr, pre.layerTC, i1, i1 + 9);
		matMul(flux, tmpAr, tGrad, 3, 3, 1);
	}
	else {
		vecToAr(tmpAr, pre.TCmat, 0, 9);
		matMul(flux, tmpAr, tGrad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::dFluxTGraddT(vector<DiffDoub0>& dFdT, vector<DiffDoub0>& dTG, double spt[], int layer, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 tmpAr1[33];
	DiffDoub0 tmpAr2[9];
	DiffDoub0 tmpAr3[30];

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	transpose(tmpAr1, dNdx, numNds, 3);
	arToVec(tmpAr1, dTG, 0, 33);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vecToAr(tmpAr2, pre.layerTC, i1, i1 + 9);
		matMul(tmpAr3, tmpAr2, tmpAr1, 3, 3, numNds);
		arToVec(tmpAr3, dFdT, 0, 30);
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

void Element::putVecToGlobMat(SparseMat& qMat, vector<DiffDoub0>& elQVec, bool forTherm, int matRow, vector<Node>& ndAr) {
	int i1;
	int ndDof = numNds*dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[i1].sortedRank;
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < totDof; i1++) {
			if (i1 < ndDof) {
				nd = nodes[dofTable[2 * i1]];
				dof = dofTable[2 * i1 + 1];
				globInd = ndAr[nd].dofIndex[dof];
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
void Element::getNdDisp(vector<DiffDoub1>& globDisp, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globDisp[i3].setVal(thisNd.displacement[i2]);
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

void Element::getNdVel(vector<DiffDoub1>& globVel, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globVel[i3].setVal(thisNd.velocity[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdAcc(vector<DiffDoub1>& globAcc, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < dofPerNd; i2++) {
			globAcc[i3].setVal(thisNd.acceleration[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdFlVel(vector<DiffDoub1>& flVel, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			flVel[i3].setVal(thisNd.flVel[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdFlVDot(vector<DiffDoub1>& flVDot, vector<Node>& ndAr) {
	int i1;
	int i2;
	int i3;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		i3 = i1;
		for (i2 = 0; i2 < 3; i2++) {
			flVDot[i3].setVal(thisNd.flVelDot[i2]);
			i3 += numNds;
		}
	}

	return;
}

void Element::getNdTemp(vector<DiffDoub1>& globTemp, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		globTemp[i1].setVal(thisNd.temperature);
	}
	return;
}

void Element::getNdTdot(vector<DiffDoub1>& globTdot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		globTdot[i1].setVal(thisNd.tempChangeRate);
	}
	return;
}

void Element::getNdFlDen(vector<DiffDoub1>& flDen, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		flDen[i1].setVal(thisNd.flDen);
	}

	return;
}

void Element::getNdFlDenDot(vector<DiffDoub1>& flDenDot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		flDenDot[i1].setVal(thisNd.flDenDot);
	}

	return;
}

void Element::getNdTurbE(vector<DiffDoub1>& turbE, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		turbE[i1].setVal(thisNd.turbE);
	}

	return;
}

void Element::getNdTurbEDot(vector<DiffDoub1>& turbEDot, vector<Node>& ndAr) {
	int i1;

	for (i1 = 0; i1 < numNds; i1++) {
		Node& thisNd = ndAr[nodes[i1]];
		turbEDot[i1].setVal(thisNd.turbEDot);
	}

	return;
}

void Element::evalN(DiffDoub1 nVec[], DiffDoub1 dNds[], double spt[]) {
	if(type == 4 || type == 400) {
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
	} else if(type == 6 || type == 600) {
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
	} else if(type == 8 || type == 81 || type == 800) {
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
	else if (type == 10 || type == 1000) {
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

void Element::getIpData(DiffDoub1 nVec[], DiffDoub1 dNdx[], DiffDoub1& detJ, vector<DiffDoub1>& locNds, double spt[]) {
	int i1;
	int i2;
	DiffDoub1 nCent[11];
	DiffDoub1 dNds[33];
	DiffDoub1 dNdsCent[33];
	DiffDoub1 jMat[9];
	DiffDoub1 tmpNds[30];
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
	vecToAr(tmpNds, locNds, 0, 30);
	matMul(jMat,tmpNds,dNds,3,numNds,3);

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
		matMul(jCent, tmpNds, dNdsCent, 3, numNds, 3);

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

void Element::getInstOri(vector<DiffDoub1>& instOriMat, vector<DiffDoub1>& locOri, vector<DiffDoub1>& globDisp, int stat) {
	// stat = 1: nonlinear geometry, element ori from nodal theta, no derivatives, 1st order diff for nodes
	// stat = 2: nonlinear geometry, 2nd order diff for DiffDoub1 version, 1st order for DiffDoub1 version
	DiffDoub1 rot[3];
	DiffDoub1 nnds;
	DiffDoub1 one;
	DiffDoub1 tmpOri[9];
	DiffDoub1 tmpInst[9];
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

	vecToAr(tmpOri, locOri, 0, 9);
	
	if(stat == 1) {
		dOridThet(tmpInst, tmpOri, rot, 0, 0);
		arToVec(tmpInst, instOriMat, 0, 9);
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*nDim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			stIndex = 144 * i1;
			dOridThet(tmpInst, tmpOri, rot, 0, 0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
			for (i2 = 1; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
			dOridThet(tmpInst,tmpOri,rot,i2,0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
			dOridThet(tmpInst,tmpOri,rot,i2,0);
			arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				stIndex = 36*i2 + 9*i6;
				dOridThet(tmpInst,tmpOri,rot,i2,i6);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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
				dOridThet(tmpInst,tmpOri,rot,i2,0);
				arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					stIndex = 144*i1 + 36*i2 + 9*i6;
					dOridThet(tmpInst,tmpOri,rot,i2,i6);
					arToVec(tmpInst, instOriMat, stIndex, stIndex + 9);
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

void Element::getInstDisp(DiffDoub1 instDisp[], vector<DiffDoub1>& globDisp, vector<DiffDoub1>& instOriMat, vector<DiffDoub1>& locOri, vector<DiffDoub1>& xGlob, bool nLGeom, int dv1, int dv2) {
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

void Element::getStressPrereq(DiffDoub1StressPrereq& pre, vector<Section>& secAr, vector<Material>& matAr, vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	int numLay;
	DiffDoub1 offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, secAr, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	getNdVel(pre.globVel, ndAr);
	getNdAcc(pre.globAcc, ndAr);
	getNdTemp(pre.globTemp, ndAr);
	getNdTdot(pre.globTdot, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		if (type != 2) {
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, secAr, dvAr);
			getLayerAngle(pre.layerAng, secAr, dvAr);
			getLayerQ(pre.layerQ, secAr, matAr, dvAr);
			getLayerD(pre.layerD, secAr, matAr, dvAr);
			getLayerThExp(pre.layerTE, secAr, matAr, dvAr);
			getLayerEinit(pre.layerE0, secAr, dvAr);
			getLayerDen(pre.layerDen, secAr, matAr, dvAr);
			getLayerCond(pre.layerTC, secAr, matAr, dvAr);
			getLayerSpecHeat(pre.layerSH, secAr, matAr, dvAr);
			getABD(pre.Cmat, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerAng, secAr);
			getShellDamp(pre.Dmat, pre.layerThk, pre.layerZ, pre.layerD, pre.layerAng, secAr);
			getShellExpLoad(pre.thermExp, pre.Einit, pre.layerThk, pre.layerZ, pre.layerQ, pre.layerTE, pre.layerE0, pre.layerAng, secAr);
			getShellMass(pre.Mmat, pre.layerThk, pre.layerZ, pre.layerDen, secAr, dvAr);
			getShellCond(pre.TCmat, pre.layerThk, pre.layerAng, pre.layerTC, secAr, dvAr);
			getShellSpecHeat(pre.SpecHeat, pre.layerThk, pre.layerSH, pre.layerDen, secAr);
		}
		else {
			getBeamStiff(pre.Cmat, secAr, matAr, dvAr);
			getBeamDamp(pre.Dmat, secAr, matAr, dvAr);
			getBeamExpLoad(pre.thermExp, pre.Einit, secAr, matAr, dvAr);
			getBeamMass(pre.Mmat, secAr, matAr, dvAr);
			getBeamCond(pre.TCmat, secAr, matAr, dvAr);
			getBeamSpecHeat(pre.SpecHeat, secAr, matAr, dvAr);
		}
	}
	else if (type == 21) {
		getFrcFldConst(pre.frcFldCoef, pre.frcFldExp, secAr, dvAr);
		getThrmFldConst(pre.thrmFldCoef, pre.refTemp, secAr, dvAr);
	}
	else if (type == 1) {
		getMassPerEl(pre.massPerEl, secAr, dvAr);
		getSpecificHeat(pre.SpecHeat, secAr, matAr, dvAr);
	}
	else {
		getSolidStiff(pre.Cmat, secAr, matAr, dvAr);
		getSolidDamp(pre.Dmat, secAr, matAr, dvAr);
		getThermalExp(pre.thermExp, pre.Einit, secAr, matAr, dvAr);
		getDensity(pre.Mmat[0], 0, secAr, matAr, dvAr);
		getConductivity(pre.TCmat, secAr, matAr, dvAr);
		getSpecificHeat(pre.SpecHeat, secAr, matAr, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getFluidPrereq(DiffDoub1FlPrereq& pre, vector<Section>& secAr, vector<Fluid>& flAr, vector<Node>& ndAr, vector<DesignVariable>& dvAr) {
	int i1;
	int i2;
	double secProp;
	
	getNdCrds(pre.globNds, ndAr, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	i2 = 3 * numNds;
	for (i1 = 0; i1 < i2; i1++) {
		pre.globNds[i1].add(pre.globDisp[i1]);
	}
	getNdVel(pre.globVel, ndAr);
	getNdFlDen(pre.flDen, ndAr);
	getNdFlVel(pre.flVel, ndAr);
	getNdTemp(pre.flTemp, ndAr);
	getNdTurbE(pre.flTurbE, ndAr);
	getNdFlDenDot(pre.flDenDot, ndAr);
	getNdFlVDot(pre.flVelDot, ndAr);
	getNdTdot(pre.flTDot, ndAr);
	getNdTurbEDot(pre.flTurbEDot, ndAr);

	Section& thisSec = secAr[sectPtr];
	Fluid& thisFl = flAr[thisSec.flPtr];

	secProp = thisSec.refTurbE;
	pre.refTurbE.setVal(secProp);
	getGenProp(pre.refTurbE, "refTurbE", dvAr);

	secProp = thisSec.gradVTurbCoef;
	pre.gradVTurbCoef.setVal(secProp);
	getGenProp(pre.refTurbE, "gradVTurbCoef", dvAr);

	secProp = thisSec.dissTurbCoef;
	pre.dissTurbCoef.setVal(secProp);
	getGenProp(pre.refTurbE, "dissTurbCoef", dvAr);

	secProp = thisFl.viscosity;
	pre.refVisc.setVal(secProp);
	getGenProp(pre.refVisc, "viscosity", dvAr);

	secProp = thisSec.denVisCoef;
	pre.denVisCoef.setVal(secProp);
	getGenProp(pre.denVisCoef, "denVisCoef", dvAr);

	secProp = thisSec.tempVisCoef;
	pre.tempVisCoef.setVal(secProp);
	getGenProp(pre.tempVisCoef, "tempVisCoef", dvAr);

	secProp = thisSec.turbVisCoef;
	pre.turbVisCoef.setVal(secProp);
	getGenProp(pre.turbVisCoef, "turbVisCoef", dvAr);

	secProp = thisSec.refEnth;
	pre.refEnth.setVal(secProp);
	getGenProp(pre.refEnth, "refEnth", dvAr);

	secProp = thisSec.enthCoef;
	pre.denEnthCoef.setVal(secProp);
	getGenProp(pre.denEnthCoef, "denEnthCoef", dvAr);

	secProp = thisSec.enthExp;
	pre.denEnthExp.setVal(secProp);
	getGenProp(pre.denEnthExp, "denEnthExp", dvAr);

	secProp = thisSec.presCoef;
	pre.denPresCoef.setVal(secProp);
	getGenProp(pre.denPresCoef, "denPresCoef", dvAr);

	secProp = thisSec.presExp;
	pre.denPresExp.setVal(secProp);
	getGenProp(pre.denPresExp, "denPresExp", dvAr);

	secProp = thisSec.refDen;
	pre.refDen.setVal(secProp);
	getGenProp(pre.refDen, "refDen", dvAr);

	secProp = thisSec.refTemp;
	pre.refTemp.setVal(secProp);
	getGenProp(pre.refTemp, "refTemp", dvAr);

	secProp = thisFl.thermCond;
	pre.thermCond.setVal(secProp);
	getGenProp(pre.thermCond, "thermCond", dvAr);

	secProp = thisFl.specHeat;
	pre.specHeat.setVal(secProp);
	getGenProp(pre.specHeat, "specHeat", dvAr);

	secProp = thisFl.idealGas;
	pre.iGConst.setVal(secProp);
	getGenProp(pre.iGConst, "iGConst", dvAr);

	return;
}

void Element::getVolume(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer, vector<Section>& secAr, vector<DesignVariable>& dvAr) {
	int i1;
	int i2;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 thk;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;
    double sTmp[3];

	if (type == 2) {
		thk.setVal(secAr[sectPtr].area); //rem: update to factor in dvars;
		for (auto& dv : designVars) {
			DesignVariable& thisDV = dvAr[dv.intDat];
			if (thisDV.category == "area") {
				tmp.setVal(dv.doubDat);
				thisDV.getValue(dvVal);
				dvVal.mult(tmp);
				thk.add(dvVal);
			}
		}
	} else if (type == 3 || type == 41) {
		thk.setVal(pre.layerThk[layer]);
	} else {
		thk.setVal(1.0);
	}

	vol.setVal(0.0);
	for (i1 = 0; i1 < numIP; i1++) {
		i2 = 3 * i1;
		vecToAr(sTmp, intPts, i2, i2 + 3);
		getIpData(nVec, dNdx, detJ, pre.locNds, sTmp);
		tmp.setVal(ipWt[i1]);
		tmp.mult(detJ);
		tmp.mult(thk);
		vol.add(tmp);
	}

	return;
}

void Element::getSectionDef(DiffDoub1 secDef[], vector<DiffDoub1>& globDisp,  vector<DiffDoub1>& instOriMat, vector<DiffDoub1>& locOri, vector<DiffDoub1>& xGlob, DiffDoub1 dNdx[], DiffDoub1 nVec[], bool nLGeom, int dv1, int dv2) {
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

void Element::getSolidStrain(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 dNdx[], vector<DiffDoub1>& locOri, int dv1, int dv2, bool nLGeom) {
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
	DiffDoub1 tmpOri[9];
	DiffDoub1 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		strnMat[i1].setVal(0.0);
	}
	if(dv1 + dv2 == -2) {
		vecToAr(tmpOri, locOri, 0, 9);
		matMul(uxL,tmpOri,ux,3,3,3);
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
	DiffDoub1 tmpAr[60];

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
		vecToAr(tmpAr, pre.globDisp, 0, 60);
		matMul(ux, tmpAr, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		for (i1 = 0; i1 < 6; i1++) {
			tmp.setVal(pre.thermExp[i1]);
			tmp.mult(ipTemp);
			adjStn[i1].setVal(strain[i1]);
			adjStn[i1].sub(tmp);
			adjStn[i1].sub(pre.Einit[i1]);
		}
		vecToAr(tmpAr, pre.Cmat, 0, 36);
		matMul(stress, tmpAr, adjStn, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(vector<DiffDoub1>& dsdU, vector<DiffDoub1>& dedU, vector<DiffDoub1>& dsdT, double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre) {
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
	DiffDoub1 tmpAr[60];
	DiffDoub1 tmpAr2[6];

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
		vecToAr(tmpAr, pre.globDisp, 0, 60);
		matMul(ux, tmpAr, dNdx, 3, nDim, 3);
		vecToAr(tmpAr, pre.Cmat, 0, 36);
		for (i1 = 0; i1 < totDof; i1++) {
			getSolidStrain(dStrain, ux, dNdx, pre.locOri, i1, -1, nLGeom);
			matMul(dStress, tmpAr, dStrain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				dedU[i3].setVal(dStrain[i2]);
				dsdU[i3].setVal(dStress[i2]);
				i3 += totDof;
			}
		}
		vecToAr(tmpAr2, pre.thermExp, 0, 6);
		matMul(CTE, tmpAr, tmpAr2, 6, 6, 1);
		for (i1 = 0; i1 < 6; i1++) {
			CTE[i1].neg();
		}
		matMul(tmpAr, CTE, nVec, 6, 1, numNds);
		arToVec(tmpAr, dsdT, 0, 60);
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
	DiffDoub1 tmpAr[36];

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

	vecToAr(tmpAr, pre.Cmat, 0, 36);
	matMul(frcMom, tmpAr, def, defDim, defDim, 1);

	for (i1 = 0; i1 < 6; i1++) {
		frcMom[i1].sub(pre.Einit[i1]);
		tmp.setVal(ptTemp);
		tmp.mult(pre.thermExp[i1]);
		frcMom[i1].sub(tmp);
	}

	return;
}

void Element::dDefFrcMomdU(vector<DiffDoub1>& dDefdU, vector<DiffDoub1>& dFrcMomdU, vector<DiffDoub1>& dFrcMomdT, double spt[], bool nLGeom, DiffDoub1StressPrereq& pre) {
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
	DiffDoub1 tmpAr[60];
	DiffDoub1 tmpAr2[6];

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

	vecToAr(tmpAr2, pre.thermExp, 0, 6);
	matMul(tmpAr, tmpAr2, nVec, 6, 1, numNds);
	arToVec(tmpAr, dFrcMomdT, 0, 60);

	return;
}

void Element::getFluxTGrad(DiffDoub1 flux[], DiffDoub1 tGrad[], double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 tmpAr[10];

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	vecToAr(tmpAr, pre.globTemp, 0, 10);
	matMul(tGrad, tmpAr, dNdx, 1, numNds, 3);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vecToAr(tmpAr, pre.layerTC, i1, i1 + 9);
		matMul(flux, tmpAr, tGrad, 3, 3, 1);
	}
	else {
		vecToAr(tmpAr, pre.TCmat, 0, 9);
		matMul(flux, tmpAr, tGrad, 3, 3, 1);
	}
	flux[0].neg();
	flux[1].neg();
	flux[2].neg();

	return;
}

void Element::dFluxTGraddT(vector<DiffDoub1>& dFdT, vector<DiffDoub1>& dTG, double spt[], int layer, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 tmpAr1[33];
	DiffDoub1 tmpAr2[9];
	DiffDoub1 tmpAr3[30];

	getIpData(nVec, dNdx, detJ, pre.locNds, spt);
	transpose(tmpAr1, dNdx, numNds, 3);
	arToVec(tmpAr1, dTG, 0, 33);
	if (type == 41 || type == 3) {
		i1 = 9 * layer;
		vecToAr(tmpAr2, pre.layerTC, i1, i1 + 9);
		matMul(tmpAr3, tmpAr2, tmpAr1, 3, 3, numNds);
		arToVec(tmpAr3, dFdT, 0, 30);
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

void Element::putVecToGlobMat(SparseMat& qMat, vector<DiffDoub1>& elQVec, bool forTherm, int matRow, vector<Node>& ndAr) {
	int i1;
	int ndDof = numNds*dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[i1].sortedRank;
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}
	else {
		for (i1 = 0; i1 < totDof; i1++) {
			if (i1 < ndDof) {
				nd = nodes[dofTable[2 * i1]];
				dof = dofTable[2 * i1 + 1];
				globInd = ndAr[nd].dofIndex[dof];
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
 
 
 
 
 
 
void Element::getElVec(vector<double>& elVec, vector<double>& globVec, bool forTherm, bool intnl, vector<Node>& ndAr) {
	int i1;
	int i2;
	int nd;
	int dof;
	int globInd;
	int ndDof;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[nd].sortedRank;
			elVec[i1] = globVec[globInd];
		}
	}
	else {
		ndDof = numNds * dofPerNd;
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = nodes[dofTable[i2]];
			dof = dofTable[i2 + 1];
			globInd = ndAr[nd].dofIndex[dof];
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

void Element::addToGlobVec(vector<double>& elVec, vector<double>& globVec, bool forTherm, bool intnl, vector<Node>& ndAr) {
	int i1;
	int i2;
	int nd;
	int dof;
	int globInd;
	int ndDof;

	if (forTherm) {
		for (i1 = 0; i1 < numNds; i1++) {
			nd = nodes[i1];
			globInd = ndAr[nd].sortedRank;
			globVec[globInd] += elVec[i1];
		}
	}
	else {
		ndDof = numNds * dofPerNd;
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = nodes[dofTable[i2]];
			dof = dofTable[i2 + 1];
			globInd = ndAr[nd].dofIndex[dof];
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

