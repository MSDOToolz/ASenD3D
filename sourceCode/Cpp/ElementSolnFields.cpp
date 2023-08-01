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
void Element::getNdDisp(Doub globDisp[], NdPt ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	Doub ndDisp[6];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]].ptr;
		nPtr->getDisp(ndDisp);
		i3 = i1;
		for (i2 = 0; i2 < 6; i2++) {
			globDisp[i3].setVal(ndDisp[i2]);
			i3+= nDim;
		}
	}
	
	if(numIntDof > 0) {
		i2 = numNds*dofPerNd;
		for (i3 = 0; i3 < numIntDof; i3++) {
			i4 = dofTable[2*i2];
			i5 = dofTable[2*i2+1];
			i6 = nDim*i5 + i4;
			globDisp[i6].setVal(internalDisp[i3]);
			i2++;
		}
	}
	
	return;
}

//double Element::getDispSqr(Doub globDisp[]) {
//	int i1;
//	int i2;
//	int i3;
//	double dispVec[3] = {0.0,0.0,0.0};
//	double dp = 0;
//	for (i1 = 0; i1 < 3; i1++) {
//		i3 = i1 * nDim;
//		for (i2 = 0; i2 < numNds; i2++) {
//			dispVec[i1] += globDisp[i3].val;
//			i3++;
//		}
//		dispVec[i1] /= numNds;
//		dp += dispVec[i1] * dispVec[i1];
//	}
//	return dp;
//}

void Element::evalN(Doub nVec[], Doub dNds[], double spt[]) {
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
		dNds[5].setVal(-0.5);
		
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
	} else if(type == 3) {
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

void Element::getIpData(Doub nVec[], Doub dNdx[], Doub& detJ, Doub locNds[], double spt[]) {
	int i1;
	int i2;
	int i3;
	
	Doub dNds[33];
	evalN(nVec, dNds, spt);
	Doub jMat[9];
	Doub jInv[9];
	Doub xVec[3];
	Doub bVec[3];
	Doub zDir;
	Doub tmp;
	
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
	
	matMul(dNdx,dNds,jInv,nDim,3,3);
	
	return;
}

void Element::getInstOri(Doub instOriMat[], Doub locOri[], Doub globDisp[], bool stat) {
	Doub rot[3];
	Doub nnds;
	Doub one;
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
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
	i2 = 3*nDim;
	for (i1 = 0; i1 < numNds; i1++) {
		rot[0].add(globDisp[i2]);
		rot[1].add(globDisp[i2+nDim]);
		rot[2].add(globDisp[i2+2*nDim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);
	
	if(stat) {
		dOridThet(&instOriMat[0],locOri,rot,0,0);
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

void Element::getInstDisp(Doub instDisp[], Doub globDisp[], Doub instOriMat[], Doub locOri[], Doub xGlob[], int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int nd;
	int ndOInd;
	int dof;
	int dofOInd;
	int nd2;
	int nd2OInd;
	int dof2;
	int dof2OInd;
	Doub tmp;
	Doub tmp2;
	Doub nnInv;
	Doub nnInv2;
	
	i2 = 6*nDim;
	for (i1 = 0; i1 < i2; i1++) {
		instDisp[i1].setVal(0.0);
	}
	
	if(dv1 + dv2 == -2) {
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < numNds; i2++) {
				i4 = i1*nDim + i2;
				i5 = i2;
				i6 = i2;
				i7 = i1*3;
				for (i3 = 0; i3 < 3; i3++) {
					tmp.setVal(globDisp[i5]);
					tmp.add(xGlob[i6]);
					tmp.mult(instOriMat[i7]);
					tmp2.setVal(locOri[i7]);
					tmp2.mult(xGlob[i6]);
					tmp.sub(tmp2);
					instDisp[i4].add(tmp);
					i5+= nDim;
					i6+= numNds;
					i7++;
				}
			}
		}
		
		i2 = numNds*dofPerNd;
		i3 = i2 + numIntDof;
		for (i1 = i2; i1 < i3; i1++) {
			i4 = dofTable[2*i1];
			i5 = dofTable[2*i1+1];
			i6 = i5*nDim + i4;
			instDisp[i6].setVal(globDisp[i6]);
		}
		
		for (i1 = 0; i1 < numNds; i1++) {
			i3 = 3*nDim + i1; // indexes of thetax, y and z for node i1
			i4 = 4*nDim + i1;
			i5 = 5*nDim + i1;
			ndOInd = 144*(i1+1);
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
	} else if((dv1+1)*(dv2+1) == 0) {
		dv1 = dv1 + dv2 + 1;
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		if(dof < 3) {
			if(nd < numNds) {
				i1 = nd; // Index in instDisp
				i2 = dof; //Index in instOriMat
				instDisp[i1].setVal(instOriMat[i2]);
				i1 = nDim + nd;
				i2 = 3 + dof;
				instDisp[i1].setVal(instOriMat[i2]);
				i1 = 2*nDim + nd;
				i2 = 6 + dof;
				instDisp[i1].setVal(instOriMat[i2]);
			} else {
				i1 = dof*nDim + nd;
				instDisp[i1].setVal(1.0);
			}
		} else { // dof is rotation
		    nnInv.setVal(1.0/numNds);
			dofOInd = 36*(dof-2);
		    for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = i1*nDim + i2;
					i5 = i2;
					i6 = i2;
					i7 = dofOInd + i1*3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i4].add(tmp);
						i5+= nDim;
						i6+= numNds;
						i7++;
					}
				}
			}
			
			for (i1 = 0; i1 < numNds; i1++) {
				i3 = 3*nDim + i1; // indexes of thetax, y and z for node i1
				i4 = 4*nDim + i1;
				i5 = 5*nDim + i1;
				ndOInd = 144*(i1+1);
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
					if(i1 == nd) {
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
	} else {
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		nd2 = dofTable[2*dv2];
		dof2 = dofTable[2*dv2+1];
		nnInv.setVal(1.0/numNds);
		nnInv2.setVal(nnInv);
		nnInv2.sqr();
		if(dof > 2 && dof2 > 2) {
			for (i1 = 0; i1 < 3; i1++) {
				dofOInd = 36*(dof-2) + 9*(dof2-2) + 3*i1;
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = nDim*i1 + i2;
					i5 = i2;
					i6 = i2;
					i7 = dofOInd;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i4].add(tmp);
						i5+= nDim;
						i6+= numNds;
						i7++;
					}
				}
			}
			
			dofOInd = 36*(dof-2);
			dof2OInd = 9*(dof2-2);
			for (i1 = 0; i1 < numNds; i1++) {
				i3 = 3*nDim + i1; // indexes of thetax, y and z for node i1
				i4 = 4*nDim + i1;
				i5 = 5*nDim + i1;
				ndOInd = 144*(i1+1);
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
					if(i1 == nd) {
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
					if(i1 == nd2) {
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
					if(i1 == nd && i1 == nd2) {
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
		} else if(dof < 3 && dof2 < 3) {
			return;
		} else {
			if(dof > 2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i1 = 36*(dof2-2) + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = nd;
			instDisp[i2].add(tmp);
			i1 = 36*(dof2-2) + 3 + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = nDim + nd;
			instDisp[i2].add(tmp);
			i1 = 36*(dof2-2) + 6 + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = 2*nDim + nd;
			instDisp[i2].add(tmp);
		}
	}
	
	return;
}

void Element::getStressPrereq(DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]) {
	int numLay;
	Doub offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, true);
		if (type != 2) {
			numLay = sectPtr->getNumLayers();
			if (numLay > pre.currentLayLen) {
				if (pre.currentLayLen > 0) {
					pre.destroy();
				}
				pre.layerZ = new Doub[numLay];
				pre.layerThk = new Doub[numLay];
				pre.layerAng = new Doub[numLay];
				pre.layerQ = new Doub[9 * numLay];
				pre.currentLayLen = numLay;
			}
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, dvAr);
			getLayerAngle(pre.layerAng, dvAr);
			getLayerQ(pre.layerQ, dvAr);
		}
	} else {
		getSolidStiff(pre.Cmat, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getVolume(Doub& vol, DoubStressPrereq& pre, int layer) {
	int i1;
	Doub nVec[11];
	Doub dNdx[33];
	Doub detJ;
	Doub thk;
	Doub tmp;

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

void Element::getSectionDef(Doub secDef[], Doub globDisp[],  Doub instOriMat[], Doub locOri[], Doub xGlob[], Doub dNdx[], Doub nVec[], int dv1, int dv2) {
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
	Doub instDisp[60];
	Doub ux[9];
	Doub rx[9];
	Doub rot[3];
	Doub tmp;
	
	getInstDisp(instDisp, globDisp, instOriMat, locOri, xGlob, dv1, dv2);
	
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

void Element::getSolidStrain(Doub strain[], Doub ux[], Doub dNdx[], Doub locOri[], int dv1, int dv2, bool nLGeom) {
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
	Doub uxL[9];
	Doub strnMat[9];
	Doub tmp;
	
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

void Element::getStressStrain(Doub stress[], Doub strain[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre) {
	Doub nVec[11];
	Doub dNdx[33];
	Doub detJ;
	Doub ux[9];
	Doub secDef[9];
	Doub sectStrn[3];
	Doub tmp;

	if (type == 41 || type == 3) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, -1, -1);
		
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
		matMul(stress, &pre.layerQ[9 * layer], strain, 3, 3, 1);

		strain[3].setVal(strain[2]);
		strain[2].setVal(0.0);
		strain[4].setVal(0.0);
		strain[5].setVal(0.0);

		stress[3].setVal(stress[2]);
		stress[2].setVal(0.0);
		stress[4].setVal(0.0);
		stress[5].setVal(0.0);

	} else if(type != 2) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		// rem: Add temperature dependence
		matMul(stress, pre.Cmat, strain, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(Doub dsdU[], Doub dedU[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	Doub nVec[11];
	Doub dNdx[33];
	Doub detJ;
	Doub ux[9];
	Doub secDef[9];
	Doub sectStrn[3];
	Doub dStrain[6];
	Doub dStress[6];
	Doub tmp;

	int totDof = numNds * dofPerNd + numIntDof;
	if (type == 41 || type == 3) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		for (i1 = 0; i1 < totDof; i1++) {
			getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, i1, -1);

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
			dedU[totDof + i1].setVal(dStrain[1]);
			dedU[3 * totDof + i1].setVal(dStrain[2]);

			dsdU[i1].setVal(dStress[0]);
			dsdU[totDof + i1].setVal(dStress[1]);
			dsdU[3 * totDof + i1].setVal(dStress[2]);
		}

	}
	else if (type != 2) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		for (i1 = 0; i1 < totDof; i1++) {
			getSolidStrain(dStrain, ux, dNdx, pre.locOri, i1, -1, nLGeom);
			// rem: Add temperature dependence
			matMul(dStress, pre.Cmat, dStrain, 6, 6, 1);
			i3 = i1;
			for (i2 = 0; i2 < 6; i2++) {
				dedU[i3].setVal(dStrain[i2]);
				dsdU[i3].setVal(dStress[i2]);
				i3 += totDof;
			}
		}
	}

	return;
}

void Element::putVecToGlobMat(SparseMat& qMat, Doub elQVec[], int matRow, NdPt ndAr[]) {
	int i1;
	int ndDof = numNds + dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;
	for (i1 = 0; i1 < totDof; i1++) {
		if (i1 < ndDof) {
			nd = nodes[dofTable[2 * i1]];
			dof = dofTable[2 * i1 + 1];
			globInd = ndAr[nd].ptr->getDofIndex(dof);
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
		else {
			dof = i1 - ndDof;
			globInd = intDofIndex + dof;
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}

	return;
}

//end dup//skip 
 
//DiffDoub versions: 
void Element::getNdDisp(DiffDoub globDisp[], NdPt ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub ndDisp[6];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]].ptr;
		nPtr->getDisp(ndDisp);
		i3 = i1;
		for (i2 = 0; i2 < 6; i2++) {
			globDisp[i3].setVal(ndDisp[i2]);
			i3+= nDim;
		}
	}
	
	if(numIntDof > 0) {
		i2 = numNds*dofPerNd;
		for (i3 = 0; i3 < numIntDof; i3++) {
			i4 = dofTable[2*i2];
			i5 = dofTable[2*i2+1];
			i6 = nDim*i5 + i4;
			globDisp[i6].setVal(internalDisp[i3]);
			i2++;
		}
	}
	
	return;
}


void Element::evalN(DiffDoub nVec[], DiffDoub dNds[], double spt[]) {
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
		dNds[5].setVal(-0.5);
		
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
	} else if(type == 3) {
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

void Element::getIpData(DiffDoub nVec[], DiffDoub dNdx[], DiffDoub& detJ, DiffDoub locNds[], double spt[]) {
	int i1;
	int i2;
	int i3;
	
	DiffDoub dNds[33];
	evalN(nVec, dNds, spt);
	DiffDoub jMat[9];
	DiffDoub jInv[9];
	DiffDoub xVec[3];
	DiffDoub bVec[3];
	DiffDoub zDir;
	DiffDoub tmp;
	
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
	
	matMul(dNdx,dNds,jInv,nDim,3,3);
	
	return;
}

void Element::getInstOri(DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub globDisp[], bool stat) {
	DiffDoub rot[3];
	DiffDoub nnds;
	DiffDoub one;
    int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
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
	i2 = 3*nDim;
	for (i1 = 0; i1 < numNds; i1++) {
		rot[0].add(globDisp[i2]);
		rot[1].add(globDisp[i2+nDim]);
		rot[2].add(globDisp[i2+2*nDim]);
		nnds.add(one);
		i2++;
	}
	rot[0].dvd(nnds);
	rot[1].dvd(nnds);
	rot[2].dvd(nnds);
	
	if(stat) {
		dOridThet(&instOriMat[0],locOri,rot,0,0);
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

void Element::getInstDisp(DiffDoub instDisp[], DiffDoub globDisp[], DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], int dv1, int dv2) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int nd;
	int ndOInd;
	int dof;
	int dofOInd;
	int nd2;
	int nd2OInd;
	int dof2;
	int dof2OInd;
	DiffDoub tmp;
	DiffDoub tmp2;
	DiffDoub nnInv;
	DiffDoub nnInv2;
	
	i2 = 6*nDim;
	for (i1 = 0; i1 < i2; i1++) {
		instDisp[i1].setVal(0.0);
	}
	
	if(dv1 + dv2 == -2) {
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < numNds; i2++) {
				i4 = i1*nDim + i2;
				i5 = i2;
				i6 = i2;
				i7 = i1*3;
				for (i3 = 0; i3 < 3; i3++) {
					tmp.setVal(globDisp[i5]);
					tmp.add(xGlob[i6]);
					tmp.mult(instOriMat[i7]);
					tmp2.setVal(locOri[i7]);
					tmp2.mult(xGlob[i6]);
					tmp.sub(tmp2);
					instDisp[i4].add(tmp);
					i5+= nDim;
					i6+= numNds;
					i7++;
				}
			}
		}
		
		i2 = numNds*dofPerNd;
		i3 = i2 + numIntDof;
		for (i1 = i2; i1 < i3; i1++) {
			i4 = dofTable[2*i1];
			i5 = dofTable[2*i1+1];
			i6 = i5*nDim + i4;
			instDisp[i6].setVal(globDisp[i6]);
		}
		
		for (i1 = 0; i1 < numNds; i1++) {
			i4 = 4*nDim + i1;
			i5 = 5*nDim + i1;
			ndOInd = 144*(i1+1);
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
	} else if((dv1+1)*(dv2+1) == 0) {
		dv1 = dv1 + dv2 + 1;
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		if(dof < 3) {
			if(nd < numNds) {
				instDisp[i1].setVal(instOriMat[i2]);
				i1 = nDim + nd;
				i2 = 3 + dof;
				instDisp[i1].setVal(instOriMat[i2]);
				i1 = 2*nDim + nd;
				i2 = 6 + dof;
				instDisp[i1].setVal(instOriMat[i2]);
			} else {
				i1 = dof*nDim + nd;
				instDisp[i1].setVal(1.0);
			}
		    nnInv.setVal(1.0/numNds);
			dofOInd = 36*(dof-2);
		    for (i1 = 0; i1 < 3; i1++) {
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = i1*nDim + i2;
					i5 = i2;
					i6 = i2;
					i7 = dofOInd + i1*3;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv);
						instDisp[i4].add(tmp);
						i5+= nDim;
						i6+= numNds;
						i7++;
					}
				}
			}
			
			for (i1 = 0; i1 < numNds; i1++) {
				i4 = 4*nDim + i1;
				i5 = 5*nDim + i1;
				ndOInd = 144*(i1+1);
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
					if(i1 == nd) {
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
	} else {
		nd = dofTable[2*dv1];
		dof = dofTable[2*dv1+1];
		nd2 = dofTable[2*dv2];
		dof2 = dofTable[2*dv2+1];
		nnInv.setVal(1.0/numNds);
		nnInv2.setVal(nnInv);
		nnInv2.sqr();
		if(dof > 2 && dof2 > 2) {
			for (i1 = 0; i1 < 3; i1++) {
				dofOInd = 36*(dof-2) + 9*(dof2-2) + 3*i1;
				for (i2 = 0; i2 < numNds; i2++) {
					i4 = nDim*i1 + i2;
					i5 = i2;
					i6 = i2;
					i7 = dofOInd;
					for (i3 = 0; i3 < 3; i3++) {
						tmp.setVal(globDisp[i5]);
						tmp.add(xGlob[i6]);
						tmp.mult(instOriMat[i7]);
						tmp.mult(nnInv2);
						instDisp[i4].add(tmp);
						i5+= nDim;
						i6+= numNds;
						i7++;
					}
				}
			}
			
			dofOInd = 36*(dof-2);
			dof2OInd = 9*(dof2-2);
			for (i1 = 0; i1 < numNds; i1++) {
				i4 = 4*nDim + i1;
				i5 = 5*nDim + i1;
				ndOInd = 144*(i1+1);
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
					if(i1 == nd) {
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
					if(i1 == nd2) {
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
					if(i1 == nd && i1 == nd2) {
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
		} else if(dof < 3 && dof2 < 3) {
			return;
		} else {
			if(dof > 2) {
				i1 = dof;
				dof = dof2;
				dof2 = i1;
				i1 = nd;
				nd = nd2;
				nd2 = i1;
			}
			i1 = 36*(dof2-2) + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = nd;
			instDisp[i2].add(tmp);
			i1 = 36*(dof2-2) + 3 + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = nDim + nd;
			instDisp[i2].add(tmp);
			i1 = 36*(dof2-2) + 6 + dof;
			tmp.setVal(instOriMat[i1]);
			tmp.mult(nnInv);
			i2 = 2*nDim + nd;
			instDisp[i2].add(tmp);
		}
	}
	
	return;
}

void Element::getStressPrereq(DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]) {
	int numLay;
	DiffDoub offset;
	getNdCrds(pre.globNds, ndAr, dvAr);
	getLocOri(pre.locOri, dvAr);
	getNdDisp(pre.globDisp, ndAr);
	if (dofPerNd == 6) {
		correctOrient(pre.locOri, pre.globNds);
		getInstOri(pre.instOri, pre.locOri, pre.globDisp, true);
		if (type != 2) {
			numLay = sectPtr->getNumLayers();
			if (numLay > pre.currentLayLen) {
				if (pre.currentLayLen > 0) {
					pre.destroy();
				}
				pre.layerZ = new DiffDoub[numLay];
				pre.layerThk = new DiffDoub[numLay];
				pre.layerAng = new DiffDoub[numLay];
				pre.layerQ = new DiffDoub[9 * numLay];
				pre.currentLayLen = numLay;
			}
			getLayerThkZ(pre.layerThk, pre.layerZ, offset, dvAr);
			getLayerAngle(pre.layerAng, dvAr);
			getLayerQ(pre.layerQ, dvAr);
		}
	} else {
		getSolidStiff(pre.Cmat, dvAr);
	}
	matMul(pre.locNds, pre.locOri, pre.globNds, 3, 3, numNds);


	return;
}

void Element::getVolume(DiffDoub& vol, DiffDoubStressPrereq& pre, int layer) {
	int i1;
	DiffDoub nVec[11];
	DiffDoub dNdx[33];
	DiffDoub detJ;
	DiffDoub thk;
	DiffDoub tmp;

	if (type == 2) {
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

void Element::getSectionDef(DiffDoub secDef[], DiffDoub globDisp[],  DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], DiffDoub dNdx[], DiffDoub nVec[], int dv1, int dv2) {
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
	DiffDoub instDisp[60];
	DiffDoub ux[9];
	DiffDoub rx[9];
	DiffDoub rot[3];
	DiffDoub tmp;
	
	getInstDisp(instDisp, globDisp, instOriMat, locOri, xGlob, dv1, dv2);
	
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

void Element::getSolidStrain(DiffDoub strain[], DiffDoub ux[], DiffDoub dNdx[], DiffDoub locOri[], int dv1, int dv2, bool nLGeom) {
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
	DiffDoub uxL[9];
	DiffDoub strnMat[9];
	DiffDoub tmp;
	
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

void Element::getStressStrain(DiffDoub stress[], DiffDoub strain[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre) {
	DiffDoub nVec[11];
	DiffDoub dNdx[33];
	DiffDoub detJ;
	DiffDoub ux[9];
	DiffDoub secDef[9];
	DiffDoub sectStrn[3];
	DiffDoub tmp;

	if (type == 41 || type == 3) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, -1, -1);
		
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
		matMul(stress, &pre.layerQ[9 * layer], strain, 3, 3, 1);

		strain[3].setVal(strain[2]);
		strain[2].setVal(0.0);
		strain[4].setVal(0.0);
		strain[5].setVal(0.0);

		stress[3].setVal(stress[2]);
		stress[2].setVal(0.0);
		stress[4].setVal(0.0);
		stress[5].setVal(0.0);

	} else if(type != 2) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
		getSolidStrain(strain, ux, dNdx, pre.locOri, -1, -1, nLGeom);
		matMul(stress, pre.Cmat, strain, 6, 6, 1);
	}
	return;
}

void Element::dStressStraindU(DiffDoub dsdU[], DiffDoub dedU[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub nVec[11];
	DiffDoub dNdx[33];
	DiffDoub detJ;
	DiffDoub ux[9];
	DiffDoub secDef[9];
	DiffDoub sectStrn[3];
	DiffDoub dStrain[6];
	DiffDoub dStress[6];
	DiffDoub tmp;

	int totDof = numNds * dofPerNd + numIntDof;
	if (type == 41 || type == 3) {
		getIpData(nVec, dNdx, detJ, pre.locNds, spt);
		for (i1 = 0; i1 < totDof; i1++) {
			getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, i1, -1);

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
			dedU[totDof + i1].setVal(dStrain[1]);
			dedU[3 * totDof + i1].setVal(dStrain[2]);

			dsdU[i1].setVal(dStress[0]);
			dsdU[totDof + i1].setVal(dStress[1]);
			dsdU[3 * totDof + i1].setVal(dStress[2]);
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
	}

	return;
}

void Element::putVecToGlobMat(SparseMat& qMat, DiffDoub elQVec[], int matRow, NdPt ndAr[]) {
	int i1;
	int ndDof = numNds + dofPerNd;
	int totDof = ndDof + numIntDof;
	int nd;
	int dof;
	int globInd;
	for (i1 = 0; i1 < totDof; i1++) {
		if (i1 < ndDof) {
			nd = nodes[dofTable[2 * i1]];
			dof = dofTable[2 * i1 + 1];
			globInd = ndAr[nd].ptr->getDofIndex(dof);
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
		else {
			dof = i1 - ndDof;
			globInd = intDofIndex + dof;
			qMat.addEntry(matRow, globInd, elQVec[i1].val);
		}
	}

	return;
}

 
//end skip 
