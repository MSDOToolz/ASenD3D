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

void Element::getRuk(Doub Rvec[], double dRdu[], bool getMatrix, bool nLGeom, NdPt ndAr[], DVPt dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int totDof;
	Doub instOriMat[720];
	Doub xGlob[30];
	Doub xLoc[30];
	Doub locOri[9];
	Doub Cmat[81];
	Doub globDisp[60];
	
	Doub nVec[11];
	Doub dNdx[33];
	Doub detJ;
	Doub dJwt;
	
	Doub strain[6];
	Doub stress[6];
	Doub ux[9];
	Doub secDef[9];
	Doub secFcMom[9];
	Doub dStndU[288];
	Doub CdSdU[288];
	
	Doub tmp;
	
	bool stat;
	if(nLGeom) {
		stat = false;
	} else {
		stat = true;
	}
	
	getNdCrds(xGlob,ndAr,dvAr);
	getLocOri(locOri,dvAr);
	getNdDisp(globDisp,ndAr);
	getStiffMat(Cmat,dvAr);
	
	if(dofPerNd == 6) {
		correctOrient(locOri,xGlob);
		getInstOri(instOriMat,locOri,globDisp,stat);
	}
	
	matMul(xLoc,locOri,xGlob,3,3,numNds);
	
	totDof = numNds*dofPerNd + numIntDof;
	
	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec,dNdx,detJ,xLoc,&intPts[i1][0]);
		dJwt.setVal(detJ);
		tmp.setVal(ipWt[i1]);
		dJwt.mult(tmp);
		if(dofPerNd == 3) {
			matMul(ux,globDisp,dNdx,3,nDim,3);
			getSolidStrain(strain,ux,dNdx,locOri,-1,-1,nLGeom);
			matMul(stress,Cmat,strain,6,6,1);
			for (i2 = 0; i2 < totDof; i2++) {
				getSolidStrain(strain,ux,dNdx,locOri,i2,-1,nLGeom);
				i4 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					tmp.setVal(stress[i3]);
					tmp.mult(strain[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						dStndU[i4].setVal(strain[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSolidStrain(strain,ux,dNdx,locOri,i2,i3,nLGeom);
						for (i4 = 0; i4 < 6; i4++) {
							dRdu[i5]+= stress[i4].val*strain[i4].val*dJwt.val;
						}
						dRdu[i6] = dRdu[i5];
						i5++;
						i6+= totDof;
					}
				}
			}
		} else {
			getSectionDef(secDef,globDisp,instOriMat,locOri,xGlob,dNdx,nVec,-1,-1);
			matMul(secFcMom,Cmat,secDef,defDim,defDim,1);
			for (i2 = 0; i2 < totDof; i2++) {
				getSectionDef(secDef,globDisp,instOriMat,locOri,xGlob,dNdx,nVec,i2,-1);
				i4 = i2;
				for (i3 = 0; i3 < defDim; i3++) {
					tmp.setVal(secFcMom[i3]);
					tmp.mult(secDef[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						dStndU[i4].setVal(secDef[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSectionDef(secDef,globDisp,instOriMat,locOri,xGlob,dNdx,nVec,i2,i3);
						for (i4 = 0; i4 < defDim; i4++) {
							dRdu[i5]+= secFcMom[i4].val*secDef[i4].val*dJwt.val;
						}
						dRdu[i6] = dRdu[i5];
						i5++;
						i6+= totDof;
					}					
				}
			}
		}
		if(getMatrix) {
			matMul(CdSdU,Cmat,dStndU,defDim,defDim,totDof);
			i3 = defDim*totDof;
			for (i2 = 0; i2 < i3; i2++) {
				CdSdU[i2].mult(dJwt);
			}
			for (i2 = 0; i2 < totDof; i2++) {
				i5 = totDof*i2;
				for (i3 = 0; i3 < totDof; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < defDim; i4++) {
						dRdu[i5]+= dStndU[i6].val*CdSdU[i7].val;
						i6+= totDof;
						i7+= totDof;
					}
					i5++;
				}
			}
		}		
	}
	
	return;
}

void Element::condenseMat(double mat[]) {
	return;
}

void Element::getRu(Doub globR[], SparseMat& globdRdu, bool getMatrix, bool dyn, bool nLGeom, NdPt ndAr[], DVPt dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int globInd;
	int globInd2;
	int ndDof;
	int totDof;
	Doub Rvec[33];
	double dRdu[1089];
	
	ndDof = numNds*dofPerNd;
	totDof = ndDof + numIntDof;
	for (i1 = 0; i1 < totDof; i1++) {
		Rvec[i1].setVal(0.0);
		if(getMatrix) {
			i3 = totDof*i1;
			for (i2 = 0; i2 < totDof; i2++) {
				dRdu[i3] = 0.0;
				i3++;
			}
		}
	}
	
	getRuk(Rvec, dRdu, getMatrix, nLGeom, ndAr[], dvAr[]);
	
	if(dyn) {
	}
	
	if(numIntDof > 0) {
		i2 = ndDof;
		for (i1 = 0; i1 < numIntDof; i1++) {
			internalRu[i1].setVal(Rvec[i2]);
			i2++;
		}
		if(getMatrix) {
			i4 = 0;
			for (i1 = 0; i1 < totDof; i1++) {
				i3 = i1*totDof + ndDof;
				for (i2 = ndDof; i2 < totDof; i2++) {
					internalMat[i4] = dRdu[i3];
					i3++;
					i4++;
				}
			}
		}
		//Condense Matrix
	}
	
	for (i1 = 0; i1 < totDof; i1++) {
		nd = dofTable[i1][0];
		dof = dofTable[i1][1];
		if(nd < numNds) {
			nd = nodes[nd];
			globInd = ndAr[nd].ptr->getDofIndex(dof);
			globR[globInd].add(Rvec[i1])
			if(getMatrix) {
				i3 = i1*totDof;
				for (i2 = 0; i2 < ndDof; i2++) {
					nd2 = nodes[dofTable[i2][0]];
					dof2 = dofTable[i2][1];
					globInd2 = ndAr[nd2].ptr->getDofIndex(dof2);
					globdRdu.addEntry(globInd,globInd2,dRdu[i3])
					i3++;
				}
			}
		} else {
			i2 = i1 - numNds*dofPerNd;
			internalRu[i2].setVal(Rvec[i1]);
		}
	}
	
	return;
}

//end dup
