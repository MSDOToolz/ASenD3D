#include <string>
#include "ElementClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "JobClass.h"
#include "matrixFunctions.h"

using namespace std;

const double r_2o3 = 0.666666666666666667;
const double r_1o3 = 0.333333333333333333;
const double r_1o6 = 0.166666666666666667;
const double r_pi = 3.14159265358979324;
const double r_pio180 = 0.0174532925199432958;
const double r_1ort3 = 0.577350269189625765;

void Element::condenseMat(double mat[]) {
	int stRow;
	int endRow;
	int endCol;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double xVec[33];
	double bVec[33];
	
	stRow = numNds*dofPerNd;
	endRow = stRow + numIntDof - 1;
	endCol = numIntDof - 1;
	qRFactor(internalMat, numIntDof, stRow, endRow, 0, endCol, 0);
	
	for (i1 = 0; i1 < stRow; i1++) {
		bVec[i1] = 0.0;
	}
	
	for (i1 = 0; i1 < stRow; i1++) {
		i3 = stRow;
		i4 = i1*numIntDof;
		for (i2 = 0; i2 < numIntDof; i2++) {
			bVec[i3] = internalMat[i4];
			i3++;
			i4++;
		}
		solveqRxEqb(xVec,internalMat,bVec,numIntDof,stRow,endRow,0,endCol,0);
		i4 = i1;
		for (i2 = 0; i2 < stRow; i2++) {
			i5 = i2*numIntDof;
			for (i3 = 0; i3 < numIntDof; i3++) {
				mat[i4]-= internalMat[i5]*xVec[i3];
				i5++;
			}
			i4+= (endRow+1);
		}
	}
	
	return;
}

void Element::updateExternal(double extVec[], int forSoln, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int stRow;
	int endRow;
	int endCol;
	int nd;
	int dof;
	int glbInd;
	double xVec[33];
	double bVec[33];

	if (numIntDof == 0) {
		return;
	}

    stRow = numNds*dofPerNd;
    endRow = stRow + numIntDof - 1;
	endCol = numIntDof - 1;
	for (i1 = 0; i1 < stRow; i1++) {
		bVec[i1] = 0.0;
	}
	
	if(forSoln == 1) {
		i2 = stRow;
		for (i1 = 0; i1 < numIntDof; i1++) {
			bVec[i2] = internalRu[i1].val;
			i2++;
		}
	} else {
		i2 = stRow;
		for (i1 = 0; i1 < numIntDof; i1++) {
			bVec[i2] = -internaldLdu[i1];
			i2++;
		}
	}
	
	solveqRxEqb(xVec,internalMat,bVec,numIntDof,stRow,endRow,0,endCol,0);
	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < stRow; i1++) {
		nd = nodes[dofTable[i4]];
		dof = dofTable[i4+1];
		glbInd = ndAr[nd]->getDofIndex(dof);
		for (i2 = 0; i2 < numIntDof; i2++) {
			extVec[glbInd]+= internalMat[i3]*xVec[i2];
			i3++;
		}
		i4 += 2;
	}		
	
	return;
}

void Element::updateInternal(double extVec[], int forSoln, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int stRow;
	int endRow;
	int endCol;
	int nd;
	int dof;
	int glbInd;
	double xVec[33];
	double bVec[33];

	if (numIntDof == 0) {
		return;
	}

    stRow = numNds*dofPerNd;
    endRow = stRow + numIntDof - 1;
	endCol = numIntDof - 1;
	for (i1 = 0; i1 < stRow; i1++) {
		bVec[i1] = 0.0;
	}
	
	if(forSoln == 1) {
		i2 = stRow;
		for (i1 = 0; i1 < numIntDof; i1++) {
			bVec[i2] = -internalRu[i1].val;
			i2++;
		}
	} else {
		i2 = stRow;
		for (i1 = 0; i1 < numIntDof; i1++) {
			bVec[i2] = internaldLdu[i1];
			i2++;
		}
	}
	
	for (i1 = 0; i1 < stRow; i1++) {
		nd = nodes[dofTable[2*i1]];
		dof = dofTable[2*i1+1];
		glbInd = ndAr[nd]->getDofIndex(dof);
		i3 = i1*numIntDof;
		i4 = stRow;
		for (i2 = 0; i2 < numIntDof; i2++) {
			bVec[i4]-= internalMat[i3]*extVec[glbInd];
			i3++;
			i4++;
		}
	}
	
	solveqRxEqb(xVec,internalMat,bVec,numIntDof,stRow,endRow,0,endCol,0);
	
	if(forSoln == 1) {
		for (i1 = 0; i1 < numIntDof; i1++) {
			internalDisp[i1] += xVec[i1];
		}
	} else {
		for (i1 = 0; i1 < numIntDof; i1++) {
			internalAdj[i1] = xVec[i1];
		}
	}
	
	return;
}

double Element::getIntAdjdRdD() {
	int i1;
	double prod = 0.0;
	for (i1 = 0; i1 < numIntDof; i1++) {
		prod += internalAdj[i1] * internalRu[i1].dval;
	}
	return prod;
}

//dup1

void Element::getRuk(DiffDoub0 Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int totDof;
	
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 dJwt;
	
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	DiffDoub0 thrmStn[6];
	DiffDoub0 ipTemp;
	DiffDoub0 CTE[6];
	DiffDoub0 CTEN[90];
	DiffDoub0 ux[9];
	DiffDoub0 secDef[9];
	DiffDoub0 secFcMom[9];
	
	DiffDoub0 tmp;
	
	totDof = numNds*dofPerNd + numIntDof;

	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < totDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < totDof; i2++) {
				dRdu[i3] = 0.0;
				i3++;
			}
			for (i2 = 0; i2 < numNds; i2++) {
				dRdT[i4] = 0.0;
				i4++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (dofPerNd == 6) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}
	}
	
	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec,dNdx,detJ,pre.locNds,&intPts[3*i1]);
		dJwt.setVal(detJ);
		tmp.setVal(ipWt[i1]);
		dJwt.mult(tmp);
		ipTemp.setVal(0.0);
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(pre.globTemp[i2]);
			tmp.mult(nVec[i2]);
			ipTemp.add(tmp);
		}
		if(dofPerNd == 3) {
			matMul(ux,pre.globDisp,dNdx,3,nDim,3);
			getSolidStrain(strain,ux,dNdx,pre.locOri,-1,-1,nLGeom);
			for (i2 = 0; i2 < 6; i2++) {
				thrmStn[i2].setVal(pre.thermExp[i2]);
				thrmStn[i2].mult(ipTemp);
				strain[i2].sub(thrmStn[i2]);
				strain[i2].sub(pre.Einit[i2]);
			}
			matMul(stress,pre.Cmat,strain,6,6,1);
			if (getMatrix) {
				matMul(CTE, pre.Cmat, pre.thermExp, 6, 6, 1);
				matMul(CTEN, CTE, nVec, 6, 1, numNds);
			}
			for (i2 = 0; i2 < totDof; i2++) {
				getSolidStrain(strain,ux,dNdx,pre.locOri,i2,-1,nLGeom);
				i4 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					tmp.setVal(stress[i3]);
					tmp.mult(strain[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						pre.BMat[i4].setVal(strain[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSolidStrain(strain,ux,dNdx,pre.locOri,i2,i3,nLGeom);
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
			getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,-1,-1);
			matMul(secFcMom,pre.Cmat,secDef,defDim,defDim,1);
			for (i2 = 0; i2 < 6; i2++) {
				tmp.setVal(pre.thermExp[i2]);
				tmp.mult(ipTemp);
				tmp.add(pre.Einit[i2]);
				secFcMom[i2].sub(tmp);
			}
			if (getMatrix) {
				matMul(CTEN, pre.thermExp, nVec, 6, 1, numNds);
			}
			for (i2 = 0; i2 < totDof; i2++) {
				getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,i2,-1);
				i4 = i2;
				for (i3 = 0; i3 < defDim; i3++) {
					tmp.setVal(secFcMom[i3]);
					tmp.mult(secDef[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						pre.BMat[i4].setVal(secDef[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,i2,i3);
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
			matMul(pre.CBMat,pre.Cmat,pre.BMat,defDim,defDim,totDof);
			i3 = defDim*totDof;
			for (i2 = 0; i2 < i3; i2++) {
				pre.CBMat[i2].mult(dJwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < totDof; i2++) {
				for (i3 = 0; i3 < totDof; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < defDim; i4++) {
						dRdu[i5]+= pre.BMat[i6].val*pre.CBMat[i7].val;
						i6+= totDof;
						i7+= totDof;
					}
					i5++;
				}
			}
			for (i2 = 0; i2 < 6 * numNds; i2++) {
				CTEN[i2].mult(dJwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < totDof; i2++) {
				for (i3 = 0; i3 < numNds; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < 6; i4++) {
						dRdT[i5] -= pre.BMat[i6].val * CTEN[i7].val;
						i6 += totDof;
						i7 += numNds;
					}
					i5++;
				}
			}
		}		
	}
	
	return;
}

void Element::getRum(DiffDoub0 Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int nd1;
	int dof1;
	int nd2;
	int dof2;
	int ndDof = numNds * dofPerNd;

	DiffDoub0 instDisp[60];

	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 dJwt;

	DiffDoub0 tmp;
	DiffDoub0 tmp61[6];
	DiffDoub0 tmp62[6];
	DiffDoub0 detJwt;

	DiffDoub0 Rtmp[30];

	DiffDoub0 saveM[36];

	i3 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < ndDof; i2++) {
				dRdA[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1) {
		i2 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			Rvec[i1].setVal(pre.globAcc[i1]);
			Rvec[i1].mult(pre.massPerEl);
			if (getMatrix) {
				dRdA[i2] = pre.massPerEl.val;
				i2 += 4;
			}
		}
		return;
	}
	else if (type == 21) {
		return;
	}

	if (dofPerNd == 6) {
		if (nLGeom) {
		    getInstOri(pre.instOri, pre.locOri, pre.globDisp, 1);
		}
		if (!actualProps) {
			for (i1 = 0; i1 < 36; i1++) {
				saveM[i1].setVal(pre.Mmat[i1]);
				pre.Mmat[i1].setVal(0.0);
			}
			for (i1 = 0; i1 < 36; i1 += 7) {
				pre.Mmat[i1].setVal(1.0);
			}
		}
	}
	else {
		if (!actualProps) {
			saveM[0].setVal(pre.Mmat[0]);
			pre.Mmat[0].setVal(1.0);
		}
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		detJwt.setVal(ipWt[i1]);
		detJwt.mult(detJ);
		if (dofPerNd == 6) {
			// Build matrix [dU^I/dU^g] * {N}
			i7 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd1 = dofTable[i7];
				dof1 = dofTable[i7 + 1];
				getInstDisp(instDisp, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, nLGeom, i2, -1);
				i6 = 0;
				i5 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					pre.BMat[i5].setVal(0.0);
					for (i4 = 0; i4 < numNds; i4++) {
						tmp.setVal(nVec[i4]);
						tmp.mult(instDisp[i6]);
						pre.BMat[i5].add(tmp);
						i6++;
					}
					i5 += ndDof;
					i6 += (nDim - numNds);
				}
				i7 += 2;
			}
			// tmp61 = BMat*{Acc}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				i4 = 0;
				tmp61[i2].setVal(0.0);
				for (i3 = 0; i3 < ndDof; i3++) {
					nd1 = dofTable[i4];
					dof1 = dofTable[i4 + 1];
					i6 = dof1 * numNds + nd1;
					tmp.setVal(pre.BMat[i5]);
					tmp.mult(pre.globAcc[i6]);
					tmp61[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			// tmp62 = [M][B]{A}
			matMul(tmp62, pre.Mmat, tmp61, 6, 6, 1);
			matMul(Rtmp, tmp62, pre.BMat, 1, 6, ndDof);
			// Update Rvec
			for (i2 = 0; i2 < ndDof; i2++) {
				tmp.setVal(Rtmp[i2]);
				tmp.mult(detJwt);
				Rvec[i2].add(tmp);
			}
			if (getMatrix) {
				matMul(pre.CBMat, pre.Mmat, pre.BMat, 6, 6, ndDof);
				i4 = 6 * ndDof;
				for (i2 = 0; i2 < i4; i2++) {
					pre.CBMat[i2].mult(detJwt);
				}
				i5 = 0;
				for (i2 = 0; i2 < ndDof; i2++) {
					i7 = 0;
					for (i3 = 0; i3 < ndDof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							dRdA[i5] += pre.CBMat[i6].val * pre.BMat[i7].val;
							i6 += ndDof;
							i7 += ndDof;
						}
						i5++;
					}
				}
			}
		}
		else {
			matMul(pre.BMat, pre.globAcc, nVec, 3, numNds, 1);
			for (i2 = 0; i2 < ndDof; i2++) {
				nd1 = dofTable[2 * i2];
				dof1 = dofTable[2 * i2 + 1];
				tmp.setVal(nVec[nd1]);
				tmp.mult(pre.Mmat[0]);
				tmp.mult(pre.BMat[dof1]);
				tmp.mult(detJwt);
				Rvec[i2].add(tmp);
				if (getMatrix) {
					for (i3 = 0; i3 < ndDof; i3++) {
						nd2 = dofTable[2 * i3];
						dof2 = dofTable[2 * i3 + 1];
						if (dof2 == dof1) {
							tmp.setVal(nVec[nd1]);
							tmp.mult(nVec[nd2]);
							tmp.mult(pre.Mmat[0]);
							tmp.mult(detJwt);
							i4 = i2 * ndDof + i3;
							dRdA[i4] += tmp.val;
						}
					}
				}
			}
		}
	}

	if (dofPerNd == 6) {
		if (!actualProps) {
			for (i1 = 0; i1 < 36; i1++) {
				pre.Mmat[i1].setVal(saveM[i1]);
			}
		}
	}
	else {
		if (!actualProps) {
			pre.Mmat[0].setVal(saveM[0]);
		}
	}

	return;
}

void Element::getRud(DiffDoub0 Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int ndDof = numNds * dofPerNd;
	DiffDoub0 tmp;
	DiffDoub0 Rtmp[33];
	//double dRtmp[1089];
	double* dRtmp = &pre.scratch[2508];
	//double dRdT[330];
	double* dRdT = &pre.scratch[3597];
	bool dNonZero;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 wtDetJ;
	DiffDoub0 strain[6];
	DiffDoub0 secDef[9];
	DiffDoub0 ux[9];
	DiffDoub0 Bvel[9];
	DiffDoub0 DBvel[9];

	i3 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < ndDof; i2++) {
				dRdV[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (cmd->rayDampCM > 0.0) {
		tmp.setVal(cmd->rayDampCM);
		for (i1 = 0; i1 < 36; i1++) {
			pre.Mmat[i1].mult(tmp);
		}
		for (i1 = 0; i1 < ndDof; i1++) {
			pre.globAcc[i1].setVal(pre.globVel[i1]);
		}
		getRum(Rtmp, dRtmp, getMatrix, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
		for (i1 = 0; i1 < ndDof; i1++) {
			Rvec[i1].add(Rtmp[i1]);
		}
		if (getMatrix) {
			i3 = 0;
			for (i1 = 0; i1 < ndDof; i1++) {
				for (i2 = 0; i2 < ndDof; i2++) {
					dRdV[i3] += dRtmp[i3];
					i3++;
				}
			}
		}
	}
	if (cmd->rayDampCK > 0.0) {
		tmp.setVal(cmd->rayDampCK);
		for (i1 = 0; i1 < 81; i1++) {
			pre.Cmat[i1].mult(tmp);
		}
		i3 = 0;
		i4 = 0;
		for (i1 = 0; i1 < dofPerNd; i1++) {
			for (i2 = 0; i2 < nDim; i2++) {
				if (i2 >= numNds) {
					pre.globDisp[i3].setVal(0.0);
				}
				else {
					pre.globDisp[i3].setVal(pre.globVel[i4]);
					i4++;
				}
				i3++;
			}
		}
		getRuk(Rtmp, dRtmp, dRdT, getMatrix, cmd->nonlinearGeom, pre, ndAr, dvAr);
		for (i1 = 0; i1 < ndDof; i1++) {
			Rvec[i1].add(Rtmp[i1]);
		}
		if (getMatrix) {
			i3 = 0;
			i4 = 0;
			for (i1 = 0; i1 < ndDof; i1++) {
				for (i2 = 0; i2 < ndDof; i2++) {
					dRdV[i3] += dRtmp[i4];
					i3++;
					i4++;
				}
				i3 += numIntDof;
			}
		}
	}

	// Material damping
	dNonZero = false;
	i1 = 0;
	while (!dNonZero && i1 < 36) {
		if (pre.Dmat[i1].val > 0.0) {
			dNonZero = true;
		}
		i1++;
	}

	if (dNonZero) {
		if (cmd->nonlinearGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}
		for (i1 = 0; i1 < numIP; i1++) {
			getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[i1 * 3]);
			wtDetJ.setVal(ipWt[i1]);
			wtDetJ.mult(detJ);
			if (dofPerNd == 3) {
				matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
				for (i2 = 0; i2 < ndDof; i2++) {
					getSolidStrain(strain, ux, dNdx, pre.locOri, i2, -1, cmd->nonlinearGeom);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						//i4 = i3 * ndDof + i2;
						pre.BMat[i4].setVal(strain[i3]);
						i4 += ndDof;
					}
				}
			}
			else {
				for (i2 = 0; i2 < ndDof; i2++) {
					getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, cmd->nonlinearGeom, i2, -1);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						pre.BMat[i4].setVal(secDef[i3]);
						i4 += ndDof;
					}
				}
			}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				Bvel[i2].setVal(0.0);
				i4 = 0;
				for (i3 = 0; i3 < ndDof; i3++) {
					nd = dofTable[i4];
					dof = dofTable[i4 + 1];
					tmp.setVal(pre.BMat[i5]);
					tmp.mult(pre.globVel[dof * numNds + nd]);
					Bvel[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			matMul(DBvel, pre.Dmat, Bvel, defDim, defDim, 1);
			matMul(Rtmp, DBvel, pre.BMat, 1, defDim, ndDof);
			for (i2 = 0; i2 < ndDof; i2++) {
				Rtmp[i2].mult(wtDetJ);
				Rvec[i2].add(Rtmp[i2]);
			}
			if (getMatrix) {
				matMul(pre.CBMat, pre.Dmat, pre.BMat, defDim, defDim, ndDof);
				i3 = ndDof * defDim;
				for (i2 = 0; i2 < i3; i2++) {
					pre.CBMat[i2].mult(wtDetJ);
				}
				i5 = 0;
				for (i2 = 0; i2 < ndDof; i2++) {
					for (i3 = 0; i3 < ndDof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							//i5 = i2 * ndDof + i3;
							//i6 = i4 * ndDof + i2;
							//i7 = i4 * ndDof + i3;
							dRdV[i5] += pre.BMat[i6].val * pre.CBMat[i7].val;
							i6 += ndDof;
							i7 += ndDof;
						}
						i5++;
					}
				}
			}
		}
	}

	return;
}


void Element::getRu(DiffDoub0 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int globInd;
	int globInd2;
	int ndDof;
	int totDof;
	DiffDoub0 tempAcc[30];
	DiffDoub0 Rvec[33];
	//double dRdu[1089];
	double* dRdu = &pre.scratch[0];
	//double dRdT[330];
	double* dRdT = &pre.scratch[1089];
	DiffDoub0 Rtmp[33];
	//double dRtmp[1089];
	double* dRtmp = &pre.scratch[1419];
	double c1;
	double c2;
	DiffDoub0 tmp;
	
	ndDof = numNds*dofPerNd;
	totDof = ndDof + numIntDof;

	for (i1 = 0; i1 < totDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix && i1 < ndDof) {
			i3 = i1 * totDof;
			for (i2 = 0; i2 < totDof; i2++) {
				dRdu[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		getRuFrcFld(globR, globdRdu, getMatrix, cmd, pre, ndAr);
		return;
	}

	if (type != 1) {
		getRuk(Rvec, dRdu, dRdT, getMatrix, cmd->nonlinearGeom, pre, ndAr, dvAr);
	}

	if (numIntDof > 0) {
		i2 = ndDof;
		for (i1 = 0; i1 < numIntDof; i1++) {
			internalRu[i1].setVal(Rvec[i2]);
			i2++;
		}
		if (getMatrix) {
			i4 = 0;
			for (i1 = 0; i1 < totDof; i1++) {
				i3 = i1 * totDof + ndDof;
				for (i2 = ndDof; i2 < totDof; i2++) {
					internalMat[i4] = dRdu[i3];
					i3++;
					i4++;
				}
			}
			//Condense Matrix
			condenseMat(dRdu);
		}
	}
	
	if(cmd->dynamic) {
		if (getMatrix) {
			if (cmd->lumpMass) {
				i2 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					nd = dofTable[i2];
					dof = dofTable[i2 + 1];
					i3 = dof * numNds + nd;
					tempAcc[i1].setVal(pre.globAcc[i3]);
					pre.globAcc[i3].setVal(1.0);
					i2 += 2;
				}
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				c1 = cmd->timeStep;
				c1 = -1.0 / (c1 * c1 * (cmd->newmarkBeta - cmd->newmarkGamma));
				i2 = 0;
				i3 = totDof + 1;
				for (i1 = 0; i1 < ndDof; i1++) {
					tempAcc[i1].mult(Rtmp[i1]);
					Rvec[i1].add(tempAcc[i1]);
					dRdu[i2] += c1 * Rtmp[i1].val;
					i2 += i3;
				}
			}
			else {
				getRum(Rtmp, dRtmp, true, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
				c1 = cmd->timeStep;
				c1 = -1.0 / (c1 * c1 * (cmd->newmarkBeta - cmd->newmarkGamma));
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					for (i2 = 0; i2 < ndDof; i2++) {
						dRdu[i3] += c1 * dRtmp[i4];
						i3++;
						i4++;
					}
					i3 += numIntDof;
				}
			}
			if (type != 1) {
				getRud(Rtmp, dRtmp, true, cmd, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
				c2 = cmd->timeStep * cmd->newmarkGamma * c1;
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					for (i2 = 0; i2 < ndDof; i2++) {
						dRdu[i3] += c2 * dRtmp[i4];
						i3++;
						i4++;
					}
					i3 += numIntDof;
				}
			}
		}
		else {
			if (cmd->lumpMass) {
				i2 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					nd = dofTable[i2];
					dof = dofTable[i2 + 1];
					i3 = dof * numNds + nd;
					tempAcc[i1].setVal(pre.globAcc[i3]);
					pre.globAcc[i3].setVal(1.0);
					i2 += 2;
				}
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					tempAcc[i1].mult(Rtmp[i1]);
					Rvec[i1].add(tempAcc[i1]);
				}
			}
			else {
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
			}
			if (type != 1) {
				getRud(Rtmp, dRtmp, false, cmd, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
			}
		}
	}
	
	i4 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i4]];
		dof = dofTable[i4+1];
		globInd = ndAr[nd]->getDofIndex(dof);
		globR[globInd].add(Rvec[i1]);
		if(getMatrix) {
			i3 = i1*totDof;
			i5 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd2 = nodes[dofTable[i5]];
				dof2 = dofTable[i5+1];
				globInd2 = ndAr[nd2]->getDofIndex(dof2);
				globdRdu.addEntry(globInd, globInd2, dRdu[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4+= 2;
	}
	
	return;
}

void Element::getRtk(DiffDoub0 Rvec[], double dRdT[], bool getMatrix, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 dNT[33];
	DiffDoub0 detJ;
	DiffDoub0 tmp;
	DiffDoub0 gradT[3];
	DiffDoub0 qVec[3];
	DiffDoub0 Rtmp[10];
	DiffDoub0 dRtmp[100];

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				dRdT[i3] = 0.0;
				i3++;
			}
		}
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		detJ.mult(tmp);
		matMul(gradT, pre.globTemp, dNdx, 1, numNds, 3);
		matMul(qVec, pre.TCmat, gradT, 3, 3, 1);
		matMul(Rtmp, dNdx, qVec, numNds, 3, 1);
		for (i2 = 0; i2 < numNds; i2++) {
			Rtmp[i2].mult(detJ);
			Rvec[i2].add(Rtmp[i2]);
		}
		if (getMatrix) {
			transpose(dNT, dNdx, numNds,3);
			matMul(pre.BMat, pre.TCmat, dNT, 3, 3, numNds);
			matMul(dRtmp, dNdx, pre.BMat, numNds, 3, numNds);
			i3 = numNds * numNds;
			for (i2 = 0; i2 < i3; i2++) {
				dRdT[i2] += dRtmp[i2].val*detJ.val;
			}
		}
	}

	return;
}

void Element::getRtm(DiffDoub0 Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoub0StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub0 nVec[11];
	DiffDoub0 dNdx[33];
	DiffDoub0 detJ;
	DiffDoub0 tmp;
	DiffDoub0 ptTdot;
	DiffDoub0 Rtmp[10];
	DiffDoub0 dRtmp[100];
	DiffDoub0 saveCp;

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				dRdTdot[i3] = 0.0;
				i3++;
			}
		}
	}

	if (!actualProps) {
		saveCp.setVal(pre.SpecHeat);
		pre.SpecHeat.setVal(1.0);
	}
	else if (dofPerNd == 3) {
		pre.SpecHeat.mult(pre.Mmat[0]);
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		detJ.mult(tmp);
		ptTdot.setVal(0.0);
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(nVec[i2]);
			tmp.mult(pre.globTdot[i2]);
			ptTdot.add(tmp);
		}
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(nVec[i2]);
			tmp.mult(pre.SpecHeat);
			tmp.mult(ptTdot);
			tmp.mult(detJ);
			Rvec[i2].add(tmp);
		}
		if (getMatrix) {
			matMul(dRtmp,nVec,nVec,numNds,1,numNds);
			i3 = numNds * numNds;
			for (i2 = 0; i2 < i3; i2++) {
				dRdTdot[i2] += dRtmp[i2].val * pre.SpecHeat.val * detJ.val;
			}
		}
	}

	if (!actualProps) {
		pre.SpecHeat.setVal(saveCp);
	}
	else if (dofPerNd == 3) {
		pre.SpecHeat.dvd(pre.Mmat[0]);
	}

	return;
}

void Element::getRt(DiffDoub0 globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int globInd1;
	int globInd2;
	double c1 = 1.0/(cmd->timeStep*cmd->newmarkGamma);
	DiffDoub0 Rvec[10];
	double dRdT[100];
	DiffDoub0 Rtmp[10];
	double dRtmp[100];

	getRtk(Rvec, dRdT, getMatrix, pre);

	if (cmd->dynamic) {
		getRtm(Rtmp, dRtmp, getMatrix, true, pre);
		i3 = 0;
		for (i1 = 0; i1 < numNds; i1++) {
			Rvec[i1].add(Rtmp[i1]);
			if (getMatrix) {
				for (i2 = 0; i2 < numNds; i2++) {
					dRdT[i3] += c1 * dRtmp[i3];
					i3++;
				}
			}
		}
	}

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		globInd1 = ndAr[nodes[i1]]->getSortedRank();
		globR[globInd1].add(Rvec[i1]);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				globInd2 = ndAr[nodes[i2]]->getSortedRank();
				globdRdT.addEntry(globInd1, globInd2, dRdT[i3]);
				i3++;
			}
		}
	}

	return;
}

void Element::getRuFrcFld(DiffDoub0 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int ndDof = 6;
	int totDof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int globInd;
	int globInd2;
	DiffDoub0 Rvec[6];
	double dRdU[36];
	DiffDoub0 dVec[3];
	DiffDoub0 dist;
	DiffDoub0 dvVec[3];
	DiffDoub0 dDvecdU[18];
	DiffDoub0 dDistdU[6];
	DiffDoub0 fN1[3];
	DiffDoub0 dfN1dU[18];
	DiffDoub0 dtoP;
	DiffDoub0 dtoP1;
	DiffDoub0 dtoP2;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	// Potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.setVal(pre.globNds[i2]);
		tmp.add(pre.globDisp[i2]);
		i2 -= 1;
		tmp.sub(pre.globNds[i2]);
		tmp.sub(pre.globDisp[i2]);
		dVec[i1].setVal(tmp);
	}

	dist.setVal(dVec[0]);
	dist.sqr();
	tmp.setVal(dVec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.setVal(dVec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	dDvecdU[0].setVal(-1.0);
	dDvecdU[3].setVal(1.0);
	dDvecdU[7].setVal(-1.0);
	dDvecdU[10].setVal(1.0);
	dDvecdU[14].setVal(-1.0);
	dDvecdU[17].setVal(1.0);

	matMul(dDistdU, dVec, dDvecdU, 1, 3, 6);
	tmp.setVal(1.0);
	tmp.dvd(dist);
	for (i1 = 0; i1 < 6; i1++) {
		dDistdU[i1].mult(tmp);
	}

	dtoP.setVal(dist);
	i1 = 1;
	while (i1 < pre.frcFldExp[0].val) {
		dtoP.mult(dist);
		i1++;
	}
	dtoP1.setVal(dtoP);
	dtoP1.mult(dist);
	dtoP2.setVal(dtoP1);
	dtoP2.mult(dist);

	tmp.setVal(pre.frcFldCoef[0]);
	tmp.dvd(dtoP1);
	for (i1 = 0; i1 < 3; i1++) {
		fN1[i1].setVal(dVec[i1]);
		fN1[i1].mult(tmp);
	}

	matMul(dfN1dU, dVec, dDistdU, 3, 1, 6);
	tmp.setVal(1.0);
	tmp.add(pre.frcFldExp[0]);
	tmp.neg();
	tmp.mult(pre.frcFldCoef[0]);
	tmp.dvd(dtoP2);
	for (i1 = 0; i1 < 18; i1++) {
		dfN1dU[i1].mult(tmp);
	}

	tmp.setVal(pre.frcFldCoef[0]);
	tmp.dvd(dtoP1);
	for (i1 = 0; i1 < 18; i1++) {
		tmp2.setVal(tmp);
		tmp2.mult(dDvecdU[i1]);
		dfN1dU[i1].add(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		Rvec[i1].setVal(fN1[i1]);
		Rvec[i1].neg();
		Rvec[i1 + 3].setVal(fN1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		dRdU[i1] = -dfN1dU[i1].val;
		dRdU[i1 + 18] = dfN1dU[i1].val;
	}

	// Damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.setVal(pre.globVel[i2]);
		i2 -= 1;
		tmp.sub(pre.globVel[i2]);
		dvVec[i1].setVal(tmp);
	}

	tmp.setVal(pre.frcFldCoef[1]);
	tmp.dvd(dtoP);

	for (i1 = 0; i1 < 3; i1++) {
		fN1[i1].setVal(tmp);
		fN1[i1].mult(dvVec[i1]);
	}

	/*matMul(dfN1dU, dvVec, dDistdU, 3, 1, 6);

	tmp2.setVal(pre.frcFldCoef[1]);
	tmp2.mult(pre.frcFldExp[1]);
	tmp2.neg();
	tmp2.dvd(dtoP1);

	for (i1 = 0; i1 < 18; i1++) {
		dfN1dU[i1].mult(tmp2);
	}*/

	tmp2.setVal(-cmd->newmarkGamma / (cmd->timeStep * (cmd->newmarkBeta - cmd->newmarkGamma)));
	tmp.mult(tmp2);

	for (i1 = 0; i1 < 18; i1++) {
		tmp2.setVal(tmp);
		tmp2.mult(dDvecdU[i1]);
		//dfN1dU[i1].add(tmp2);
		dfN1dU[i1].setVal(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		Rvec[i1].sub(fN1[i1]);
		Rvec[i1 + 3].add(fN1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		dRdU[i1] -= dfN1dU[i1].val;
		dRdU[i1 + 18] += dfN1dU[i1].val;
	}

	i4 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i4]];
		dof = dofTable[i4 + 1];
		globInd = ndAr[nd]->getDofIndex(dof);
		globR[globInd].add(Rvec[i1]);
		if (getMatrix) {
			i3 = i1 * totDof;
			i5 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd2 = nodes[dofTable[i5]];
				dof2 = dofTable[i5 + 1];
				globInd2 = ndAr[nd2]->getDofIndex(dof2);
				globdRdu.addEntry(globInd, globInd2, dRdU[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4 += 2;
	}

	return;
}

void Element::getAppLoad(DiffDoub0 AppLd[], Load* ldPt, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int globInd;
	int nd;
	int dof;
	int ndDof = numNds*dofPerNd;
	int numLay = sectPtr->getNumLayers();
	string ldType = ldPt->getType();
	double ldLoad[6];
	double ldCent[3];
	double ldAxis[3];
	double ldAngVel;
	double ldNorm[3];
	double dRdA[2];
	Face* fcPt;
	int* fcLocNd;
	int fcNumNds;
	bool ndInFace;
	DiffDoub0 elAppLd[30];
	DiffDoub0 totNdF[6];
	DiffDoub0 inpMag;
	DiffDoub0 vecMag;
	DiffDoub0 elVol;
	DiffDoub0 scFact;
	DiffDoub0 elCent[3];
	DiffDoub0 centToEl[3];
	DiffDoub0 axToEl[3];
	DiffDoub0 angVel2;
	DiffDoub0 nnInv;
	DiffDoub0 fcArea;
	DiffDoub0 fcNorm[3];
	DiffDoub0 trac[3];
	DiffDoub0 dp;
	DiffDoub0 tmp;

	for (i1 = 0; i1 < ndDof; i1++) {
		pre.globAcc[i1].setVal(0.0);
	}

	if (ldType == "bodyForce") {
		ldPt->getLoad(ldLoad);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			pre.globAcc[i3].setVal(ldLoad[dof]);
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, false, nLGeom, pre, ndAr, dvAr);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			dof = dofTable[i2 + 1];
			totNdF[dof].add(elAppLd[i1]);
			i2 += 2;
		}
		for (i1 = 0; i1 < dofPerNd; i1++) {
			tmp.setVal(ldLoad[i1]);
			tmp.sqr();
			inpMag.add(tmp);
			tmp.setVal(totNdF[i1]);
			tmp.sqr();
			vecMag.add(tmp);
		}
		inpMag.sqt();
		vecMag.sqt();
		if (type == 41 || type == 3) {
			elVol.setVal(0.0);
			for (i1 = 0; i1 < numLay; i1++) {
				getVolume(tmp, pre, i1);
				elVol.add(tmp);
			}
		}
		else {
			getVolume(elVol, pre, 0);
		}
		scFact.setVal(inpMag);
		scFact.mult(elVol);
		scFact.dvd(vecMag);
		for (i1 = 0; i1 < ndDof; i1++) {
			elAppLd[i1].mult(scFact);
		}
	}
	else if (ldType == "gravitational") {
		ldPt->getLoad(ldLoad);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			if (dof < 3) {
				pre.globAcc[i3].setVal(ldLoad[dof]);
			}
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, true, nLGeom, pre, ndAr, dvAr);
	}
	else if (ldType == "centrifugal") {
		ldPt->getCenter(ldCent);
		ldPt->getAxis(ldAxis);
		ldAngVel = ldPt->getAngVel();
		angVel2.setVal(ldAngVel);
		angVel2.sqr();
		i3 = 0;
		nnInv.setVal(1.0 / numNds);
		dp.setVal(0.0);
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < numNds; i2++) {
				elCent[i1].add(pre.globNds[i3]);
				i3++;
			}
			elCent[i1].mult(nnInv);
			centToEl[i1].setVal(elCent[i1]);
			tmp.setVal(ldCent[i1]);
			centToEl[i1].sub(tmp);
			tmp.setVal(ldAxis[i1]);
			tmp.mult(centToEl[i1]);
			dp.add(tmp);
		}
		for (i1 = 0; i1 < 3; i1++) {
			axToEl[i1].setVal(centToEl[i1]);
			tmp.setVal(ldAxis[i1]);
			tmp.mult(dp);
			axToEl[i1].sub(tmp);
			axToEl[i1].mult(angVel2);
		}
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			if (dof < 3) {
				pre.globAcc[i3].setVal(axToEl[dof]);
			}
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, true, nLGeom, pre, ndAr, dvAr);
	}
	else if (ldType == "surfaceTraction" || ldType == "surfacePressure") {
		ldPt->getNormDir(ldNorm);
		ldPt->getLoad(ldLoad);
		fcPt = faces.getFirst();
		while (fcPt) {
			if (fcPt->onSurface()) {
				fcPt->getAreaNormal(fcArea, fcNorm, ndAr, dvAr);
				dp.setVal(0.0);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.setVal(ldNorm[i1]);
					tmp.mult(fcNorm[i1]);
					dp.add(tmp);
				}
				tmp.setVal(r_pio180* ldPt->getNormTol());
				tmp.cs();
				if (dp.val > tmp.val) {
					if (ldType == "surfacePressure") {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].setVal(ldLoad[0]);
							trac[i1].mult(fcNorm[i1]);
							trac[i1].neg();
						}
					}
					else {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].setVal(ldLoad[i1]);
						}
					}
					fcNumNds = fcPt->getNumNds();
					fcLocNd = fcPt->getLocNds();
					for (i1 = 0; i1 < fcNumNds; i1++) {
						i4 = fcLocNd[i1];
						for (i3 = 0; i3 < 3; i3++) {
							//i4 = i3 * numNds + fcLocNd[i1];
							pre.globAcc[i4].setVal(trac[i3]);
							i4 += numNds;
						}
					}
					getRum(elAppLd, dRdA, false, false, nLGeom, pre, ndAr, dvAr);
					i2 = 0;
					for (i1 = 0; i1 < ndDof; i1++) {
						nd = dofTable[i2];
						dof = dofTable[i2 + 1];
						ndInFace = false;
						for (i3 = 0; i3 < fcNumNds; i3++) {
							if (fcLocNd[i3] == nd) {
								ndInFace = true;
							}
						}
						if (!ndInFace) {
							elAppLd[i1].setVal(0.0);
						}
						totNdF[dof].add(elAppLd[i1]);
						i2 += 2;
					}
					inpMag.setVal(0.0);
					vecMag.setVal(0.0);
					for (i1 = 0; i1 < 3; i1++) {
						tmp.setVal(totNdF[i1]);
						tmp.sqr();
						vecMag.add(tmp);
						tmp.setVal(trac[i1]);
						tmp.sqr();
						inpMag.add(tmp);
					}
					inpMag.sqt();
					vecMag.sqt();
					tmp.setVal(fcArea);
					tmp.mult(inpMag);
					tmp.dvd(vecMag);
					for (i1 = 0; i1 < ndDof; i1++) {
						elAppLd[i1].mult(tmp);
					}
				}
			}
			fcPt = fcPt->getNext();
		}
	}

	i2 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i2]];
		dof = dofTable[i2 + 1];
		globInd = ndAr[nd]->getDofIndex(dof);
		AppLd[globInd].add(elAppLd[i1]);
		i2 += 2;
	}

	return;
}

void Element::getAppThermLoad(DiffDoub0 AppLd[], Load* ldPt, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int globInd;
	int numLay = sectPtr->getNumLayers();
	string ldType = ldPt->getType();
	double ldLoad[6];
	double ldNorm[3];
	DiffDoub0 elAppLd[10];
	double dRdT[2];
	Face* fcPt;
	int fcNumNds;
	int* fcLocNds;
	bool ndInFace;
	DiffDoub0 fcArea;
	DiffDoub0 fcNorm[3];
	DiffDoub0 totHG;
	DiffDoub0 elVol;
	DiffDoub0 dp;
	DiffDoub0 tmp;

	for (i1 = 0; i1 < numNds; i1++) {
		pre.globTdot[i1].setVal(0.0);
	}

	ldPt->getLoad(ldLoad);
	if (ldType == "bodyHeatGen") {
		for (i1 = 0; i1 < numNds; i1++) {
			pre.globTdot[i1].setVal(ldLoad[0]);
		}
		getRtm(elAppLd, dRdT,false,false,pre);
		totHG.setVal(0.0);
		for (i1 = 0; i1 < numNds; i1++) {
			totHG.add(elAppLd[i1]);
		}
		if (numLay > 0) {
			elVol.setVal(0.0);
			for (i1 = 0; i1 < numLay; i1++) {
				getVolume(tmp, pre, i1);
				elVol.add(tmp);
			}
		}
		else {
			getVolume(elVol, pre, 0);
		}
		tmp.setVal(ldLoad[0]);
		tmp.mult(elVol);
		tmp.dvd(totHG);
		for (i1 = 0; i1 < numNds; i1++) {
			elAppLd[i1].mult(tmp);
		}
	}
	else if (ldType == "surfaceFlux") {
		fcPt = faces.getFirst();
		while (fcPt) {
			if (fcPt->onSurface()) {
				ldPt->getNormDir(ldNorm);
				fcPt->getAreaNormal(fcArea, fcNorm, ndAr, dvAr);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.setVal(ldNorm[i1]);
					tmp.mult(fcNorm[i1]);
					dp.add(tmp);
				}
				tmp.setVal(r_pio180 * ldPt->getNormTol());
				tmp.cs();
				if (dp.val > tmp.val) {
					fcNumNds = fcPt->getNumNds();
					fcLocNds = fcPt->getLocNds();
					for (i1 = 0; i1 < fcNumNds; i1++) {
						i2 = fcLocNds[i1];
						pre.globTdot[i2].setVal(ldLoad[0]);
					}
					getRtm(elAppLd, dRdT, false, false, pre);
					totHG.setVal(0.0);
					for (i1 = 0; i1 < numNds; i1++) {
						ndInFace = false;
						for (i2 = 0; i2 < fcNumNds; i2++) {
							if (fcLocNds[i2] == i1) {
								ndInFace = true;
							}
						}
						if (!ndInFace) {
							elAppLd[i1].setVal(0.0);
						}
						totHG.add(elAppLd[i1]);
					}
					tmp.setVal(ldLoad[0]);
					tmp.mult(fcArea);
					tmp.dvd(totHG);
					for (i1 = 0; i1 < numNds; i1++) {
						elAppLd[i1].mult(tmp);
					}
				}
			}
			fcPt = fcPt->getNext();
		}
	}

	for (i1 = 0; i1 < numNds; i1++) {
		globInd = ndAr[nodes[i1]]->getSortedRank();
		AppLd[globInd].add(elAppLd[i1]);
	}

	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void Element::getRuk(DiffDoub1 Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int totDof;
	
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 dJwt;
	
	DiffDoub1 strain[6];
	DiffDoub1 stress[6];
	DiffDoub1 thrmStn[6];
	DiffDoub1 ipTemp;
	DiffDoub1 CTE[6];
	DiffDoub1 CTEN[90];
	DiffDoub1 ux[9];
	DiffDoub1 secDef[9];
	DiffDoub1 secFcMom[9];
	
	DiffDoub1 tmp;
	
	totDof = numNds*dofPerNd + numIntDof;

	i3 = 0;
	i4 = 0;
	for (i1 = 0; i1 < totDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < totDof; i2++) {
				dRdu[i3] = 0.0;
				i3++;
			}
			for (i2 = 0; i2 < numNds; i2++) {
				dRdT[i4] = 0.0;
				i4++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (dofPerNd == 6) {
		if (nLGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}
	}
	
	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec,dNdx,detJ,pre.locNds,&intPts[3*i1]);
		dJwt.setVal(detJ);
		tmp.setVal(ipWt[i1]);
		dJwt.mult(tmp);
		ipTemp.setVal(0.0);
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(pre.globTemp[i2]);
			tmp.mult(nVec[i2]);
			ipTemp.add(tmp);
		}
		if(dofPerNd == 3) {
			matMul(ux,pre.globDisp,dNdx,3,nDim,3);
			getSolidStrain(strain,ux,dNdx,pre.locOri,-1,-1,nLGeom);
			for (i2 = 0; i2 < 6; i2++) {
				thrmStn[i2].setVal(pre.thermExp[i2]);
				thrmStn[i2].mult(ipTemp);
				strain[i2].sub(thrmStn[i2]);
				strain[i2].sub(pre.Einit[i2]);
			}
			matMul(stress,pre.Cmat,strain,6,6,1);
			if (getMatrix) {
				matMul(CTE, pre.Cmat, pre.thermExp, 6, 6, 1);
				matMul(CTEN, CTE, nVec, 6, 1, numNds);
			}
			for (i2 = 0; i2 < totDof; i2++) {
				getSolidStrain(strain,ux,dNdx,pre.locOri,i2,-1,nLGeom);
				i4 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					tmp.setVal(stress[i3]);
					tmp.mult(strain[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						pre.BMat[i4].setVal(strain[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSolidStrain(strain,ux,dNdx,pre.locOri,i2,i3,nLGeom);
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
			getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,-1,-1);
			matMul(secFcMom,pre.Cmat,secDef,defDim,defDim,1);
			for (i2 = 0; i2 < 6; i2++) {
				tmp.setVal(pre.thermExp[i2]);
				tmp.mult(ipTemp);
				tmp.add(pre.Einit[i2]);
				secFcMom[i2].sub(tmp);
			}
			if (getMatrix) {
				matMul(CTEN, pre.thermExp, nVec, 6, 1, numNds);
			}
			for (i2 = 0; i2 < totDof; i2++) {
				getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,i2,-1);
				i4 = i2;
				for (i3 = 0; i3 < defDim; i3++) {
					tmp.setVal(secFcMom[i3]);
					tmp.mult(secDef[i3]);
					tmp.mult(dJwt);
					Rvec[i2].add(tmp);
					if(getMatrix) {
						pre.BMat[i4].setVal(secDef[i3]);
						i4+= totDof;
					}
				}
				if(nLGeom && getMatrix) {
					i5 = (totDof + 1)*i2;
					i6 = i5;
					for (i3 = i2; i3 < totDof; i3++) {
						getSectionDef(secDef,pre.globDisp,pre.instOri,pre.locOri,pre.globNds,dNdx,nVec,nLGeom,i2,i3);
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
			matMul(pre.CBMat,pre.Cmat,pre.BMat,defDim,defDim,totDof);
			i3 = defDim*totDof;
			for (i2 = 0; i2 < i3; i2++) {
				pre.CBMat[i2].mult(dJwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < totDof; i2++) {
				for (i3 = 0; i3 < totDof; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < defDim; i4++) {
						dRdu[i5]+= pre.BMat[i6].val*pre.CBMat[i7].val;
						i6+= totDof;
						i7+= totDof;
					}
					i5++;
				}
			}
			for (i2 = 0; i2 < 6 * numNds; i2++) {
				CTEN[i2].mult(dJwt);
			}
			i5 = 0;
			for (i2 = 0; i2 < totDof; i2++) {
				for (i3 = 0; i3 < numNds; i3++) {
					i6 = i2;
					i7 = i3;
					for (i4 = 0; i4 < 6; i4++) {
						dRdT[i5] -= pre.BMat[i6].val * CTEN[i7].val;
						i6 += totDof;
						i7 += numNds;
					}
					i5++;
				}
			}
		}		
	}
	
	return;
}

void Element::getRum(DiffDoub1 Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int nd1;
	int dof1;
	int nd2;
	int dof2;
	int ndDof = numNds * dofPerNd;

	DiffDoub1 instDisp[60];

	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 dJwt;

	DiffDoub1 tmp;
	DiffDoub1 tmp61[6];
	DiffDoub1 tmp62[6];
	DiffDoub1 detJwt;

	DiffDoub1 Rtmp[30];

	DiffDoub1 saveM[36];

	i3 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < ndDof; i2++) {
				dRdA[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1) {
		i2 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			Rvec[i1].setVal(pre.globAcc[i1]);
			Rvec[i1].mult(pre.massPerEl);
			if (getMatrix) {
				dRdA[i2] = pre.massPerEl.val;
				i2 += 4;
			}
		}
		return;
	}
	else if (type == 21) {
		return;
	}

	if (dofPerNd == 6) {
		if (nLGeom) {
		    getInstOri(pre.instOri, pre.locOri, pre.globDisp, 1);
		}
		if (!actualProps) {
			for (i1 = 0; i1 < 36; i1++) {
				saveM[i1].setVal(pre.Mmat[i1]);
				pre.Mmat[i1].setVal(0.0);
			}
			for (i1 = 0; i1 < 36; i1 += 7) {
				pre.Mmat[i1].setVal(1.0);
			}
		}
	}
	else {
		if (!actualProps) {
			saveM[0].setVal(pre.Mmat[0]);
			pre.Mmat[0].setVal(1.0);
		}
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		detJwt.setVal(ipWt[i1]);
		detJwt.mult(detJ);
		if (dofPerNd == 6) {
			// Build matrix [dU^I/dU^g] * {N}
			i7 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd1 = dofTable[i7];
				dof1 = dofTable[i7 + 1];
				getInstDisp(instDisp, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, nLGeom, i2, -1);
				i6 = 0;
				i5 = i2;
				for (i3 = 0; i3 < 6; i3++) {
					pre.BMat[i5].setVal(0.0);
					for (i4 = 0; i4 < numNds; i4++) {
						tmp.setVal(nVec[i4]);
						tmp.mult(instDisp[i6]);
						pre.BMat[i5].add(tmp);
						i6++;
					}
					i5 += ndDof;
					i6 += (nDim - numNds);
				}
				i7 += 2;
			}
			// tmp61 = BMat*{Acc}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				i4 = 0;
				tmp61[i2].setVal(0.0);
				for (i3 = 0; i3 < ndDof; i3++) {
					nd1 = dofTable[i4];
					dof1 = dofTable[i4 + 1];
					i6 = dof1 * numNds + nd1;
					tmp.setVal(pre.BMat[i5]);
					tmp.mult(pre.globAcc[i6]);
					tmp61[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			// tmp62 = [M][B]{A}
			matMul(tmp62, pre.Mmat, tmp61, 6, 6, 1);
			matMul(Rtmp, tmp62, pre.BMat, 1, 6, ndDof);
			// Update Rvec
			for (i2 = 0; i2 < ndDof; i2++) {
				tmp.setVal(Rtmp[i2]);
				tmp.mult(detJwt);
				Rvec[i2].add(tmp);
			}
			if (getMatrix) {
				matMul(pre.CBMat, pre.Mmat, pre.BMat, 6, 6, ndDof);
				i4 = 6 * ndDof;
				for (i2 = 0; i2 < i4; i2++) {
					pre.CBMat[i2].mult(detJwt);
				}
				i5 = 0;
				for (i2 = 0; i2 < ndDof; i2++) {
					i7 = 0;
					for (i3 = 0; i3 < ndDof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							dRdA[i5] += pre.CBMat[i6].val * pre.BMat[i7].val;
							i6 += ndDof;
							i7 += ndDof;
						}
						i5++;
					}
				}
			}
		}
		else {
			matMul(pre.BMat, pre.globAcc, nVec, 3, numNds, 1);
			for (i2 = 0; i2 < ndDof; i2++) {
				nd1 = dofTable[2 * i2];
				dof1 = dofTable[2 * i2 + 1];
				tmp.setVal(nVec[nd1]);
				tmp.mult(pre.Mmat[0]);
				tmp.mult(pre.BMat[dof1]);
				tmp.mult(detJwt);
				Rvec[i2].add(tmp);
				if (getMatrix) {
					for (i3 = 0; i3 < ndDof; i3++) {
						nd2 = dofTable[2 * i3];
						dof2 = dofTable[2 * i3 + 1];
						if (dof2 == dof1) {
							tmp.setVal(nVec[nd1]);
							tmp.mult(nVec[nd2]);
							tmp.mult(pre.Mmat[0]);
							tmp.mult(detJwt);
							i4 = i2 * ndDof + i3;
							dRdA[i4] += tmp.val;
						}
					}
				}
			}
		}
	}

	if (dofPerNd == 6) {
		if (!actualProps) {
			for (i1 = 0; i1 < 36; i1++) {
				pre.Mmat[i1].setVal(saveM[i1]);
			}
		}
	}
	else {
		if (!actualProps) {
			pre.Mmat[0].setVal(saveM[0]);
		}
	}

	return;
}

void Element::getRud(DiffDoub1 Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int nd;
	int dof;
	int ndDof = numNds * dofPerNd;
	DiffDoub1 tmp;
	DiffDoub1 Rtmp[33];
	//double dRtmp[1089];
	double* dRtmp = &pre.scratch[2508];
	//double dRdT[330];
	double* dRdT = &pre.scratch[3597];
	bool dNonZero;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 wtDetJ;
	DiffDoub1 strain[6];
	DiffDoub1 secDef[9];
	DiffDoub1 ux[9];
	DiffDoub1 Bvel[9];
	DiffDoub1 DBvel[9];

	i3 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < ndDof; i2++) {
				dRdV[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 1 || type == 21) {
		return;
	}

	if (cmd->rayDampCM > 0.0) {
		tmp.setVal(cmd->rayDampCM);
		for (i1 = 0; i1 < 36; i1++) {
			pre.Mmat[i1].mult(tmp);
		}
		for (i1 = 0; i1 < ndDof; i1++) {
			pre.globAcc[i1].setVal(pre.globVel[i1]);
		}
		getRum(Rtmp, dRtmp, getMatrix, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
		for (i1 = 0; i1 < ndDof; i1++) {
			Rvec[i1].add(Rtmp[i1]);
		}
		if (getMatrix) {
			i3 = 0;
			for (i1 = 0; i1 < ndDof; i1++) {
				for (i2 = 0; i2 < ndDof; i2++) {
					dRdV[i3] += dRtmp[i3];
					i3++;
				}
			}
		}
	}
	if (cmd->rayDampCK > 0.0) {
		tmp.setVal(cmd->rayDampCK);
		for (i1 = 0; i1 < 81; i1++) {
			pre.Cmat[i1].mult(tmp);
		}
		i3 = 0;
		i4 = 0;
		for (i1 = 0; i1 < dofPerNd; i1++) {
			for (i2 = 0; i2 < nDim; i2++) {
				if (i2 >= numNds) {
					pre.globDisp[i3].setVal(0.0);
				}
				else {
					pre.globDisp[i3].setVal(pre.globVel[i4]);
					i4++;
				}
				i3++;
			}
		}
		getRuk(Rtmp, dRtmp, dRdT, getMatrix, cmd->nonlinearGeom, pre, ndAr, dvAr);
		for (i1 = 0; i1 < ndDof; i1++) {
			Rvec[i1].add(Rtmp[i1]);
		}
		if (getMatrix) {
			i3 = 0;
			i4 = 0;
			for (i1 = 0; i1 < ndDof; i1++) {
				for (i2 = 0; i2 < ndDof; i2++) {
					dRdV[i3] += dRtmp[i4];
					i3++;
					i4++;
				}
				i3 += numIntDof;
			}
		}
	}

	// Material damping
	dNonZero = false;
	i1 = 0;
	while (!dNonZero && i1 < 36) {
		if (pre.Dmat[i1].val > 0.0) {
			dNonZero = true;
		}
		i1++;
	}

	if (dNonZero) {
		if (cmd->nonlinearGeom) {
			getInstOri(pre.instOri, pre.locOri, pre.globDisp, 2);
		}
		for (i1 = 0; i1 < numIP; i1++) {
			getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[i1 * 3]);
			wtDetJ.setVal(ipWt[i1]);
			wtDetJ.mult(detJ);
			if (dofPerNd == 3) {
				matMul(ux, pre.globDisp, dNdx, 3, nDim, 3);
				for (i2 = 0; i2 < ndDof; i2++) {
					getSolidStrain(strain, ux, dNdx, pre.locOri, i2, -1, cmd->nonlinearGeom);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						//i4 = i3 * ndDof + i2;
						pre.BMat[i4].setVal(strain[i3]);
						i4 += ndDof;
					}
				}
			}
			else {
				for (i2 = 0; i2 < ndDof; i2++) {
					getSectionDef(secDef, pre.globDisp, pre.instOri, pre.locOri, pre.globNds, dNdx, nVec, cmd->nonlinearGeom, i2, -1);
					i4 = i2;
					for (i3 = 0; i3 < 6; i3++) {
						pre.BMat[i4].setVal(secDef[i3]);
						i4 += ndDof;
					}
				}
			}
			i5 = 0;
			for (i2 = 0; i2 < 6; i2++) {
				Bvel[i2].setVal(0.0);
				i4 = 0;
				for (i3 = 0; i3 < ndDof; i3++) {
					nd = dofTable[i4];
					dof = dofTable[i4 + 1];
					tmp.setVal(pre.BMat[i5]);
					tmp.mult(pre.globVel[dof * numNds + nd]);
					Bvel[i2].add(tmp);
					i4 += 2;
					i5++;
				}
			}
			matMul(DBvel, pre.Dmat, Bvel, defDim, defDim, 1);
			matMul(Rtmp, DBvel, pre.BMat, 1, defDim, ndDof);
			for (i2 = 0; i2 < ndDof; i2++) {
				Rtmp[i2].mult(wtDetJ);
				Rvec[i2].add(Rtmp[i2]);
			}
			if (getMatrix) {
				matMul(pre.CBMat, pre.Dmat, pre.BMat, defDim, defDim, ndDof);
				i3 = ndDof * defDim;
				for (i2 = 0; i2 < i3; i2++) {
					pre.CBMat[i2].mult(wtDetJ);
				}
				i5 = 0;
				for (i2 = 0; i2 < ndDof; i2++) {
					for (i3 = 0; i3 < ndDof; i3++) {
						i6 = i2;
						i7 = i3;
						for (i4 = 0; i4 < 6; i4++) {
							//i5 = i2 * ndDof + i3;
							//i6 = i4 * ndDof + i2;
							//i7 = i4 * ndDof + i3;
							dRdV[i5] += pre.BMat[i6].val * pre.CBMat[i7].val;
							i6 += ndDof;
							i7 += ndDof;
						}
						i5++;
					}
				}
			}
		}
	}

	return;
}


void Element::getRu(DiffDoub1 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int globInd;
	int globInd2;
	int ndDof;
	int totDof;
	DiffDoub1 tempAcc[30];
	DiffDoub1 Rvec[33];
	//double dRdu[1089];
	double* dRdu = &pre.scratch[0];
	//double dRdT[330];
	double* dRdT = &pre.scratch[1089];
	DiffDoub1 Rtmp[33];
	//double dRtmp[1089];
	double* dRtmp = &pre.scratch[1419];
	double c1;
	double c2;
	DiffDoub1 tmp;
	
	ndDof = numNds*dofPerNd;
	totDof = ndDof + numIntDof;

	for (i1 = 0; i1 < totDof; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix && i1 < ndDof) {
			i3 = i1 * totDof;
			for (i2 = 0; i2 < totDof; i2++) {
				dRdu[i3] = 0.0;
				i3++;
			}
		}
	}

	if (type == 21) {
		getRuFrcFld(globR, globdRdu, getMatrix, cmd, pre, ndAr);
		return;
	}

	if (type != 1) {
		getRuk(Rvec, dRdu, dRdT, getMatrix, cmd->nonlinearGeom, pre, ndAr, dvAr);
	}

	if (numIntDof > 0) {
		i2 = ndDof;
		for (i1 = 0; i1 < numIntDof; i1++) {
			internalRu[i1].setVal(Rvec[i2]);
			i2++;
		}
		if (getMatrix) {
			i4 = 0;
			for (i1 = 0; i1 < totDof; i1++) {
				i3 = i1 * totDof + ndDof;
				for (i2 = ndDof; i2 < totDof; i2++) {
					internalMat[i4] = dRdu[i3];
					i3++;
					i4++;
				}
			}
			//Condense Matrix
			condenseMat(dRdu);
		}
	}
	
	if(cmd->dynamic) {
		if (getMatrix) {
			if (cmd->lumpMass) {
				i2 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					nd = dofTable[i2];
					dof = dofTable[i2 + 1];
					i3 = dof * numNds + nd;
					tempAcc[i1].setVal(pre.globAcc[i3]);
					pre.globAcc[i3].setVal(1.0);
					i2 += 2;
				}
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				c1 = cmd->timeStep;
				c1 = -1.0 / (c1 * c1 * (cmd->newmarkBeta - cmd->newmarkGamma));
				i2 = 0;
				i3 = totDof + 1;
				for (i1 = 0; i1 < ndDof; i1++) {
					tempAcc[i1].mult(Rtmp[i1]);
					Rvec[i1].add(tempAcc[i1]);
					dRdu[i2] += c1 * Rtmp[i1].val;
					i2 += i3;
				}
			}
			else {
				getRum(Rtmp, dRtmp, true, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
				c1 = cmd->timeStep;
				c1 = -1.0 / (c1 * c1 * (cmd->newmarkBeta - cmd->newmarkGamma));
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					for (i2 = 0; i2 < ndDof; i2++) {
						dRdu[i3] += c1 * dRtmp[i4];
						i3++;
						i4++;
					}
					i3 += numIntDof;
				}
			}
			if (type != 1) {
				getRud(Rtmp, dRtmp, true, cmd, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
				c2 = cmd->timeStep * cmd->newmarkGamma * c1;
				i3 = 0;
				i4 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					for (i2 = 0; i2 < ndDof; i2++) {
						dRdu[i3] += c2 * dRtmp[i4];
						i3++;
						i4++;
					}
					i3 += numIntDof;
				}
			}
		}
		else {
			if (cmd->lumpMass) {
				i2 = 0;
				for (i1 = 0; i1 < ndDof; i1++) {
					nd = dofTable[i2];
					dof = dofTable[i2 + 1];
					i3 = dof * numNds + nd;
					tempAcc[i1].setVal(pre.globAcc[i3]);
					pre.globAcc[i3].setVal(1.0);
					i2 += 2;
				}
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					tempAcc[i1].mult(Rtmp[i1]);
					Rvec[i1].add(tempAcc[i1]);
				}
			}
			else {
				getRum(Rtmp, dRtmp, false, true, cmd->nonlinearGeom, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
			}
			if (type != 1) {
				getRud(Rtmp, dRtmp, false, cmd, pre, ndAr, dvAr);
				for (i1 = 0; i1 < ndDof; i1++) {
					Rvec[i1].add(Rtmp[i1]);
				}
			}
		}
	}
	
	i4 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i4]];
		dof = dofTable[i4+1];
		globInd = ndAr[nd]->getDofIndex(dof);
		globR[globInd].add(Rvec[i1]);
		if(getMatrix) {
			i3 = i1*totDof;
			i5 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd2 = nodes[dofTable[i5]];
				dof2 = dofTable[i5+1];
				globInd2 = ndAr[nd2]->getDofIndex(dof2);
				globdRdu.addEntry(globInd, globInd2, dRdu[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4+= 2;
	}
	
	return;
}

void Element::getRtk(DiffDoub1 Rvec[], double dRdT[], bool getMatrix, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 dNT[33];
	DiffDoub1 detJ;
	DiffDoub1 tmp;
	DiffDoub1 gradT[3];
	DiffDoub1 qVec[3];
	DiffDoub1 Rtmp[10];
	DiffDoub1 dRtmp[100];

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				dRdT[i3] = 0.0;
				i3++;
			}
		}
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		detJ.mult(tmp);
		matMul(gradT, pre.globTemp, dNdx, 1, numNds, 3);
		matMul(qVec, pre.TCmat, gradT, 3, 3, 1);
		matMul(Rtmp, dNdx, qVec, numNds, 3, 1);
		for (i2 = 0; i2 < numNds; i2++) {
			Rtmp[i2].mult(detJ);
			Rvec[i2].add(Rtmp[i2]);
		}
		if (getMatrix) {
			transpose(dNT, dNdx, numNds,3);
			matMul(pre.BMat, pre.TCmat, dNT, 3, 3, numNds);
			matMul(dRtmp, dNdx, pre.BMat, numNds, 3, numNds);
			i3 = numNds * numNds;
			for (i2 = 0; i2 < i3; i2++) {
				dRdT[i2] += dRtmp[i2].val*detJ.val;
			}
		}
	}

	return;
}

void Element::getRtm(DiffDoub1 Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoub1StressPrereq& pre) {
	int i1;
	int i2;
	int i3;
	DiffDoub1 nVec[11];
	DiffDoub1 dNdx[33];
	DiffDoub1 detJ;
	DiffDoub1 tmp;
	DiffDoub1 ptTdot;
	DiffDoub1 Rtmp[10];
	DiffDoub1 dRtmp[100];
	DiffDoub1 saveCp;

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		Rvec[i1].setVal(0.0);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				dRdTdot[i3] = 0.0;
				i3++;
			}
		}
	}

	if (!actualProps) {
		saveCp.setVal(pre.SpecHeat);
		pre.SpecHeat.setVal(1.0);
	}
	else if (dofPerNd == 3) {
		pre.SpecHeat.mult(pre.Mmat[0]);
	}

	for (i1 = 0; i1 < numIP; i1++) {
		getIpData(nVec, dNdx, detJ, pre.locNds, &intPts[3 * i1]);
		tmp.setVal(ipWt[i1]);
		detJ.mult(tmp);
		ptTdot.setVal(0.0);
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(nVec[i2]);
			tmp.mult(pre.globTdot[i2]);
			ptTdot.add(tmp);
		}
		for (i2 = 0; i2 < numNds; i2++) {
			tmp.setVal(nVec[i2]);
			tmp.mult(pre.SpecHeat);
			tmp.mult(ptTdot);
			tmp.mult(detJ);
			Rvec[i2].add(tmp);
		}
		if (getMatrix) {
			matMul(dRtmp,nVec,nVec,numNds,1,numNds);
			i3 = numNds * numNds;
			for (i2 = 0; i2 < i3; i2++) {
				dRdTdot[i2] += dRtmp[i2].val * pre.SpecHeat.val * detJ.val;
			}
		}
	}

	if (!actualProps) {
		pre.SpecHeat.setVal(saveCp);
	}
	else if (dofPerNd == 3) {
		pre.SpecHeat.dvd(pre.Mmat[0]);
	}

	return;
}

void Element::getRt(DiffDoub1 globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int globInd1;
	int globInd2;
	double c1 = 1.0/(cmd->timeStep*cmd->newmarkGamma);
	DiffDoub1 Rvec[10];
	double dRdT[100];
	DiffDoub1 Rtmp[10];
	double dRtmp[100];

	getRtk(Rvec, dRdT, getMatrix, pre);

	if (cmd->dynamic) {
		getRtm(Rtmp, dRtmp, getMatrix, true, pre);
		i3 = 0;
		for (i1 = 0; i1 < numNds; i1++) {
			Rvec[i1].add(Rtmp[i1]);
			if (getMatrix) {
				for (i2 = 0; i2 < numNds; i2++) {
					dRdT[i3] += c1 * dRtmp[i3];
					i3++;
				}
			}
		}
	}

	i3 = 0;
	for (i1 = 0; i1 < numNds; i1++) {
		globInd1 = ndAr[nodes[i1]]->getSortedRank();
		globR[globInd1].add(Rvec[i1]);
		if (getMatrix) {
			for (i2 = 0; i2 < numNds; i2++) {
				globInd2 = ndAr[nodes[i2]]->getSortedRank();
				globdRdT.addEntry(globInd1, globInd2, dRdT[i3]);
				i3++;
			}
		}
	}

	return;
}

void Element::getRuFrcFld(DiffDoub1 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int ndDof = 6;
	int totDof = 6;
	int nd;
	int nd2;
	int dof;
	int dof2;
	int globInd;
	int globInd2;
	DiffDoub1 Rvec[6];
	double dRdU[36];
	DiffDoub1 dVec[3];
	DiffDoub1 dist;
	DiffDoub1 dvVec[3];
	DiffDoub1 dDvecdU[18];
	DiffDoub1 dDistdU[6];
	DiffDoub1 fN1[3];
	DiffDoub1 dfN1dU[18];
	DiffDoub1 dtoP;
	DiffDoub1 dtoP1;
	DiffDoub1 dtoP2;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	// Potential force
	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.setVal(pre.globNds[i2]);
		tmp.add(pre.globDisp[i2]);
		i2 -= 1;
		tmp.sub(pre.globNds[i2]);
		tmp.sub(pre.globDisp[i2]);
		dVec[i1].setVal(tmp);
	}

	dist.setVal(dVec[0]);
	dist.sqr();
	tmp.setVal(dVec[1]);
	tmp.sqr();
	dist.add(tmp);
	tmp.setVal(dVec[2]);
	tmp.sqr();
	dist.add(tmp);
	dist.sqt();

	dDvecdU[0].setVal(-1.0);
	dDvecdU[3].setVal(1.0);
	dDvecdU[7].setVal(-1.0);
	dDvecdU[10].setVal(1.0);
	dDvecdU[14].setVal(-1.0);
	dDvecdU[17].setVal(1.0);

	matMul(dDistdU, dVec, dDvecdU, 1, 3, 6);
	tmp.setVal(1.0);
	tmp.dvd(dist);
	for (i1 = 0; i1 < 6; i1++) {
		dDistdU[i1].mult(tmp);
	}

	dtoP.setVal(dist);
	i1 = 1;
	while (i1 < pre.frcFldExp[0].val) {
		dtoP.mult(dist);
		i1++;
	}
	dtoP1.setVal(dtoP);
	dtoP1.mult(dist);
	dtoP2.setVal(dtoP1);
	dtoP2.mult(dist);

	tmp.setVal(pre.frcFldCoef[0]);
	tmp.dvd(dtoP1);
	for (i1 = 0; i1 < 3; i1++) {
		fN1[i1].setVal(dVec[i1]);
		fN1[i1].mult(tmp);
	}

	matMul(dfN1dU, dVec, dDistdU, 3, 1, 6);
	tmp.setVal(1.0);
	tmp.add(pre.frcFldExp[0]);
	tmp.neg();
	tmp.mult(pre.frcFldCoef[0]);
	tmp.dvd(dtoP2);
	for (i1 = 0; i1 < 18; i1++) {
		dfN1dU[i1].mult(tmp);
	}

	tmp.setVal(pre.frcFldCoef[0]);
	tmp.dvd(dtoP1);
	for (i1 = 0; i1 < 18; i1++) {
		tmp2.setVal(tmp);
		tmp2.mult(dDvecdU[i1]);
		dfN1dU[i1].add(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		Rvec[i1].setVal(fN1[i1]);
		Rvec[i1].neg();
		Rvec[i1 + 3].setVal(fN1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		dRdU[i1] = -dfN1dU[i1].val;
		dRdU[i1 + 18] = dfN1dU[i1].val;
	}

	// Damping force

	for (i1 = 0; i1 < 3; i1++) {
		i2 = i1 * 2 + 1;
		tmp.setVal(pre.globVel[i2]);
		i2 -= 1;
		tmp.sub(pre.globVel[i2]);
		dvVec[i1].setVal(tmp);
	}

	tmp.setVal(pre.frcFldCoef[1]);
	tmp.dvd(dtoP);

	for (i1 = 0; i1 < 3; i1++) {
		fN1[i1].setVal(tmp);
		fN1[i1].mult(dvVec[i1]);
	}

	/*matMul(dfN1dU, dvVec, dDistdU, 3, 1, 6);

	tmp2.setVal(pre.frcFldCoef[1]);
	tmp2.mult(pre.frcFldExp[1]);
	tmp2.neg();
	tmp2.dvd(dtoP1);

	for (i1 = 0; i1 < 18; i1++) {
		dfN1dU[i1].mult(tmp2);
	}*/

	tmp2.setVal(-cmd->newmarkGamma / (cmd->timeStep * (cmd->newmarkBeta - cmd->newmarkGamma)));
	tmp.mult(tmp2);

	for (i1 = 0; i1 < 18; i1++) {
		tmp2.setVal(tmp);
		tmp2.mult(dDvecdU[i1]);
		//dfN1dU[i1].add(tmp2);
		dfN1dU[i1].setVal(tmp2);
	}

	for (i1 = 0; i1 < 3; i1++) {
		Rvec[i1].sub(fN1[i1]);
		Rvec[i1 + 3].add(fN1[i1]);
	}

	for (i1 = 0; i1 < 18; i1++) {
		dRdU[i1] -= dfN1dU[i1].val;
		dRdU[i1 + 18] += dfN1dU[i1].val;
	}

	i4 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i4]];
		dof = dofTable[i4 + 1];
		globInd = ndAr[nd]->getDofIndex(dof);
		globR[globInd].add(Rvec[i1]);
		if (getMatrix) {
			i3 = i1 * totDof;
			i5 = 0;
			for (i2 = 0; i2 < ndDof; i2++) {
				nd2 = nodes[dofTable[i5]];
				dof2 = dofTable[i5 + 1];
				globInd2 = ndAr[nd2]->getDofIndex(dof2);
				globdRdu.addEntry(globInd, globInd2, dRdU[i3]);
				i3++;
				i5 += 2;
			}
		}
		i4 += 2;
	}

	return;
}

void Element::getAppLoad(DiffDoub1 AppLd[], Load* ldPt, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int globInd;
	int nd;
	int dof;
	int ndDof = numNds*dofPerNd;
	int numLay = sectPtr->getNumLayers();
	string ldType = ldPt->getType();
	double ldLoad[6];
	double ldCent[3];
	double ldAxis[3];
	double ldAngVel;
	double ldNorm[3];
	double dRdA[2];
	Face* fcPt;
	int* fcLocNd;
	int fcNumNds;
	bool ndInFace;
	DiffDoub1 elAppLd[30];
	DiffDoub1 totNdF[6];
	DiffDoub1 inpMag;
	DiffDoub1 vecMag;
	DiffDoub1 elVol;
	DiffDoub1 scFact;
	DiffDoub1 elCent[3];
	DiffDoub1 centToEl[3];
	DiffDoub1 axToEl[3];
	DiffDoub1 angVel2;
	DiffDoub1 nnInv;
	DiffDoub1 fcArea;
	DiffDoub1 fcNorm[3];
	DiffDoub1 trac[3];
	DiffDoub1 dp;
	DiffDoub1 tmp;

	for (i1 = 0; i1 < ndDof; i1++) {
		pre.globAcc[i1].setVal(0.0);
	}

	if (ldType == "bodyForce") {
		ldPt->getLoad(ldLoad);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			pre.globAcc[i3].setVal(ldLoad[dof]);
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, false, nLGeom, pre, ndAr, dvAr);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			dof = dofTable[i2 + 1];
			totNdF[dof].add(elAppLd[i1]);
			i2 += 2;
		}
		for (i1 = 0; i1 < dofPerNd; i1++) {
			tmp.setVal(ldLoad[i1]);
			tmp.sqr();
			inpMag.add(tmp);
			tmp.setVal(totNdF[i1]);
			tmp.sqr();
			vecMag.add(tmp);
		}
		inpMag.sqt();
		vecMag.sqt();
		if (type == 41 || type == 3) {
			elVol.setVal(0.0);
			for (i1 = 0; i1 < numLay; i1++) {
				getVolume(tmp, pre, i1);
				elVol.add(tmp);
			}
		}
		else {
			getVolume(elVol, pre, 0);
		}
		scFact.setVal(inpMag);
		scFact.mult(elVol);
		scFact.dvd(vecMag);
		for (i1 = 0; i1 < ndDof; i1++) {
			elAppLd[i1].mult(scFact);
		}
	}
	else if (ldType == "gravitational") {
		ldPt->getLoad(ldLoad);
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			if (dof < 3) {
				pre.globAcc[i3].setVal(ldLoad[dof]);
			}
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, true, nLGeom, pre, ndAr, dvAr);
	}
	else if (ldType == "centrifugal") {
		ldPt->getCenter(ldCent);
		ldPt->getAxis(ldAxis);
		ldAngVel = ldPt->getAngVel();
		angVel2.setVal(ldAngVel);
		angVel2.sqr();
		i3 = 0;
		nnInv.setVal(1.0 / numNds);
		dp.setVal(0.0);
		for (i1 = 0; i1 < 3; i1++) {
			for (i2 = 0; i2 < numNds; i2++) {
				elCent[i1].add(pre.globNds[i3]);
				i3++;
			}
			elCent[i1].mult(nnInv);
			centToEl[i1].setVal(elCent[i1]);
			tmp.setVal(ldCent[i1]);
			centToEl[i1].sub(tmp);
			tmp.setVal(ldAxis[i1]);
			tmp.mult(centToEl[i1]);
			dp.add(tmp);
		}
		for (i1 = 0; i1 < 3; i1++) {
			axToEl[i1].setVal(centToEl[i1]);
			tmp.setVal(ldAxis[i1]);
			tmp.mult(dp);
			axToEl[i1].sub(tmp);
			axToEl[i1].mult(angVel2);
		}
		i2 = 0;
		for (i1 = 0; i1 < ndDof; i1++) {
			nd = dofTable[i2];
			dof = dofTable[i2 + 1];
			i3 = dof * numNds + nd;
			if (dof < 3) {
				pre.globAcc[i3].setVal(axToEl[dof]);
			}
			i2 += 2;
		}
		getRum(elAppLd, dRdA, false, true, nLGeom, pre, ndAr, dvAr);
	}
	else if (ldType == "surfaceTraction" || ldType == "surfacePressure") {
		ldPt->getNormDir(ldNorm);
		ldPt->getLoad(ldLoad);
		fcPt = faces.getFirst();
		while (fcPt) {
			if (fcPt->onSurface()) {
				fcPt->getAreaNormal(fcArea, fcNorm, ndAr, dvAr);
				dp.setVal(0.0);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.setVal(ldNorm[i1]);
					tmp.mult(fcNorm[i1]);
					dp.add(tmp);
				}
				tmp.setVal(r_pio180* ldPt->getNormTol());
				tmp.cs();
				if (dp.val > tmp.val) {
					if (ldType == "surfacePressure") {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].setVal(ldLoad[0]);
							trac[i1].mult(fcNorm[i1]);
							trac[i1].neg();
						}
					}
					else {
						for (i1 = 0; i1 < 3; i1++) {
							trac[i1].setVal(ldLoad[i1]);
						}
					}
					fcNumNds = fcPt->getNumNds();
					fcLocNd = fcPt->getLocNds();
					for (i1 = 0; i1 < fcNumNds; i1++) {
						i4 = fcLocNd[i1];
						for (i3 = 0; i3 < 3; i3++) {
							//i4 = i3 * numNds + fcLocNd[i1];
							pre.globAcc[i4].setVal(trac[i3]);
							i4 += numNds;
						}
					}
					getRum(elAppLd, dRdA, false, false, nLGeom, pre, ndAr, dvAr);
					i2 = 0;
					for (i1 = 0; i1 < ndDof; i1++) {
						nd = dofTable[i2];
						dof = dofTable[i2 + 1];
						ndInFace = false;
						for (i3 = 0; i3 < fcNumNds; i3++) {
							if (fcLocNd[i3] == nd) {
								ndInFace = true;
							}
						}
						if (!ndInFace) {
							elAppLd[i1].setVal(0.0);
						}
						totNdF[dof].add(elAppLd[i1]);
						i2 += 2;
					}
					inpMag.setVal(0.0);
					vecMag.setVal(0.0);
					for (i1 = 0; i1 < 3; i1++) {
						tmp.setVal(totNdF[i1]);
						tmp.sqr();
						vecMag.add(tmp);
						tmp.setVal(trac[i1]);
						tmp.sqr();
						inpMag.add(tmp);
					}
					inpMag.sqt();
					vecMag.sqt();
					tmp.setVal(fcArea);
					tmp.mult(inpMag);
					tmp.dvd(vecMag);
					for (i1 = 0; i1 < ndDof; i1++) {
						elAppLd[i1].mult(tmp);
					}
				}
			}
			fcPt = fcPt->getNext();
		}
	}

	i2 = 0;
	for (i1 = 0; i1 < ndDof; i1++) {
		nd = nodes[dofTable[i2]];
		dof = dofTable[i2 + 1];
		globInd = ndAr[nd]->getDofIndex(dof);
		AppLd[globInd].add(elAppLd[i1]);
		i2 += 2;
	}

	return;
}

void Element::getAppThermLoad(DiffDoub1 AppLd[], Load* ldPt, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int globInd;
	int numLay = sectPtr->getNumLayers();
	string ldType = ldPt->getType();
	double ldLoad[6];
	double ldNorm[3];
	DiffDoub1 elAppLd[10];
	double dRdT[2];
	Face* fcPt;
	int fcNumNds;
	int* fcLocNds;
	bool ndInFace;
	DiffDoub1 fcArea;
	DiffDoub1 fcNorm[3];
	DiffDoub1 totHG;
	DiffDoub1 elVol;
	DiffDoub1 dp;
	DiffDoub1 tmp;

	for (i1 = 0; i1 < numNds; i1++) {
		pre.globTdot[i1].setVal(0.0);
	}

	ldPt->getLoad(ldLoad);
	if (ldType == "bodyHeatGen") {
		for (i1 = 0; i1 < numNds; i1++) {
			pre.globTdot[i1].setVal(ldLoad[0]);
		}
		getRtm(elAppLd, dRdT,false,false,pre);
		totHG.setVal(0.0);
		for (i1 = 0; i1 < numNds; i1++) {
			totHG.add(elAppLd[i1]);
		}
		if (numLay > 0) {
			elVol.setVal(0.0);
			for (i1 = 0; i1 < numLay; i1++) {
				getVolume(tmp, pre, i1);
				elVol.add(tmp);
			}
		}
		else {
			getVolume(elVol, pre, 0);
		}
		tmp.setVal(ldLoad[0]);
		tmp.mult(elVol);
		tmp.dvd(totHG);
		for (i1 = 0; i1 < numNds; i1++) {
			elAppLd[i1].mult(tmp);
		}
	}
	else if (ldType == "surfaceFlux") {
		fcPt = faces.getFirst();
		while (fcPt) {
			if (fcPt->onSurface()) {
				ldPt->getNormDir(ldNorm);
				fcPt->getAreaNormal(fcArea, fcNorm, ndAr, dvAr);
				for (i1 = 0; i1 < 3; i1++) {
					tmp.setVal(ldNorm[i1]);
					tmp.mult(fcNorm[i1]);
					dp.add(tmp);
				}
				tmp.setVal(r_pio180 * ldPt->getNormTol());
				tmp.cs();
				if (dp.val > tmp.val) {
					fcNumNds = fcPt->getNumNds();
					fcLocNds = fcPt->getLocNds();
					for (i1 = 0; i1 < fcNumNds; i1++) {
						i2 = fcLocNds[i1];
						pre.globTdot[i2].setVal(ldLoad[0]);
					}
					getRtm(elAppLd, dRdT, false, false, pre);
					totHG.setVal(0.0);
					for (i1 = 0; i1 < numNds; i1++) {
						ndInFace = false;
						for (i2 = 0; i2 < fcNumNds; i2++) {
							if (fcLocNds[i2] == i1) {
								ndInFace = true;
							}
						}
						if (!ndInFace) {
							elAppLd[i1].setVal(0.0);
						}
						totHG.add(elAppLd[i1]);
					}
					tmp.setVal(ldLoad[0]);
					tmp.mult(fcArea);
					tmp.dvd(totHG);
					for (i1 = 0; i1 < numNds; i1++) {
						elAppLd[i1].mult(tmp);
					}
				}
			}
			fcPt = fcPt->getNext();
		}
	}

	for (i1 = 0; i1 < numNds; i1++) {
		globInd = ndAr[nodes[i1]]->getSortedRank();
		AppLd[globInd].add(elAppLd[i1]);
	}

	return;
}

//end dup
 
//end skip 
 
 
