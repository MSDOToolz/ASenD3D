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


Element::Element(int newType, int newLabel, int newNodes[]) {
	type = newType;
	label = newLabel;
	int i1;
	int i2;
	int i3;
	int i4;
	if(type == 4) {
		numNds = 4;
		dofPerNd = 3;
		numIntDof = 0;
		nDim = 4;
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
	
	for (i1 = 0; i1 < numNds; i1++) {
		nodes[i1] = newNodes[i1];
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
	} else {
		internalDisp = new double[numIntDof];
	}
	
	sectPtr = NULL;
	nextEl = NULL;
	
	return;
}

void Element::setNext(Element *newEl) {
	nextEl = newEl;
}

//dup1

void Element::getStiffMat(Doub& Cmat, DVPt& dvAr) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material *matPt;
	double *stiffMat;
	double *elastic;
	Doub stiffMatDV[21];
	Doub elasticDV[9];
	Doub Smat[36];
	Doub dvVal;
	Doub coef;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *thisDVpt;
	string dvCat;
	
	if(dofPerNd == 3) {
		matPt = sectPtr->getMaterial();
		stiffMat = matPt->getStiffMat();
		if(stiffMat[0] > 0.0) {
			for (i1 = 0; i1 < 21; i1++) {
				stiffMatDV[i1].setVal(stiffMat[i1]);
			}
			thisDV = designVars.getFirst();
			thisCoef = dvCoef.getFirst();
			while(thisDV) {
				dvInd = thisDV->value;
				thisDVpt = dvAr[dvInd].ptr;
				if(thisDVpt->getCategory() == "stiffnessMat") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					stiffMatDV[dvComp-1].add(dvVal);
				}
				thisDV = thisDV->next;
				thisCoef = thisCoef->next;
			}
			i3 = 0;
			for (i1 = 0; i1 < 6; i1++) {
				i4 = 7*i1;
				i5 = 7*i1;
				for (i2 = i1; i2 < 6; i2++) {
					i4 = 6*i1 + i2;
					i5 = 6*i2 + i1;
					Cmat[i4].setVal(stiffMatDV[i3]);
					Cmat[i5].setVal(stiffMatDV[i3]);
					i4++;
					i5+= 6;
					i3++;
				}
			}
		} else {
			elastic = matPt->getElastic();
			for (i1 = 0; i1 < 9; i1++) {
				elasticDV[i1].setVal(elastic[i1]);
			}
			thisDV = designVars.getFirst();
			thisCoef = dvCoef.getFirst();
			while(thisDV) {
				dvInd = thisDV->value;
				thisDVpt = dvAr[dvInd].ptr;
				dvCat = thisDVpt->getCategory;
				if(dvCat == "modulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					elasticDV[dvComp-1].add(dvVal);
				} else if(dvCat == "poissonRatio") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					elasticDV[dvComp+2].add(dvVal);
				} else if(dvCat == "shearModulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					elasticDV[dvComp+5].add(dvVal);
				}
				thisDV = thisDV->next;
				thisCoef = thisCoef->next;
			}
			for (i1 = 0; i1 < 36; i1++) {
				Smat[i1].setVal(0.0);
			}
			
			Smat[0].setVal(1.0);
			Smat[0].dvd(elasticDV[0]);
			Smat[1].setVal(elasticDV[3]);
			Smat[1].neg();
			Smat[1].dvd(elasticDV[0]);
			Smat[2].setVal(elasticDV[4]);
			Smat[2].neg();
			Smat[2].dvd(elasticDV[0]);
			Smat[6].setVal(Smat[1]);
			Smat[7].setVal(1.0);
			Smat[7].dvd(elasticDV[1]);
			Smat[8].setVal(elasticDV[5]);
			Smat[8].neg();
			Smat[8].dvd(elasticDV[1]);
			Smat[12].setVal(Smat[2]);
			Smat[13].setVal(Smat[8]);
			Smat[14].setVal(1.0);
			Smat[14].dvd(elasticDV[2]);
			Smat[21].setVal(1.0);
			Smat[21].dvd(elasticDV[6]);
			Smat[28].setVal(1.0);
			Smat[28].dvd(elasticDV[7]);
			Smat[35].setVal(1.0);
			Smat[35].dvd(elasticDV[8]);
			
			getDetInv(&coef,Cmat,&Smat,6,0);
		}
	} else {
	}
	
	
	return;
}

void Element::getNdCrds(Doub& xGlob, NdPt& ndAr, DVPt& dvAr) {
	int i1;
	int i2;
	Doub ndCrd[3];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]].ptr;
		nPtr->getCrd(ndCrd,dvAr);
		xGlob[i1].setVal(ndCrd[0]);
		xGlob[i1+numNds].setVal(ndCrd[1]);
		xGlob[i1+2*numNds].setVal(ndCrd[2]);
	}
	
	return;
}

void Element::getNdDisp(Doub& globDisp, NdPt& ndAr) {
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
			i4 = dofTable[i2][0];
			i5 = dofTable[i2][1];
			i6 = nDim*i5 + i4;
			globDisp[i6].setVal(internalDisp[i3]);
			i2++;
		}
	}
	
	return;
}

void Element::evalN(Doub& nVec, Doub& dNds, double spt[]) {
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
		
	    nVec[3].setVal(0.5*(1.0-spt[0]-spt[1])*(1.0+spt[2])
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
		nVec[0].setVal(0.125*(1.0-spt[0])*(1.0-spt[1])*(1.0-spt[2])):
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

void Element::getIpData(Doub& nVec, Doub& dNdx, Doub& detJ, Doub& locNds, double spt[]) {
	int i1;
	int i2;
	int i3;
	
	Doub dNds[33];
	evalN(nVec, &dNds, spt);
	Doub jMat[9];
	Doub jInv[9];
	Doub zDir;
	Doub tmp;
	
	matMul(&jMat,locNds,&dNds,3,numNds,3);
	if(type == 41 || type == 3) {
		zDir.setVal(jMat[0]);
		zDir.mult(jMat[4]);
		tmp.setVal(jMat[3]);
		tmp.mult(jMat[1]);
		zDir.sub(tmp);
		if(zDir.val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0)
		}
	} else if(type == 2) {
		jMat[4].setVal(1.0);
		if(jMat[0].val > 0.0) {
			jMat[8].setVal(1.0);
		} else {
			jMat[8].setVal(-1.0);
		}
	}
	
	getDetInv(detJ,&jInv,&jMat,3,0);
	
	matMul(dNdx,&dNds,&jInv,nDim,3,3);
	
	return;
}

void Element::getInstOri(Doub& instOriMat, Doub& locOri, Doub& globDisp, bool stat) {
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
		instOriMat[i1].setValue(0.0);
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
		dOridThet(&instOriMat[0],locOri,&rot,0,0);
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*Dim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,&rot,i2,0);
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
			dOridThet(&instOriMat[stIndex],locOri,&rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
		}
		for (i1 = 1; i1 <= numNds; i1++) {
			i2 = 3*Dim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,&rot,i2,0);
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
			dOridThet(&instOriMat[stIndex],locOri,&rot,i2,0);
			i3 = stIndex;
			i4 = 9*i2;
			for (i5 = 0; i5 < 9; i5++) {
				instOriMat[i4].setVal(instOriMat[i3]);
				i3++;
				i4++;
			}
			for (i6 = i2; i6 < 4; i6++) {
				stIndex = 36*i2 + 9*i6;
				dOridThet(&instOriMat[stIndex],locOri,&rot,i2,i6);
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
			i2 = 3*Dim + i1 - 1;
			rot[0].setVal(globDisp[i2]);
			rot[1].setVal(globDisp[i2+nDim]);
			rot[2].setVal(globDisp[i2+2*nDim]);
			for (i2 = 0; i2 < 4; i2++) {
				stIndex = 144*i1 + 36*i2;
				dOridThet(&instOriMat[stIndex],locOri,&rot,i2,0);
				i3 = stIndex;
				i4 = 144*i1 + 9*i2;
				for (i5 = 0; i5 < 9; i5++) {
					instOriMat[i4].setVal(instOriMat[i3]);
					i3++;
					i4++;
				}
				for (i6 = i2; i6 < 4; i6++) {
					stIndex = 144*i1 + 36*i2 + 9*i6;
					dOridThet(&instOriMat[stIndex],locOri,&rot,i2,i6);
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

void Element::getInstDisp(Doub& instDisp, Doub& globDisp, Doub& instOriMat, Doub& locOri, Doub& xGlob, int dv1, int dv2) {
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
			i4 = dofTable[i1][0];
			i5 = dofTable[i1][1];
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
		nd = dofTable[dv1][0];
		dof = dofTable[dv1][1];
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
		nd = dofTable[dv1][0];
		dof = dofTable[dv1][1];
		nd2 = dofTable[dv2][0];
		dof2 = dofTable[dv2][1];
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

void Element::getSectionDef(Doub& secDef, Doub& globDisp,  Doub& instOriMat, Doub& locOri, Doub& xGlob, Doub& dNdx, Doub& nVec, int dv1, int dv2) {
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
	
	getInstDisp(&instDisp, globDisp, instOriMat, locOri, xGlob, int dv1, int dv2);
	
	if(dv1 + dv2 == -2) {
		matMul(&ux,&instDisp,dNdx,3,nDim,3);
		i1 = 3*nDim;
		matMul(&rx,&instDisp[i1],dNdx,3,nDim,3);
		matMul(&rot,&instDisp[i1],nVec,3,nDim,1);
	} else if((dv1+1)*(dv2+1) == 0) {
		dv1 = dv1 + dv2 + 1;
		nd = dofTable[dv1][0];
		dof = dofTable[dv1][1];
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
		nd = dofTable[dv1][0];
		dof = dofTable[dv1][1];
		nd2 = dofTable[dv2][0];
		dof2 = dofTable[dv2][1];
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

void Element::getSolidStrain(Doub& strain, Doub& ux, Doub& dNdx, Doub& locOri, int dv1, int dv2, bool nLGeom) {
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
		matMul(&uxL,locOri,ux,3,3,3);
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
			nd = dofTable[dv1][0];
			dof = dofTable[dv1][1];
		} else {
			nd = dofTable[dv2][0];
			dof = dofTable[dv2][1];
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
			nd = dofTable[dv1][0];
			dof = dofTable[dv1][1];
			nd2 = dofTable[dv2][0];
			dof2 = dofTable[dv2][1];
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

//end dup

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

void ElementList::addElement(int newType, int newLabel; int newNodes[]) {
	Element *newEl = new Element(newType,newLabel,newNodes);
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