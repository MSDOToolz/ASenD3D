#include <string>
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

void Element::getLayerThkZ(DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0& zOffset, DesignVariable* dvAr[]) {
	//zOffset = 1: upper Z surface is reference plane
	//zOffset = -1: lower Z surface is reference plane
	int layi = 0;
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DesignVariable* thisDVpt;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DiffDoub0 totThk;
	DiffDoub0 tmp;
	DiffDoub0 zCrd;
	DiffDoub0 zNext;
	DiffDoub0 zMid;

	totThk.setVal(0.0);
	Layer* thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		layThk[layi].setVal(thisLay->getThickness());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVpt = dvAr[thisDV->value];
			if (thisDVpt->getCategory() == "thickness" && thisDVpt->getLayer() == layi) {
				thisDVpt->getValue(dvVal);
				coef.setVal(thisCoef->value);
				dvVal.mult(coef);
				layThk[layi].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		totThk.add(layThk[layi]);
		thisLay = thisLay->getNext();
		layi++;
	}

	zOffset.setVal(sectPtr->getZOffset());
	thisDV = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisDV) {
		thisDVpt = dvAr[thisDV->value];
		if (thisDVpt->getCategory() == "zOffset") {
			thisDVpt->getValue(dvVal);
			coef.setVal(thisCoef->value);
			dvVal.mult(coef);
			zOffset.add(dvVal);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}

	tmp.setVal(1.0);
	tmp.add(zOffset);
	zCrd.setVal(-0.5);
	zCrd.mult(totThk);
	zCrd.mult(tmp);

	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		zNext.setVal(zCrd);
		zNext.add(layThk[layi]);
		tmp.setVal(zCrd);
		tmp.add(zNext);
		zMid.setVal(0.5);
		zMid.mult(tmp);
		layZ[layi].setVal(zMid);
		thisLay = thisLay->getNext();
		layi++;
		zCrd.setVal(zNext);
	}

	return;
}

void Element::getLayerQ(DiffDoub0 layQ[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	int dvInd;
	int dvComp;
	Material* matPt;
	Layer* thisLay;
	double* Emod;
	double* pr;
	double* shMod;
	DiffDoub0 modulusDV[3];
	DiffDoub0 poissonDV[3];
	DiffDoub0 shearModDV[3];
	DiffDoub0 Smat[9];
	DiffDoub0 Qmat[9];
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DiffDoub0 xVec[3];
	DiffDoub0 bVec[3];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		matPt = thisLay->getMatPt();
		Emod = matPt->getModulus();
		pr = matPt->getPoissonRatio();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			poissonDV[i1].setVal(pr[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "modulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					modulusDV[dvComp - 1].add(dvVal);
				}
				else if (dvCat == "poissonRatio") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					poissonDV[dvComp - 1].add(dvVal);
				}
				else if (dvCat == "shearModulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					shearModDV[dvComp - 1].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 9; i1++) {
			Smat[i1].setVal(0.0);
		}
		Smat[0].setVal(1.0);
		Smat[0].dvd(modulusDV[0]);
		Smat[1].setVal(poissonDV[0]);
		Smat[1].neg();
		Smat[1].dvd(modulusDV[0]);
		Smat[3].setVal(Smat[1]);
		Smat[4].setVal(1.0);
		Smat[4].dvd(modulusDV[1]);
		Smat[8].setVal(1.0);
		Smat[8].dvd(shearModDV[0]);
		getDetInv(coef, Qmat, Smat, 3, 0, xVec, bVec);

		for (i1 = 0; i1 < 9; i1++) {
			layQ[i2].setVal(Qmat[i1]);
			i2++;
		}

		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerD(DiffDoub0 layD[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int layi;
	int dvInd;
	int dvComp;
	Material* matPt;
	Layer* thisLay;
	double* dampMat;
	DiffDoub0 dampMatDV[36];
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		matPt = thisLay->getMatPt();
		dampMat = matPt->getDamping();
		for (i1 = 0; i1 < 36; i1++) {
			dampMatDV[i1].setVal(dampMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "dampingMat") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					dampMatDV[dvComp].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i3 = 1; i3 < 6; i3++) {
			i5 = 6 * i3; // Lower tri term
			i6 = i3; // Upper tri term
			for (i4 = 0; i4 < i3; i4++) {
				dampMatDV[i5].setVal(dampMatDV[i6]);
				i5++;
				i6 += 6;
			}
		}

		layD[i2].setVal(dampMatDV[0]);
		layD[i2 + 1].setVal(dampMatDV[1]);
		layD[i2 + 2].setVal(dampMatDV[3]);
		layD[i2 + 3].setVal(dampMatDV[6]);
		layD[i2 + 4].setVal(dampMatDV[7]);
		layD[i2 + 5].setVal(dampMatDV[9]);
		layD[i2 + 6].setVal(dampMatDV[18]);
		layD[i2 + 7].setVal(dampMatDV[19]);
		layD[i2 + 8].setVal(dampMatDV[21]);

		thisLay = thisLay->getNext();
		layi++;
		i2 += 9;
	}
	return;
}

void Element::getLayerAngle(DiffDoub0 layAng[], DesignVariable* dvAr[]) {
	int i2;
	int layi;
	int dvInd;
	Layer* thisLay;
	DiffDoub0 angle;
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		angle.setVal(thisLay->getAngle());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "angle") {
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					angle.add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		layAng[layi].setVal(angle);

		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerThExp(DiffDoub0 layThExp[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	double* matExp;
	DiffDoub0 tExpDV[6];
	Material* thisMat;
	Layer* thisLay = sectPtr->getFirstLayer();
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;
	
	layi = 0;
	i2 = 0;
	while (thisLay) {
		thisMat = thisLay->getMatPt();
		matExp = thisMat->getThermExp();
		for (i1 = 0; i1 < 6; i1++) {
			tExpDV[i1].setVal(matExp[i1]);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermalExp" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				tExpDV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layThExp[i2].setVal(tExpDV[0]);
		i2++;
		layThExp[i2].setVal(tExpDV[1]);
		i2++;
		layThExp[i2].setVal(tExpDV[3]);
		i2++;
		thisLay = thisLay->getNext();
		layi++;
	}
	return;
}

void Element::getLayerEinit(DiffDoub0 layEinit[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	DiffDoub0 E0DV[6];
	Material* thisMat;
	Layer* thisLay = sectPtr->getFirstLayer();
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;

	layi = 0;
	i2 = 0;
	while (thisLay) {
		thisMat = thisLay->getMatPt();
		for (i1 = 0; i1 < 6; i1++) {
			E0DV[i1].setVal(0.0);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "initialStrain" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				E0DV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layEinit[i2].setVal(E0DV[0]);
		i2++;
		layEinit[i2].setVal(E0DV[1]);
		i2++;
		layEinit[i2].setVal(E0DV[3]);
		i2++;
		thisLay = thisLay->getNext();
		layi++;
	}
	return;
}

void Element::getLayerDen(DiffDoub0 layerDen[], DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double matDen;
	DiffDoub0 denDV;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;

	layi = 0;
	while (thisLay) {
		matDen = thisLay->getMatPt()->getDensity();
		denDV.setVal(matDen);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "density" && thisDV->getLayer() == layi) {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				denDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layerDen[layi].setVal(denDV);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerCond(DiffDoub0 layCond[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double* matCond;
	DiffDoub0 condDV[6];
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;
	int dvComp;

	i1 = 0;
	layi = 0;
	while (thisLay) {
		matCond = thisLay->getMatPt()->getConductivity();
		for (i2 = 0; i2 < 6; i2++) {
			condDV[i1].setVal(matCond[i2]);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermalCond" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layCond[i1].setVal(condDV[0]);
		layCond[i1 + 1].setVal(condDV[3]);
		layCond[i1 + 2].setVal(condDV[4]);
		layCond[i1 + 3].setVal(condDV[3]);
		layCond[i1 + 4].setVal(condDV[1]);
		layCond[i1 + 5].setVal(condDV[5]);
		layCond[i1 + 6].setVal(condDV[4]);
		layCond[i1 + 7].setVal(condDV[5]);
		layCond[i1 + 8].setVal(condDV[2]);
		thisLay = thisLay->getNext();
		i1 += 9;
		layi++;
	}

	return;
}

void Element::getLayerSpecHeat(DiffDoub0 laySH[], DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double matSH;
	DiffDoub0 shDV;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;

	layi = 0;
	while (thisLay) {
		matSH = thisLay->getMatPt()->getSpecificHeat();
		shDV.setVal(matSH);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "specHeat" && thisDV->getLayer() == layi) {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				shDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		laySH[layi].setVal(shDV);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::transformStrain(DiffDoub0 stnNew[], DiffDoub0 stnOrig[], DiffDoub0& angle) {
	DiffDoub0 angleRad;
	DiffDoub0 a11;
	DiffDoub0 a12;
	DiffDoub0 a21;
	DiffDoub0 a22;
	DiffDoub0 tmp;
	DiffDoub0 Te[9];


	angleRad.setVal(r_pio180);
	angleRad.mult(angle);
	a11.setVal(angleRad);
	a11.cs();
	a21.setVal(angleRad);
	a21.sn();
	a12.setVal(a21);
	a12.neg();
	a22.setVal(a11);

	Te[0].setVal(a11);
	Te[0].sqr();
	Te[1].setVal(a12);
	Te[1].sqr();
	Te[2].setVal(a11);
	Te[2].mult(a12);
	Te[3].setVal(a21);
	Te[3].sqr();
	Te[4].setVal(a22);
	Te[4].sqr();
	Te[5].setVal(a22);
	Te[5].mult(a21);
	Te[6].setVal(2.0);
	Te[6].mult(a11);
	Te[6].mult(a21);
	Te[7].setVal(2.0);
	Te[7].mult(a12);
	Te[7].mult(a22);
	Te[8].setVal(a11);
	Te[8].mult(a22);
	tmp.setVal(a12);
	tmp.mult(a21);
	Te[8].add(tmp);

	matMul(stnNew, Te, stnOrig, 3, 3, 1);

	return;
}

void Element::transformQ(DiffDoub0 qNew[], DiffDoub0 qOrig[], DiffDoub0& angle) {
	int i1;
	DiffDoub0 angleRad;
	DiffDoub0 a11;
	DiffDoub0 a12;
	DiffDoub0 a21;
	DiffDoub0 a22;
	DiffDoub0 tmp;
	DiffDoub0 coef;
	DiffDoub0 Ts[9];
	DiffDoub0 Te[9];
	DiffDoub0 TeInv[9];
	DiffDoub0 xVec[3];
	DiffDoub0 bVec[3];


	angleRad.setVal(r_pio180);
	angleRad.mult(angle);
	a11.setVal(angleRad);
	a11.cs();
	a21.setVal(angleRad);
	a21.sn();
	a12.setVal(a21);
	a12.neg();
	a22.setVal(a11);

	Ts[0].setVal(a11);
	Ts[0].sqr();
	Ts[1].setVal(a12);
	Ts[1].sqr();
	Ts[2].setVal(2.0);
	Ts[2].mult(a11);
	Ts[2].mult(a12);
	Ts[3].setVal(a21);
	Ts[3].sqr();
	Ts[4].setVal(a22);
	Ts[4].sqr();
	Ts[5].setVal(2.0);
	Ts[5].mult(a22);
	Ts[5].mult(a21);
	Ts[6].setVal(a11);
	Ts[6].mult(a21);
	Ts[7].setVal(a12);
	Ts[7].mult(a22);
	Ts[8].setVal(a11);
	Ts[8].mult(a22);
	tmp.setVal(a12);
	tmp.mult(a21);
	Ts[8].add(tmp);

	for (i1 = 0; i1 < 9; i1++) {
		Te[i1].setVal(Ts[i1]);
	}
	tmp.setVal(0.5);
	Te[2].mult(tmp);
	Te[5].mult(tmp);
	tmp.setVal(2.0);
	Te[6].mult(tmp);
	Te[7].mult(tmp);

	getDetInv(coef, TeInv, Te, 3, 0, xVec, bVec);

	matMul(Te, qOrig, TeInv, 3, 3, 3);
	matMul(qNew, Ts, Te, 3, 3, 3);

	return;
}

void Element::getSolidStiff(DiffDoub0 Cmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* stiffMat;
	double* Emod;
	double* pr;
	double* shMod;
	DiffDoub0 modulusDV[3];
	DiffDoub0 poissonDV[3];
	DiffDoub0 shearModDV[3];
	DiffDoub0 Smat[36];
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DiffDoub0 xVec[6];
	DiffDoub0 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;


	matPt = sectPtr->getMatPtr();
	stiffMat = matPt->getStiffMat();
	if (stiffMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Cmat[i1].setVal(stiffMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "stiffnessMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Cmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i4].setVal(Cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		Emod = matPt->getModulus();
		pr = matPt->getPoissonRatio();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			poissonDV[i1].setVal(pr[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			dvCat = thisDVpt->getCategory();
			if (dvCat == "modulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				modulusDV[dvComp - 1].add(dvVal);
			}
			else if (dvCat == "poissonRatio") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				poissonDV[dvComp - 1].add(dvVal);
			}
			else if (dvCat == "shearModulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				shearModDV[dvComp - 1].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 0; i1 < 36; i1++) {
			Smat[i1].setVal(0.0);
		}

		Smat[0].setVal(1.0);
		Smat[0].dvd(modulusDV[0]);
		Smat[1].setVal(poissonDV[0]);
		Smat[1].neg();
		Smat[1].dvd(modulusDV[0]);
		Smat[2].setVal(poissonDV[1]);
		Smat[2].neg();
		Smat[2].dvd(modulusDV[0]);
		Smat[6].setVal(Smat[1]);
		Smat[7].setVal(1.0);
		Smat[7].dvd(modulusDV[1]);
		Smat[8].setVal(poissonDV[2]);
		Smat[8].neg();
		Smat[8].dvd(modulusDV[1]);
		Smat[12].setVal(Smat[2]);
		Smat[13].setVal(Smat[8]);
		Smat[14].setVal(1.0);
		Smat[14].dvd(modulusDV[2]);
		Smat[21].setVal(1.0);
		Smat[21].dvd(shearModDV[0]);
		Smat[28].setVal(1.0);
		Smat[28].dvd(shearModDV[1]);
		Smat[35].setVal(1.0);
		Smat[35].dvd(shearModDV[2]);

		getDetInv(coef, Cmat, Smat, 6, 0, xVec, bVec);
	}

	return;
}

void Element::getABD(DiffDoub0 Cmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layQ[], DiffDoub0 layAng[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int numLay;
	DiffDoub0 zMax;
	DiffDoub0 zMin;
	DiffDoub0 thk;
	DiffDoub0 Qmat[9];
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	for (i1 = 0; i1 < 81; i1++) {
		Cmat[i1].setVal(0.0);
	}

	numLay = sectPtr->getNumLayers();
	i2 = 0;
	for (i1 = 0; i1 < numLay; i1++) {
		thk.setVal(0.5);
		thk.mult(layThk[i1]);
		zMin.setVal(layZ[i1]);
		zMin.sub(thk);
		zMax.setVal(layZ[i1]);
		zMax.add(thk);
		transformQ(Qmat, &layQ[i2], layAng[i1]);
		
		// A matrix portion
		tmp.setVal(zMax);
		tmp.sub(zMin);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// B matrix portion
		tmp.setVal(zMax);
		tmp.sqr();
		tmp2.setVal(zMin);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp2.setVal(0.5);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3 + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				//i6 = i3*3 + i4;
				//i5 = i3*9 + (i4 + 3)
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// D matrix portion				
		tmp.setVal(zMax);
		tmp.sqr();
		tmp.mult(zMax);
		tmp2.setVal(zMin);
		tmp2.sqr();
		tmp2.mult(zMin);
		tmp.sub(tmp2);
		tmp2.setVal(r_1o3);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * (i3 + 3) + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}
		i2 += 9;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			//i3 = i1*9 + i2
			//i4 = i2*9 + i1
			Cmat[i3].setVal(Cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	tmp.setVal(1.0);
	tmp.mult(Cmat[20]);
	Cmat[60].setVal(tmp);
	Cmat[70].setVal(tmp);
	Cmat[80].setVal(tmp);

	return;
}

void Element::getBeamStiff(DiffDoub0 Cmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* stiffMat;
	double* Emod;
	double* shMod;
	double* areaMom;
	DiffDoub0 modulusDV[3];
	DiffDoub0 shearModDV[3];
	DiffDoub0 areaDV;
	DiffDoub0 IDV[5];
	DiffDoub0 JDV;
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DiffDoub0 xVec[6];
	DiffDoub0 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	matPt = sectPtr->getMatPtr();
	stiffMat = sectPtr->getStiffMat();
	if (stiffMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Cmat[i1].setVal(stiffMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "stiffnessMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Cmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i4].setVal(Cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		Emod = matPt->getModulus();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		JDV.setVal(sectPtr->getPolarMoment());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "modulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				modulusDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "shearModulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				shearModDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "area") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (thisDVpt->getCategory() == "areaMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				IDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "polarMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				Cmat[i3].setVal(0.0);
				i3++;
			}
		}

		Cmat[0].setVal(modulusDV[0]);
		Cmat[0].mult(areaDV);
		Cmat[4].setVal(modulusDV[0]);
		Cmat[4].mult(IDV[0]);
		Cmat[5].setVal(modulusDV[0]);
		Cmat[5].neg();
		Cmat[5].mult(IDV[1]);
		Cmat[7].setVal(shearModDV[0]);
		Cmat[7].mult(areaDV);
		Cmat[9].setVal(shearModDV[0]);
		Cmat[9].neg();
		Cmat[9].mult(IDV[0]);
		Cmat[14].setVal(shearModDV[1]);
		Cmat[14].mult(areaDV);
		Cmat[15].setVal(shearModDV[2]);
		Cmat[15].mult(IDV[1]);
		Cmat[21].setVal(shearModDV[0]);
		Cmat[21].mult(JDV);
		Cmat[28].setVal(modulusDV[0]);
		Cmat[28].mult(IDV[2]);
		Cmat[29].setVal(modulusDV[0]);
		Cmat[29].neg();
		Cmat[29].mult(IDV[4]);
		Cmat[35].setVal(modulusDV[0]);
		Cmat[35].mult(IDV[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i3].setVal(Cmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::getThermalExp(DiffDoub0 thExp[], DiffDoub0 Einit[], DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub0 dvVal;
	DiffDoub0 tmp;
	
	double* matTExp = sectPtr->getMatPtr()->getThermExp();
	for (i1 = 0; i1 < 6; i1++) {
		thExp[i1].setVal(matTExp[i1]);
		Einit[i1].setVal(0.0);
	}

	thisDVEnt = designVars->getFirst();
	thisCEnt = dvCoef->getFirst();
	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "thermalExp") {
			dvComp = thisDV->getComponent() - 1;
			tmp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(tmp);
			thExp[dvComp].add(dvVal);
		}
		else if (thisDV->getCategory() == "initialStrain") {
			dvComp = thisDV->getComponent() - 1;
			tmp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(tmp);
			Einit[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	return;
}

void Element::getShellExpLoad(DiffDoub0 expLd[], DiffDoub0 E0Ld[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layQ[], DiffDoub0 layThExp[], DiffDoub0 layEinit[], DiffDoub0 layAng[]) {
	int i1;
	int numLay;
	int layi;
	int qi;
	int exi;
	DiffDoub0 sectQ[9];
	DiffDoub0 sectTE[3];
	DiffDoub0 sectE0[3];
	DiffDoub0 QTeProd[3];
	DiffDoub0 QE0Prod[3];
	DiffDoub0 zMin;
	DiffDoub0 zMax;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	for (i1 = 0; i1 < 6; i1++) {
		expLd[i1].setVal(0.0);
		E0Ld[i1].setVal(0.0);
	}

	numLay = sectPtr->getNumLayers();
	qi = 0;
	exi = 0;
	for (layi = 0; layi < numLay; layi++) {
		transformQ(sectQ, &layQ[qi], layAng[layi]);
		transformStrain(sectTE, &layThExp[exi], layAng[layi]);
		transformStrain(sectE0, &layEinit[exi], layAng[layi]);
		matMul(QTeProd, sectQ, sectTE, 3, 3, 1);
		matMul(QE0Prod, sectQ, sectE0, 3, 3, 1);
		
		for (i1 = 0; i1 < 3; i1++) {
			tmp.setVal(QTeProd[i1]);
			tmp.mult(layThk[layi]);
			expLd[i1].add(tmp);
			tmp.setVal(QE0Prod[i1]);
			tmp.mult(layThk[layi]);
			E0Ld[i1].add(tmp);
		}

		tmp.setVal(0.5);
		tmp.mult(layThk[layi]);
		zMin.setVal(layZ[layi]);
		zMin.sub(tmp);
		zMin.sqr();
		zMax.setVal(layZ[layi]);
		zMax.add(tmp);
		zMax.sqr();
		tmp.setVal(0.5);
		tmp2.setVal(zMax);
		tmp2.sub(zMin);
		tmp.mult(tmp2); // tmp = 0.5*(zMax^2 - zMin^2)
		for (i1 = 0; i1 < 3; i1++) {
			tmp2.setVal(QTeProd[i1]);
			tmp2.mult(tmp);
			expLd[i1 + 3].add(tmp2);
			tmp2.setVal(QE0Prod[i1]);
			tmp2.mult(tmp);
			E0Ld[i1 + 3].add(tmp2);
		}
		qi += 9;
		exi += 3;
	}
	
	return;
}

void Element::getBeamExpLoad(DiffDoub0 expLd[], DiffDoub0 E0Ld[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	DiffDoub0 dvVal;
	string cat;
	string catList;
	DiffDoub0 tmp;
	Material* matPt;
	int dvComp;
	double* matMod;
	DiffDoub0 modDV[3];
	double* matG;
	DiffDoub0 shrModDV[3];
	double* matTE;
	DiffDoub0 teCoefDV[6];
	DiffDoub0 E0DV[6];
	DiffDoub0 areaDV;
	double* sectI;
	DiffDoub0 IDV[5];
	DiffDoub0 Qmat[9];
	DiffDoub0 QTE[3];
	DiffDoub0 QE0[3];
	DiffDoub0 dedgu[18];
	double* secExpLd = sectPtr->getExpLoad();
	
	if (secExpLd[0] > 0.0) {
		for (i1 = 0; i1 < 6; i1++) {
			expLd[i1].setVal(secExpLd[i1]);
			E0Ld[i1].setVal(0.0);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			if (cat == "thermalExp") {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				expLd[dvComp].add(dvVal);
			}
			else if (cat == "initialStrain") {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				E0Ld[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	}
	else {
		matPt = sectPtr->getMatPtr();
		matMod = matPt->getModulus();
		matG = matPt->getShearMod();
		matTE = matPt->getThermExp();
		for (i1 = 0; i1 < 3; i1++) {
			modDV[i1].setVal(matMod[i1]);
			shrModDV[i1].setVal(matG[i1]);
			teCoefDV[i1].setVal(matTE[i1]);
			teCoefDV[i1 + 3].setVal(matTE[i1]);
			E0DV[i1].setVal(0.0);
			E0DV[i1 + 3].setVal(0.0);
		}
		areaDV.setVal(sectPtr->getArea());
		sectI = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(sectI[i1]);
		}
		catList = "modulus shearModulus thermalExp initialStrain area areaMoment";
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			i2 = catList.find(cat);
			if (i2 > -1) {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				if (cat == "modulus") {
					modDV[dvComp].add(dvVal);
				}
				else if (cat == "shearModulus") {
					shrModDV[dvComp].add(dvVal);
				}
				else if (cat == "thermalExp") {
					teCoefDV[dvComp].add(dvVal);
				}
				else if (cat == "initialStrain") {
					E0DV[dvComp].add(dvVal);
				}
				else if (cat == "area") {
					areaDV.add(dvVal);
				}
				else if (cat == "areaMoment") {
					IDV[dvComp].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 8; i1++) {
			Qmat[i1].setVal(0.0);
		}
		Qmat[0].setVal(modDV[0]);
		Qmat[4].setVal(shrModDV[0]);
		Qmat[8].setVal(shrModDV[1]);
		teCoefDV[1].setVal(teCoefDV[3]);
		teCoefDV[2].setVal(teCoefDV[4]);
		matMul(QTE, Qmat, teCoefDV, 3, 3, 1);
		E0DV[1].setVal(E0DV[3]);
		E0DV[2].setVal(E0DV[4]);
		matMul(QE0, Qmat, E0DV, 3, 3, 1);
		for (i1 = 0; i1 < 18; i1++) {
			dedgu[i1].setVal(0.0);
		}
		dedgu[0].setVal(areaDV);
		dedgu[4].setVal(areaDV);
		dedgu[8].setVal(areaDV);
		dedgu[10].setVal(IDV[0]);
		dedgu[10].neg();
		dedgu[11].setVal(IDV[1]);
		dedgu[12].setVal(IDV[0]);
		dedgu[15].setVal(IDV[1]);
		dedgu[15].neg();

		matMul(expLd, dedgu, QTE, 6, 3, 1);
		matMul(E0Ld, dedgu, QE0, 6, 3, 1);
	}

	return;
}

void Element::getDensity(DiffDoub0& den, int layer, DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay;
	Material* thisMat;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	string cat;
	int dvLay;
	DiffDoub0 coef;
	DiffDoub0 dvVal;

	if (type == 3 || type == 41) {
		thisLay = sectPtr->getFirstLayer();
		layi = 0;
		while (layi < layer && thisLay) {
			thisLay = thisLay->getNext();
			layi++;
		}
		thisMat = thisLay->getMatPt();
		den.setVal(thisMat->getDensity());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			dvLay = thisDVPt->getLayer();
			if (cat == "density" && dvLay == layer) {
				coef.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(coef);
				den.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	} else {
		thisMat = sectPtr->getMatPtr();
		den.setVal(thisMat->getDensity());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			dvLay = thisDVPt->getLayer();
			if (cat == "density") {
				coef.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(coef);
				den.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	}

	return;
}

void Element::getShellMass(DiffDoub0 Mmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layDen[], DesignVariable* dvAr[]) {
	int i1;
	int layi;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 zMin;
	DiffDoub0 zMin2;
	DiffDoub0 zMax;
	DiffDoub0 zMax2;

	for (i1 = 0; i1 < 36; i1++) {
		Mmat[i1].setVal(0.0);
	}

	i1 = sectPtr->getNumLayers();
	for (layi = 0; layi < i1; layi++) {
		tmp.setVal(layDen[layi]);
		tmp.mult(layThk[layi]);
		Mmat[0].add(tmp);
		Mmat[7].add(tmp);
		Mmat[14].add(tmp);

		tmp.setVal(0.5);
		tmp.mult(layThk[layi]);
		zMin.setVal(layZ[layi]);
		zMin.sub(tmp);
		zMin2.setVal(zMin);
		zMin2.sqr();
		zMax.setVal(layZ[layi]);
		zMax.add(tmp);
		zMax2.setVal(zMax);
		zMax2.sqr();
		tmp.setVal(zMax2);
		tmp.sub(zMin2);  // tmp == zMax^2 - zMin^2
		tmp2.setVal(0.5);
		tmp2.mult(layDen[layi]);
		tmp2.mult(tmp); // tmp2 = 0.5*rho*(zMax^2 - zMin^2)
		Mmat[4].add(tmp2);
		Mmat[24].add(tmp2);
		Mmat[9].sub(tmp2);
		Mmat[19].sub(tmp2);

		zMax2.mult(zMax);
		zMin2.mult(zMin);
		tmp.setVal(zMax2);
		tmp.sub(zMin2); // tmp == zMax^3 - zMin^3
		tmp2.setVal(r_1o3);
		tmp2.mult(layDen[layi]);
		tmp2.mult(tmp);
		Mmat[21].add(tmp2);
		Mmat[28].add(tmp2);
		Mmat[35].add(tmp2);
	}

	return;
}

void Element::getBeamMass(DiffDoub0 Mmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dvComp;

	DesignVariable* thisDV;
	DiffDoub0 dvVal;
	DoubListEnt* coefEnt;
	DiffDoub0 coef;
	IntListEnt* dvEnt;

	string dCat;
	Material* matPt;
	double* areaMom;
	DiffDoub0 denDV;
	DiffDoub0 areaDV;
	DiffDoub0 IDV[5];
	DiffDoub0 JDV;
	DiffDoub0 tmp;

	double* massMat = sectPtr->getMassMat();
	if (massMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(massMat[i1]);
		}
		dvEnt = designVars->getFirst();
		coefEnt = dvCoef->getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value];
			if (thisDV->getCategory() == "massMat") {
				dvComp = thisDV->getComponent();
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				Mmat[dvComp].add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1; // Lower tri term
			i4 = i1; // Upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				Mmat[i3].setVal(Mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(0.0);
		}
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		JDV.setVal(sectPtr->getPolarMoment());
		matPt = sectPtr->getMatPtr();
		denDV.setVal(matPt->getDensity());
		dvEnt = designVars->getFirst();
		coefEnt = dvCoef->getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value];
			dCat = thisDV->getCategory();
			if (dCat == "density") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				denDV.add(dvVal);
			}
			else if (dCat == "area") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (dCat == "areaMoment") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				dvComp = thisDV->getComponent();
				IDV[dvComp - 1].add(dvVal);
			}
			else if (dCat == "polarMoment") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		tmp.setVal(denDV);
		tmp.mult(areaDV);
		Mmat[0].setVal(tmp);
		Mmat[7].setVal(tmp);
		Mmat[14].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[0]);
		Mmat[4].setVal(tmp);
		Mmat[9].setVal(tmp);
		Mmat[9].neg();
		tmp.setVal(denDV);
		tmp.mult(IDV[1]);
		Mmat[5].setVal(tmp);
		Mmat[5].neg();
		Mmat[15].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(JDV);
		Mmat[21].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[2]);
		Mmat[28].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[4]);
		Mmat[29].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[3]);
		Mmat[35].setVal(tmp);

		for (i1 = 3; i1 < 6; i1++) {
			i3 = 6 * i1; // Lower tri term
			i4 = i1; // Upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				Mmat[i3].setVal(Mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::getSolidDamp(DiffDoub0 Dmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub0 temp;
	DiffDoub0 dvVal;
	int dvComp;
	double* matDamp = sectPtr->getMatPtr()->getDamping();

	for (i1 = 0; i1 < 36; i1++) {
		Dmat[i1].setVal(matDamp[i1]);
	}

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "dampingMat") {
			thisDV->getValue(dvVal);
			dvComp = thisDV->getComponent();
			temp.setVal(thisCEnt->value);
			dvVal.mult(temp);
			Dmat[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 6 * i1; // Lower tri term
		i4 = i1; // Upper tri term
		for (i2 = 0; i2 < i1; i2++) {
			Dmat[i3].setVal(Dmat[i4]);
			i3++;
			i4 += 6;
		}
	}

	return;
}

void Element::getShellDamp(DiffDoub0 Dmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layD[], DiffDoub0 layAng[]) {
	getABD(Dmat, layThk, layZ, layD, layAng);
	Dmat[60].setVal(0.0);
	Dmat[70].setVal(0.0);
	Dmat[80].setVal(0.0);
	return;
}

void Element::getBeamDamp(DiffDoub0 Dmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* dampMat;
	double* areaMom;
	DiffDoub0 DmatDV[36];
	DiffDoub0 areaDV;
	DiffDoub0 IDV[5];
	DiffDoub0 JDV;
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DiffDoub0 xVec[6];
	DiffDoub0 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	matPt = sectPtr->getMatPtr();
	dampMat = sectPtr->getDampMat();
	if (dampMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Dmat[i1].setVal(dampMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "dampingMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Dmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Dmat[i4].setVal(Dmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		dampMat = matPt->getDamping();
		for (i1 = 0; i1 < 36; i1++) {
			DmatDV[i1].setVal(dampMat[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		JDV.setVal(sectPtr->getPolarMoment());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "dampingMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				DmatDV[dvComp].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "area") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (thisDVpt->getCategory() == "areaMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				IDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "polarMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				Dmat[i3].setVal(0.0);
				i3++;
			}
		}

		Dmat[0].setVal(DmatDV[0]);
		Dmat[0].mult(areaDV);
		Dmat[4].setVal(DmatDV[0]);
		Dmat[4].mult(IDV[0]);
		Dmat[5].setVal(DmatDV[0]);
		Dmat[5].neg();
		Dmat[5].mult(IDV[1]);
		Dmat[7].setVal(DmatDV[21]);
		Dmat[7].mult(areaDV);
		Dmat[9].setVal(DmatDV[21]);
		Dmat[9].neg();
		Dmat[9].mult(IDV[0]);
		Dmat[14].setVal(DmatDV[28]); // 28
		Dmat[14].mult(areaDV);
		Dmat[15].setVal(DmatDV[35]); // 35
		Dmat[15].mult(IDV[1]);
		Dmat[21].setVal(DmatDV[21]); // 21
		Dmat[21].mult(JDV);
		Dmat[28].setVal(DmatDV[0]);
		Dmat[28].mult(IDV[2]);
		Dmat[29].setVal(DmatDV[0]);
		Dmat[29].neg();
		Dmat[29].mult(IDV[4]);
		Dmat[35].setVal(DmatDV[0]);
		Dmat[35].mult(IDV[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Dmat[i3].setVal(Dmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::getConductivity(DiffDoub0 tCond[], DesignVariable* dvAr[]) {
	int i1;
	Material* matPt = sectPtr->getMatPtr();
	double* secCond = matPt->getConductivity();
	DiffDoub0 condDV[6];
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub0 temp;
	DiffDoub0 dvVal;
	int dvComp;

	for (i1 = 0; i1 < 6; i1++) {
		condDV[i1].setVal(secCond[i1]);
	}

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "thermalCond") {
			dvComp = thisDV->getComponent() - 1;
			temp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(temp);
			condDV[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	tCond[0] = condDV[0];
	tCond[1] = condDV[3];
	tCond[2] = condDV[4];
	tCond[3] = condDV[3];
	tCond[4] = condDV[1];
	tCond[5] = condDV[5];
	tCond[6] = condDV[4];
	tCond[7] = condDV[5];
	tCond[8] = condDV[2];

	return;
}

void Element::getShellCond(DiffDoub0 tCond[], DiffDoub0 layThk[], DiffDoub0 layAng[], DiffDoub0 layCond[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int numLay = sectPtr->getNumLayers();
	double* matCond;
	DiffDoub0 condDV[6];
	DiffDoub0 layerMat[9];
	DiffDoub0 alMat[9];
	DiffDoub0 alT[9];
	DiffDoub0 tmp[9];

	for (i1 = 0; i1 < 9; i1++) {
		tCond[i1].setVal(0.0);
		alMat[i1].setVal(0.0);
	}
	alMat[8].setVal(1.0);

	for (i1 = 0; i1 < numLay; i1++) {
		alMat[0].setVal(r_pio180);
		alMat[0].mult(layAng[i1]);
		alMat[0].cs();
		alMat[4].setVal(alMat[0]);
		alMat[3].setVal(r_pio180);
		alMat[3].mult(layAng[i1]);
		alMat[3].sn();
		alMat[1].setVal(alMat[3]);
		alMat[1].neg();
		transpose(alT, alMat, 3, 3);
		matMul(tmp, &layCond[9 * i1], alT, 3, 3, 3);
		matMul(layerMat, alMat, tmp, 3, 3, 3);
		for (i2 = 0; i2 < 9; i2++) {
			layerMat[i2].mult(layThk[i1]);
			tCond[i2].add(layerMat[i2]);
		}
	}

	return;
}

void Element::getBeamCond(DiffDoub0 tCond[], DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	string cat;
	DiffDoub0 condDV;
	DiffDoub0 areaDV;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;
	double secCon = sectPtr->getConductivity();
	double* matCon = sectPtr->getMatPtr()->getConductivity();
	
	if (secCon > 0.0) {
		condDV.setVal(secCon);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermCond") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
	}
	else {
		condDV.setVal(matCon[0]);
		areaDV.setVal(sectPtr->getArea());
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			cat = thisDV->getCategory();
			if (cat == "thermCond") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV.add(dvVal);
			}
			else if (cat == "area") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				areaDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		condDV.mult(areaDV);
	}

	for (i1 = 1; i1 < 8; i1++) {
		tCond[i1].setVal(0.0);
	}

	tCond[0].setVal(condDV);
	tCond[4].setVal(condDV);
	tCond[8].setVal(condDV);

	return;
}

void Element::getSpecificHeat(DiffDoub0& specHeat, DesignVariable* dvAr[]) {
	int i1;
	Material* matPt = sectPtr->getMatPtr();
	double secSpecHeat = matPt->getSpecificHeat();
	DiffDoub0 specHeatDV;
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub0 temp;
	DiffDoub0 dvVal;
	int dvComp;

	specHeat.setVal(secSpecHeat);

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "specHeat") {
			temp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(temp);
			specHeat.add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}
	
	return;
}

void Element::getShellSpecHeat(DiffDoub0& specHeat, DiffDoub0 layThk[], DiffDoub0 laySH[], DiffDoub0 layDen[]) {
	int i1;
	int numLay;
	DiffDoub0 tmp;

	specHeat.setVal(0.0);
	numLay = sectPtr->getNumLayers();
	for (i1 = 0; i1 < numLay; i1++) {
		tmp.setVal(layDen[i1]);
		tmp.mult(laySH[i1]);
		tmp.mult(layThk[i1]);
		specHeat.add(tmp);
	}

	return;
}

void Element::getBeamSpecHeat(DiffDoub0& specHeat, DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	string cat;
	DiffDoub0 densityDV;
	DiffDoub0 areaDV;
	DiffDoub0 tmp;
	DiffDoub0 dvVal;
	double secSH = sectPtr->getSpecificHeat();
	double matSH = sectPtr->getMatPtr()->getSpecificHeat();
	double matDen = sectPtr->getMatPtr()->getDensity();

	if (secSH > 0.0) {
		specHeat.setVal(secSH);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "specHeat") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				specHeat.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
	}
	else {
		specHeat.setVal(matSH);
		densityDV.setVal(matDen);
		areaDV.setVal(sectPtr->getArea());
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			cat = thisDV->getCategory();
			if (cat == "specHeatDV") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				specHeat.add(dvVal);
			}
			else if (cat == "density") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				densityDV.add(dvVal);
			}
			else if (cat == "area") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				areaDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		specHeat.mult(densityDV);
		specHeat.mult(areaDV);
	}

	return;
}

void Element::getNdCrds(DiffDoub0 xGlob[], Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	DiffDoub0 ndCrd[3];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getCrd(ndCrd,dvAr);
		xGlob[i1].setVal(ndCrd[0]);
		xGlob[i1+numNds].setVal(ndCrd[1]);
		xGlob[i1+2*numNds].setVal(ndCrd[2]);
	}
	
	return;
}

void Element::getLocOri(DiffDoub0 locOri[], DesignVariable* dvAr[]) {
	int i1;
	int dvInd;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *thisDVpt;
	string dvCat;
	DiffDoub0 rot[3];
	DiffDoub0 dvVal;
	DiffDoub0 coef;
	DiffDoub0 oriCopy[9];
	
	double *sectnOri = sectPtr->getOrientation();
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(sectnOri[i1]);
	}
	
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
    rot[2].setVal(0.0);
	thisDV = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while(thisDV) {
		dvInd = thisDV->value;
		thisDVpt = dvAr[dvInd];
		dvCat = thisDVpt->getCategory();
		if(dvCat == "orientation") {
			i1 = thisDVpt->getComponent() - 1;
			coef.setVal(thisCoef->value);
			thisDVpt->getValue(dvVal);
			dvVal.mult(coef);
			rot[i1].add(dvVal);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	
	rotateOrient(locOri, oriCopy, rot);
	
	return;
}

void Element::correctOrient(DiffDoub0 locOri[], DiffDoub0 xGlob[]) {
	int i1;
	DiffDoub0 v1[3];
	DiffDoub0 v2[3];
	DiffDoub0 v3[3];
	DiffDoub0 rot[3];
	DiffDoub0 oriCopy[9];
	
	DiffDoub0 dp;
	DiffDoub0 magv3;
	DiffDoub0 magCp;
	DiffDoub0 theta;
	DiffDoub0 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(locOri[i1]);
	}
	
	if(type == 3 || type == 41) {
		if(type == 3) {
			v1[0].setVal(xGlob[1]);
			v1[0].sub(xGlob[0]);
			v1[1].setVal(xGlob[4]);
			v1[1].sub(xGlob[3]);
			v1[2].setVal(xGlob[7]);
			v1[2].sub(xGlob[6]);
			
			v2[0].setVal(xGlob[2]);
			v2[0].sub(xGlob[0]);
			v2[1].setVal(xGlob[5]);
			v2[1].sub(xGlob[3]);
			v2[2].setVal(xGlob[8]);
			v2[2].sub(xGlob[6]);		
		} else {
			v1[0].setVal(xGlob[2]);
			v1[0].sub(xGlob[0]);
			v1[1].setVal(xGlob[6]);
			v1[1].sub(xGlob[4]);
			v1[2].setVal(xGlob[10]);
			v1[2].sub(xGlob[8]);
			
			v2[0].setVal(xGlob[3]);
			v2[0].sub(xGlob[1]);
			v2[1].setVal(xGlob[7]);
			v2[1].sub(xGlob[5]);
			v2[2].setVal(xGlob[11]);
			v2[2].sub(xGlob[9]);	
		}
		crossProd(v3, v1, v2);
		
		dp.setVal(v3[0]);
		dp.mult(locOri[6]);
		tmp.setVal(v3[1]);
		tmp.mult(locOri[7]);
		dp.add(tmp);
		tmp.setVal(v3[2]);
		tmp.mult(locOri[8]);
		dp.add(tmp);
		
		if(dp.val < 0.0) {
			v3[0].neg();
			v3[1].neg();
			v3[2].neg();
		}
		
		magv3.setVal(v3[0]);
		magv3.sqr();
		tmp.setVal(v3[1]);
		tmp.sqr();
		magv3.add(tmp);
		tmp.setVal(v3[2]);
		tmp.sqr();
		magv3.add(tmp);
		magv3.sqt();
		
		crossProd(rot,&locOri[6],v3);
		magCp.setVal(rot[0]);
		magCp.sqr();
		tmp.setVal(rot[1]);
		tmp.sqr();
		magCp.add(tmp);
		tmp.setVal(rot[2]);
		tmp.sqr();
		magCp.add(tmp);
		magCp.sqt();
		if (magCp.val < 1e-12) {
			return;
		}
		
		theta.setVal(magCp);
		theta.dvd(magv3);
		theta.asn();
		
		tmp.setVal(theta);
		tmp.dvd(magCp);
		
		rot[0].mult(tmp);
		rot[1].mult(tmp);
		rot[2].mult(tmp);
		
		rotateOrient(locOri, oriCopy, rot);
	}
	
	return;
}

void Element::getFrcFldConst(DiffDoub0 coef[], DiffDoub0 exp[], DesignVariable* dvAr[]) {
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	DoubListEnt* thisCoef;
	DiffDoub0 dvVal;
	DiffDoub0 tmp;
	string cat;

	coef[0].setVal(sectPtr->getPotCoef());
	coef[1].setVal(sectPtr->getDampCoef());
	exp[0].setVal(sectPtr->getPotExp());
	exp[1].setVal(sectPtr->getDampExp());

	thisEnt = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisEnt) {
		thisDV = dvAr[thisEnt->value];
		cat = thisDV->getCategory();
		if (cat == "potFldCoef") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			coef[0].add(dvVal);
		}
		else if (cat == "dampFldCoef") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			coef[1].add(dvVal);
		}
		thisEnt = thisEnt->next;
		thisCoef = thisCoef->next;
	}

	return;
}

void Element::getMassPerEl(DiffDoub0& massPerEl, DesignVariable* dvAr[]) {
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	DoubListEnt* thisCoef;
	DiffDoub0 dvVal;
	DiffDoub0 tmp;
	string cat;

	massPerEl.setVal(sectPtr->getMassPerEl());

	thisEnt = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisEnt) {
		thisDV = dvAr[thisEnt->value];
		cat = thisDV->getCategory();
		if (cat == "massPerEl") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			massPerEl.add(dvVal);
		}
		thisEnt = thisEnt->next;
		thisCoef = thisCoef->next;
	}

	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void Element::getLayerThkZ(DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1& zOffset, DesignVariable* dvAr[]) {
	//zOffset = 1: upper Z surface is reference plane
	//zOffset = -1: lower Z surface is reference plane
	int layi = 0;
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DesignVariable* thisDVpt;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DiffDoub1 totThk;
	DiffDoub1 tmp;
	DiffDoub1 zCrd;
	DiffDoub1 zNext;
	DiffDoub1 zMid;

	totThk.setVal(0.0);
	Layer* thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		layThk[layi].setVal(thisLay->getThickness());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVpt = dvAr[thisDV->value];
			if (thisDVpt->getCategory() == "thickness" && thisDVpt->getLayer() == layi) {
				thisDVpt->getValue(dvVal);
				coef.setVal(thisCoef->value);
				dvVal.mult(coef);
				layThk[layi].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		totThk.add(layThk[layi]);
		thisLay = thisLay->getNext();
		layi++;
	}

	zOffset.setVal(sectPtr->getZOffset());
	thisDV = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisDV) {
		thisDVpt = dvAr[thisDV->value];
		if (thisDVpt->getCategory() == "zOffset") {
			thisDVpt->getValue(dvVal);
			coef.setVal(thisCoef->value);
			dvVal.mult(coef);
			zOffset.add(dvVal);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}

	tmp.setVal(1.0);
	tmp.add(zOffset);
	zCrd.setVal(-0.5);
	zCrd.mult(totThk);
	zCrd.mult(tmp);

	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		zNext.setVal(zCrd);
		zNext.add(layThk[layi]);
		tmp.setVal(zCrd);
		tmp.add(zNext);
		zMid.setVal(0.5);
		zMid.mult(tmp);
		layZ[layi].setVal(zMid);
		thisLay = thisLay->getNext();
		layi++;
		zCrd.setVal(zNext);
	}

	return;
}

void Element::getLayerQ(DiffDoub1 layQ[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	int dvInd;
	int dvComp;
	Material* matPt;
	Layer* thisLay;
	double* Emod;
	double* pr;
	double* shMod;
	DiffDoub1 modulusDV[3];
	DiffDoub1 poissonDV[3];
	DiffDoub1 shearModDV[3];
	DiffDoub1 Smat[9];
	DiffDoub1 Qmat[9];
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DiffDoub1 xVec[3];
	DiffDoub1 bVec[3];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		matPt = thisLay->getMatPt();
		Emod = matPt->getModulus();
		pr = matPt->getPoissonRatio();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			poissonDV[i1].setVal(pr[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "modulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					modulusDV[dvComp - 1].add(dvVal);
				}
				else if (dvCat == "poissonRatio") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					poissonDV[dvComp - 1].add(dvVal);
				}
				else if (dvCat == "shearModulus") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					shearModDV[dvComp - 1].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 9; i1++) {
			Smat[i1].setVal(0.0);
		}
		Smat[0].setVal(1.0);
		Smat[0].dvd(modulusDV[0]);
		Smat[1].setVal(poissonDV[0]);
		Smat[1].neg();
		Smat[1].dvd(modulusDV[0]);
		Smat[3].setVal(Smat[1]);
		Smat[4].setVal(1.0);
		Smat[4].dvd(modulusDV[1]);
		Smat[8].setVal(1.0);
		Smat[8].dvd(shearModDV[0]);
		getDetInv(coef, Qmat, Smat, 3, 0, xVec, bVec);

		for (i1 = 0; i1 < 9; i1++) {
			layQ[i2].setVal(Qmat[i1]);
			i2++;
		}

		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerD(DiffDoub1 layD[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int layi;
	int dvInd;
	int dvComp;
	Material* matPt;
	Layer* thisLay;
	double* dampMat;
	DiffDoub1 dampMatDV[36];
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		matPt = thisLay->getMatPt();
		dampMat = matPt->getDamping();
		for (i1 = 0; i1 < 36; i1++) {
			dampMatDV[i1].setVal(dampMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "dampingMat") {
					dvComp = thisDVpt->getComponent();
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					dampMatDV[dvComp].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i3 = 1; i3 < 6; i3++) {
			i5 = 6 * i3; // Lower tri term
			i6 = i3; // Upper tri term
			for (i4 = 0; i4 < i3; i4++) {
				dampMatDV[i5].setVal(dampMatDV[i6]);
				i5++;
				i6 += 6;
			}
		}

		layD[i2].setVal(dampMatDV[0]);
		layD[i2 + 1].setVal(dampMatDV[1]);
		layD[i2 + 2].setVal(dampMatDV[3]);
		layD[i2 + 3].setVal(dampMatDV[6]);
		layD[i2 + 4].setVal(dampMatDV[7]);
		layD[i2 + 5].setVal(dampMatDV[9]);
		layD[i2 + 6].setVal(dampMatDV[18]);
		layD[i2 + 7].setVal(dampMatDV[19]);
		layD[i2 + 8].setVal(dampMatDV[21]);

		thisLay = thisLay->getNext();
		layi++;
		i2 += 9;
	}
	return;
}

void Element::getLayerAngle(DiffDoub1 layAng[], DesignVariable* dvAr[]) {
	int i2;
	int layi;
	int dvInd;
	Layer* thisLay;
	DiffDoub1 angle;
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		angle.setVal(thisLay->getAngle());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getLayer() == layi) {
				dvCat = thisDVpt->getCategory();
				if (dvCat == "angle") {
					coef.setVal(thisCoef->value);
					thisDVpt->getValue(dvVal);
					dvVal.mult(coef);
					angle.add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		layAng[layi].setVal(angle);

		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerThExp(DiffDoub1 layThExp[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	double* matExp;
	DiffDoub1 tExpDV[6];
	Material* thisMat;
	Layer* thisLay = sectPtr->getFirstLayer();
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;
	
	layi = 0;
	i2 = 0;
	while (thisLay) {
		thisMat = thisLay->getMatPt();
		matExp = thisMat->getThermExp();
		for (i1 = 0; i1 < 6; i1++) {
			tExpDV[i1].setVal(matExp[i1]);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermalExp" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				tExpDV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layThExp[i2].setVal(tExpDV[0]);
		i2++;
		layThExp[i2].setVal(tExpDV[1]);
		i2++;
		layThExp[i2].setVal(tExpDV[3]);
		i2++;
		thisLay = thisLay->getNext();
		layi++;
	}
	return;
}

void Element::getLayerEinit(DiffDoub1 layEinit[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	DiffDoub1 E0DV[6];
	Material* thisMat;
	Layer* thisLay = sectPtr->getFirstLayer();
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;

	layi = 0;
	i2 = 0;
	while (thisLay) {
		thisMat = thisLay->getMatPt();
		for (i1 = 0; i1 < 6; i1++) {
			E0DV[i1].setVal(0.0);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "initialStrain" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				E0DV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layEinit[i2].setVal(E0DV[0]);
		i2++;
		layEinit[i2].setVal(E0DV[1]);
		i2++;
		layEinit[i2].setVal(E0DV[3]);
		i2++;
		thisLay = thisLay->getNext();
		layi++;
	}
	return;
}

void Element::getLayerDen(DiffDoub1 layerDen[], DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double matDen;
	DiffDoub1 denDV;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;

	layi = 0;
	while (thisLay) {
		matDen = thisLay->getMatPt()->getDensity();
		denDV.setVal(matDen);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "density" && thisDV->getLayer() == layi) {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				denDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layerDen[layi].setVal(denDV);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getLayerCond(DiffDoub1 layCond[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double* matCond;
	DiffDoub1 condDV[6];
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;
	int dvComp;

	i1 = 0;
	layi = 0;
	while (thisLay) {
		matCond = thisLay->getMatPt()->getConductivity();
		for (i2 = 0; i2 < 6; i2++) {
			condDV[i1].setVal(matCond[i2]);
		}
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermalCond" && thisDV->getLayer() == layi) {
				dvComp = thisDV->getComponent() - 1;
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV[dvComp].add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		layCond[i1].setVal(condDV[0]);
		layCond[i1 + 1].setVal(condDV[3]);
		layCond[i1 + 2].setVal(condDV[4]);
		layCond[i1 + 3].setVal(condDV[3]);
		layCond[i1 + 4].setVal(condDV[1]);
		layCond[i1 + 5].setVal(condDV[5]);
		layCond[i1 + 6].setVal(condDV[4]);
		layCond[i1 + 7].setVal(condDV[5]);
		layCond[i1 + 8].setVal(condDV[2]);
		thisLay = thisLay->getNext();
		i1 += 9;
		layi++;
	}

	return;
}

void Element::getLayerSpecHeat(DiffDoub1 laySH[], DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay = sectPtr->getFirstLayer();
	double matSH;
	DiffDoub1 shDV;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;

	layi = 0;
	while (thisLay) {
		matSH = thisLay->getMatPt()->getSpecificHeat();
		shDV.setVal(matSH);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "specHeat" && thisDV->getLayer() == layi) {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				shDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		laySH[layi].setVal(shDV);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::transformStrain(DiffDoub1 stnNew[], DiffDoub1 stnOrig[], DiffDoub1& angle) {
	DiffDoub1 angleRad;
	DiffDoub1 a11;
	DiffDoub1 a12;
	DiffDoub1 a21;
	DiffDoub1 a22;
	DiffDoub1 tmp;
	DiffDoub1 Te[9];


	angleRad.setVal(r_pio180);
	angleRad.mult(angle);
	a11.setVal(angleRad);
	a11.cs();
	a21.setVal(angleRad);
	a21.sn();
	a12.setVal(a21);
	a12.neg();
	a22.setVal(a11);

	Te[0].setVal(a11);
	Te[0].sqr();
	Te[1].setVal(a12);
	Te[1].sqr();
	Te[2].setVal(a11);
	Te[2].mult(a12);
	Te[3].setVal(a21);
	Te[3].sqr();
	Te[4].setVal(a22);
	Te[4].sqr();
	Te[5].setVal(a22);
	Te[5].mult(a21);
	Te[6].setVal(2.0);
	Te[6].mult(a11);
	Te[6].mult(a21);
	Te[7].setVal(2.0);
	Te[7].mult(a12);
	Te[7].mult(a22);
	Te[8].setVal(a11);
	Te[8].mult(a22);
	tmp.setVal(a12);
	tmp.mult(a21);
	Te[8].add(tmp);

	matMul(stnNew, Te, stnOrig, 3, 3, 1);

	return;
}

void Element::transformQ(DiffDoub1 qNew[], DiffDoub1 qOrig[], DiffDoub1& angle) {
	int i1;
	DiffDoub1 angleRad;
	DiffDoub1 a11;
	DiffDoub1 a12;
	DiffDoub1 a21;
	DiffDoub1 a22;
	DiffDoub1 tmp;
	DiffDoub1 coef;
	DiffDoub1 Ts[9];
	DiffDoub1 Te[9];
	DiffDoub1 TeInv[9];
	DiffDoub1 xVec[3];
	DiffDoub1 bVec[3];


	angleRad.setVal(r_pio180);
	angleRad.mult(angle);
	a11.setVal(angleRad);
	a11.cs();
	a21.setVal(angleRad);
	a21.sn();
	a12.setVal(a21);
	a12.neg();
	a22.setVal(a11);

	Ts[0].setVal(a11);
	Ts[0].sqr();
	Ts[1].setVal(a12);
	Ts[1].sqr();
	Ts[2].setVal(2.0);
	Ts[2].mult(a11);
	Ts[2].mult(a12);
	Ts[3].setVal(a21);
	Ts[3].sqr();
	Ts[4].setVal(a22);
	Ts[4].sqr();
	Ts[5].setVal(2.0);
	Ts[5].mult(a22);
	Ts[5].mult(a21);
	Ts[6].setVal(a11);
	Ts[6].mult(a21);
	Ts[7].setVal(a12);
	Ts[7].mult(a22);
	Ts[8].setVal(a11);
	Ts[8].mult(a22);
	tmp.setVal(a12);
	tmp.mult(a21);
	Ts[8].add(tmp);

	for (i1 = 0; i1 < 9; i1++) {
		Te[i1].setVal(Ts[i1]);
	}
	tmp.setVal(0.5);
	Te[2].mult(tmp);
	Te[5].mult(tmp);
	tmp.setVal(2.0);
	Te[6].mult(tmp);
	Te[7].mult(tmp);

	getDetInv(coef, TeInv, Te, 3, 0, xVec, bVec);

	matMul(Te, qOrig, TeInv, 3, 3, 3);
	matMul(qNew, Ts, Te, 3, 3, 3);

	return;
}

void Element::getSolidStiff(DiffDoub1 Cmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* stiffMat;
	double* Emod;
	double* pr;
	double* shMod;
	DiffDoub1 modulusDV[3];
	DiffDoub1 poissonDV[3];
	DiffDoub1 shearModDV[3];
	DiffDoub1 Smat[36];
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DiffDoub1 xVec[6];
	DiffDoub1 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;


	matPt = sectPtr->getMatPtr();
	stiffMat = matPt->getStiffMat();
	if (stiffMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Cmat[i1].setVal(stiffMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "stiffnessMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Cmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i4].setVal(Cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		Emod = matPt->getModulus();
		pr = matPt->getPoissonRatio();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			poissonDV[i1].setVal(pr[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			dvCat = thisDVpt->getCategory();
			if (dvCat == "modulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				modulusDV[dvComp - 1].add(dvVal);
			}
			else if (dvCat == "poissonRatio") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				poissonDV[dvComp - 1].add(dvVal);
			}
			else if (dvCat == "shearModulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				shearModDV[dvComp - 1].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 0; i1 < 36; i1++) {
			Smat[i1].setVal(0.0);
		}

		Smat[0].setVal(1.0);
		Smat[0].dvd(modulusDV[0]);
		Smat[1].setVal(poissonDV[0]);
		Smat[1].neg();
		Smat[1].dvd(modulusDV[0]);
		Smat[2].setVal(poissonDV[1]);
		Smat[2].neg();
		Smat[2].dvd(modulusDV[0]);
		Smat[6].setVal(Smat[1]);
		Smat[7].setVal(1.0);
		Smat[7].dvd(modulusDV[1]);
		Smat[8].setVal(poissonDV[2]);
		Smat[8].neg();
		Smat[8].dvd(modulusDV[1]);
		Smat[12].setVal(Smat[2]);
		Smat[13].setVal(Smat[8]);
		Smat[14].setVal(1.0);
		Smat[14].dvd(modulusDV[2]);
		Smat[21].setVal(1.0);
		Smat[21].dvd(shearModDV[0]);
		Smat[28].setVal(1.0);
		Smat[28].dvd(shearModDV[1]);
		Smat[35].setVal(1.0);
		Smat[35].dvd(shearModDV[2]);

		getDetInv(coef, Cmat, Smat, 6, 0, xVec, bVec);
	}

	return;
}

void Element::getABD(DiffDoub1 Cmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layQ[], DiffDoub1 layAng[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int numLay;
	DiffDoub1 zMax;
	DiffDoub1 zMin;
	DiffDoub1 thk;
	DiffDoub1 Qmat[9];
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	for (i1 = 0; i1 < 81; i1++) {
		Cmat[i1].setVal(0.0);
	}

	numLay = sectPtr->getNumLayers();
	i2 = 0;
	for (i1 = 0; i1 < numLay; i1++) {
		thk.setVal(0.5);
		thk.mult(layThk[i1]);
		zMin.setVal(layZ[i1]);
		zMin.sub(thk);
		zMax.setVal(layZ[i1]);
		zMax.add(thk);
		transformQ(Qmat, &layQ[i2], layAng[i1]);
		
		// A matrix portion
		tmp.setVal(zMax);
		tmp.sub(zMin);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// B matrix portion
		tmp.setVal(zMax);
		tmp.sqr();
		tmp2.setVal(zMin);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp2.setVal(0.5);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3 + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				//i6 = i3*3 + i4;
				//i5 = i3*9 + (i4 + 3)
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// D matrix portion				
		tmp.setVal(zMax);
		tmp.sqr();
		tmp.mult(zMax);
		tmp2.setVal(zMin);
		tmp2.sqr();
		tmp2.mult(zMin);
		tmp.sub(tmp2);
		tmp2.setVal(r_1o3);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * (i3 + 3) + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i6]);
				Cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}
		i2 += 9;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			//i3 = i1*9 + i2
			//i4 = i2*9 + i1
			Cmat[i3].setVal(Cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	tmp.setVal(1.0);
	tmp.mult(Cmat[20]);
	Cmat[60].setVal(tmp);
	Cmat[70].setVal(tmp);
	Cmat[80].setVal(tmp);

	return;
}

void Element::getBeamStiff(DiffDoub1 Cmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* stiffMat;
	double* Emod;
	double* shMod;
	double* areaMom;
	DiffDoub1 modulusDV[3];
	DiffDoub1 shearModDV[3];
	DiffDoub1 areaDV;
	DiffDoub1 IDV[5];
	DiffDoub1 JDV;
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DiffDoub1 xVec[6];
	DiffDoub1 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	matPt = sectPtr->getMatPtr();
	stiffMat = sectPtr->getStiffMat();
	if (stiffMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Cmat[i1].setVal(stiffMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "stiffnessMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Cmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i4].setVal(Cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		Emod = matPt->getModulus();
		shMod = matPt->getShearMod();
		for (i1 = 0; i1 < 3; i1++) {
			modulusDV[i1].setVal(Emod[i1]);
			shearModDV[i1].setVal(shMod[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		JDV.setVal(sectPtr->getPolarMoment());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "modulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				modulusDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "shearModulus") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				shearModDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "area") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (thisDVpt->getCategory() == "areaMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				IDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "polarMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				Cmat[i3].setVal(0.0);
				i3++;
			}
		}

		Cmat[0].setVal(modulusDV[0]);
		Cmat[0].mult(areaDV);
		Cmat[4].setVal(modulusDV[0]);
		Cmat[4].mult(IDV[0]);
		Cmat[5].setVal(modulusDV[0]);
		Cmat[5].neg();
		Cmat[5].mult(IDV[1]);
		Cmat[7].setVal(shearModDV[0]);
		Cmat[7].mult(areaDV);
		Cmat[9].setVal(shearModDV[0]);
		Cmat[9].neg();
		Cmat[9].mult(IDV[0]);
		Cmat[14].setVal(shearModDV[1]);
		Cmat[14].mult(areaDV);
		Cmat[15].setVal(shearModDV[2]);
		Cmat[15].mult(IDV[1]);
		Cmat[21].setVal(shearModDV[0]);
		Cmat[21].mult(JDV);
		Cmat[28].setVal(modulusDV[0]);
		Cmat[28].mult(IDV[2]);
		Cmat[29].setVal(modulusDV[0]);
		Cmat[29].neg();
		Cmat[29].mult(IDV[4]);
		Cmat[35].setVal(modulusDV[0]);
		Cmat[35].mult(IDV[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Cmat[i3].setVal(Cmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::getThermalExp(DiffDoub1 thExp[], DiffDoub1 Einit[], DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	int dvComp;
	DiffDoub1 dvVal;
	DiffDoub1 tmp;
	
	double* matTExp = sectPtr->getMatPtr()->getThermExp();
	for (i1 = 0; i1 < 6; i1++) {
		thExp[i1].setVal(matTExp[i1]);
		Einit[i1].setVal(0.0);
	}

	thisDVEnt = designVars->getFirst();
	thisCEnt = dvCoef->getFirst();
	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "thermalExp") {
			dvComp = thisDV->getComponent() - 1;
			tmp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(tmp);
			thExp[dvComp].add(dvVal);
		}
		else if (thisDV->getCategory() == "initialStrain") {
			dvComp = thisDV->getComponent() - 1;
			tmp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(tmp);
			Einit[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	return;
}

void Element::getShellExpLoad(DiffDoub1 expLd[], DiffDoub1 E0Ld[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layQ[], DiffDoub1 layThExp[], DiffDoub1 layEinit[], DiffDoub1 layAng[]) {
	int i1;
	int numLay;
	int layi;
	int qi;
	int exi;
	DiffDoub1 sectQ[9];
	DiffDoub1 sectTE[3];
	DiffDoub1 sectE0[3];
	DiffDoub1 QTeProd[3];
	DiffDoub1 QE0Prod[3];
	DiffDoub1 zMin;
	DiffDoub1 zMax;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	for (i1 = 0; i1 < 6; i1++) {
		expLd[i1].setVal(0.0);
		E0Ld[i1].setVal(0.0);
	}

	numLay = sectPtr->getNumLayers();
	qi = 0;
	exi = 0;
	for (layi = 0; layi < numLay; layi++) {
		transformQ(sectQ, &layQ[qi], layAng[layi]);
		transformStrain(sectTE, &layThExp[exi], layAng[layi]);
		transformStrain(sectE0, &layEinit[exi], layAng[layi]);
		matMul(QTeProd, sectQ, sectTE, 3, 3, 1);
		matMul(QE0Prod, sectQ, sectE0, 3, 3, 1);
		
		for (i1 = 0; i1 < 3; i1++) {
			tmp.setVal(QTeProd[i1]);
			tmp.mult(layThk[layi]);
			expLd[i1].add(tmp);
			tmp.setVal(QE0Prod[i1]);
			tmp.mult(layThk[layi]);
			E0Ld[i1].add(tmp);
		}

		tmp.setVal(0.5);
		tmp.mult(layThk[layi]);
		zMin.setVal(layZ[layi]);
		zMin.sub(tmp);
		zMin.sqr();
		zMax.setVal(layZ[layi]);
		zMax.add(tmp);
		zMax.sqr();
		tmp.setVal(0.5);
		tmp2.setVal(zMax);
		tmp2.sub(zMin);
		tmp.mult(tmp2); // tmp = 0.5*(zMax^2 - zMin^2)
		for (i1 = 0; i1 < 3; i1++) {
			tmp2.setVal(QTeProd[i1]);
			tmp2.mult(tmp);
			expLd[i1 + 3].add(tmp2);
			tmp2.setVal(QE0Prod[i1]);
			tmp2.mult(tmp);
			E0Ld[i1 + 3].add(tmp2);
		}
		qi += 9;
		exi += 3;
	}
	
	return;
}

void Element::getBeamExpLoad(DiffDoub1 expLd[], DiffDoub1 E0Ld[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	DiffDoub1 dvVal;
	string cat;
	string catList;
	DiffDoub1 tmp;
	Material* matPt;
	int dvComp;
	double* matMod;
	DiffDoub1 modDV[3];
	double* matG;
	DiffDoub1 shrModDV[3];
	double* matTE;
	DiffDoub1 teCoefDV[6];
	DiffDoub1 E0DV[6];
	DiffDoub1 areaDV;
	double* sectI;
	DiffDoub1 IDV[5];
	DiffDoub1 Qmat[9];
	DiffDoub1 QTE[3];
	DiffDoub1 QE0[3];
	DiffDoub1 dedgu[18];
	double* secExpLd = sectPtr->getExpLoad();
	
	if (secExpLd[0] > 0.0) {
		for (i1 = 0; i1 < 6; i1++) {
			expLd[i1].setVal(secExpLd[i1]);
			E0Ld[i1].setVal(0.0);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			if (cat == "thermalExp") {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				expLd[dvComp].add(dvVal);
			}
			else if (cat == "initialStrain") {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				E0Ld[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	}
	else {
		matPt = sectPtr->getMatPtr();
		matMod = matPt->getModulus();
		matG = matPt->getShearMod();
		matTE = matPt->getThermExp();
		for (i1 = 0; i1 < 3; i1++) {
			modDV[i1].setVal(matMod[i1]);
			shrModDV[i1].setVal(matG[i1]);
			teCoefDV[i1].setVal(matTE[i1]);
			teCoefDV[i1 + 3].setVal(matTE[i1]);
			E0DV[i1].setVal(0.0);
			E0DV[i1 + 3].setVal(0.0);
		}
		areaDV.setVal(sectPtr->getArea());
		sectI = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(sectI[i1]);
		}
		catList = "modulus shearModulus thermalExp initialStrain area areaMoment";
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			i2 = catList.find(cat);
			if (i2 > -1) {
				dvComp = thisDVPt->getComponent() - 1;
				tmp.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(tmp);
				if (cat == "modulus") {
					modDV[dvComp].add(dvVal);
				}
				else if (cat == "shearModulus") {
					shrModDV[dvComp].add(dvVal);
				}
				else if (cat == "thermalExp") {
					teCoefDV[dvComp].add(dvVal);
				}
				else if (cat == "initialStrain") {
					E0DV[dvComp].add(dvVal);
				}
				else if (cat == "area") {
					areaDV.add(dvVal);
				}
				else if (cat == "areaMoment") {
					IDV[dvComp].add(dvVal);
				}
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 8; i1++) {
			Qmat[i1].setVal(0.0);
		}
		Qmat[0].setVal(modDV[0]);
		Qmat[4].setVal(shrModDV[0]);
		Qmat[8].setVal(shrModDV[1]);
		teCoefDV[1].setVal(teCoefDV[3]);
		teCoefDV[2].setVal(teCoefDV[4]);
		matMul(QTE, Qmat, teCoefDV, 3, 3, 1);
		E0DV[1].setVal(E0DV[3]);
		E0DV[2].setVal(E0DV[4]);
		matMul(QE0, Qmat, E0DV, 3, 3, 1);
		for (i1 = 0; i1 < 18; i1++) {
			dedgu[i1].setVal(0.0);
		}
		dedgu[0].setVal(areaDV);
		dedgu[4].setVal(areaDV);
		dedgu[8].setVal(areaDV);
		dedgu[10].setVal(IDV[0]);
		dedgu[10].neg();
		dedgu[11].setVal(IDV[1]);
		dedgu[12].setVal(IDV[0]);
		dedgu[15].setVal(IDV[1]);
		dedgu[15].neg();

		matMul(expLd, dedgu, QTE, 6, 3, 1);
		matMul(E0Ld, dedgu, QE0, 6, 3, 1);
	}

	return;
}

void Element::getDensity(DiffDoub1& den, int layer, DesignVariable* dvAr[]) {
	int layi;
	Layer* thisLay;
	Material* thisMat;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	string cat;
	int dvLay;
	DiffDoub1 coef;
	DiffDoub1 dvVal;

	if (type == 3 || type == 41) {
		thisLay = sectPtr->getFirstLayer();
		layi = 0;
		while (layi < layer && thisLay) {
			thisLay = thisLay->getNext();
			layi++;
		}
		thisMat = thisLay->getMatPt();
		den.setVal(thisMat->getDensity());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			dvLay = thisDVPt->getLayer();
			if (cat == "density" && dvLay == layer) {
				coef.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(coef);
				den.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	} else {
		thisMat = sectPtr->getMatPtr();
		den.setVal(thisMat->getDensity());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value];
			cat = thisDVPt->getCategory();
			dvLay = thisDVPt->getLayer();
			if (cat == "density") {
				coef.setVal(thisCoef->value);
				thisDVPt->getValue(dvVal);
				dvVal.mult(coef);
				den.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
	}

	return;
}

void Element::getShellMass(DiffDoub1 Mmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layDen[], DesignVariable* dvAr[]) {
	int i1;
	int layi;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 zMin;
	DiffDoub1 zMin2;
	DiffDoub1 zMax;
	DiffDoub1 zMax2;

	for (i1 = 0; i1 < 36; i1++) {
		Mmat[i1].setVal(0.0);
	}

	i1 = sectPtr->getNumLayers();
	for (layi = 0; layi < i1; layi++) {
		tmp.setVal(layDen[layi]);
		tmp.mult(layThk[layi]);
		Mmat[0].add(tmp);
		Mmat[7].add(tmp);
		Mmat[14].add(tmp);

		tmp.setVal(0.5);
		tmp.mult(layThk[layi]);
		zMin.setVal(layZ[layi]);
		zMin.sub(tmp);
		zMin2.setVal(zMin);
		zMin2.sqr();
		zMax.setVal(layZ[layi]);
		zMax.add(tmp);
		zMax2.setVal(zMax);
		zMax2.sqr();
		tmp.setVal(zMax2);
		tmp.sub(zMin2);  // tmp == zMax^2 - zMin^2
		tmp2.setVal(0.5);
		tmp2.mult(layDen[layi]);
		tmp2.mult(tmp); // tmp2 = 0.5*rho*(zMax^2 - zMin^2)
		Mmat[4].add(tmp2);
		Mmat[24].add(tmp2);
		Mmat[9].sub(tmp2);
		Mmat[19].sub(tmp2);

		zMax2.mult(zMax);
		zMin2.mult(zMin);
		tmp.setVal(zMax2);
		tmp.sub(zMin2); // tmp == zMax^3 - zMin^3
		tmp2.setVal(r_1o3);
		tmp2.mult(layDen[layi]);
		tmp2.mult(tmp);
		Mmat[21].add(tmp2);
		Mmat[28].add(tmp2);
		Mmat[35].add(tmp2);
	}

	return;
}

void Element::getBeamMass(DiffDoub1 Mmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dvComp;

	DesignVariable* thisDV;
	DiffDoub1 dvVal;
	DoubListEnt* coefEnt;
	DiffDoub1 coef;
	IntListEnt* dvEnt;

	string dCat;
	Material* matPt;
	double* areaMom;
	DiffDoub1 denDV;
	DiffDoub1 areaDV;
	DiffDoub1 IDV[5];
	DiffDoub1 JDV;
	DiffDoub1 tmp;

	double* massMat = sectPtr->getMassMat();
	if (massMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(massMat[i1]);
		}
		dvEnt = designVars->getFirst();
		coefEnt = dvCoef->getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value];
			if (thisDV->getCategory() == "massMat") {
				dvComp = thisDV->getComponent();
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				Mmat[dvComp].add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1; // Lower tri term
			i4 = i1; // Upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				Mmat[i3].setVal(Mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(0.0);
		}
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		JDV.setVal(sectPtr->getPolarMoment());
		matPt = sectPtr->getMatPtr();
		denDV.setVal(matPt->getDensity());
		dvEnt = designVars->getFirst();
		coefEnt = dvCoef->getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value];
			dCat = thisDV->getCategory();
			if (dCat == "density") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				denDV.add(dvVal);
			}
			else if (dCat == "area") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (dCat == "areaMoment") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				dvComp = thisDV->getComponent();
				IDV[dvComp - 1].add(dvVal);
			}
			else if (dCat == "polarMoment") {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		tmp.setVal(denDV);
		tmp.mult(areaDV);
		Mmat[0].setVal(tmp);
		Mmat[7].setVal(tmp);
		Mmat[14].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[0]);
		Mmat[4].setVal(tmp);
		Mmat[9].setVal(tmp);
		Mmat[9].neg();
		tmp.setVal(denDV);
		tmp.mult(IDV[1]);
		Mmat[5].setVal(tmp);
		Mmat[5].neg();
		Mmat[15].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(JDV);
		Mmat[21].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[2]);
		Mmat[28].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[4]);
		Mmat[29].setVal(tmp);
		tmp.setVal(denDV);
		tmp.mult(IDV[3]);
		Mmat[35].setVal(tmp);

		for (i1 = 3; i1 < 6; i1++) {
			i3 = 6 * i1; // Lower tri term
			i4 = i1; // Upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				Mmat[i3].setVal(Mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::getSolidDamp(DiffDoub1 Dmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub1 temp;
	DiffDoub1 dvVal;
	int dvComp;
	double* matDamp = sectPtr->getMatPtr()->getDamping();

	for (i1 = 0; i1 < 36; i1++) {
		Dmat[i1].setVal(matDamp[i1]);
	}

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "dampingMat") {
			thisDV->getValue(dvVal);
			dvComp = thisDV->getComponent();
			temp.setVal(thisCEnt->value);
			dvVal.mult(temp);
			Dmat[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 6 * i1; // Lower tri term
		i4 = i1; // Upper tri term
		for (i2 = 0; i2 < i1; i2++) {
			Dmat[i3].setVal(Dmat[i4]);
			i3++;
			i4 += 6;
		}
	}

	return;
}

void Element::getShellDamp(DiffDoub1 Dmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layD[], DiffDoub1 layAng[]) {
	getABD(Dmat, layThk, layZ, layD, layAng);
	Dmat[60].setVal(0.0);
	Dmat[70].setVal(0.0);
	Dmat[80].setVal(0.0);
	return;
}

void Element::getBeamDamp(DiffDoub1 Dmat[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dvInd;
	int dvComp;
	Material* matPt;
	double* dampMat;
	double* areaMom;
	DiffDoub1 DmatDV[36];
	DiffDoub1 areaDV;
	DiffDoub1 IDV[5];
	DiffDoub1 JDV;
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DiffDoub1 xVec[6];
	DiffDoub1 bVec[6];
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	matPt = sectPtr->getMatPtr();
	dampMat = sectPtr->getDampMat();
	if (dampMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Dmat[i1].setVal(dampMat[i1]);
		}
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "dampingMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				Dmat[dvComp].add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Dmat[i4].setVal(Dmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		dampMat = matPt->getDamping();
		for (i1 = 0; i1 < 36; i1++) {
			DmatDV[i1].setVal(dampMat[i1]);
		}
		areaDV.setVal(sectPtr->getArea());
		areaMom = sectPtr->getAreaMoment();
		for (i1 = 0; i1 < 5; i1++) {
			IDV[i1].setVal(areaMom[i1]);
		}
		JDV.setVal(sectPtr->getPolarMoment());
		thisDV = designVars->getFirst();
		thisCoef = dvCoef->getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd];
			if (thisDVpt->getCategory() == "dampingMat") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				DmatDV[dvComp].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "area") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				areaDV.add(dvVal);
			}
			else if (thisDVpt->getCategory() == "areaMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				IDV[dvComp - 1].add(dvVal);
			}
			else if (thisDVpt->getCategory() == "polarMoment") {
				dvComp = thisDVpt->getComponent();
				coef.setVal(thisCoef->value);
				thisDVpt->getValue(dvVal);
				dvVal.mult(coef);
				JDV.add(dvVal);
			}
			thisDV = thisDV->next;
			thisCoef = thisCoef->next;
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				Dmat[i3].setVal(0.0);
				i3++;
			}
		}

		Dmat[0].setVal(DmatDV[0]);
		Dmat[0].mult(areaDV);
		Dmat[4].setVal(DmatDV[0]);
		Dmat[4].mult(IDV[0]);
		Dmat[5].setVal(DmatDV[0]);
		Dmat[5].neg();
		Dmat[5].mult(IDV[1]);
		Dmat[7].setVal(DmatDV[21]);
		Dmat[7].mult(areaDV);
		Dmat[9].setVal(DmatDV[21]);
		Dmat[9].neg();
		Dmat[9].mult(IDV[0]);
		Dmat[14].setVal(DmatDV[28]); // 28
		Dmat[14].mult(areaDV);
		Dmat[15].setVal(DmatDV[35]); // 35
		Dmat[15].mult(IDV[1]);
		Dmat[21].setVal(DmatDV[21]); // 21
		Dmat[21].mult(JDV);
		Dmat[28].setVal(DmatDV[0]);
		Dmat[28].mult(IDV[2]);
		Dmat[29].setVal(DmatDV[0]);
		Dmat[29].neg();
		Dmat[29].mult(IDV[4]);
		Dmat[35].setVal(DmatDV[0]);
		Dmat[35].mult(IDV[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				Dmat[i3].setVal(Dmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::getConductivity(DiffDoub1 tCond[], DesignVariable* dvAr[]) {
	int i1;
	Material* matPt = sectPtr->getMatPtr();
	double* secCond = matPt->getConductivity();
	DiffDoub1 condDV[6];
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub1 temp;
	DiffDoub1 dvVal;
	int dvComp;

	for (i1 = 0; i1 < 6; i1++) {
		condDV[i1].setVal(secCond[i1]);
	}

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "thermalCond") {
			dvComp = thisDV->getComponent() - 1;
			temp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(temp);
			condDV[dvComp].add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}

	tCond[0] = condDV[0];
	tCond[1] = condDV[3];
	tCond[2] = condDV[4];
	tCond[3] = condDV[3];
	tCond[4] = condDV[1];
	tCond[5] = condDV[5];
	tCond[6] = condDV[4];
	tCond[7] = condDV[5];
	tCond[8] = condDV[2];

	return;
}

void Element::getShellCond(DiffDoub1 tCond[], DiffDoub1 layThk[], DiffDoub1 layAng[], DiffDoub1 layCond[], DesignVariable* dvAr[]) {
	int i1;
	int i2;
	int numLay = sectPtr->getNumLayers();
	double* matCond;
	DiffDoub1 condDV[6];
	DiffDoub1 layerMat[9];
	DiffDoub1 alMat[9];
	DiffDoub1 alT[9];
	DiffDoub1 tmp[9];

	for (i1 = 0; i1 < 9; i1++) {
		tCond[i1].setVal(0.0);
		alMat[i1].setVal(0.0);
	}
	alMat[8].setVal(1.0);

	for (i1 = 0; i1 < numLay; i1++) {
		alMat[0].setVal(r_pio180);
		alMat[0].mult(layAng[i1]);
		alMat[0].cs();
		alMat[4].setVal(alMat[0]);
		alMat[3].setVal(r_pio180);
		alMat[3].mult(layAng[i1]);
		alMat[3].sn();
		alMat[1].setVal(alMat[3]);
		alMat[1].neg();
		transpose(alT, alMat, 3, 3);
		matMul(tmp, &layCond[9 * i1], alT, 3, 3, 3);
		matMul(layerMat, alMat, tmp, 3, 3, 3);
		for (i2 = 0; i2 < 9; i2++) {
			layerMat[i2].mult(layThk[i1]);
			tCond[i2].add(layerMat[i2]);
		}
	}

	return;
}

void Element::getBeamCond(DiffDoub1 tCond[], DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	string cat;
	DiffDoub1 condDV;
	DiffDoub1 areaDV;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;
	double secCon = sectPtr->getConductivity();
	double* matCon = sectPtr->getMatPtr()->getConductivity();
	
	if (secCon > 0.0) {
		condDV.setVal(secCon);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "thermCond") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
	}
	else {
		condDV.setVal(matCon[0]);
		areaDV.setVal(sectPtr->getArea());
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			cat = thisDV->getCategory();
			if (cat == "thermCond") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				condDV.add(dvVal);
			}
			else if (cat == "area") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				areaDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		condDV.mult(areaDV);
	}

	for (i1 = 1; i1 < 8; i1++) {
		tCond[i1].setVal(0.0);
	}

	tCond[0].setVal(condDV);
	tCond[4].setVal(condDV);
	tCond[8].setVal(condDV);

	return;
}

void Element::getSpecificHeat(DiffDoub1& specHeat, DesignVariable* dvAr[]) {
	int i1;
	Material* matPt = sectPtr->getMatPtr();
	double secSpecHeat = matPt->getSpecificHeat();
	DiffDoub1 specHeatDV;
	DesignVariable* thisDV;
	IntListEnt* thisDVEnt = designVars->getFirst();
	DoubListEnt* thisCEnt = dvCoef->getFirst();
	DiffDoub1 temp;
	DiffDoub1 dvVal;
	int dvComp;

	specHeat.setVal(secSpecHeat);

	while (thisDVEnt) {
		thisDV = dvAr[thisDVEnt->value];
		if (thisDV->getCategory() == "specHeat") {
			temp.setVal(thisCEnt->value);
			thisDV->getValue(dvVal);
			dvVal.mult(temp);
			specHeat.add(dvVal);
		}
		thisDVEnt = thisDVEnt->next;
		thisCEnt = thisCEnt->next;
	}
	
	return;
}

void Element::getShellSpecHeat(DiffDoub1& specHeat, DiffDoub1 layThk[], DiffDoub1 laySH[], DiffDoub1 layDen[]) {
	int i1;
	int numLay;
	DiffDoub1 tmp;

	specHeat.setVal(0.0);
	numLay = sectPtr->getNumLayers();
	for (i1 = 0; i1 < numLay; i1++) {
		tmp.setVal(layDen[i1]);
		tmp.mult(laySH[i1]);
		tmp.mult(layThk[i1]);
		specHeat.add(tmp);
	}

	return;
}

void Element::getBeamSpecHeat(DiffDoub1& specHeat, DesignVariable* dvAr[]) {
	int i1;
	IntListEnt* thisDVEnt;
	DoubListEnt* thisCEnt;
	DesignVariable* thisDV;
	string cat;
	DiffDoub1 densityDV;
	DiffDoub1 areaDV;
	DiffDoub1 tmp;
	DiffDoub1 dvVal;
	double secSH = sectPtr->getSpecificHeat();
	double matSH = sectPtr->getMatPtr()->getSpecificHeat();
	double matDen = sectPtr->getMatPtr()->getDensity();

	if (secSH > 0.0) {
		specHeat.setVal(secSH);
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			if (thisDV->getCategory() == "specHeat") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				specHeat.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
	}
	else {
		specHeat.setVal(matSH);
		densityDV.setVal(matDen);
		areaDV.setVal(sectPtr->getArea());
		thisDVEnt = designVars->getFirst();
		thisCEnt = dvCoef->getFirst();
		while (thisDVEnt) {
			thisDV = dvAr[thisDVEnt->value];
			cat = thisDV->getCategory();
			if (cat == "specHeatDV") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				specHeat.add(dvVal);
			}
			else if (cat == "density") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				densityDV.add(dvVal);
			}
			else if (cat == "area") {
				tmp.setVal(thisCEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(tmp);
				areaDV.add(dvVal);
			}
			thisDVEnt = thisDVEnt->next;
			thisCEnt = thisCEnt->next;
		}
		specHeat.mult(densityDV);
		specHeat.mult(areaDV);
	}

	return;
}

void Element::getNdCrds(DiffDoub1 xGlob[], Node* ndAr[], DesignVariable* dvAr[]) {
	int i1;
	DiffDoub1 ndCrd[3];
	Node *nPtr;
	
	for (i1 = 0; i1 < numNds; i1++) {
		nPtr = ndAr[nodes[i1]];
		nPtr->getCrd(ndCrd,dvAr);
		xGlob[i1].setVal(ndCrd[0]);
		xGlob[i1+numNds].setVal(ndCrd[1]);
		xGlob[i1+2*numNds].setVal(ndCrd[2]);
	}
	
	return;
}

void Element::getLocOri(DiffDoub1 locOri[], DesignVariable* dvAr[]) {
	int i1;
	int dvInd;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *thisDVpt;
	string dvCat;
	DiffDoub1 rot[3];
	DiffDoub1 dvVal;
	DiffDoub1 coef;
	DiffDoub1 oriCopy[9];
	
	double *sectnOri = sectPtr->getOrientation();
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(sectnOri[i1]);
	}
	
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
    rot[2].setVal(0.0);
	thisDV = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while(thisDV) {
		dvInd = thisDV->value;
		thisDVpt = dvAr[dvInd];
		dvCat = thisDVpt->getCategory();
		if(dvCat == "orientation") {
			i1 = thisDVpt->getComponent() - 1;
			coef.setVal(thisCoef->value);
			thisDVpt->getValue(dvVal);
			dvVal.mult(coef);
			rot[i1].add(dvVal);
		}
		thisDV = thisDV->next;
		thisCoef = thisCoef->next;
	}
	
	rotateOrient(locOri, oriCopy, rot);
	
	return;
}

void Element::correctOrient(DiffDoub1 locOri[], DiffDoub1 xGlob[]) {
	int i1;
	DiffDoub1 v1[3];
	DiffDoub1 v2[3];
	DiffDoub1 v3[3];
	DiffDoub1 rot[3];
	DiffDoub1 oriCopy[9];
	
	DiffDoub1 dp;
	DiffDoub1 magv3;
	DiffDoub1 magCp;
	DiffDoub1 theta;
	DiffDoub1 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(locOri[i1]);
	}
	
	if(type == 3 || type == 41) {
		if(type == 3) {
			v1[0].setVal(xGlob[1]);
			v1[0].sub(xGlob[0]);
			v1[1].setVal(xGlob[4]);
			v1[1].sub(xGlob[3]);
			v1[2].setVal(xGlob[7]);
			v1[2].sub(xGlob[6]);
			
			v2[0].setVal(xGlob[2]);
			v2[0].sub(xGlob[0]);
			v2[1].setVal(xGlob[5]);
			v2[1].sub(xGlob[3]);
			v2[2].setVal(xGlob[8]);
			v2[2].sub(xGlob[6]);		
		} else {
			v1[0].setVal(xGlob[2]);
			v1[0].sub(xGlob[0]);
			v1[1].setVal(xGlob[6]);
			v1[1].sub(xGlob[4]);
			v1[2].setVal(xGlob[10]);
			v1[2].sub(xGlob[8]);
			
			v2[0].setVal(xGlob[3]);
			v2[0].sub(xGlob[1]);
			v2[1].setVal(xGlob[7]);
			v2[1].sub(xGlob[5]);
			v2[2].setVal(xGlob[11]);
			v2[2].sub(xGlob[9]);	
		}
		crossProd(v3, v1, v2);
		
		dp.setVal(v3[0]);
		dp.mult(locOri[6]);
		tmp.setVal(v3[1]);
		tmp.mult(locOri[7]);
		dp.add(tmp);
		tmp.setVal(v3[2]);
		tmp.mult(locOri[8]);
		dp.add(tmp);
		
		if(dp.val < 0.0) {
			v3[0].neg();
			v3[1].neg();
			v3[2].neg();
		}
		
		magv3.setVal(v3[0]);
		magv3.sqr();
		tmp.setVal(v3[1]);
		tmp.sqr();
		magv3.add(tmp);
		tmp.setVal(v3[2]);
		tmp.sqr();
		magv3.add(tmp);
		magv3.sqt();
		
		crossProd(rot,&locOri[6],v3);
		magCp.setVal(rot[0]);
		magCp.sqr();
		tmp.setVal(rot[1]);
		tmp.sqr();
		magCp.add(tmp);
		tmp.setVal(rot[2]);
		tmp.sqr();
		magCp.add(tmp);
		magCp.sqt();
		if (magCp.val < 1e-12) {
			return;
		}
		
		theta.setVal(magCp);
		theta.dvd(magv3);
		theta.asn();
		
		tmp.setVal(theta);
		tmp.dvd(magCp);
		
		rot[0].mult(tmp);
		rot[1].mult(tmp);
		rot[2].mult(tmp);
		
		rotateOrient(locOri, oriCopy, rot);
	}
	
	return;
}

void Element::getFrcFldConst(DiffDoub1 coef[], DiffDoub1 exp[], DesignVariable* dvAr[]) {
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	DoubListEnt* thisCoef;
	DiffDoub1 dvVal;
	DiffDoub1 tmp;
	string cat;

	coef[0].setVal(sectPtr->getPotCoef());
	coef[1].setVal(sectPtr->getDampCoef());
	exp[0].setVal(sectPtr->getPotExp());
	exp[1].setVal(sectPtr->getDampExp());

	thisEnt = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisEnt) {
		thisDV = dvAr[thisEnt->value];
		cat = thisDV->getCategory();
		if (cat == "potFldCoef") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			coef[0].add(dvVal);
		}
		else if (cat == "dampFldCoef") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			coef[1].add(dvVal);
		}
		thisEnt = thisEnt->next;
		thisCoef = thisCoef->next;
	}

	return;
}

void Element::getMassPerEl(DiffDoub1& massPerEl, DesignVariable* dvAr[]) {
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	DoubListEnt* thisCoef;
	DiffDoub1 dvVal;
	DiffDoub1 tmp;
	string cat;

	massPerEl.setVal(sectPtr->getMassPerEl());

	thisEnt = designVars->getFirst();
	thisCoef = dvCoef->getFirst();
	while (thisEnt) {
		thisDV = dvAr[thisEnt->value];
		cat = thisDV->getCategory();
		if (cat == "massPerEl") {
			thisDV->getValue(dvVal);
			tmp.setVal(thisCoef->value);
			dvVal.mult(tmp);
			massPerEl.add(dvVal);
		}
		thisEnt = thisEnt->next;
		thisCoef = thisCoef->next;
	}

	return;
}

//end dup
 
//end skip 
 
 
 
 
 
 
 
