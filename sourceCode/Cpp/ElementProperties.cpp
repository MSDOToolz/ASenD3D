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

void Element::getLayerThkZ(Doub layThk[], Doub layZ[], Doub& zOffset, DVPt dvAr[]) {
	int layi = 0;
	Doub dvVal;
	Doub coef;
	DesignVariable* thisDVpt;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	Doub totThk;
	Doub tmp;
	Doub zCrd;
	Doub zNext;
	Doub zMid;

	totThk.setVal(0.0);
	Layer* thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		layThk[layi].setVal(thisLay->getThickness());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVpt = dvAr[thisDV->value].ptr;
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
	thisDV = designVars.getFirst();
	thisCoef = dvCoef.getFirst();
	while (thisDV) {
		thisDVpt = dvAr[thisDV->value].ptr;
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
	}

	return;
}

void Element::getLayerQ(Doub layQ[], DVPt dvAr[]) {
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
	Doub modulusDV[3];
	Doub poissonDV[3];
	Doub shearModDV[3];
	Doub Smat[9];
	Doub Qmat[9];
	Doub dvVal;
	Doub coef;
	Doub xVec[3];
	Doub bVec[3];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getLayerAngle(Doub layAng[], DVPt dvAr[]) {
	int i2;
	int layi;
	int dvInd;
	Layer* thisLay;
	Doub angle;
	Doub dvVal;
	Doub coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		angle.setVal(thisLay->getAngle());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::transformStrain(Doub stnNew[], Doub stnOrig[], Doub& angle) {
	Doub angleRad;
	Doub a11;
	Doub a12;
	Doub a21;
	Doub a22;
	Doub tmp;
	Doub Te[9];


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

void Element::transformQ(Doub qNew[], Doub qOrig[], Doub& angle) {
	int i1;
	Doub angleRad;
	Doub a11;
	Doub a12;
	Doub a21;
	Doub a22;
	Doub tmp;
	Doub coef;
	Doub Ts[9];
	Doub Te[9];
	Doub TeInv[9];
	Doub xVec[3];
	Doub bVec[3];


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

void Element::getSolidStiff(Doub Cmat[], DVPt dvAr[]) {
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
	Doub modulusDV[3];
	Doub poissonDV[3];
	Doub shearModDV[3];
	Doub Smat[36];
	Doub dvVal;
	Doub coef;
	Doub xVec[6];
	Doub bVec[6];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getABD(Doub Cmat[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layAng[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int numLay;
	Doub zMax;
	Doub zMin;
	Doub thk;
	Doub Qmat[9];
	Doub tmp;
	Doub tmp2;

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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * i1;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * i1 + 3;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * (i1 + 3) + 3;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
			}
		}
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			Cmat[i3].setVal(Cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	Cmat[60].setVal(Cmat[20]);
	Cmat[70].setVal(Cmat[20]);
	Cmat[80].setVal(Cmat[20]);

	return;
}

void Element::getBeamStiff(Doub Cmat[], DVPt dvAr[]) {
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
	Doub modulusDV[3];
	Doub shearModDV[3];
	Doub areaDV;
	Doub IDV[5];
	Doub JDV;
	Doub dvVal;
	Doub coef;
	Doub xVec[6];
	Doub bVec[6];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getDensity(Doub& den, int layer, DVPt dvAr[]) {
	int layi;
	Layer* thisLay;
	Material* thisMat;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	string cat;
	int dvLay;
	Doub coef;
	Doub dvVal;

	if (type == 3 || type == 4) {
		thisLay = sectPtr->getFirstLayer();
		layi = 0;
		while (layi < layer && thisLay) {
			thisLay = thisLay->getNext();
			layi++;
		}
		thisMat = thisLay->getMatPt();
		den.setVal(thisMat->getDensity());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value].ptr;
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

void Element::getShellMass(Doub Mmat[], Doub layThk[], Doub layZ[], DVPt dvAr[]) {
	int i1;
	int layi;
	Doub layDen;
	DesignVariable* thisDV;
	Doub dvVal;
	DoubListEnt* coefEnt;
	Doub coef;
	IntListEnt* dvEnt;

	for (i1 = 0; i1 < 36; i1++) {
		Mmat[i1].setVal(0.0);
	}

	Layer* thisLay = sectPtr->getFirstLayer();
	layi = 0;
	while (thisLay) {
		layDen.setVal(thisLay->getMatPt()->getDensity());
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
			if (thisDV->getCategory() == "density" && thisDV->getLayer() == layi) {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				layDen.add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		layDen.mult(layThk[layi]);
		Mmat[0].add(layDen);
		Mmat[7].add(layDen);
		Mmat[14].add(layDen);
		layDen.mult(layZ[layi]);
		Mmat[4].add(layDen);
		Mmat[24].add(layDen);
		Mmat[9].sub(layDen);
		Mmat[19].sub(layDen);
		layDen.mult(layZ[layi]);
		Mmat[21].add(layDen);
		Mmat[28].add(layDen);
		Mmat[35].add(layDen);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getBeamMass(Doub Mmat[], DVPt dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dvComp;

	DesignVariable* thisDV;
	Doub dvVal;
	DoubListEnt* coefEnt;
	Doub coef;
	IntListEnt* dvEnt;

	string dCat;
	Material* matPt;
	double* areaMom;
	Doub denDV;
	Doub areaDV;
	Doub IDV[5];
	Doub JDV;
	Doub tmp;

	double* massMat = sectPtr->getMassMat();
	if (massMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(massMat[i1]);
		}
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
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
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
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

void Element::getNdCrds(Doub xGlob[], NdPt ndAr[], DVPt dvAr[]) {
	int i1;
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

void Element::getLocOri(Doub locOri[], DVPt dvAr[]) {
	int i1;
	int dvInd;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *thisDVpt;
	string dvCat;
	Doub rot[3];
	Doub dvVal;
	Doub coef;
	Doub oriCopy[9];
	
	double *sectnOri = sectPtr->getOrientation();
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(sectnOri[i1]);
	}
	
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
    rot[2].setVal(0.0);
	thisDV = designVars.getFirst();
	thisCoef = dvCoef.getFirst();
	while(thisDV) {
		dvInd = thisDV->value;
		thisDVpt = dvAr[dvInd].ptr;
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

void Element::correctOrient(Doub locOri[], Doub xGlob[]) {
	int i1;
	Doub v1[3];
	Doub v2[3];
	Doub v3[3];
	Doub rot[3];
	Doub oriCopy[9];
	
	Doub dp;
	Doub magv3;
	Doub magCp;
	Doub theta;
	Doub tmp;
	
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

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

void Element::getLayerThkZ(DiffDoub layThk[], DiffDoub layZ[], DiffDoub& zOffset, DVPt dvAr[]) {
	int layi = 0;
	DiffDoub dvVal;
	DiffDoub coef;
	DesignVariable* thisDVpt;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DiffDoub totThk;
	DiffDoub tmp;
	DiffDoub zCrd;
	DiffDoub zNext;
	DiffDoub zMid;

	totThk.setVal(0.0);
	Layer* thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		layThk[layi].setVal(thisLay->getThickness());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVpt = dvAr[thisDV->value].ptr;
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
	thisDV = designVars.getFirst();
	thisCoef = dvCoef.getFirst();
	while (thisDV) {
		thisDVpt = dvAr[thisDV->value].ptr;
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
	}

	return;
}

void Element::getLayerQ(DiffDoub layQ[], DVPt dvAr[]) {
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
	DiffDoub modulusDV[3];
	DiffDoub poissonDV[3];
	DiffDoub shearModDV[3];
	DiffDoub Smat[9];
	DiffDoub Qmat[9];
	DiffDoub dvVal;
	DiffDoub coef;
	DiffDoub xVec[3];
	DiffDoub bVec[3];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getLayerAngle(DiffDoub layAng[], DVPt dvAr[]) {
	int i2;
	int layi;
	int dvInd;
	Layer* thisLay;
	DiffDoub angle;
	DiffDoub dvVal;
	DiffDoub coef;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVpt;
	string dvCat;

	i2 = 0;
	layi = 0;
	thisLay = sectPtr->getFirstLayer();
	while (thisLay) {
		angle.setVal(thisLay->getAngle());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::transformStrain(DiffDoub stnNew[], DiffDoub stnOrig[], DiffDoub& angle) {
	DiffDoub angleRad;
	DiffDoub a11;
	DiffDoub a12;
	DiffDoub a21;
	DiffDoub a22;
	DiffDoub tmp;
	DiffDoub Te[9];


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

void Element::transformQ(DiffDoub qNew[], DiffDoub qOrig[], DiffDoub& angle) {
	int i1;
	DiffDoub angleRad;
	DiffDoub a11;
	DiffDoub a12;
	DiffDoub a21;
	DiffDoub a22;
	DiffDoub tmp;
	DiffDoub coef;
	DiffDoub Ts[9];
	DiffDoub Te[9];
	DiffDoub TeInv[9];
	DiffDoub xVec[3];
	DiffDoub bVec[3];


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

void Element::getSolidStiff(DiffDoub Cmat[], DVPt dvAr[]) {
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
	DiffDoub modulusDV[3];
	DiffDoub poissonDV[3];
	DiffDoub shearModDV[3];
	DiffDoub Smat[36];
	DiffDoub dvVal;
	DiffDoub coef;
	DiffDoub xVec[6];
	DiffDoub bVec[6];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getABD(DiffDoub Cmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layAng[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int numLay;
	DiffDoub zMax;
	DiffDoub zMin;
	DiffDoub thk;
	DiffDoub Qmat[9];
	DiffDoub tmp;
	DiffDoub tmp2;

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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * i1;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * i1 + 3;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
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
		for (i1 = 0; i1 < 3; i1++) {
			i3 = 9 * (i1 + 3) + 3;
			i4 = 3 * i1;
			for (i2 = 0; i2 < 3; i2++) {
				tmp2.setVal(tmp);
				tmp2.mult(Qmat[i4]);
				Cmat[i3].add(tmp2);
				i3++;
				i4++;
			}
		}
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			Cmat[i3].setVal(Cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	Cmat[60].setVal(Cmat[20]);
	Cmat[70].setVal(Cmat[20]);
	Cmat[80].setVal(Cmat[20]);

	return;
}

void Element::getBeamStiff(DiffDoub Cmat[], DVPt dvAr[]) {
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
	DiffDoub modulusDV[3];
	DiffDoub shearModDV[3];
	DiffDoub areaDV;
	DiffDoub IDV[5];
	DiffDoub JDV;
	DiffDoub dvVal;
	DiffDoub coef;
	DiffDoub xVec[6];
	DiffDoub bVec[6];
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			dvInd = thisDV->value;
			thisDVpt = dvAr[dvInd].ptr;
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

void Element::getDensity(DiffDoub& den, int layer, DVPt dvAr[]) {
	int layi;
	Layer* thisLay;
	Material* thisMat;
	IntListEnt* thisDV;
	DoubListEnt* thisCoef;
	DesignVariable* thisDVPt;
	string cat;
	int dvLay;
	DiffDoub coef;
	DiffDoub dvVal;

	if (type == 3 || type == 4) {
		thisLay = sectPtr->getFirstLayer();
		layi = 0;
		while (layi < layer && thisLay) {
			thisLay = thisLay->getNext();
			layi++;
		}
		thisMat = thisLay->getMatPt();
		den.setVal(thisMat->getDensity());
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value].ptr;
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
		thisDV = designVars.getFirst();
		thisCoef = dvCoef.getFirst();
		while (thisDV) {
			thisDVPt = dvAr[thisDV->value].ptr;
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

void Element::getShellMass(DiffDoub Mmat[], DiffDoub layThk[], DiffDoub layZ[], DVPt dvAr[]) {
	int i1;
	int layi;
	DiffDoub layDen;
	DesignVariable* thisDV;
	DiffDoub dvVal;
	DoubListEnt* coefEnt;
	DiffDoub coef;
	IntListEnt* dvEnt;

	for (i1 = 0; i1 < 36; i1++) {
		Mmat[i1].setVal(0.0);
	}

	Layer* thisLay = sectPtr->getFirstLayer();
	layi = 0;
	while (thisLay) {
		layDen.setVal(thisLay->getMatPt()->getDensity());
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
			if (thisDV->getCategory() == "density" && thisDV->getLayer() == layi) {
				coef.setVal(coefEnt->value);
				thisDV->getValue(dvVal);
				dvVal.mult(coef);
				layDen.add(dvVal);
			}
			dvEnt = dvEnt->next;
			coefEnt = coefEnt->next;
		}
		layDen.mult(layThk[layi]);
		Mmat[0].add(layDen);
		Mmat[7].add(layDen);
		Mmat[14].add(layDen);
		layDen.mult(layZ[layi]);
		Mmat[4].add(layDen);
		Mmat[24].add(layDen);
		Mmat[9].sub(layDen);
		Mmat[19].sub(layDen);
		layDen.mult(layZ[layi]);
		Mmat[21].add(layDen);
		Mmat[28].add(layDen);
		Mmat[35].add(layDen);
		thisLay = thisLay->getNext();
		layi++;
	}

	return;
}

void Element::getBeamMass(DiffDoub Mmat[], DVPt dvAr[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dvComp;

	DesignVariable* thisDV;
	DiffDoub dvVal;
	DoubListEnt* coefEnt;
	DiffDoub coef;
	IntListEnt* dvEnt;

	string dCat;
	Material* matPt;
	double* areaMom;
	DiffDoub denDV;
	DiffDoub areaDV;
	DiffDoub IDV[5];
	DiffDoub JDV;
	DiffDoub tmp;

	double* massMat = sectPtr->getMassMat();
	if (massMat[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			Mmat[i1].setVal(massMat[i1]);
		}
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
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
		dvEnt = designVars.getFirst();
		coefEnt = dvCoef.getFirst();
		while (dvEnt) {
			thisDV = dvAr[dvEnt->value].ptr;
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

void Element::getNdCrds(DiffDoub xGlob[], NdPt ndAr[], DVPt dvAr[]) {
	int i1;
	DiffDoub ndCrd[3];
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

void Element::getLocOri(DiffDoub locOri[], DVPt dvAr[]) {
	int i1;
	int dvInd;
	IntListEnt *thisDV;
	DoubListEnt *thisCoef;
	DesignVariable *thisDVpt;
	string dvCat;
	DiffDoub rot[3];
	DiffDoub dvVal;
	DiffDoub coef;
	DiffDoub oriCopy[9];
	
	double *sectnOri = sectPtr->getOrientation();
	for (i1 = 0; i1 < 9; i1++) {
		oriCopy[i1].setVal(sectnOri[i1]);
	}
	
	rot[0].setVal(0.0);
	rot[1].setVal(0.0);
    rot[2].setVal(0.0);
	thisDV = designVars.getFirst();
	thisCoef = dvCoef.getFirst();
	while(thisDV) {
		dvInd = thisDV->value;
		thisDVpt = dvAr[dvInd].ptr;
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

void Element::correctOrient(DiffDoub locOri[], DiffDoub xGlob[]) {
	int i1;
	DiffDoub v1[3];
	DiffDoub v2[3];
	DiffDoub v3[3];
	DiffDoub rot[3];
	DiffDoub oriCopy[9];
	
	DiffDoub dp;
	DiffDoub magv3;
	DiffDoub magCp;
	DiffDoub theta;
	DiffDoub tmp;
	
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

//end dup
 
//end skip 
 
