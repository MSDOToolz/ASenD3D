#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "ModelClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "JobClass.h"

using namespace std;

void Model::writeTimeStepSoln(int tStep) {
	int i1;
	int dofPerNd;
	int numIntDof;
	double ndDat[6];
	double elDat[9];
	char* buf;
	string fullFile = solveCmd->fileName + "/solnTStep" + to_string(tStep) + ".out";
	ofstream outFile;
	outFile.open(fullFile, std::ofstream::binary);

	buf = reinterpret_cast<char*>(&ndDat[0]);
	Node* thisNd = nodes->getFirst();
	while (thisNd) {
		if (solveCmd->thermal) {
			ndDat[0] = thisNd->getPrevTemp();
			outFile.write(buf, 8);
			ndDat[0] = thisNd->getPrevTdot();
			outFile.write(buf, 8);
		}
		if (solveCmd->elastic) {
			dofPerNd = thisNd->getNumDof();
			i1 = 8 * dofPerNd;
			thisNd->getPrevDisp(ndDat);
			outFile.write(buf, i1);
			thisNd->getPrevVel(ndDat);
			outFile.write(buf, i1);
			thisNd->getPrevAcc(ndDat);
			outFile.write(buf, i1);
		}
		thisNd = thisNd->getNext();
	}

	if (solveCmd->elastic) {
		buf = reinterpret_cast<char*>(&elDat[0]);
		Element* thisEl = elements->getFirst();
		while (thisEl) {
			numIntDof = thisEl->getNumIntDof();
			if (numIntDof > 0) {
				i1 = 8 * numIntDof;
				thisEl->getIntPrevDisp(elDat);
				outFile.write(buf, i1);
			}
			thisEl = thisEl->getNext();
		}
	}

	outFile.close();
	return;
}

void Model::writeNodeResults(string fileName, string nodeSet, StringList& fields, int timeStep) {
	int i1;
	int globInd;
	Set* setPt;
	string errStr;
	string thisField;
	IntListEnt *thisEnt;
	Node *ndPt;
	int ndLabel;
	double ndDat[6];
	int ndof;
	double time;
	
	if(timeStep >= 0) {
		// Read the results from the time step file and store them in nodes
		try {
			readTimeStepSoln(timeStep);
			ndPt = nodes->getFirst();
			while (ndPt) {
				if (solveCmd->thermal) {
					ndPt->backstepTemp();
				}
				if (solveCmd->elastic) {
					ndPt->backstepDisp();
				}
				ndPt = ndPt->getNext();
			}
		}
		catch (...) {
			cout << "Error: time step " << timeStep << " requested for node results is out of range, or solution history was not saved" << endl;
			cout << "Aborting writeNodeResults" << endl;
			return;
		}
	}
	
	try {
		i1 = nsMap.at(nodeSet);
		setPt = nsArray[i1];
	}
	catch (...) {
		errStr = "Warning: there is no node set named " + nodeSet + ". Defaulting to all nodes in writeNodeResults";
		cout << errStr << endl;
		i1 = nsMap.at("all");
		setPt = nsArray[i1];
	}
	
	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);
	
	outFile << "nodeResults:\n";
	time = solveCmd->timeStep * timeStep;
	outFile << "    time: " << time << "\n";
	outFile << "    nodeSet: " << nodeSet << "\n";
	
	StringListEnt *strPt = fields.getFirst();
	while(strPt) {
		thisField = strPt->value;
		outFile << "    " << thisField << ":\n";
		// -----------------
		// Calculate reaction force if necessary
		if (thisField == "reactionForce") {
			buildElasticSolnLoad(tempV1, false);
			thisEnt = setPt->getFirstEntry();
			while (thisEnt) {
				ndLabel = thisEnt->value;
				ndPt = nodeArray[ndLabel];
				outFile << "        - [" << ndLabel;
				ndof = ndPt->getNumDof();
				for (i1 = 0; i1 < ndof; i1++) {
					globInd = ndPt->getDofIndex(i1);
					outFile << ", " << -tempV1[globInd];
				}
				outFile << "]\n";
				thisEnt = thisEnt->next;
			}
		}
		else if (thisField == "reactionHeatGen") {
			buildThermalSolnLoad(tempV1, false);
			thisEnt = setPt->getFirstEntry();
			while (thisEnt) {
				ndLabel = thisEnt->value;
				ndPt = nodeArray[ndLabel];
				outFile << "        - [" << ndLabel;
				globInd = ndPt->getSortedRank();
				outFile << ", " << -tempV1[globInd] << "\n";
				thisEnt = thisEnt->next;
			}
		}
		else {
			thisEnt = setPt->getFirstEntry();
			while (thisEnt) {
				ndLabel = thisEnt->value;
				ndPt = nodeArray[ndLabel];
				outFile << "        - [" << ndLabel;
				if (thisField == "displacement") {
					ndPt->getDisp(ndDat);
					ndof = ndPt->getNumDof();
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndDat[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "velocity") {
					ndPt->getVel(ndDat);
					ndof = ndPt->getNumDof();
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndDat[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "acceleration") {
					ndPt->getAcc(ndDat);
					ndof = ndPt->getNumDof();
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndDat[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "temperature") {
					ndDat[0] = ndPt->getTemperature();
					outFile << ", " << ndDat[0] << "]\n";
				}
				else if (thisField == "tdot") {
					ndDat[0] = ndPt->getTdot();
					outFile << ", " << ndDat[0] << "]\n";
				}
				thisEnt = thisEnt->next;
			}
		}
		strPt = strPt->next;
	}
	
	outFile.close();
	
	return;
}

void Model::writeElementResults(string fileName, string elSet, StringList& fields, string position, int timeStep) {
	int i1;
	int i2;
	int i3;
	string errStr;
	string thisField;
	IntListEnt* thisEnt;
	Element* elPt;
	Set* setPt;
	int type;
	int elLabel;
	int numIP;
	int numLay;
	double* intPts;
	double cent[3] = { 0.0,0.0,0.0 };
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	double SEDen;
	DiffDoub0 def[9];
	DiffDoub0 frcMom[9];
	DiffDoub0 tGrad[3];
	DiffDoub0 flux[3];
	string fieldList;
	double time;

	if (timeStep >= 0) {
		// Read the results from the time step file and store them in nodes
		try {
			readTimeStepSoln(timeStep);
		}
		catch (...) {
			cout << "Error: time step " << timeStep << " requested for element results is out of range," << endl;
			cout << "or solution history was not saved in the solve options. Aborting writeElementResults." << endl;
			return;
		}
	}

	try {
		i1 = esMap.at(elSet);
		setPt = esArray[i1];
	}
	catch (...) {
		errStr = "Warning: there is no element set named " + elSet + ". Defaulting to all elements in writeNodeResults";
		cout << errStr << endl;
		i1 = nsMap.at("all");
		setPt = esArray[i1];
	}

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "elementResults:\n";
	time = solveCmd->timeStep * timeStep;
	outFile << "    time: " << time << "\n";
	outFile << "    elSet: " << elSet << "\n";

	StringListEnt* strPt = fields.getFirst();
	while (strPt) {
		thisField = strPt->value;
		outFile << "    " << thisField << ":\n";
		fieldList = "stress strain";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ##  - [element label, integration pt, layer, S11, S22, S33, S12, S13, S23]\n";
		}
		fieldList = "strainEnergyDen";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ##  - [element label, integration pt, layer, strain energy]\n";
		}
		fieldList = "sectionDef sectionFrcMom";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ## for shells:  - [element label, integration pt, S11, S22, S12, K11, K22, K12]\n";
			outFile << "    ## for beams:  - [element label, integration pt, S11, S12, S13, K11, K12, K13]\n";
		}
		if (thisField == "heatFlux") {
			outFile << "    ##  - [element label, integration pt, layer, f1, f2, f3]\n";
		}
		if (thisField == "tempGradient") {
			outFile << "    ##  - [element label, integration pt, layer, dT/dx, dT/dy, dT/dz]\n";
		}

		thisEnt = setPt->getFirstEntry();
		while (thisEnt) {
			elLabel = thisEnt->value;
			elPt = elementArray[elLabel];
			type = elPt->getType();
			
			fieldList = "stress strain strainEnergyDen";
			i2 = fieldList.find(thisField);
			if (i2 > -1) {
				elPt->getStressPrereq(*d0Pre, nodeArray, dVarArray);
				if (position == "intPts") {
					numIP = elPt->getNumIP();
					intPts = elPt->getIP();
				}
				else {
					numIP = 1;
					intPts = elPt->getSCent();
				}
				for (i1 = 0; i1 < numIP; i1++) {
					if (type == 3 || type == 41) {
						numLay = elPt->getNumLayers();
					} else {
						numLay = 1;
					}
					for (i2 = 0; i2 < numLay; i2++) {
						elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, *d0Pre);
						outFile << "        - [" << elLabel << ", ";
						if (thisField == "strain") {
							outFile << i1 << ", " << i2 << ", " << strain[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								outFile << ", " << strain[i3].val;
							}
							outFile << "]\n";
						} else if(thisField == "stress") {
							outFile << i1 << ", " << i2 << ", " << stress[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								outFile << ", " << stress[i3].val;
							}
							outFile << "]\n";
						} else {
							SEDen = 0.0;
							for (i3 = 0; i3 < 6; i3++) {
								SEDen += stress[i3].val * strain[i3].val;
							}
							SEDen *= 0.5;
							outFile << i1 << ", " << i2 << ", " << SEDen << "]\n";
						}
					}
				}
			}

			fieldList = "sectionDef sectionFrcMom";
			i2 = fieldList.find(thisField);
			if (i2 > -1 && elPt->getDofPerNd() == 6) {
				elPt->getStressPrereq(*d0Pre, nodeArray, dVarArray);
				if (position == "intPts") {
					numIP = elPt->getNumIP();
					intPts = elPt->getIP();
				}
				else {
					numIP = 1;
					intPts = elPt->getSCent();
				}
				for (i1 = 0; i1 < numIP; i1++) {
					type = elPt->getType();
					//elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, stPre);
					elPt->getDefFrcMom(def, frcMom, &intPts[3 * i1], solveCmd->nonlinearGeom, *d0Pre);
					outFile << "        - [" << elLabel << ", ";
					if (thisField == "sectionDef") {
						outFile << i1 << ", " << def[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							outFile << ", " << def[i3].val;
						}
						outFile << "]\n";
					}
					else if (thisField == "sectionFrcMom") {
						outFile << i1 << ", " << frcMom[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							outFile << ", " << frcMom[i3].val;
						}
						outFile << "]\n";
					}
				}
			}

			fieldList = "tempGradient heatFlux";
			i2 = fieldList.find(thisField);
			if (i2 > -1) {
				elPt->getStressPrereq(*d0Pre, nodeArray, dVarArray);
				if (position == "intPts") {
					numIP = elPt->getNumIP();
					intPts = elPt->getIP();
				}
				else {
					numIP = 1;
					intPts = elPt->getSCent();
				}
				for (i1 = 0; i1 < numIP; i1++) {
					type = elPt->getType();
					if (type == 3 || type == 41) {
						numLay = elPt->getNumLayers();
					}
					else {
						numLay = 1;
					}
					for (i2 = 0; i2 < numLay; i2++) {
						//elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, stPre);
						elPt->getFluxTGrad(flux, tGrad, &intPts[3 * i1], i2, *d0Pre);
						outFile << "        - [" << elLabel << ", " << i1 << ", " << i2 << ", ";
						if (thisField == "tempGradient") {
							outFile << tGrad[0].val << ", " << tGrad[1].val << ", " << tGrad[2].val << "]\n";
						}
						else if (thisField == "heatFlux") {
							outFile << flux[0].val << ", " << flux[1].val << ", " << flux[2].val << "]\n";
						}
					}
				}
			}

			thisEnt = thisEnt->next;
		}
		strPt = strPt->next;
	}

	outFile.close();

	return;
}

void Model::writeModalResults(string fileName, bool writeModes) {
	int i1;
	int i2;
	int i3;
	int nd;
	int dofPerNd;
	int globInd;
	int nMds = modalCmd->numModes;
	Node* thisNd;

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "modalResults:\n";
	outFile << "    eigenValues:\n";
	for (i1 = 0; i1 < nMds; i1++) {
		outFile << "      - " << eigVals[i1] << "\n";
	}
	if (modalCmd->type == "buckling") {
		outFile << "    loadFactors:\n";
	}
	else {
		outFile << "    frequencies:\n";
	}
	for (i1 = 0; i1 < nMds; i1++) {
		outFile << "      - " << loadFact[i1] << "\n";
	}

	if (modalCmd->writeModes) {
		outFile << "    modes:\n";
		for (i1 = 0; i1 < nMds; i1++) {
			outFile << "      - mode: " << i1 << "\n";
			outFile << "        displacement:\n";
			thisNd = nodes->getFirst();
			while (thisNd) {
				nd = thisNd->getLabel();
				outFile << "          - [" << nd;
				dofPerNd = thisNd->getNumDof();
				for (i2 = 0; i2 < dofPerNd; i2++) {
					globInd = thisNd->getDofIndex(i2);
					i3 = i1 * elMatDim + globInd;
					outFile << ", " << eigVecs[i3];
				}
				outFile << "]\n";
				thisNd = thisNd->getNext();
			}
		}
	}

	outFile.close();

	return;
}

void Model::writeObjective(string fileName, StringList& includeFields, bool writeGrad) {
	int i1;
	int numDV;
	double totObj;
	double thisVal;
	ObjectiveTerm* thisTerm;
	StringListEnt* thisFld;
	string fldStr;
	double* actTime;

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "objective:\n";
	outFile << "    terms:\n";
	thisTerm = objective->getFirst();
	totObj = 0.0;
	while (thisTerm) {
		thisVal = thisTerm->getValue();
		outFile << "        - value: " << thisVal << "\n";
		thisFld = includeFields.getFirst();
		while (thisFld) {
			fldStr = thisFld->value;
			if (fldStr == "category") {
				outFile << "          category: " << thisTerm->getCategory() << "\n";
			}
			else if (fldStr == "operator") {
				outFile << "          operator: " << thisTerm->getOperator() << "\n";
			}
			else if (fldStr == "component") {
				outFile << "          component: " << thisTerm->getComponent() << "\n";
			}
			else if (fldStr == "layer") {
				outFile << "          layer: " << thisTerm->getLayer() << "\n";
			}
			else if (fldStr == "coefficient") {
				outFile << "          coefficient: " << thisTerm->getCoef() << "\n";
			}
			else if (fldStr == "exponent") {
				outFile << "          exponent: " << thisTerm->getExpnt() << "\n";
			}
			else if (fldStr == "elementSet") {
				outFile << "          elementSet: " << thisTerm->getElsetName() << "\n";
			}
			else if (fldStr == "nodeSet") {
				outFile << "          nodeSet: " << thisTerm->getNdsetName() << "\n";
			}
			else if (fldStr == "activeTime") {
				actTime = thisTerm->getActTime();
				outFile << "          activeTime: [" << actTime[0] << ", " << actTime[1] << "]\n";
			}
			thisFld = thisFld->next;
		}
		totObj += thisVal;
		thisTerm = thisTerm->getNext();
	}
	outFile << "    totalValue: " << totObj << "\n";
	if (writeGrad) {
		outFile << "objectiveGradient:\n";
		numDV = designVars->getLength();
		for (i1 = 0; i1 < numDV; i1++) {
			outFile << "    - [" << i1 << ", " << dLdD[i1] << "]\n";
		}
	}

	outFile.close();

	return;
}