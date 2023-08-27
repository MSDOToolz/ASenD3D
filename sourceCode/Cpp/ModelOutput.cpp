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
	string fullFile = solveCmd->fileName + "solnTStep" + to_string(tStep) + ".out";
	ofstream outFile;
	outFile.open(fullFile);
	outFile << setprecision(15);

	Node* thisNd = nodes.getFirst();
	while (thisNd) {
		if (solveCmd->thermal) {

		}
		if (solveCmd->elastic) {
			dofPerNd = thisNd->getNumDof();
			thisNd->getPrevDisp(ndDat);
			for (i1 = 0; i1 < dofPerNd; i1++) {
				outFile << ndDat[i1] << "\n";
			}
			thisNd->getPrevVel(ndDat);
			for (i1 = 0; i1 < dofPerNd; i1++) {
				outFile << ndDat[i1] << "\n";
			}
			thisNd->getPrevAcc(ndDat);
			for (i1 = 0; i1 < dofPerNd; i1++) {
				outFile << ndDat[i1] << "\n";
			}
		}
		thisNd = thisNd->getNext();
	}

	Element* thisEl = elements.getFirst();
	while (thisEl) {
		numIntDof = thisEl->getNumIntDof();
		if (numIntDof < 0) {
			thisEl->getIntPrevDisp(elDat);
			for (i1 = 0; i1 < numIntDof; i1++) {
				outFile << elDat[i1] << "\n";
			}
		}
		thisEl = thisEl->getNext();
	}

	outFile.close();
	return;
}

void Model::writeNodeResults(string fileName, string nodeSet, StringList& fields, int timeStep) {
	int i1;
	Set* setPt;
	string errStr;
	string thisField;
	IntListEnt *thisEnt;
	Node *ndPt;
	int ndLabel;
	double ndDat[6];
	int ndof;
	
	if(timeStep >= 0) {
		// Read the results from the time step file and store them in nodes
		try {
			readTimeStepSoln(timeStep);
		}
		catch (...) {
			cout << "Error: time step " << timeStep << " requested for node results is out of range, or solution history was not saved" << endl;
			cout << "Aborting writeNodeResults" << endl;
			return;
		}
	}
	
	try {
		i1 = nsMap.at(nodeSet);
		setPt = nsArray[i1].ptr;
	}
	catch (...) {
		errStr = "Warning: there is no node set named " + nodeSet + ". Defaulting to all nodes in writeNodeResults";
		cout << errStr << endl;
		i1 = nsMap.at("all");
		setPt = nsArray[i1].ptr;
	}
	
	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);
	
	outFile << "nodeResults:\n";
	outFile << "    nodeSet: " << nodeSet << "\n";
	
	StringListEnt *strPt = fields.getFirst();
	while(strPt) {
		thisField = strPt->value;
		outFile << "    " << thisField << ":\n";
		// -----------------
		// Calculate reaction force if necessary
        // ------------------		
		thisEnt = setPt->getFirstEntry();
		while(thisEnt) {
			ndLabel = thisEnt->value;
			ndPt = nodeArray[ndLabel].ptr;
			outFile << "        - [" << ndLabel << ", ";
			if(thisField == "displacement") {
				ndPt->getDisp(ndDat);
				outFile << ndDat[0];
				ndof = ndPt->getNumDof();
				for (i1 = 1; i1 < ndof; i1++) {
					outFile << ", " << ndDat[i1];
				}
				outFile << "]\n";
			} else if(thisField == "velocity") {
				ndPt->getVel(ndDat);
				outFile << ndDat[0];
				ndof = ndPt->getNumDof();
				for (i1 = 1; i1 < ndof; i1++) {
					outFile << ", " << ndDat[i1];
				}
				outFile << "]\n";
			} else if(thisField == "acceleration") {
				ndPt->getAcc(ndDat);
				outFile << ndDat[0];
				ndof = ndPt->getNumDof();
				for (i1 = 1; i1 < ndof; i1++) {
					outFile << ", " << ndDat[i1];
				}
				outFile << "]\n";
			} else if(thisField == "temperature") {
				ndDat[0] = ndPt->getTemperature();
				outFile << ndDat[0] << "]\n";
			} else if(thisField == "tdot") {
				ndDat[0] = ndPt->getTdot();
				outFile << ndDat[0] << "]\n";
			}
			thisEnt = thisEnt->next;
		}
		strPt = strPt->next;
	}
	
	outFile.close();
	
	return;
}

void Model::writeElementResults(string fileName, string elSet, StringList& fields, int timeStep) {
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
	Doub strain[6];
	Doub stress[6];
	double SEDen;
	string fieldList;
	DoubStressPrereq stPre;

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
		setPt = esArray[i1].ptr;
	}
	catch (...) {
		errStr = "Warning: there is no element set named " + elSet + ". Defaulting to all elements in writeNodeResults";
		cout << errStr << endl;
		i1 = nsMap.at("all");
		setPt = esArray[i1].ptr;
	}

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "elementResults:\n";
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
		thisEnt = setPt->getFirstEntry();
		while (thisEnt) {
			elLabel = thisEnt->value;
			elPt = elementArray[elLabel].ptr;
			fieldList = "stress strain strainEnergyDen";
			i2 = fieldList.find(thisField);
			if (i2 > -1) {
				elPt->getStressPrereq(stPre, !solveCmd->nonlinearGeom, nodeArray, dVarArray);
				numIP = elPt->getNumIP();
				intPts = elPt->getIP();
				for (i1 = 0; i1 < numIP; i1++) {
					type = elPt->getType();
					if (type == 3 || type == 41) {
						numLay = elPt->getNumLayers();
					} else {
						numLay = 1;
					}
					for (i2 = 0; i2 < numLay; i2++) {
						elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, stPre);
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
			thisEnt = thisEnt->next;
		}
		strPt = strPt->next;
	}

	outFile.close();

	if (stPre.layerZ) {
		stPre.destroy();
	}

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
	thisTerm = objective.getFirst();
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
		numDV = designVars.getLength();
		for (i1 = 0; i1 < numDV; i1++) {
			outFile << "    - [" << i1 << ", " << dLdD[i1] << "]\n";
		}
	}

	outFile.close();

	return;
}