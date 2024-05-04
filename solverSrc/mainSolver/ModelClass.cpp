#include <fstream>
#include <iostream>
#include <string>
#include "ModelClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

using namespace std;

Model::Model() {
	nodeArray = nullptr;
	elementArray = nullptr;
	nsArray = nullptr;
	esArray = nullptr;
	dVarArray = nullptr;
	anPrepRun = false;
	timeStepsSaved = 0;
	elMatDim = 0;
	totGlobDof = 0;
	elasticScaled = false;
	thermScaled = false;
	solveCmd = nullptr;
	modalCmd = nullptr;

	eigVecs = nullptr;
	eigVals = nullptr;
	diagMass = nullptr;
	loadFact = nullptr;

	tempV1 = nullptr;
	tempV2 = nullptr;
	tempV3 = nullptr;
	tempV4 = nullptr;
	tempV5 = nullptr;
	tempV6 = nullptr;

	tempD1 = nullptr;

	dLdU = nullptr;
	dLdV = nullptr;
	dLdA = nullptr;
	dLdT = nullptr;
	dLdTdot = nullptr;
	uAdj = nullptr;
	vAdj = nullptr;
	aAdj = nullptr;
	tAdj = nullptr;
	tdotAdj = nullptr;

	dRudD = nullptr;
	dRtdD = nullptr;
	dLdD = nullptr;
	elInD = nullptr;
	return;
}

NodeList* Model::getNodes() {
	return &nodes;
}

ElementList* Model::getElements() {
	return &elements;
}

SetList* Model::getNodeSets() {
	return &nodeSets;
}

SetList* Model::getElementSets() {
	return &elementSets;
}

SectionList* Model::getSections() {
	return &sections;
}

MaterialList* Model::getMaterials() {
	return &materials;
}

ConstraintList* Model::getElasticConstraints() {
	return &elasticConst;
}

DesVarList* Model::getDesignVars() {
	return &designVars;
}

void Model::executeJob() {
	int i1;
	int i2;
	int numTsteps;
	double time;
	IntListEnt* thisEnt;
	JobCommand* thisCmd = job.getFirst();
	Node* thisNd;
	Element* thisEl;
	DoubListEnt* thisLdTm;
	string cmdStr;
	string fileName;
	string exten;
	string fullFname;

	while(thisCmd) {
		cmdStr = thisCmd->cmdString;
		fileName = thisCmd->fileName;
		if(cmdStr == "readModelInput") {
			cout << "reading model input: " << fileName << endl;
			readModelInput(fileName);
			readConstraintInput(fileName);
			readLoadInput(fileName);
			readInitialState(fileName);
		} else if(cmdStr == "readConstraints") {
			cout << "reading constraints: " << fileName << endl;
			readConstraintInput(fileName);
		} else if(cmdStr == "readLoads") {
			cout << "reading loads: " << fileName << endl;
			readLoadInput(fileName);
		} else if(cmdStr == "readInitialState") {
			cout << "reading initial state: " << fileName << endl;
			readInitialState(fileName);
		} else if(cmdStr == "readDesignVarInput") {
			cout << "reading design variable input: " << fileName << endl;
			readDesVarInput(fileName);
		} else if(cmdStr == "readObjectiveInput") {
			cout << "reading objective input: " << fileName << endl;
			readObjectiveInput(fileName);
		}
		else if (cmdStr == "readDesignVarValues") {
			readDesVarValues(fileName);
		}
		else if (cmdStr == "readNodeResults") {
			readNodeResults(thisCmd->fileName);
		}
		else if (cmdStr == "solvePrep") {
			solveCmd = thisCmd;
			analysisPrep();
		}
		else if (cmdStr == "solve") {
			cout << "running main analysis " << endl;
			solveCmd = thisCmd;
			solve(thisCmd);
		}
		else if (cmdStr == "zeroSolution") {
			zeroSolution(thisCmd->fields);
		}
		else if (cmdStr == "modalAnalysis") {
			cout << "running modal analysis " << fileName << endl;
			modalCmd = thisCmd;
			eigenSolve(thisCmd);
		}
		else if (cmdStr == "setSolnToMode") {
			setSolnToMode(thisCmd->solnField, thisCmd->mode, thisCmd->maxAmplitude);
		}
		else if (cmdStr == "calcObjective") {
			cout << "calculating objective function " << fileName << endl;
			if (solveCmd->dynamic) {
				if (timeStepsSaved == 0) {
					cout << "Warning: dynamic analysis was run but time history was not saved." << endl;
					cout << "         The objective can only be computed based on the final state." << endl;
					objective.clearValues();
					objective.calculateTerms(solveCmd->simPeriod, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray,d0Pre);
				}
				else {
					objective.clearValues();
					readTimeStepSoln(timeStepsSaved);
					i1 = timeStepsSaved - 1;
					time = solveCmd->timeStep * timeStepsSaved;
					while (i1 >= 0) {
						thisNd = nodes.getFirst();
						while (thisNd) {
							if (solveCmd->elastic) {
								thisNd->backstepDisp();
							}
							if (solveCmd->thermal) {
								thisNd->backstepTemp();
							}
							thisNd = thisNd->getNext();
						}
						if (solveCmd->elastic) {
							thisEl = elements.getFirst();
							while (thisEl) {
								thisEl->backstepIntDisp();
								thisEl = thisEl->getNext();
							}
						}
						readTimeStepSoln(i1);
						objective.calculateTerms(time, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray, d0Pre);
						i1--;
						time -= solveCmd->timeStep;
					}
				}
			}
			else {
				objective.clearValues();
				thisLdTm = solveCmd->staticLoadTime.getFirst();
				i1 = 0;
				while (thisLdTm) {
					readTimeStepSoln(i1);
					thisNd = nodes.getFirst();
					while (thisNd) {
						if (solveCmd->elastic) {
							thisNd->backstepDisp();
						}
						if (solveCmd->thermal) {
							thisNd->backstepTemp();
						}
						thisNd = thisNd->getNext();
					}
					if (solveCmd->elastic) {
						thisEl = elements.getFirst();
						while (thisEl) {
							thisEl->backstepIntDisp();
							thisEl = thisEl->getNext();
						}
					}
					objective.calculateTerms(thisLdTm->value, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray, d0Pre);
					thisLdTm = thisLdTm->next;
					i1++;
				}
			}
		}
		else if (cmdStr == "calcObjGradient") {
			cout << "calculating objective gradient " << fileName << endl;
			getObjGradient();
		} else if (cmdStr == "writeNodeResults" || cmdStr == "writeElementResults") {
			if (cmdStr == "writeNodeResults") {
				cout << "writing node results" << endl;
			}
			else {
				cout << "writing element results" << endl;
			}
			numTsteps = thisCmd->timeSteps.getLength();
			if (numTsteps > 0) {
				i2 = thisCmd->fileName.find(".");
				if (i2 > -1) {
					exten = thisCmd->fileName.substr(i2);
					fileName = thisCmd->fileName.substr(0, i2);
				}
				else {
					exten = "";
					fileName = thisCmd->fileName;
				}
				thisEnt = thisCmd->timeSteps.getFirst();
				while (thisEnt) {
					fullFname = fileName + "_timestep" + to_string(thisEnt->value) + exten;
					if (cmdStr == "writeNodeResults") {
						writeNodeResults(fullFname, thisCmd->nodeSet, thisCmd->fields, thisEnt->value);
					}
					else {
						writeElementResults(fullFname, thisCmd->elementSet, thisCmd->fields, thisCmd->position, thisEnt->value);
					}
					thisEnt = thisEnt->next;
				}
			}
			else {
				if (cmdStr == "writeNodeResults") {
					writeNodeResults(thisCmd->fileName, thisCmd->nodeSet, thisCmd->fields, -1);
				}
				else {
					writeElementResults(thisCmd->fileName, thisCmd->elementSet, thisCmd->fields, thisCmd->position, -1);
				}
			}
		}
		else if (cmdStr == "writeModalResults") {
			cout << "writing modal results" << endl;
			writeModalResults(thisCmd->fileName, thisCmd->writeModes);
		}
		else if (cmdStr == "writeObjective") {
			cout << "writing objective results" << endl;
			writeObjective(thisCmd->fileName, thisCmd->objInclude, thisCmd->writeGradient);
		}
		
		thisCmd = thisCmd->next;
	}
	
	return;
}

Model::~Model() {
	delete[] nodeArray;
	delete[] elementArray;
	delete[] nsArray;
	delete[] esArray;
	delete[] dVarArray;

	if (eigVecs) {
		delete[] eigVecs;
		delete[] eigVals;
		delete[] loadFact;
		delete[] diagMass;
	}

	if (tempV1) {
		delete[] tempV1;
		delete[] tempV2;
		delete[] tempV3;
		delete[] tempV4;
		delete[] tempV5;
		delete[] tempV6;
		delete[] tempD1;
		delete[] dLdU;
		delete[] dLdV;
		delete[] dLdA;
		delete[] dLdT;
		delete[] dLdTdot;
		delete[] uAdj;
		delete[] vAdj;
		delete[] aAdj;
		delete[] tAdj;
		delete[] tdotAdj;
		delete[] dRudD;
		delete[] dRtdD;
		delete[] elInD;
	}

	if (dLdD) {
		delete[] dLdD;
	}

	return;
}