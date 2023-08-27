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
	anPrepRun = false;
	timeStepsSaved = 0;
	elasticScaled = false;
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

ConstraintList* Model::getConstraints() {
	return &constraints;
}

DesVarList* Model::getDesignVars() {
	return &designVars;
}

void Model::executeJob() {
	int i2;
	int numTsteps;
	IntListEnt* thisEnt;
	JobCommand *thisCmd = job.getFirst();
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
			readConstraintInput(fileName);
		} else if(cmdStr == "readLoads") {
			readLoadInput(fileName);
		} else if(cmdStr == "readInitialState") {
			readInitialState(fileName);
		} else if(cmdStr == "readDesignVarInput") {
			cout << "reading design variable input: " << fileName << endl;
			readDesVarInput(fileName);
		} else if(cmdStr == "readObjectiveInput") {
			cout << "reading objective input: " << fileName << endl;
			readObjectiveInput(fileName);
		} else if (cmdStr == "solve") {
			solve(thisCmd);
			solveCmd = thisCmd;
		} else if (cmdStr == "calcObjective") {
			if (solveCmd->dynamic) {
				if (timeStepsSaved == 0) {
					cout << "Warning: dynamic analysis was run but time history was not saved." << endl;
					cout << "         The objective can only be computed based on the final state." << endl;
					objective.clearValues();
					objective.calculateTerms(thisCmd->simPeriod, thisCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
				}
				else {
					// rem: finish for time integration
				}
			}
			else {
				objective.clearValues();
				objective.calculateTerms(thisCmd->staticLoadTime, thisCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
			}
		}
		else if (cmdStr == "calcObjGradient") {
			getObjGradient();
		} else if (cmdStr == "writeNodeResults" || cmdStr == "writeElementResults") {
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
						writeElementResults(fullFname, thisCmd->elementSet, thisCmd->fields, thisEnt->value);
					}
					thisEnt = thisEnt->next;
				}
			}
			else {
				if (cmdStr == "writeNodeResults") {
					writeNodeResults(thisCmd->fileName, thisCmd->nodeSet, thisCmd->fields, -1);
				}
				else {
					writeElementResults(thisCmd->fileName, thisCmd->elementSet, thisCmd->fields, -1);
				}
			}
		} else if (cmdStr == "writeObjective") {
			writeObjective(thisCmd->fileName, thisCmd->objInclude, thisCmd->writeGradient);
		}
		
		thisCmd = thisCmd->next;
	}
	
	return;
}

void Model::destroy() {
	nodes.destroy();
	delete[] nodeArray;
	elements.destroy();
	delete[] elementArray;
	nodeSets.destroy();
	delete[] nsArray;
	elementSets.destroy();
	delete[] esArray;
	sections.destroy();
	materials.destroy();
	constraints.destroy();
	loads.destroy();
	designVars.destroy();
	delete[] dVarArray;
	objective.destroy();
	job.destroy();

	elasticMat.destroy();
	elasticLT.destroy();

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
	}

	if (dLdD) {
		delete[] dLdD;
	}

	return;
}