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
	elMatDim = 0;
	totGlobDof = 0;
	elasticScaled = false;
	thermScaled = false;
	solveCmd = -1;
	modalCmd = -1;

	return;
}

void Model::executeJob() {
	int i1;
	int i2;
	int numTsteps;
	double time;
	string cmdStr;
	string fileName;
	string exten;
	string fullFname;

	i1 = 0;
	for (auto& thisCmd : job) {
		cmdStr = thisCmd.cmdString;
		fileName = thisCmd.fileName;
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
			readNodeResults(thisCmd.fileName);
		}
		else if (cmdStr == "solvePrep") {
			solveCmd = i1;
			analysisPrep();
		}
		else if (cmdStr == "solve") {
			cout << "running main analysis " << endl;
			solveCmd = i1;
			solve();
		}
		else if (cmdStr == "zeroSolution") {
			zeroSolution(thisCmd.fields);
		}
		else if (cmdStr == "modalAnalysis") {
			cout << "running modal analysis " << fileName << endl;
			modalCmd = i1;
			eigenSolve();
		}
		else if (cmdStr == "setSolnToMode") {
			setSolnToMode(thisCmd.solnField, thisCmd.mode, thisCmd.maxAmplitude);
		}
		else if (cmdStr == "calcObjective") {
			cout << "calculating objective function " << fileName << endl;
			getObjective();
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
			numTsteps = thisCmd.timeSteps.size();
			if (numTsteps > 0) {
				i2 = thisCmd.fileName.find(".");
				if (i2 > -1) {
					exten = thisCmd.fileName.substr(i2);
					fileName = thisCmd.fileName.substr(0, i2);
				}
				else {
					exten = "";
					fileName = thisCmd.fileName;
				}
				for (auto& ts : thisCmd.timeSteps) {
					fullFname = fileName + "_timestep" + to_string(ts) + exten;
					if (cmdStr == "writeNodeResults") {
						writeNodeResults(fullFname, thisCmd.nodeSet, thisCmd.fields, ts);
					}
					else {
						writeElementResults(fullFname, thisCmd.elementSet, thisCmd.fields, thisCmd.position, ts);
					}
				}
			}
			else {
				if (cmdStr == "writeNodeResults") {
					writeNodeResults(thisCmd.fileName, thisCmd.nodeSet, thisCmd.fields, -1);
				}
				else {
					writeElementResults(thisCmd.fileName, thisCmd.elementSet, thisCmd.fields, thisCmd.position, -1);
				}
			}
		}
		else if (cmdStr == "writeModalResults") {
			cout << "writing modal results" << endl;
			writeModalResults(thisCmd.fileName, thisCmd.writeModes);
		}
		else if (cmdStr == "writeObjective") {
			cout << "writing objective results" << endl;
			writeObjective(thisCmd.fileName, thisCmd.objInclude, thisCmd.writeGradient);
		}
		
		i1++;
	}
	
	return;
}