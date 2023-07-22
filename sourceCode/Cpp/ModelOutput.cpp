#include <fstream>
#include <iostream>
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

using namespace std;

void Model::writeNodeResults(string fileName, string nodeSet, StringList& fields, int timeStep) {
	int i1;
	string errStr;
	string thisField;
	IntListEnt *thisEnt;
	Node *ndPt;
	int ndLabel;
	double ndDat[6];
	int ndof;
	
	if(timestep >= 0) {
		// Read the results from the time step file and store them in nodes
	}
	
	Set *setPt = nodeSets.getFirst();
	bool found = false;
	while(!found && setPt) {
		if(setPt->getName() == nodeSet) {
			found = true;
		} else {
			setPt = setPt->getNext();
		}
	}
	
	if(!found) {
		errStr = "Warning: there is no node set named " + nodeSet + ", aborting writeNodeResults";
		cout << errStr << endl;
		return;
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