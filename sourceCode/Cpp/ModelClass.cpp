#include <fstream>
#include <iostream>
#include <string>
#include "ModelClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

using namespace std;

Model::Model() {
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
	JobCommand *thisCmd = job.getFirst();
	string cmdStr;
	string fileName;
	while(thisCmd) {
		cmdStr = thisCmd->cmdString;
		fileName = thisCmd->fileName;
		if(cmdStr == "readModelInput") {
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
			readDesVarInput(fileName);
		} else if(cmdStr == "readObjectiveInput") {
			readObjectiveInput(fileName);
		}
		
		thisCmd = thisCmd->next;
	}
	
	return;
}

void Model::destroy() {
	return;
}