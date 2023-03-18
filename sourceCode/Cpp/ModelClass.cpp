#include "ModelClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"


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