#ifndef MODEL
#define MODEL
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

class Model {
	private:
	    NodeList nodes;
		NdPt *nodeArray;
		ElementList elements;
		ElPt *elementArray;
		SetList nodeSets;
		SetList elementSets;
		SectionList sections;
		MaterialList materials;
		ConstraintList constraints;
		DesVarList designVars;
		DVPt *dVarArray;
		Objective objective;
	
	public:
	    Model();
		
		NodeList* getNodes();
		
		ElementList* getElements();
		
		SetList* getNodeSets();
		
		SetList* getElementSets();
		
		SectionList* getSections();
		
		MaterialList* getMaterials();
		
		ConstraintList* getConstraints();
		
		DesVarList* getDesignVars();
};

#endif