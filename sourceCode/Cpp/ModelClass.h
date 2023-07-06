#ifndef MODEL
#define MODEL
#include <fstream>
#include <iostream>
#include "JobClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "FaceClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "DiffDoubClass.h"

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
		LoadList loads;
		DesVarList designVars;
		DVPt *dVarArray;
		Objective objective;
		Job job;
		
		int elMatDim;
		SparseMat elasticMat;
	
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
		
		void executeJob();

        void destroy();	// Finish	
		
		// Input
		
		void readInputLine(std::ifstream& inFile,std::string& fileLine,std::string headings[],std::string data[], int& dataLen);
		
		void readJob(std::string fileName);
		
		void readModelInput(std::string fileName);
		
		void readConstraintInput(std::string fileName);
		
		void readLoadInput(std::string fileName);
		
		void readInitialState(std::string fileName);
		
		void readDesVarInput(std::string fileName);
		
		void readObjectiveInput(std::string fileName);
		
		void readDesVarValues(std::string fileName);
		
		// Analysis
		
		void reorderNodes(int blockDim);
		
		void updateReference();
		
		void findSurfaceFaces();
		
		void analysisPrep(int blockDim);
		
		void buildElasticAppLoad(double appLd[], double time);
		
		void buildElasticSolnLoad();
		
		void solveStep();
		
		void solve();
		
		void eigenSolve();
		
		void getObjGradient();
		
		// Output
		
		//
};

#endif