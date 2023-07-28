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
#include "LoadClass.h"
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
		bool anPrepRun;
		int timeStepsSaved;
		bool nonlinearGeom;
		SparseMat elasticMat;
		LowerTriMat elasticLT;
		
		double *tempV1;
		double *tempV2;
		Doub *tempD1;
	
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
		
		void readInputLine(std::ifstream& inFile,std::string& fileLine,std::string headings[],int hdLdSpace[], std::string data[], int& dataLen);
		
		void readJob(std::string fileName);
		
		void readModelInput(std::string fileName);
		
		void readConstraintInput(std::string fileName);
		
		void readLoadInput(std::string fileName);
		
		void readInitialState(std::string fileName);
		
		void readDesVarInput(std::string fileName);
		
		void readObjectiveInput(std::string fileName);
		
		void readDesVarValues(std::string fileName);

		void readTimeStepSoln();
		
		// Analysis
		
		void reorderNodes(int blockDim);
		
		void updateReference();
		
		void findSurfaceFaces();
		
		void analysisPrep(int blockDim);
		
		void buildElasticAppLoad(double appLd[], double time);
		
		void buildElasticSolnLoad(double solnLd[], bool buildMat, bool dyn, bool nLGeom);
		
		void solveStep(JobCommand *cmd, double time, double appLdFact);
		
		void solve(JobCommand *cmd);
		
		void eigenSolve();
		
		void getObjGradient();
		
		// Output

		void writeTimeStepSoln();
		
		void writeNodeResults(std::string fileName, std::string nodeSet, StringList& fields, int timeStep);

		void writeElementResults(std::string fileName, std::string elSet, StringList& fields, int timeStep, NdPt ndAr[], DVPt dvAr[]);
		
		//
};

#endif