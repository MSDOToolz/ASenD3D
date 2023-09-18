#ifndef MODEL
#define MODEL
#include <fstream>
#include <iostream>
#include <map>
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
		SetPt* nsArray;
		std::map<std::string, int> nsMap;
		SetList elementSets;
		SetPt* esArray;
		std::map<std::string, int> esMap;
		SectionList sections;
		MaterialList materials;
		ConstraintList elasticConst;
		ConstraintList thermalConst;
		LoadList elasticLoads;
		LoadList thermalLoads;
		DesVarList designVars;
		DVPt *dVarArray;
		Objective objective;
		Job job;
		
		int elMatDim;
		int totGlobDof;
		bool anPrepRun;
		int timeStepsSaved;
		JobCommand* solveCmd;
		
		SparseMat elasticMat;
		LowerTriMat elasticLT;
		bool elasticScaled;
		SparseMat thermMat;
		LowerTriMat thermLT;
		bool thermScaled;
		
		double* tempV1;
		double* tempV2;
		double* tempV3;
		double* tempV4;
		double* tempV5;
		double* tempV6;

		Doub *tempD1;

		double* dLdU;
		double* dLdV;
		double* dLdA;
		double* dLdT;
		double* dLdTdot;
		double* uAdj;
		double* vAdj;
		double* aAdj;
		double* tAdj;
		double* tdotAdj;

		DiffDoub* dRudD;
		DiffDoub* dRtdD;

		double* dLdD;
	
	public:
	    Model();
		
		NodeList* getNodes();
		
		ElementList* getElements();
		
		SetList* getNodeSets();
		
		SetList* getElementSets();
		
		SectionList* getSections();
		
		MaterialList* getMaterials();
		
		ConstraintList* getElasticConstraints();
		
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

		void readTimeStepSoln(int tStep);
		
		// Analysis
		
		void reorderNodes(int blockDim);

		void buildConstraintMats();
		
		void updateReference();
		
		void findSurfaceFaces();
		
		void analysisPrep(int blockDim);
		
		void buildElasticAppLoad(double appLd[], double time);

		void buildThermalAppLoad(double appLd[], double time);
		
		void buildElasticSolnLoad(double solnLd[], bool buildMat, bool dyn, bool nLGeom);

		void scaleElasticConst();

		void scaleThermalConst();

		void buildElasticConstLoad(double constLd[]);
		
		void solveStep(JobCommand *cmd, double time, double appLdFact);
		
		void solve(JobCommand *cmd);
		
		void eigenSolve();

		void backupElastic();

		void restoreElastic();

		void setSolutionVec();

		void augmentdLdU();

		void solveForAdjoint();

		void dRelasticdD(int dVarNum);
		
		void getObjGradient();
		
		// Output

		void writeTimeStepSoln(int tStep);
		
		void writeNodeResults(std::string fileName, std::string nodeSet, StringList& fields, int timeStep);

		void writeElementResults(std::string fileName, std::string elSet, StringList& fields, int timeStep);

		void writeObjective(std::string fileName, StringList& includeFields, bool writeGrad);
		
		//
};

#endif